# Load libraries
library(shiny)
library(shinythemes)
library(shinyalert)
library(spsComps)
library(DT)
library(jsonlite)

HAS_SHINYWIDGETS <- requireNamespace("shinyWidgets", quietly = TRUE)
if (HAS_SHINYWIDGETS) {
  library(shinyWidgets)
}


# Shiny caps a request at 5 MB unless told otherwise, and never says so in a way
# that survives into a log -- the browser simply reports the upload as too large.
# A paired FASTQ run in this lab is 450-562 MB, so "Method of input: Upload" has
# been unusable for real data since the app was written: all 308 runs in
# docs/run-log.tsv came from storage1 or storage2, and not one from an upload.
#
# The staging area is the container's own disk (1.9 TB free), and the files are
# handed straight to rsync afterwards, so the ceiling only needs to be larger
# than a plausible run rather than carefully tuned.
# An uploaded file's name is chosen entirely by whoever uploads it, and it used
# to travel through goUpload.bash into a command run over ssh on the login node.
# A name containing a single quote closed the quoting there and the rest ran as
# user. Accept only what a FASTQ is actually called; the script re-checks,
# because neither layer should be the only thing standing there.
# The email box is free text and its value is handed to the mail command as a
# recipient. Unvalidated, that turns the server into something anyone can send
# mail from, wearing this institution's name. Empty is allowed -- the field is
# optional and an absent address simply means no notification is sent.
# A barcode file has a closed format: one line per well, a well label and a DNA
# barcode separated by a tab. Nothing else is a barcode file.
#
# That closedness is what makes this safe to lock down completely. The uploaded
# bytes are stored on the cluster and rsynced to lab storage, so accepting
# arbitrary content would make this a way to park anything anywhere in the lab's
# infrastructure, whatever the uploader called the file. Two things already limit
# the damage -- the file is never executed, and the original filename is
# discarded because the copy is always written as BC1.txt or BC2.txt, so hiding
# an extension achieves nothing -- but neither stops arbitrary bytes being
# stored and distributed.
#
# Validating the content rather than the name closes it: a file that must consist
# solely of well labels and ACGT cannot carry a payload. This also catches the
# far more common case of somebody uploading a spreadsheet, a Word file, or a
# FASTA by mistake.
BARCODE_MAX_BYTES <- 262144L    # a real set is ~2 KB; this is generous
BARCODE_MAX_LINES <- 4096L      # a 384-well set is 384

validate_barcode_file <- function(path, label = "Barcode file") {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(sprintf("%s is missing.", label))
  }
  size <- file.info(path)$size
  if (is.na(size) || size == 0) {
    return(sprintf("%s is empty.", label))
  }
  if (size > BARCODE_MAX_BYTES) {
    return(sprintf("%s is %.1f MB. A barcode file is a few kilobytes; this is not one.",
                   label, size / 1048576))
  }

  # Read as raw first: a text-only reader would happily accept a binary file
  # that merely lacks newlines, and NUL bytes are the clearest possible sign
  # that this is not a text file at all.
  raw_bytes <- readBin(path, what = "raw", n = min(size, BARCODE_MAX_BYTES))
  if (any(raw_bytes == as.raw(0))) {
    return(sprintf("%s is a binary file, not a barcode list.", label))
  }
  # Allow tab, newline, carriage return; reject other control characters.
  ctrl <- raw_bytes[raw_bytes < as.raw(32)]
  if (length(setdiff(as.integer(ctrl), c(9L, 10L, 13L)))) {
    return(sprintf("%s contains control characters and is not a barcode list.", label))
  }

  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) NULL)
  if (is.null(lines)) {
    return(sprintf("%s could not be read as text.", label))
  }
  lines <- trimws(lines[nzchar(trimws(lines))])
  if (!length(lines)) {
    return(sprintf("%s has no barcode lines.", label))
  }
  if (length(lines) > BARCODE_MAX_LINES) {
    return(sprintf("%s has %d lines. A barcode set has at most %d.",
                   label, length(lines), BARCODE_MAX_LINES))
  }

  ok <- grepl("^[A-Za-z0-9._-]{1,32}[ \t]+[ACGTNacgtn]{4,64}$", lines)
  if (!all(ok)) {
    bad <- head(which(!ok), 3)
    return(sprintf(
      "%s is not a barcode list. Line %s does not look like a well label and a DNA barcode: \"%s\"",
      label, paste(bad, collapse = ", "), substr(lines[bad[1]], 1, 40)))
  }

  dup <- duplicated(vapply(strsplit(lines, "[ \t]+"), `[`, character(1), 1))
  if (any(dup)) {
    return(sprintf("%s repeats a well label, which would make demultiplexing ambiguous.", label))
  }
  NULL   # valid
}

# Shiny inputs are NULL before their widget has rendered, which is not the same
# as empty. Without this, every readiness check would need its own is.null guard.
`%||%` <- function(a, b) if (is.null(a)) b else a

is_valid_email <- function(x) {
  if (is.null(x)) return(TRUE)
  x <- trimws(as.character(x))
  if (!length(x) || is.na(x[1]) || !nzchar(x[1])) return(TRUE)
  length(x) == 1 && nchar(x) <= 254 &&
    grepl("^[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\\.[A-Za-z]{2,}$", x)
}

is_valid_fastq_filename <- function(x) {
  is.character(x) && length(x) == 1 && !is.na(x) &&
    nzchar(x) && nchar(x) <= 255 &&
    grepl("^[A-Za-z0-9._-]+$", x) &&
    !grepl("^[.]", x) &&
    grepl("[.](fastq|fq)([.]gz)?$", x, ignore.case = TRUE)
}

TBRC_MAX_UPLOAD_MB <- {
  configured <- suppressWarnings(as.numeric(Sys.getenv("TBRC_MAX_UPLOAD_MB", unset = "")))
  if (is.na(configured) || configured <= 0) 8192 else configured
}
options(shiny.maxRequestSize = TBRC_MAX_UPLOAD_MB * 1024^2)

SERVER_CONFIG_FILE <- "server_config.json"
RUN_STATS_FILE <- "run_stats.tsv"
APP_DIR <- normalizePath(
  if (file.exists(SERVER_CONFIG_FILE)) "." else "frontend",
  winslash = "/",
  mustWork = TRUE
)
SERVER_CONFIG_FILE <- file.path(APP_DIR, SERVER_CONFIG_FILE)

# Code and state are separate concerns. In a container the app directory is
# mounted read-only and these three writable locations are volumes, so allow each
# to be relocated. Unset, they sit inside the app directory exactly as before.
app_state_dir <- function(env_name, default_subdir) {
  configured <- Sys.getenv(env_name, unset = "")
  path <- if (nzchar(configured)) configured else file.path(APP_DIR, default_subdir)
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

STATE_DIR <- app_state_dir("TBRC_STATE_DIR", ".")
RUNS_DIR  <- app_state_dir("TBRC_RUNS_DIR", "runs")
BC_DIR    <- app_state_dir("TBRC_BC_DIR", "bc")

RUN_STATS_FILE <- file.path(STATE_DIR, RUN_STATS_FILE)

format_duration_compact <- function(total_seconds) {
  total_seconds <- as.numeric(total_seconds)
  if (is.na(total_seconds) || total_seconds < 0) {
    return("not available")
  }

  hours <- total_seconds %/% 3600
  minutes <- (total_seconds %% 3600) %/% 60

  if (hours > 0) {
    sprintf("%d hours and %d minutes", hours, minutes)
  } else {
    sprintf("%d minutes", minutes)
  }
}

read_run_stats <- function() {
  if (!file.exists(RUN_STATS_FILE)) {
    return(data.frame())
  }

  tryCatch(
    read.table(RUN_STATS_FILE, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(...) data.frame()
  )
}

append_run_stat <- function(run_stat) {
  if (!is.data.frame(run_stat) || nrow(run_stat) != 1) {
    stop("append_run_stat expects a one-row data.frame", call. = FALSE)
  }

  # Appending a row with more columns than the file's header silently shifts
  # every field after the new one, and the damage is invisible until something
  # reads the file back. When the schema has grown, rewrite the whole file with
  # the union of columns instead. One row per run, so this stays cheap.
  if (file.exists(RUN_STATS_FILE)) {
    existing <- tryCatch(
      utils::read.table(RUN_STATS_FILE, sep = "\t", header = TRUE, quote = "\"",
                        stringsAsFactors = FALSE, check.names = FALSE,
                        colClasses = "character", fill = TRUE),
      error = function(...) NULL
    )
    if (!is.null(existing) && !identical(names(existing), names(run_stat))) {
      for (column in setdiff(names(run_stat), names(existing))) existing[[column]] <- ""
      for (column in setdiff(names(existing), names(run_stat))) run_stat[[column]] <- ""
      merged <- rbind(existing[, names(run_stat), drop = FALSE],
                      as.data.frame(lapply(run_stat, as.character),
                                    stringsAsFactors = FALSE, check.names = FALSE))
      write.table(merged, file = RUN_STATS_FILE, sep = "\t", quote = TRUE,
                  row.names = FALSE, col.names = TRUE, append = FALSE)
      return(invisible(NULL))
    }
  }

  write.table(
    run_stat,
    file = RUN_STATS_FILE,
    sep = "\t",
    quote = TRUE,
    row.names = FALSE,
    col.names = !file.exists(RUN_STATS_FILE),
    append = file.exists(RUN_STATS_FILE)
  )
}

# ------------------------------------------------------------------ costing --
# What a library cost: the sequencing run off the core facility's rate card,
# plus what the lab spent on its own bench preparing the barcode plates.
#
# Both halves are estimates and are meant to be edited rather than trusted. The
# rate card is transcribed and goes stale; the bench figure is list prices for
# reagents bought at negotiated ones, and it ignores everything that is not
# consumed per plate -- instrument time, salary, the plates that failed. It is
# useful for "is this a $2,000 experiment or a $6,000 one", not for a budget.

cost_config <- function(config) {
  costs <- config$costs
  if (is.null(costs)) list() else costs
}

cost_sequencing_choices <- function(config) {
  seq <- cost_config(config)$sequencing
  if (is.null(seq) || !length(seq)) return(character(0))
  labels <- vapply(seq, function(entry) {
    label <- entry$label
    if (is.null(label) || !nzchar(label)) "" else as.character(label)
  }, character(1))
  stats::setNames(names(seq), ifelse(nzchar(labels), labels, names(seq)))
}

# Three pricing shapes now, because the price lists have three.
#
#   flat        one charge per flowcell, from the core facility rate card
#   per_sample  a fixed charge for each sample submitted
#   custom_gb   Plasmidsaurus Custom Sequencing: a base for the first Gb, a
#               smaller charge per additional Gb, a per-sample multiplexing
#               charge, and the whole thing repeating for every shipped batch.
#               Batching is the part people get wrong -- two shipments of one
#               library each cost a full base price twice.
cost_sequencing_price <- function(config, key, user_type = "ru_user",
                                  gigabases = 1, samples = 1, batches = 1) {
  entry <- cost_config(config)$sequencing[[key]]
  if (is.null(entry)) return(NA_real_)
  num <- function(field, fallback = NA_real_) {
    value <- suppressWarnings(as.numeric(entry[[field]]))
    if (length(value) != 1 || is.na(value)) fallback else value
  }
  mode <- as.character(entry$pricing %||% "flat")

  batches <- max(1, floor(suppressWarnings(as.numeric(batches))[1] %||% 1))
  samples <- max(1, floor(suppressWarnings(as.numeric(samples))[1] %||% 1))
  gigabases <- suppressWarnings(as.numeric(gigabases))[1]
  if (!is.finite(gigabases) || gigabases < 1) gigabases <- 1

  if (identical(mode, "per_sample")) {
    per <- num("price_per_sample")
    return(if (is.na(per)) NA_real_ else per * samples)
  }

  if (identical(mode, "custom_gb")) {
    base <- num("base_price_first_gb")
    extra <- num("price_per_additional_gb", 0)
    per_sample <- num("price_per_sample", 0)
    if (is.na(base)) return(NA_real_)
    # Data volume is what actually arrives, so it is split across batches; the
    # base price and the per-sample charge are not.
    gb_per_batch <- gigabases / batches
    additional <- max(0, ceiling(gb_per_batch) - 1)
    return(batches * (base + extra * additional) + per_sample * samples)
  }

  num(user_type)
}

# One row per consumable so the breakdown can say where the money went, rather
# than presenting a single number nobody can check.
cost_prep_lines <- function(config) {
  items <- cost_config(config)$library_prep_per_plate$consumables
  if (is.null(items) || !length(items)) return(NULL)
  rows <- lapply(items, function(item) {
    unit <- suppressWarnings(as.numeric(item$unit_cost))
    per_plate <- suppressWarnings(as.numeric(item$units_per_plate))
    if (length(unit) != 1 || is.na(unit)) unit <- 0
    if (length(per_plate) != 1 || is.na(per_plate)) per_plate <- 0
    data.frame(
      item = as.character(item$item %||% "(unnamed)"),
      note = as.character(item$note %||% ""),
      cost = unit * per_plate,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

cost_prep_per_plate <- function(config) {
  lines <- cost_prep_lines(config)
  if (is.null(lines) || !nrow(lines)) return(NA_real_)
  sum(lines$cost, na.rm = TRUE)
}

estimate_run_cost <- function(config, sequencing_key = NULL, plates = 1,
                              user_type = NULL, gigabases = 1, samples = 1,
                              batches = 1) {
  costs <- cost_config(config)
  if (is.null(user_type) || !nzchar(user_type)) {
    user_type <- as.character(costs$default_user_type %||% "ru_user")
  }
  plates <- suppressWarnings(as.numeric(plates))
  if (length(plates) != 1 || !is.finite(plates) || plates < 1) plates <- 1

  prep_each <- cost_prep_per_plate(config)
  prep_total <- if (is.na(prep_each)) NA_real_ else prep_each * plates

  if (is.null(sequencing_key) || !nzchar(sequencing_key)) {
    sequencing_key <- as.character(costs$default_sequencing %||% "")
  }
  entry <- costs$sequencing[[sequencing_key]]
  seq_total <- cost_sequencing_price(config, sequencing_key, user_type,
                                     gigabases = gigabases, samples = samples,
                                     batches = batches)

  list(
    sequencing = seq_total,
    sequencing_label = as.character(entry$label %||% sequencing_key),
    pricing = as.character(entry$pricing %||% "flat"),
    prep_per_plate = prep_each,
    prep_total = prep_total,
    plates = plates,
    user_type = user_type,
    total = if (is.na(seq_total) && is.na(prep_total)) NA_real_
            else sum(c(seq_total, prep_total), na.rm = TRUE),
    complete = !is.na(seq_total) && !is.na(prep_total)
  )
}

# What one barcode plate costs once the flowcell is shared. This is the number
# that answers "is this worth doing instead of Sanger": the prep is spent per
# plate whatever happens, but the sequencing charge divides across the plates
# pooled onto the run.
cost_per_plate <- function(config, cost, pooled = NULL) {
  pooling <- cost_config(config)$pooling
  if (is.null(pooled)) {
    pooled <- suppressWarnings(as.numeric(pooling$plates_pooled_per_flowcell))
  }
  if (length(pooled) != 1 || !is.finite(pooled) || pooled < 1) pooled <- 1
  if (!identical(cost$pricing, "flat")) {
    # Per-sample and per-Gb services are already charged per library, so there
    # is nothing to divide.
    share <- if (is.finite(cost$sequencing)) cost$sequencing / max(1, cost$plates) else NA_real_
  } else {
    share <- if (is.finite(cost$sequencing)) cost$sequencing / pooled else NA_real_
  }
  prep <- cost$prep_per_plate
  list(
    sequencing_share = share,
    prep = prep,
    total = if (is.na(share) && is.na(prep)) NA_real_ else sum(c(share, prep), na.rm = TRUE),
    pooled = pooled
  )
}

cost_sanger_equivalent <- function(config, wells = 96) {
  per <- suppressWarnings(as.numeric(cost_config(config)$pooling$sanger_per_reaction))
  if (length(per) != 1 || !is.finite(per)) return(NA_real_)
  per * wells
}

# Every distinct library the log knows about, costed once.
#
# Counted by name rather than by row: re-running a library did not buy another
# flowcell. Plate counts come from what the backend counted in the delivered
# wells; libraries that predate that measurement fall back to one plate, which
# makes the total a floor and the card says so.
cost_all_runs <- function(config, log, sequencing_key = NULL, plates_by_run = NULL) {
  if (is.null(log) || !nrow(log)) return(NULL)
  names_seen <- trimws(tolower(as.character(log$run_name)))
  names_seen <- unique(names_seen[nzchar(names_seen)])
  if (!length(names_seen)) return(NULL)

  plates <- rep(NA_real_, length(names_seen))
  names(plates) <- names_seen
  if (!is.null(plates_by_run) && length(plates_by_run)) {
    known <- names(plates_by_run)[names(plates_by_run) %in% names_seen]
    plates[known] <- suppressWarnings(as.numeric(plates_by_run[known]))
  }
  measured <- sum(is.finite(plates) & plates > 0)
  plates[!is.finite(plates) | plates <= 0] <- 1

  prep_each <- cost_prep_per_plate(config)
  pooling <- cost_config(config)$pooling
  pooled <- suppressWarnings(as.numeric(pooling$plates_pooled_per_flowcell))
  if (length(pooled) != 1 || !is.finite(pooled) || pooled < 1) pooled <- 1

  one <- estimate_run_cost(config, sequencing_key = sequencing_key, plates = 1)
  total_plates <- sum(plates)

  # A flowcell is bought once and shared by the plates pooled onto it, so its
  # cost scales with plates. Plasmidsaurus Custom is billed per project, so it
  # scales with libraries. Multiplying the wrong one by the wrong count is how
  # this number would quietly become nonsense.
  if (identical(one$pricing, "flat")) {
    share <- if (is.finite(one$sequencing)) one$sequencing / pooled else NA_real_
    seq_total <- if (is.finite(share)) share * total_plates else NA_real_
  } else {
    share <- if (is.finite(one$sequencing)) one$sequencing else NA_real_
    seq_total <- if (is.finite(share)) share * length(names_seen) else NA_real_
  }
  prep_total <- if (is.finite(prep_each)) prep_each * total_plates else NA_real_

  list(
    distinct_runs = length(names_seen),
    rows = nrow(log),
    total_plates = total_plates,
    measured_runs = measured,
    pricing = one$pricing,
    per_plate = if (identical(one$pricing, "flat"))
      sum(c(share, prep_each), na.rm = TRUE) else NA_real_,
    sequencing = seq_total,
    prep = prep_total,
    pooled = pooled,
    total = if (is.na(seq_total) && is.na(prep_total)) NA_real_
            else sum(c(seq_total, prep_total), na.rm = TRUE)
  )
}

format_money <- function(value, currency = "USD") {
  if (length(value) != 1 || !is.finite(value)) return("not available")
  prefix <- if (identical(currency, "USD")) "$" else paste0(currency, " ")
  paste0(prefix, formatC(round(value), format = "d", big.mark = ","))
}

get_metric_value <- function(config, field_name, fallback) {
  value <- suppressWarnings(as.numeric(config$metrics[[field_name]]))
  if (length(value) != 1 || is.na(value)) fallback else value
}

compute_run_metrics <- function(
  run_stats,
  input_bytes = NA_real_,
  manual_minutes_per_dataset = 6,
  manual_minutes_per_100mb = 29,
  manual_hours_per_run = 0,
  run_index_offset = 0
) {
  if (!is.data.frame(run_stats) || nrow(run_stats) == 0) {
    return(list(
      dataset_count = 0,
      run_index = run_index_offset + 1,
      total_time_saved_sec = 0,
      estimated_runtime_sec = NA_real_,
      has_runtime_model = FALSE
    ))
  }

  completed_mask <- rep(TRUE, nrow(run_stats))
  if ("status" %in% names(run_stats)) {
    completed_mask <- tolower(trimws(as.character(run_stats$status))) == "completed"
  }

  completed_runs <- run_stats[completed_mask, , drop = FALSE]
  dataset_count <- nrow(completed_runs)

  # Counts datasets successfully processed, not submissions attempted. The
  # historical offset is derived from result archives on the storage servers, and
  # a run that failed never produced one -- so attempts are unknowable before this
  # ledger existed and completed runs are the only consistent basis.
  run_index <- dataset_count + run_index_offset + 1
  valid_saved_size <- suppressWarnings(as.numeric(completed_runs$input_bytes))
  saved_size_keep <- is.finite(valid_saved_size) & valid_saved_size > 0
  # Two components: one that scales with how much data was processed, and a flat
  # cost per run for the setting-up, copy-pasting and waiting that does not.
  manual_seconds_per_byte <- (manual_minutes_per_100mb * 60) / (100 * 1024 * 1024)
  total_time_saved_sec <- sum(valid_saved_size[saved_size_keep] * manual_seconds_per_byte, na.rm = TRUE) +
    dataset_count * manual_hours_per_run * 3600

  if (!any(saved_size_keep)) {
    total_time_saved_sec <- dataset_count * manual_minutes_per_dataset * 60
  }

  if (dataset_count == 0) {
    return(list(
      dataset_count = 0,
      run_index = run_index,
      total_time_saved_sec = 0,
      estimated_runtime_sec = NA_real_,
      has_runtime_model = FALSE
    ))
  }

  valid_runtime <- suppressWarnings(as.numeric(completed_runs$runtime_sec))
  valid_size <- suppressWarnings(as.numeric(completed_runs$input_bytes))
  keep <- is.finite(valid_runtime) & valid_runtime > 0 & is.finite(valid_size) & valid_size > 0

  if (!any(keep)) {
    estimated_runtime_sec <- mean(valid_runtime[is.finite(valid_runtime) & valid_runtime > 0], na.rm = TRUE)
    if (is.nan(estimated_runtime_sec)) {
      estimated_runtime_sec <- NA_real_
    }
    return(list(
      dataset_count = dataset_count,
      run_index = run_index,
      total_time_saved_sec = total_time_saved_sec,
      estimated_runtime_sec = estimated_runtime_sec,
      has_runtime_model = FALSE
    ))
  }

  seconds_per_mb <- sum(valid_runtime[keep], na.rm = TRUE) / sum(valid_size[keep] / (1024 * 1024), na.rm = TRUE)
  target_input_bytes <- suppressWarnings(as.numeric(input_bytes))
  if (is.na(target_input_bytes) || !is.finite(target_input_bytes) || target_input_bytes <= 0) {
    typical_size_mb <- median(valid_size[keep] / (1024 * 1024), na.rm = TRUE)
  } else {
    typical_size_mb <- target_input_bytes / (1024 * 1024)
  }
  estimated_runtime_sec <- seconds_per_mb * typical_size_mb

  list(
    dataset_count = dataset_count,
    run_index = run_index,
    total_time_saved_sec = total_time_saved_sec,
    estimated_runtime_sec = estimated_runtime_sec,
    has_runtime_model = TRUE
  )
}

get_storage_entry <- function(config, server_name) {
  entry <- config$storage_server[[server_name]]

  if (is.null(entry)) {
    return(NULL)
  }

  if (is.character(entry)) {
    return(list(
      rawdata_path = unname(entry),
      results_path = sub("rawdata/?$", "results/", unname(entry))
    ))
  }

  entry
}

get_storage_value <- function(config, server_name, field_name) {
  entry <- get_storage_entry(config, server_name)

  if (is.null(entry)) {
    return(NA_character_)
  }

  value <- entry[[field_name]]
  if (is.null(value) || identical(value, "")) {
    return(NA_character_)
  }

  unname(value)
}

get_cluster_value <- function(config, field_name) {
  value <- config$cluster[[field_name]]

  if (is.null(value) || identical(value, "")) {
    return(NA_character_)
  }

  unname(as.character(value))
}

get_cluster_value_or <- function(config, field_name, fallback) {
  value <- get_cluster_value(config, field_name)
  if (is.na(value) || identical(value, "")) {
    return(fallback)
  }
  value
}

# User name, folder name and barcode folder name all become path segments on up
# to three machines and are interpolated into remote shell commands. Keep them to
# a conservative character set; the shell scripts re-check independently.
VALID_PATH_SEGMENT <- "^[A-Za-z0-9._-]+$"

is_valid_path_segment <- function(value) {
  if (is.null(value) || length(value) < 1) {
    return(FALSE)
  }

  value <- trimws(as.character(value[[1]]))

  nzchar(value) &&
    grepl(VALID_PATH_SEGMENT, value) &&
    !value %in% c(".", "..") &&
    !startsWith(value, "-")
}

invalid_segment_message <- function(label) {
  sprintf(
    "%s may only contain letters, numbers, dot, underscore and hyphen \u2014 no spaces, slashes or quotes.",
    label
  )
}

# Smallest number of positions in which any two barcodes in a file differ.
# The demultiplexer resolves a tie by file order without warning, so a set can
# only tolerate m sequencing errors while its minimum distance is at least 2m+1.
barcode_min_distance <- function(path) {
  rows <- tryCatch(read.table(path, stringsAsFactors = FALSE), error = function(...) NULL)
  if (is.null(rows) || ncol(rows) < 2 || nrow(rows) < 2) {
    return(NA_integer_)
  }

  codes <- toupper(trimws(as.character(rows[[2]])))
  codes <- codes[nzchar(codes)]
  if (length(codes) < 2) {
    return(NA_integer_)
  }

  best <- Inf
  for (i in seq_len(length(codes) - 1)) {
    for (j in seq(i + 1, length(codes))) {
      a <- codes[[i]]; b <- codes[[j]]
      if (nchar(a) != nchar(b)) next
      d <- sum(strsplit(a, "")[[1]] != strsplit(b, "")[[1]])
      if (d < best) best <- d
    }
  }

  if (is.infinite(best)) NA_integer_ else as.integer(best)
}

# The 5' handle differs by chain, and getting it wrong is silent: cutadapt
# --discard-untrimmed simply deletes everything. Look it up from the barcode set
# instead of relying on the operator to remember.
barcode_default <- function(config, barcode_folder, field) {
  entry <- config$barcode_defaults[[barcode_folder]]
  if (is.null(entry)) {
    return(NA_character_)
  }
  value <- entry[[field]]
  if (is.null(value) || !nzchar(as.character(value))) NA_character_ else as.character(value)
}

app_script <- function(script_name) {
  file.path(APP_DIR, "scripts", script_name)
}

app_path <- function(...) {
  file.path(APP_DIR, ...)
}

app_radio_buttons <- function(inputId, label, choices = NULL, selected = NULL, ...) {
  if (HAS_SHINYWIDGETS) {
    shinyWidgets::prettyRadioButtons(inputId, label, choices = choices, selected = selected, ...)
  } else {
    shiny::radioButtons(inputId, label, choices = choices, selected = selected)
  }
}

app_update_radio_buttons <- function(session, inputId, choices = NULL, selected = NULL) {
  if (HAS_SHINYWIDGETS) {
    shinyWidgets::updatePrettyRadioButtons(session, inputId = inputId, choices = choices, selected = selected)
  } else {
    shiny::updateRadioButtons(session, inputId = inputId, choices = choices, selected = selected)
  }
}

app_picker_input <- function(inputId, label, choices = NULL, selected = NULL, options = NULL, ...) {
  if (HAS_SHINYWIDGETS) {
    shinyWidgets::pickerInput(inputId, label, choices = choices, selected = selected, options = options, ...)
  } else {
    shiny::selectInput(inputId, label, choices = choices, selected = selected, selectize = TRUE)
  }
}

app_update_picker_input <- function(session, inputId, choices = NULL, selected = NULL) {
  if (HAS_SHINYWIDGETS) {
    shinyWidgets::updatePickerInput(session, inputId = inputId, choices = choices, selected = selected)
  } else {
    shiny::updateSelectInput(session, inputId = inputId, choices = choices, selected = selected)
  }
}

app_switch_input <- function(inputId, label, value = FALSE, ...) {
  shiny::checkboxInput(inputId, label, value)
}

app_action_button <- function(inputId, label, ...) {
  if (HAS_SHINYWIDGETS) {
    shinyWidgets::actionBttn(inputId, label, ...)
  } else {
    shiny::actionButton(inputId, label)
  }
}

# "Run on this server" only works where a sibling backend/ directory and the
# workflow tools actually exist. On a Shiny host that received only frontend/,
# neither is true -- which is why the toggle produced 7 raw shell errors for
# Sandra and CMu before it was gated.
LOCAL_BACKEND_DIR <- normalizePath(file.path(APP_DIR, "..", "backend"), winslash = "/", mustWork = FALSE)

local_backend_status <- function() {
  if (!dir.exists(LOCAL_BACKEND_DIR)) {
    return(list(available = FALSE, reason = "no backend/ directory beside the app"))
  }

  required_tools <- c("snakemake", "pandaseq", "cutadapt", "fastx_trimmer")
  missing_tools <- required_tools[!nzchar(Sys.which(required_tools))]

  if (length(missing_tools) > 0) {
    return(list(
      available = FALSE,
      reason = paste("missing on this host:", paste(missing_tools, collapse = ", "))
    ))
  }

  list(available = TRUE, reason = "")
}

LOCAL_BACKEND <- local_backend_status()

# Read JSON file and extract server names
servers_config <- fromJSON(SERVER_CONFIG_FILE, simplifyVector = FALSE)

# Extract server names
server_names <- names(servers_config$storage_server)

igblast_species_choices <- c(
  "Human" = "human",
  "Mouse" = "mouse"
)

# ig_light is no longer offered. ig_all carries the same 168 IGKV and 19 IGLV
# germlines, and a light-only panel on a plate that also holds heavy chains
# forces every heavy sequence onto an IGKV gene (CMP15_kappa: 533 of 1,425).
# The backend still accepts it from an older sample.json, and the receptor
# check widens it to ig_all when the plate warrants.
igblast_panel_choices <- c(
  "Immunoglobulin heavy chain only (IgH)" = "igh",
  "All immunoglobulins (IgH + IgK + IgL)" = "ig_all",
  "TCR alpha only (TRA)" = "tra",
  "TCR beta only (TRB)" = "trb",
  "All TCR (TRA + TRB + TRD + TRG)" = "tcr_all",
  "All TCR and immunoglobulins" = "all_receptors",
  "Mouse B6 immunoglobulin germline (IgH V/D/J + C)" = "ig_b6"
)

# The run log tab lives in its own file; sourcing it here keeps ui.R declarative.
local({
  runlog_module <- app_path("tab_runlog.R")
  if (file.exists(runlog_module)) source(runlog_module, local = FALSE)
  changelog_module <- app_path("tab_changelog.R")
  if (file.exists(changelog_module)) source(changelog_module, local = FALSE)
})
