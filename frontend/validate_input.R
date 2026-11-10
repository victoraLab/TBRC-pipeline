# The one place that decides whether a submission is allowed to run.
#
# Shiny does not enforce what a widget offers: a browser can send any value for
# any input id, in any shape. Until this check existed, several of those values
# (the trim sequence, the pipeline, the barcode set, the user name) were pasted
# into a shell, first on this host as root and later on the cluster. So every
# value the Submit form sends goes through validate_submission() before anything
# else runs. The check refuses a bad value; it does not hand back what to use.
# `clean` is not read by submission_row() or anything else that writes sample.json,
# which build their values from input$ directly, so each of them has to trim for
# itself (the user name through clean_user_id(), the email at its two places in
# server.R). The cluster refuses a value that still has padding on it.
#
# The rules are the ones in docs/superpowers/specs/2026-11-08-tbrc-shell-
# injection-fix-design.md, section 1. The cluster repeats them when it loads
# sample.json; tests/validation_cases.json is shared by both sets of tests.
#
# Sourced from global.R, so ui.R can use the choice constants below. It relies on
# is_valid_email(), is_path_segment_text(), invalid_segment_message() and
# igblast_panel_choices from global.R, all looked up when the check runs.

# ------------------------------------------------------------------ choices
# One definition each, used by the widget in ui.R and by the check below, so
# what the form offers and what the server accepts cannot drift apart. The names
# are the labels people see; the values are what the form sends.

METHOD_CHOICES <- c(
  "vsearch cluster assembly (tiago_vsearch)" = "tiago_vsearch",
  "Classic clonality assembly (classic_luka)" = "classic_luka",
  "Nanopore long-read assembly (tiago_nanopore)" = "tiago_nanopore",
  "MIBS / sulkseq multiwell incidence (tiago_sulkseq)" = "tiago_sulkseq"
)

ARCHIVE_SCOPE_CHOICES <- c(
  "Light (recommended)" = "light",
  "Full" = "full"
)

ARCHIVE_FORMAT_CHOICES <- c(
  "ZIP (Recommended)" = "zip",
  "tar.gz" = "tar.gz"
)

IGBLAST_SPECIES_CHOICES <- c(
  "Human" = "human",
  "Mouse" = "mouse"
)

RECEPTOR_CHECK_CHOICES <- c(
  "Correct it, and say so in the report" = "auto",
  "Report it, but annotate as chosen" = "warn",
  "Do not check" = "off"
)

LUKA_LIGHT_CHOICES <- c(
  "Standard preprocessing" = "standard",
  "Light-chain rescue preprocessing" = "luka_light"
)

MIBS_GENE_GROUPING_CHOICES <- c(
  "Strict (a shared gene call is required)" = "strict",
  "Transitive (a chain of ambiguous calls may link)" = "transitive"
)

MIBS_CLONE_SCOPE_CHOICES <- c(
  "Within each sample" = "sample",
  "Across the whole run, then attributed to a sample (as published)" = "global"
)

# What the form puts in the folder box while "Upload" is chosen, so the box does
# not look like it wants an answer. ui and server use this one string.
UPLOAD_FOLDER_PLACEHOLDER <- "Not used when uploading..."

# The folder radio's value when a run comes from uploaded files.
UPLOAD_SERVER_NAME <- "upload"

# Same bounds as the numeric boxes in ui.R and as the Snakefile. numericInput's
# min and max are advisory, so they have to be enforced here.
MIBS_INPUT_RANGES <- list(
  mibs_min_reads_per_well    = c(1, 10000),
  mibs_min_wells_per_clone   = c(1, 96),
  mibs_junction_max_distance = c(0, 30),
  mibs_min_reads_per_cluster = c(1, 1000)
)

N_READ_RANGE <- c(1, 500)

TRIM_SEQUENCE_PATTERN <- "^[ACGTRYSWKMBDHVN]{10,64}\\z"

# The names refusals use. They are the labels on the form, so the person who
# reads "Trim sequence must be ..." can find the box.
SUBMISSION_FIELD_LABELS <- c(
  user_id = "User name",
  email = "User email",
  type_input = "Method of input",
  server_name = "Server to send/get results",
  input_folder = "Folder in server",
  primer = "Barcode set",
  sequencing_kit = "Sequencing run",
  keep_intermediate = "Keep intermediate files on the compute host",
  archive_scope = "Results archive contents",
  run_group = "Results subfolder on the storage server",
  trim_primers = "Trim primers",
  run_local = "Run on this server (not HPC)",
  run_igblast = "Run IGBlast annotation",
  run_clonality = "Run clonality analysis from IGBlast",
  archive_format = "Final results archive format",
  igblast_species = "IGBlast species",
  igblast_panel = "IGBlast receptor scope",
  receptor_check = "If the sequences disagree with that scope",
  trim_sequence = "Trim sequence",
  n_read = "Contigs or isoforms to retrieve",
  method = "Core assembly pipeline",
  luka_light = "Barcode preprocessing mode",
  mibs_min_reads_per_well = "Minimum reads for a clone in a well",
  mibs_min_wells_per_clone = "Minimum wells a clone must appear in",
  mibs_junction_max_distance = "Junction Hamming distance for one clone",
  mibs_min_reads_per_cluster = "Minimum reads to keep a sequence in a well",
  mibs_gene_grouping = "How gene calls decide clone membership",
  mibs_clone_scope = "Where clones are assigned",
  mibs_productive_only = "Use productive rearrangements only"
)

# The input ids validate_submission() looks at.
SUBMISSION_FIELD_IDS <- names(SUBMISSION_FIELD_LABELS)

# What each choice field may hold. server_names are the keys of
# server_config.json; kit_names are the keys of its costs$sequencing table.
submission_choices <- function(server_names, kit_names) {
  list(
    method = unname(METHOD_CHOICES),
    archive_scope = unname(ARCHIVE_SCOPE_CHOICES),
    archive_format = unname(ARCHIVE_FORMAT_CHOICES),
    igblast_species = unname(IGBLAST_SPECIES_CHOICES),
    igblast_panel = unname(igblast_panel_choices),
    luka_light = unname(LUKA_LIGHT_CHOICES),
    receptor_check = unname(RECEPTOR_CHECK_CHOICES),
    mibs_gene_grouping = unname(MIBS_GENE_GROUPING_CHOICES),
    mibs_clone_scope = unname(MIBS_CLONE_SCOPE_CHOICES),
    keep_intermediate = c("TRUE", "FALSE"),
    trim_primers = c("TRUE", "FALSE"),
    type_input = c("Upload", "Server"),
    server_name = as.character(unname(server_names)),
    sequencing_kit = c("", as.character(unname(kit_names)))
  )
}

# ------------------------------------------------------------------ cleaning

# The form's trim sequence box takes any case and stray spaces; the pipeline
# wants upper case with none.
clean_trim_sequence <- function(value) toupper(trimws(value))

# ------------------------------------------------------------------ helpers

# One piece of text: a single non-NA string with nothing outside plain ASCII.
# Case-mapping would otherwise turn some non-ASCII letters into ASCII ones
# (toupper() makes U+017F an S), and a pattern in a UTF-8 locale can match more
# than it reads as. Anything that is not plain ASCII is not something a
# barcode set, user, folder or sequence can be called here.
vs_text <- function(value) {
  if (!is.character(value) || length(value) != 1L || is.na(value)) return(NULL)
  non_ascii <- tryCatch(
    grepl("[^\\x01-\\x7F]", value, perl = TRUE, useBytes = TRUE),
    error = function(e) TRUE, warning = function(w) TRUE
  )
  if (non_ascii) return(NULL)
  value
}

# TRUE or FALSE; never NA, never an error.
vs_matches <- function(pattern, value) {
  isTRUE(tryCatch(grepl(pattern, value, perl = TRUE),
                  error = function(e) FALSE, warning = function(w) FALSE))
}

# ------------------------------------------------------------------ the check

# `fields` is a named list of the raw input values, `choices` is what
# submission_choices() returns. The result is
#   ok      TRUE only if every field passed
#   errors  one plain sentence per refused field, each naming the field on the
#           form and never repeating the value that was sent
#   clean   the values to use from here on (user and trim sequence trimmed, the
#           sequence upper-cased, truth values as real TRUE/FALSE). Read it
#           only when ok is TRUE.
# It never raises: a form it cannot make sense of is a refusal.
validate_submission <- function(fields, choices) {
  errors <- character(0)
  clean <- list()

  label <- function(id) SUBMISSION_FIELD_LABELS[[id]]
  refuse <- function(text) errors <<- c(errors, text)
  value_of <- function(id) fields[[id]]

  segment <- function(id, trim = FALSE, optional = FALSE) {
    value <- value_of(id)
    if (is.null(value) && optional) { clean[[id]] <<- ""; return(invisible()) }
    if (is.null(value)) { refuse(sprintf("%s is required.", label(id))); return(invisible()) }
    text <- vs_text(value)
    if (!is.null(text) && trim) text <- trimws(text)
    if (!is.null(text) && optional && !nzchar(text)) { clean[[id]] <<- ""; return(invisible()) }
    if (is.null(text) || !is_path_segment_text(text)) {
      refuse(invalid_segment_message(label(id)))
    } else {
      clean[[id]] <<- text
    }
    invisible()
  }

  pick <- function(id, allowed, optional = FALSE, default = "",
                   otherwise = "%s must be one of the options the form offers.") {
    value <- value_of(id)
    if (is.null(value)) {
      if (optional) clean[[id]] <<- default else refuse(sprintf("%s is required.", label(id)))
      return(invisible())
    }
    text <- vs_text(value)
    if (!is.null(text) && text %in% allowed) {
      clean[[id]] <<- text
    } else {
      refuse(sprintf(otherwise, label(id)))
    }
    invisible()
  }

  flag <- function(id, optional = FALSE) {
    value <- value_of(id)
    if (is.null(value)) {
      if (optional) clean[[id]] <<- FALSE else refuse(sprintf("%s is required.", label(id)))
      return(invisible())
    }
    if (is.logical(value) && length(value) == 1L && !is.na(value)) {
      clean[[id]] <<- as.logical(value)
    } else {
      refuse(sprintf("%s must be on or off.", label(id)))
    }
    invisible()
  }

  whole_number <- function(id, limits, optional = FALSE) {
    value <- value_of(id)
    empty <- is.null(value) ||
      (length(value) == 1L && !is.character(value) && !is.list(value) && is.na(value))
    if (empty) {
      if (!optional) refuse(sprintf("%s must be a whole number from %d to %d.",
                                    label(id), as.integer(limits[1]), as.integer(limits[2])))
      return(invisible())
    }
    if (is.numeric(value) && length(value) == 1L && is.finite(value) &&
        value == floor(value) && value >= limits[1] && value <= limits[2]) {
      clean[[id]] <<- as.integer(value)
    } else {
      refuse(sprintf("%s must be a whole number from %d to %d.",
                     label(id), as.integer(limits[1]), as.integer(limits[2])))
    }
    invisible()
  }

  checked <- tryCatch({
    if (!is.list(fields)) stop("not a list")

    pick("type_input", choices$type_input)
    upload <- identical(clean$type_input, "Upload")

    segment("user_id", trim = TRUE)

    email <- value_of("email")
    if (is.null(email)) {
      clean$email <- ""
    } else if (!is.null(vs_text(email)) && is_valid_email(email)) {
      clean$email <- trimws(email)
    } else {
      refuse(sprintf("%s must be a valid email address, or empty to skip the notification.",
                     label("email")))
    }

    # Upload mode does not use the folder: the form fills the box with a note.
    # Anything else the box holds is still checked as a folder name.
    folder <- value_of("input_folder")
    if (upload && (is.null(folder) || identical(folder, UPLOAD_FOLDER_PLACEHOLDER))) {
      clean$input_folder <- ""
    } else {
      segment("input_folder", trim = TRUE, optional = upload)
    }

    # The radio keeps offering "upload" after the switch back to Server, so a
    # leftover choice of it is a plausible slip and gets its own wording.
    pick("server_name", c(choices$server_name, if (upload) UPLOAD_SERVER_NAME),
         otherwise = paste("%s must be one of the storage servers",
                           "(\"upload\" is only for uploaded files)."))
    segment("run_group", trim = TRUE, optional = TRUE)

    # Not membership in bc_names.txt: LongMS_IgA exists only on the cluster and
    # has real runs.
    segment("primer")
    pick("sequencing_kit", choices$sequencing_kit, optional = TRUE)
    pick("keep_intermediate", choices$keep_intermediate)
    pick("archive_scope", choices$archive_scope)
    pick("trim_primers", choices$trim_primers)
    pick("archive_format", choices$archive_format)
    pick("igblast_species", choices$igblast_species)
    pick("igblast_panel", choices$igblast_panel)
    pick("receptor_check", choices$receptor_check, optional = TRUE, default = "auto")
    pick("method", choices$method)
    pick("luka_light", choices$luka_light)

    flag("run_igblast")
    flag("run_clonality")
    flag("run_local", optional = TRUE)

    # The sequence is checked after cleaning, and only the cleaned one is ever
    # used. It is empty only when the trim is off.
    trim_on <- identical(clean$trim_primers, "TRUE")
    sequence <- value_of("trim_sequence")
    text <- vs_text(sequence)
    cleaned <- if (is.null(text)) NULL else clean_trim_sequence(text)
    if (is.null(sequence) || (!is.null(cleaned) && !nzchar(cleaned))) {
      if (trim_on) {
        refuse(sprintf("%s is required while Trim primers is on.", label("trim_sequence")))
      } else {
        clean$trim_sequence <- ""
      }
    } else if (!is.null(cleaned) && vs_matches(TRIM_SEQUENCE_PATTERN, cleaned)) {
      clean$trim_sequence <- cleaned
    } else {
      refuse(sprintf("%s must be 10 to 64 letters from A C G T R Y S W K M B D H V N.",
                     label("trim_sequence")))
    }

    whole_number("n_read", N_READ_RANGE)

    # These boxes sit in a panel that is hidden for every other pipeline, so an
    # emptied one matters only to the run that reads it.
    sulkseq <- identical(clean$method, "tiago_sulkseq")
    for (id in names(MIBS_INPUT_RANGES)) {
      whole_number(id, MIBS_INPUT_RANGES[[id]], optional = !sulkseq)
    }
    pick("mibs_gene_grouping", choices$mibs_gene_grouping, optional = !sulkseq,
         default = "strict")
    pick("mibs_clone_scope", choices$mibs_clone_scope, optional = !sulkseq,
         default = "sample")
    flag("mibs_productive_only", optional = TRUE)
    TRUE
  }, error = function(e) FALSE)

  if (!checked) {
    refuse("The form could not be read, so nothing was run. Reload the page and try again.")
  }

  list(ok = length(errors) == 0L, errors = errors, clean = clean)
}
