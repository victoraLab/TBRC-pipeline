# Keep server.R resilient to partial Shiny deploys where global.R might lag one
# version behind. The fallbacks below mirror the live helpers closely enough to
# keep the app usable until both files are in sync again.
if (!exists("APP_DIR")) {
  APP_DIR <- normalizePath(
    if (file.exists("server_config.json")) "." else "frontend",
    winslash = "/",
    mustWork = TRUE
  )
}

if (!exists("RUN_STATS_FILE")) {
  RUN_STATS_FILE <- file.path(APP_DIR, "run_stats.tsv")
}

if (!exists("format_duration_compact", mode = "function")) {
  format_duration_compact <- function(total_seconds) {
    total_seconds <- suppressWarnings(as.numeric(total_seconds))
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
}

if (!exists("read_run_stats", mode = "function")) {
  read_run_stats <- function() {
    if (!file.exists(RUN_STATS_FILE)) {
      return(data.frame())
    }

    tryCatch(
      read.table(RUN_STATS_FILE, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(...) data.frame()
    )
  }
}

if (!exists("append_run_stat", mode = "function")) {
  append_run_stat <- function(run_stat) {
    if (!is.data.frame(run_stat) || nrow(run_stat) != 1) {
      stop("append_run_stat expects a one-row data.frame", call. = FALSE)
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
}

if (!exists("compute_run_metrics", mode = "function")) {
  compute_run_metrics <- function(
    run_stats,
    input_bytes = NA_real_,
    manual_minutes_per_dataset = 6,
    manual_minutes_per_100mb = 29
  ) {
    if (!is.data.frame(run_stats) || nrow(run_stats) == 0) {
      return(list(
        dataset_count = 0,
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
    valid_saved_size <- suppressWarnings(as.numeric(completed_runs$input_bytes))
    saved_size_keep <- is.finite(valid_saved_size) & valid_saved_size > 0
    manual_seconds_per_byte <- (manual_minutes_per_100mb * 60) / (100 * 1024 * 1024)
    total_time_saved_sec <- sum(valid_saved_size[saved_size_keep] * manual_seconds_per_byte, na.rm = TRUE)

    if (!any(saved_size_keep)) {
      total_time_saved_sec <- dataset_count * manual_minutes_per_dataset * 60
    }

    if (dataset_count == 0) {
      return(list(
        dataset_count = 0,
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

    list(
      dataset_count = dataset_count,
      total_time_saved_sec = total_time_saved_sec,
      estimated_runtime_sec = seconds_per_mb * typical_size_mb,
      has_runtime_model = TRUE
    )
  }
}

server <- function(input, output, session) {
  storage_server_transfer_mode <- reactive({
    selected_server <- input$server_name

    if (selected_server %in% server_names) {
      configured_mode <- get_storage_value(servers_config, selected_server, "transfer_mode")

      if (!is.na(configured_mode) && nzchar(configured_mode)) {
        return(configured_mode)
      }
    }

    "push"
  })

  known_users <- reactive({
    user_lookup <- character()

    sample_files <- character()
    runs_dir <- RUNS_DIR
    if (dir.exists(runs_dir)) {
      # Accept both the current JSON handoff and legacy TSV metadata so older
      # runs still contribute user/email history.
      sample_files <- list.files(runs_dir, pattern = "^sample\\.(tsv|json)$", recursive = TRUE, full.names = TRUE)
    }

    for (sample_file in sample_files) {
      sample_data <- NULL

      if (grepl("\\.json$", sample_file, ignore.case = TRUE)) {
        sample_data <- tryCatch(
          jsonlite::fromJSON(sample_file, simplifyVector = TRUE),
          error = function(...) NULL
        )

        if (is.null(sample_data) || !all(c("user_id", "user_email") %in% names(sample_data))) {
          next
        }

        current_user <- trimws(as.character(sample_data$user_id[[1]]))
        current_email <- trimws(as.character(sample_data$user_email[[1]]))
      } else {
        sample_data <- tryCatch(
          read.table(sample_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
          error = function(...) NULL
        )

        if (is.null(sample_data) || !all(c("user_id", "user_email") %in% names(sample_data))) {
          next
        }

        if (nrow(sample_data) < 1) {
          next
        }

        current_user <- trimws(sample_data$user_id[[1]])
        current_email <- trimws(sample_data$user_email[[1]])
      }

      if (nzchar(current_user) && nzchar(current_email)) {
        user_lookup[current_user] <- current_email
      }
    }

    legacy_config <- app_path("..", "backend", "config", "config.yaml")
    if (file.exists(legacy_config)) {
      config_lines <- readLines(legacy_config, warn = FALSE)
      current_user <- NULL

      for (line in config_lines) {
        trimmed_line <- trimws(line)

        if (grepl("^user_id\\s*:\\s*", trimmed_line)) {
          current_user <- sub("^user_id\\s*:\\s*", "", trimmed_line)
        } else if (!is.null(current_user) && grepl("^user_email\\s*:\\s*", trimmed_line)) {
          current_email <- sub("^user_email\\s*:\\s*", "", trimmed_line)
          current_user <- trimws(current_user)
          current_email <- trimws(current_email)

          if (nzchar(current_user) && nzchar(current_email)) {
            user_lookup[current_user] <- current_email
          }

          current_user <- NULL
        }
      }
    }

    user_lookup
  })

  historical_run_count <- reactive({
    runs_dir <- RUNS_DIR
    if (!dir.exists(runs_dir)) {
      return(0L)
    }

    sample_files <- list.files(runs_dir, pattern = "^sample\\.(json|tsv)$", recursive = TRUE, full.names = TRUE)
    length(unique(dirname(sample_files)))
  })
  
  storage_server_rawdata_path <- reactive({
    selected_server <- input$server_name
    if (selected_server %in% server_names) {
      return(get_storage_value(servers_config, selected_server, "rawdata_path"))
    }

    NA_character_
  })

  storage_server_results_path <- reactive({
    selected_server <- input$server_name
    if (selected_server %in% server_names) {
      return(get_storage_value(servers_config, selected_server, "results_path"))
    }

    NA_character_
  })

  current_run_name <- reactive({
    if (input$type_input == "Upload" && !is.null(input$input_files$name)) {
      uploaded_names <- input$input_files$name
      if (length(uploaded_names) >= 1) {
        r1_name <- uploaded_names[[1]]
        r2_name <- if (length(uploaded_names) >= 2) uploaded_names[[2]] else uploaded_names[[1]]
        base1 <- sub("\\.fastq(\\.gz)?$", "", r1_name, ignore.case = TRUE)
        base2 <- sub("\\.fastq(\\.gz)?$", "", r2_name, ignore.case = TRUE)
        for (pattern in c("([_\\.-])R1([_\\.-]|$)", "([_\\.-])R2([_\\.-]|$)", "([_\\.-])001([_\\.-]|$)", "([_\\.-])S\\d+([_\\.-]|$)", "([_\\.-])L\\d{3}([_\\.-]|$)")) {
          base1 <- gsub(pattern, "\\1", base1, ignore.case = TRUE)
          base2 <- gsub(pattern, "\\1", base2, ignore.case = TRUE)
        }
        base1 <- gsub("[_\\.-]+$", "", base1)
        base2 <- gsub("[_\\.-]+$", "", base2)
        if (nzchar(base1) && identical(base1, base2)) {
          return(base1)
        }
        return(sub("\\.gz$", "", uploaded_names[[1]], ignore.case = TRUE))
      }
    }

    basename(trimws(input$input_folder))
  })

  classify_uploaded_read <- function(file_name) {
    if (grepl("(^|[_\\.-])R1([_\\.-]|$)", file_name, ignore.case = TRUE)) {
      return("R1")
    }
    if (grepl("(^|[_\\.-])R2([_\\.-]|$)", file_name, ignore.case = TRUE)) {
      return("R2")
    }
    NA_character_
  }

  # Clonality reads the IgBLAST table. Ticking it alone produced a run where the
  # step was silently skipped -- final.script.sh logs "Clonality skipped because
  # IGBlast output was not found" to a file nobody opens, and the delivered
  # archive simply has no clonality/ folder. Turn annotation on rather than let
  # the request evaporate.
  observeEvent(input$run_clonality, {
    if (isTRUE(input$run_clonality) && !isTRUE(input$run_igblast)) {
      updateCheckboxInput(session, "run_igblast", value = TRUE)
      showNotification(
        "Clonality needs the IgBLAST annotation, so that has been switched on too.",
        type = "message", duration = 8
      )
    }
  }, ignoreInit = TRUE)

  # And the other direction: turning annotation off strands a clonality request.
  observeEvent(input$run_igblast, {
    if (!isTRUE(input$run_igblast) && isTRUE(input$run_clonality)) {
      updateCheckboxInput(session, "run_clonality", value = FALSE)
      showNotification(
        "Clonality was switched off as well, because it has nothing to read without IgBLAST.",
        type = "message", duration = 8
      )
    }
  }, ignoreInit = TRUE)

  # "Contigs or isoforms to retrieve" means something different for MIBS. The
  # other pipelines take the top few sequences from a well holding one cell; a
  # MIBS well holds a pool, and this is how many of them to keep -- so it should
  # start at the number of cells sorted per well. 100 is the published protocol,
  # but people sort 10, 50 or 500, so it stays a slider rather than becoming a
  # constant.
  SULK_DEFAULT_CELLS <- 100L
  STANDARD_DEFAULT_CONTIGS <- 3L

  observeEvent(input$method, {
    if (identical(input$method, "tiago_sulkseq")) {
      if (!identical(as.integer(input$n_read), SULK_DEFAULT_CELLS)) {
        updateSliderInput(session, "n_read", value = SULK_DEFAULT_CELLS)
        showNotification(
          paste0("Sequences per well set to ", SULK_DEFAULT_CELLS,
                 ", the cells per well in the published MIBS protocol. ",
                 "Change it if this library sorted a different number."),
          type = "message", duration = 10
        )
      }
    } else if (identical(as.integer(input$n_read), SULK_DEFAULT_CELLS)) {
      # Only undo our own default, never a number someone chose.
      updateSliderInput(session, "n_read", value = STANDARD_DEFAULT_CONTIGS)
    }
  }, ignoreInit = TRUE)

  # Same shape as the slider above, and for the same reason: a default that is
  # wrong for the chosen pipeline is worse than no default. final.script.sh
  # exits 1 when the trim sequence matches nothing, and a nanopore well's 5'
  # handle is already consumed when the read is assigned to its well, so the
  # trim can only fail or discard the library. The field stays visible because
  # the hazard is real and hiding it would not remove it.
  observeEvent(input$method, {
    if (identical(input$method, "tiago_nanopore") &&
        isTRUE(as.logical(input$trim_primers))) {
      app_update_radio_buttons(session, inputId = "trim_primers",
                               choices = c(TRUE, FALSE), selected = "FALSE")
      showNotification(
        paste("The final primer trim was switched off: a nanopore well's 5' handle is",
              "already consumed when the read is assigned to its well, and trimming",
              "again discards the library."),
        type = "message", duration = 10
      )
    }
    # Never switched back on when the pipeline changes away. Undoing our own
    # default is one thing; overriding a choice somebody made is another.
  }, ignoreInit = TRUE)

  output$nReadHint <- renderText({
    if (identical(input$method, "tiago_sulkseq")) {
      paste0("For MIBS this is the number of cells sorted into each well: every ",
             "sequence above the read floor is kept, up to this many. The published ",
             "protocol uses 100.")
    } else {
      paste0("How many sequences to keep per well, ranked by read count. A ",
             "single-cell well needs only a few.")
    }
  })

  resolved_igblast_data <- reactive({
    get_cluster_value_or(
      servers_config,
      "igblast_data",
      file.path(get_cluster_value(servers_config, "pipeline_root"), "igblast", "internal_data")
    )
  })

  resolved_igblast_auxiliary_data <- reactive({
    # optional_file, not refs. refs is empty on the cluster, so the fallback
    # silently resolved to a directory with no aux file in it -- and IgBLAST
    # without an auxiliary file cannot read the J reading frame, which makes
    # every sequence non-productive. That has happened once already.
    get_cluster_value_or(
      servers_config,
      "igblast_auxiliary_data",
      file.path(get_cluster_value(servers_config, "pipeline_root"), "igblast", "optional_file")
    )
  })
  
  # One-row metadata table used for writing sample.tsv and validation.
  submission_row <- reactive({
    # Processing filename and path
    run_name <- current_run_name()
    source_path <- paste0(input$user_id, "/", run_name)
    
    
      
    # Creating a dataframe with the details of the submission
    data.frame(
      user_id = input$user_id,
      user_email = input$email,
      # Published so the backend does not have to re-derive it from FASTQ
      # filenames. When the two disagreed, results landed at
      # results/<user>/<folder>/<other-name>/ and were inconsistent run to run.
      sample_name = run_name,
      input_method = input$type_input,
      server_name = input$server_name,
      server_path = storage_server_rawdata_path(),
      server_results_path = storage_server_results_path(),
      server_ssh_port = get_storage_value(servers_config, input$server_name, "ssh_port"),
      server_transfer_mode = storage_server_transfer_mode(),
      primer_set = input$primer,
      trim_primers = input$trim_primers,
      sequence_number = input$n_read * 2,
      pipeline = input$method,
      luka_light_protocol = input$luka_light,
      run_igblast = input$run_igblast,
      run_clonality = input$run_clonality,
      archive_format = input$archive_format,
      igblast_species = input$igblast_species,
      igblast_panel = input$igblast_panel,
      # MIBS clone calling, per run. pipelines.yaml still holds the defaults and
      # still governs a run that omits these; when they are present they win.
      # Written for every pipeline, not only MIBS, because the required-field
      # subscript below indexes this frame by name and an absent column is an
      # "undefined columns selected" error rather than a missing value.
      mibs_junction_max_distance = as.integer(input$mibs_junction_max_distance %||% 4L),
      mibs_gene_grouping         = input$mibs_gene_grouping %||% "strict",
      mibs_clone_scope           = input$mibs_clone_scope %||% "sample",
      mibs_min_reads_per_well    = as.integer(input$mibs_min_reads_per_well %||% 6L),
      mibs_min_wells_per_clone   = as.integer(input$mibs_min_wells_per_clone %||% 2L),
      mibs_productive_only       = isTRUE(input$mibs_productive_only),
      mibs_min_reads_per_cluster = as.integer(input$mibs_min_reads_per_cluster %||% 2L),
      # Whether a design was attached, passed validation, and belongs to a
      # pipeline that reads it. The backend refuses to run a MIBS library that
      # says yes here and arrives without the file, because pooling every well
      # into one sample produces a clone count that looks entirely normal.
      # Only tiago_sulkseq is asked: stage_platemap() refuses to copy the file
      # for the other three, so recording true for them would put a claim in
      # sample.json and in the run log that nothing on disk backs.
      platemap_attached          = identical(input$method, "tiago_sulkseq") &&
                                     isTRUE(platemap_state()$ok),
      # Recorded with the submission for later statistics; neither is read by
      # the workflow. Demultiplexing uses the barcode files, not this count.
      sequencing_kit = input$sequencing_kit %||% "",
      receptor_check = input$receptor_check %||% "auto",
      igblast_bin = get_cluster_value_or(servers_config, "igblast_bin", "igblastn"),
      igblast_organism = input$igblast_species,
      igblast_db_v = get_cluster_value(servers_config, "igblast_db_v"),
      igblast_db_d = get_cluster_value(servers_config, "igblast_db_d"),
      igblast_db_j = get_cluster_value(servers_config, "igblast_db_j"),
      igblast_data = resolved_igblast_data(),
      igblast_auxiliary_data = resolved_igblast_auxiliary_data(),
      trim_seq = input$trim_sequence,
      keep_raw = input$keep_intermediate,
      # Separate from keep_raw: keep_raw decides what survives on disk in the
      # run directory, archive_scope decides what is put into the delivered
      # archive.
      archive_scope = input$archive_scope,
      run_group = trimws(input$run_group %||% ""),
      source = source_path,
      dest = paste0("results/", source_path),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  # Reactive function to create a transposed submission form for display in the UI.
  # What the submission preview shows.
  #
  # This used to be as.data.frame(t(submission_row())), which transposes a
  # one-row frame into a single unnamed column -- so the header read "V1" and the
  # field names were rownames rather than data. It also arrived in whatever order
  # the row happened to be built in, which is an implementation detail rather
  # than anything a reader would choose.
  #
  # It is a settings list, not a dataset: two named columns, in a deliberate
  # order, with labels people recognise rather than the internal field names.
  # Anything not in the map still appears, at the end, so a new field can never
  # silently vanish from the preview.
  PREVIEW_FIELDS <- c(
    sample_name          = "Run name",
    user_id              = "User",
    user_email           = "Email",
    pipeline             = "Assembly pipeline",
    primer_set           = "Barcode set",
    input_method         = "Input method",
    server_name          = "Storage server",
    server_path          = "Source folder",
    source               = "Source path",
    dest                 = "Results path",
    server_transfer_mode = "Transfer mode",
    trim_primers         = "Trim primers",
    trim_seq             = "Trim sequence",
    luka_light_protocol  = "Barcode preprocessing",
    sequence_number      = "Sequences per well",
    mibs_min_reads_per_well    = "MIBS reads per well",
    mibs_min_wells_per_clone   = "MIBS wells per clone",
    mibs_junction_max_distance = "MIBS junction distance",
    mibs_gene_grouping         = "MIBS gene grouping",
    mibs_clone_scope           = "MIBS clone scope",
    mibs_min_reads_per_cluster = "MIBS reads per cluster",
    mibs_productive_only       = "MIBS productive only",
    platemap_attached          = "Experimental groups",
    run_igblast          = "Run IgBLAST",
    igblast_panel        = "IgBLAST panel",
    igblast_species      = "IgBLAST species",
    receptor_check       = "Receptor check",
    run_clonality        = "Run clonality",
    archive_format       = "Archive format",
    archive_scope        = "Archive scope",
    run_group            = "Results subfolder",
    keep_raw             = "Keep intermediates"
  )

  # What is still missing before a run can be submitted.
  #
  # The submit button used to be always clickable, so the only way to discover a
  # missing field was to press Run and read an alert -- one field at a time,
  # since the checks return on the first failure. This computes the whole list at
  # once and shows it, so the button becomes a statement about readiness rather
  # than a thing you try.
  #
  # Each entry is (done, label, what to do about it). The wording matters: "Pick
  # a barcode set" is actionable, "primer_set is required" is not.
  # What a run cannot go without, by pipeline.
  #
  # This was one flat list applied to every pipeline, which meant a nanopore run
  # was blocked on a 5' trim handle no nanopore rule reads, and a MIBS run was
  # not blocked on the clone thresholds that decide its entire answer. The
  # names are sample.json COLUMN names, not input ids: trim_seq not
  # trim_sequence, user_email not email.
  #
  # The check that consumes this tests for NA and for "". A radio or a picker
  # can never be either, so listing one is a statement of intent rather than a
  # working guard; the fields that can actually arrive empty are the text boxes,
  # the derived paths, and the MIBS numeric boxes.
  REQUIRED_FIELDS_BASE <- c(
    "user_id",
    "user_email",
    "input_method",
    "primer_set",
    # Read unconditionally at Snakefile:151 and turned into TOPN with int() at
    # :237, for every pipeline including nanopore, which never uses the result.
    # Blank here is a crash there, so it stays in the base set.
    "sequence_number",
    "pipeline",
    "run_igblast",
    "run_clonality",
    "archive_format",
    "igblast_species",
    # All four read it: classic_luka and vsearch hand it to final.script.sh,
    # nanopore and sulkseq each choose their germline panel from it.
    "igblast_panel",
    "keep_raw",
    "archive_scope",
    "source",
    "dest"
  )

  REQUIRED_FIELDS_BY_PIPELINE <- list(
    # trim_seq and trim_primers are the 5' handle trim final.script.sh applies
    # at the end of the run. luka_light_protocol is read by genfasta.smk for
    # classic_luka and by illumina_core.smk for the two vsearch-based backends.
    classic_luka  = c("trim_primers", "trim_seq", "luka_light_protocol"),
    tiago_vsearch = c("trim_primers", "trim_seq", "luka_light_protocol"),
    # Same front half as vsearch, plus the numbers that define a MIBS clone. An
    # empty box here must stop the run: the alternative is a null in sample.json
    # that the Snakefile discards, so the run silently uses the config default
    # rather than the number the submitter was in the middle of typing.
    # mibs_productive_only is deliberately absent -- FALSE is a real answer and
    # this check cannot tell it from unset.
    tiago_sulkseq = c("trim_primers", "trim_seq", "luka_light_protocol",
                      "mibs_junction_max_distance", "mibs_gene_grouping",
                      "mibs_clone_scope", "mibs_min_reads_per_well",
                      "mibs_min_wells_per_clone", "mibs_min_reads_per_cluster"),
    # No nanopore RULE reads METH, TRIM or TRIM_PRIMERS; a nanopore well's
    # handle is consumed by the FR1 anchor during demultiplexing. Requiring a
    # trim sequence for a run that does not trim is how the flat list blocked
    # submissions on a field with nothing to fill in.
    #
    # Nothing is required here unconditionally, but nothing is unconditionally
    # safe either: rules/final.smk is included for every pipeline and hands TRIM
    # and TRIM_PRIMERS to final.script.sh, which runs cutadapt with an empty
    # adapter if the trim is on and the sequence is blank. trim_seq is therefore
    # added below whenever the trim is switched on, for every pipeline.
    tiago_nanopore = character(0)
  )

  REQUIRED_FIELDS_SERVER <- c(
    "server_name",
    "server_path",
    "server_results_path",
    "server_ssh_port",
    "server_transfer_mode"
  )

  required_fields_for <- function(pipeline, input_method, known, trim_primers = NULL) {
    fields <- c(REQUIRED_FIELDS_BASE,
                REQUIRED_FIELDS_BY_PIPELINE[[pipeline %||% ""]] %||% character(0))
    # The one field whose need is a choice rather than a pipeline.
    # final.script.sh runs cutadapt -g "$TRIMMING_SEQ" --discard-untrimmed
    # whenever the trim is on, for every pipeline including nanopore, and an
    # empty adapter discards the library after the whole assembly has run. So
    # switching the trim on is what makes the sequence required, and switching
    # it off is what makes it optional -- not which backend was chosen.
    if (isTRUE(as.logical(trim_primers %||% NA))) fields <- c(fields, "trim_seq")
    if (identical(input_method, "Server")) fields <- c(fields, REQUIRED_FIELDS_SERVER)
    # The consumer subscripts submission_row() by these names, and an unknown
    # name there is "undefined columns selected" -- a stack trace instead of the
    # alert this check exists to raise.
    intersect(unique(fields), known)
  }

  # numericInput's min and max are advisory: R accepts whatever is typed and
  # Shiny sends it. The Snakefile then range-checks each MIBS number and drops
  # one it cannot use, so a typed 0 leaves the config default in place and the
  # run answers a question nobody asked -- the exact outcome the empty-box check
  # was written to prevent. Same bounds as the widgets and as the Snakefile, in
  # one place, and a non-integer counts as out of range because int(float()) on
  # the other side would silently truncate it.
  MIBS_RANGES <- list(
    mibs_min_reads_per_well    = c(1, 10000),
    mibs_min_wells_per_clone   = c(1, 96),
    mibs_junction_max_distance = c(0, 30),
    mibs_min_reads_per_cluster = c(1, 1000)
  )

  mibs_bad_numbers <- function() {
    bad <- character(0)
    for (id in names(MIBS_RANGES)) {
      value <- suppressWarnings(as.numeric(input[[id]]))
      limits <- MIBS_RANGES[[id]]
      if (length(value) != 1L || is.na(value) || value != floor(value) ||
          value < limits[[1]] || value > limits[[2]]) {
        bad <- c(bad, id)
      }
    }
    bad
  }

  submission_checklist <- reactive({
    is_upload <- identical(input$type_input, "Upload")
    items <- list()

    add_item <- function(done, label, todo) {
      items[[length(items) + 1]] <<- list(done = isTRUE(done), label = label, todo = todo)
    }

    add_item(!is.null(input$user_id) && is_valid_path_segment(input$user_id),
             "Your name",
             "Choose or type your user name")

    # The run name is DERIVED, not typed: in Server mode it comes from the
    # source folder, in Upload mode from the uploaded file names. An earlier
    # version of this checklist watched input$folder_name, which belongs to the
    # barcode uploader -- so typing a perfectly good source folder left the list
    # demanding a run name that no visible field on the form corresponds to.
    run_nm <- tryCatch(current_run_name(), error = function(e) "")
    add_item(!is.null(run_nm) && nzchar(trimws(run_nm %||% "")),
             "Run name",
             if (is_upload) "Upload the FASTQ files; the run is named from them"
             else "Type the source folder; the run takes its name from it")

    if (is_upload) {
      n_files <- length(input$input_files$name %||% character(0))
      add_item(n_files == 2, "FASTQ files",
               if (n_files == 0) "Upload the R1 and R2 FASTQ files"
               else sprintf("Upload exactly two files (you have %d)", n_files))
    } else {
      add_item(!is.null(input$input_folder) && nzchar(trimws(input$input_folder %||% "")),
               "Source folder",
               "Type the folder name on the storage server")
      add_item(!is.null(input$server_name) && nzchar(input$server_name %||% ""),
               "Storage server",
               "Choose which server holds the data")
    }

    add_item(!is.null(input$primer) && nzchar(input$primer %||% "") &&
               !identical(input$primer, "placeholder"),
             "Barcode set",
             "Choose the barcode set used for this library")

    add_item(is_valid_email(input$email),
             "Email (optional)",
             "The email address does not look valid")

    # Only MIBS reads either of these, and the run button should say so before
    # anyone presses it. An empty numeric box is the trap: it looks like a
    # filled-in form and reaches sample.json as null.
    if (identical(input$method, "tiago_sulkseq")) {
      bad_mibs <- mibs_bad_numbers()
      add_item(!length(bad_mibs),
               "MIBS clone options",
               sprintf(paste("%s: fill in a whole number inside the range shown.",
                             "An empty or out-of-range box is not the same as the default"),
                       paste(unname(PREVIEW_FIELDS[bad_mibs]), collapse = ", ")))

      # An absent file leaves this done: the design is optional to the software,
      # and only an attached and invalid one blocks the run.
      pm <- platemap_state()
      add_item(is.null(input$platemap_file) || isTRUE(pm$ok),
               "Experimental groups",
               if (length(pm$errors)) pm$errors[[1]]
               else "The experimental groups file could not be read")
    }

    items
  })

  submission_ready <- reactive({
    all(vapply(submission_checklist(), function(x) x$done, logical(1)))
  })

  form <- reactive({
    row <- submission_row()
    if (is.null(row) || !ncol(row)) {
      return(data.frame(Setting = character(0), Value = character(0),
                        stringsAsFactors = FALSE))
    }
    vals <- vapply(names(row), function(k) {
      v <- row[[k]]
      v <- if (length(v) > 1) paste(v, collapse = ", ") else as.character(v)
      if (is.na(v) || !nzchar(trimws(v))) "\u2014" else trimws(v)
    }, character(1))

    # The MIBS knobs and the design flag are written for every pipeline so the
    # required-field subscript can never hit an undefined column. Showing "MIBS
    # clone scope" on a classic_luka submission is noise about a setting nothing
    # will read, and only tiago_sulkseq reads the design.
    if (!identical(row$pipeline, "tiago_sulkseq")) {
      vals <- vals[!grepl("^mibs_|^platemap_attached$", names(vals))]
    }

    known <- intersect(names(PREVIEW_FIELDS), names(vals))
    extra <- setdiff(names(vals), known)
    ordered <- c(known, sort(extra))
    labels <- ifelse(ordered %in% names(PREVIEW_FIELDS),
                     PREVIEW_FIELDS[ordered],
                     gsub("_", " ", ordered))

    data.frame(Setting = unname(labels),
               Value = unname(vals[ordered]),
               stringsAsFactors = FALSE)
  })
  
  # Reactive function to read barcode names file
  bc_names <- reactive({
    df <- read.table(file.path(BC_DIR, "bc_names.txt"), stringsAsFactors = FALSE)
    bcnames <- df[["bcname"]]
    names(bcnames) <- rownames(df)
    return(bcnames)
  })

  
  
  # Reactive to list folders in the "bc" directory
  folder_list <- reactive({
    list.dirs(BC_DIR, recursive = FALSE, full.names = FALSE)
  })
  
  # Global variable to hold logs
  log <- reactiveVal("")
  addlog <- function(message) {
    log(paste0(isolate(log()), "\n", message))
  }
  progress_detail <- reactiveVal("Waiting for submission.")
  set_progress_detail <- function(...) {
    progress_detail(paste(..., collapse = "\n"))
  }
  notification_log_path <- file.path(STATE_DIR, "notification_log.tsv")

  append_notification_log <- function(recipient, subject, status, context, detail = "") {
    log_entry <- data.frame(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      recipient = as.character(recipient),
      subject = as.character(subject),
      status = as.character(status),
      context = as.character(context),
      detail = as.character(detail),
      stringsAsFactors = FALSE
    )

    write.table(
      log_entry,
      file = notification_log_path,
      sep = "\t",
      quote = TRUE,
      row.names = FALSE,
      col.names = !file.exists(notification_log_path),
      append = file.exists(notification_log_path)
    )
  }
  # ---- experimental design ------------------------------------------------
  # A template the submitter can fill in and hand back. Plate rows for the
  # common case where one plate is one animal, then a commented well-level
  # example for a plate that is split.
  output$platemapTemplate <- downloadHandler(
    filename = function() "tbrc_experimental_groups_template.csv",
    contentType = "text/csv",
    content = function(file) {
      lines <- c(
        "# TBRC experimental groups. Two columns: the plate or well, then the label.",
        "# Save as CSV or tab-separated text and upload it with the run.",
        "#",
        "# Short form -- one row per plate. Repeat a label for plates that share it.",
        "# Long form  -- one row per well (P01A01), for a plate split between groups.",
        "# A third column, if you add one, is the broader experimental group.",
        "#",
        "# Uncomment the rows you need and fill in the label. Lines starting",
        "# with # are ignored, by this form and by the pipeline alike.",
        # Not plate_or_well: load_platemap() recognises that row as a header
        # and then finds no key column in it, so every row below was skipped
        # and the map came back empty with no error. A template its own parser
        # could not read.
        "well,sample,group",
        # Commented, not blank. The template used to ship twelve rows reading
        # "P01,," and a plate with no label is a row the pipeline pools into a
        # sample named after the run, so the form rejects it -- which made the
        # downloaded template fail its own validator unless every plate in it
        # was either filled in or deleted.
        paste0("# P", sprintf("%02d", 1:4), ",GF_M", 1:4, ",GF"),
        "# P01A01,GF_M1,GF",
        "# P01A02,SPF_M4,SPF"
      )
      writeLines(lines, file)
    }
  )

  # canonical_well_key() is a port of canonical_well() in
  # backend/workflow/scripts/sulkseq_clones.py, so the form judges a key exactly
  # the way the pipeline will read it, plus zero-padding of the plate number,
  # which the Python side gains in the same change. A hand-typed P1 and the
  # demultiplexer's P01 are the same plate; unpadded they were two, and the
  # wells of the second one became a sample nobody had sorted.
  #
  # P is not a row letter, so the three plate-bearing spellings cannot be
  # confused with the two plate-less ones and one ordered pass is unambiguous.
  #
  # Matched on the upper-cased token because a plate map is typed by hand and
  # p1 is what gets typed. The Python side upper-cases too, so a lower-case key
  # the pipeline reads correctly was being refused here: the form was stricter
  # than the thing it is meant to predict. The original token is what comes back
  # when nothing matches, so an error message quotes what the submitter wrote.
  canonical_well_key <- function(token) {
    token <- trimws(as.character(token %||% ""))
    if (!length(token) || is.na(token[1]) || !nzchar(token[1])) return("")
    original <- token[1]
    token <- toupper(original)
    m <- regmatches(token, regexec("^P([0-9]{1,2})([A-Ha-h])([0-9]{1,2})$", token))[[1]]
    if (length(m) == 4L) {
      return(sprintf("P%02d%s%02d", as.integer(m[2]), toupper(m[3]), as.integer(m[4])))
    }
    m <- regmatches(token, regexec("^([0-9]{1,2})P([0-9]{1,2})([A-Ha-h])$", token))[[1]]
    if (length(m) == 4L) {
      return(sprintf("P%02d%s%02d", as.integer(m[3]), toupper(m[4]), as.integer(m[2])))
    }
    m <- regmatches(token, regexec("^([A-Ha-h])([0-9]{1,2})P([0-9]{1,2})$", token))[[1]]
    if (length(m) == 4L) {
      return(sprintf("P%02d%s%02d", as.integer(m[4]), toupper(m[2]), as.integer(m[3])))
    }
    m <- regmatches(token, regexec("^([A-Ha-h])([0-9]{1,2})$", token))[[1]]
    if (length(m) == 3L) return(sprintf("%s%02d", toupper(m[2]), as.integer(m[3])))
    m <- regmatches(token, regexec("^([0-9]{1,2})([A-Ha-h])$", token))[[1]]
    if (length(m) == 3L) return(sprintf("%s%02d", toupper(m[3]), as.integer(m[2])))
    m <- regmatches(token, regexec("^P([0-9]{1,2})$", token))[[1]]
    if (length(m) == 2L) return(sprintf("P%02d", as.integer(m[2])))
    original
  }

  # Every well the chosen barcode set can produce, built the way demux_wells.sh
  # builds them: each BC1 name concatenated with each BC2 name, then put through
  # the same canonicaliser the pipeline uses. The two files split the well name
  # differently in every set -- B names columns in BC1 and plate+row in BC2, T
  # is the other way round, B2, CBC and B_Juhee split it as row+column and plate
  # -- and the concatenation is the one rule all six obey, because each output
  # file is named after the two barcodes that matched it.
  #
  # The cross product is the point rather than an implementation detail.
  # Classifying the tokens and multiplying the three vocabularies said LukaP has
  # sixteen wells; LukaP has four (A03, B04, E01 and G10, on one plate), so a
  # design naming P01A10 passed a check that only asked whether row A and column
  # 10 exist somewhere in the set.
  barcode_geometry <- function(barcode_set) {
    empty <- list(plates = character(0), rows = character(0), cols = character(0),
                  well_names = character(0), wells = 0L, ok = FALSE)
    set_name <- trimws(as.character(barcode_set %||% ""))
    if (!length(set_name) || is.na(set_name[1]) || !nzchar(set_name[1])) return(empty)
    barcode_names <- function(bc_file) {
      bc_path <- file.path(BC_DIR, set_name[1], bc_file)
      if (!file.exists(bc_path)) return(character(0))
      # warn = FALSE because the shipped files are not all newline-terminated;
      # readLines still returns the last barcode, it only complains about it.
      lines <- tryCatch(readLines(bc_path, warn = FALSE), error = function(...) character(0))
      lines <- toupper(sub("[ \t].*$", "", trimws(lines)))
      unique(lines[nzchar(lines)])
    }
    bc1 <- barcode_names("BC1.txt")
    bc2 <- barcode_names("BC2.txt")
    if (!length(bc1) || !length(bc2)) return(empty)
    wells <- unique(vapply(as.vector(outer(bc1, bc2, paste0)),
                           canonical_well_key, character(1), USE.NAMES = FALSE))
    plated <- grepl("^P[0-9]{2}[A-H][0-9]{2}$", wells)
    plateless <- grepl("^[A-H][0-9]{2}$", wells)
    wells <- sort(wells[plated | plateless])
    if (!length(wells)) return(empty)
    has_plate <- grepl("^P", wells)
    list(plates = sort(unique(substr(wells[has_plate], 1L, 3L))),
         rows = sort(unique(ifelse(has_plate, substr(wells, 4L, 4L),
                                   substr(wells, 1L, 1L)))),
         cols = sort(unique(ifelse(has_plate, substr(wells, 5L, 6L),
                                   substr(wells, 2L, 3L)))),
         well_names = wells, wells = length(wells), ok = TRUE)
  }

  # Read the uploaded design the way the pipeline will, and say what it will
  # produce before the run is dispatched rather than after it has finished. A
  # design that names a plate the library does not have costs a whole run to
  # discover, and the run it costs looks successful.
  validate_platemap <- function(path, barcode_set) {
    errors <- character(0)
    warnings <- character(0)
    shape <- "unknown"
    plate_keys <- character(0)
    well_keys <- character(0)
    samples <- character(0)
    groups <- character(0)
    duplicates <- character(0)
    unknown_keys <- character(0)
    n_rows <- 0L
    geometry <- barcode_geometry(barcode_set)
    set_label <- trimws(as.character(barcode_set %||% ""))
    if (!length(set_label) || is.na(set_label[1]) || !nzchar(set_label[1])) {
      set_label <- "(none chosen)"
    } else {
      set_label <- set_label[1]
    }

    # wells_covered counts the wells the design NAMES. A plate row stands for
    # every well the run recovers on that plate, and that number does not exist
    # until the library has been demultiplexed -- multiplying the plate count by
    # the barcode set's rows and columns claimed 192 wells for a two-plate
    # design whose samples are 16 to 32 wells each, and the preview, the
    # geometry line and the staging notification all repeated it.
    state <- function() {
      list(ok = length(errors) == 0L, errors = errors, warnings = warnings,
           rows = as.integer(n_rows), shape = shape, plate_keys = plate_keys,
           well_keys = well_keys, wells_covered = as.integer(length(well_keys)),
           design_plates = sort(unique(c(
             plate_keys, substr(well_keys[grepl("^P", well_keys)], 1L, 3L)))),
           samples = samples, groups = groups, duplicates = duplicates,
           unknown_keys = unknown_keys)
    }

    if (is.null(path) || !length(path) || is.na(path[1]) || !nzchar(path[1]) ||
        !file.exists(path[1])) {
      errors <- c(errors, "The experimental groups file could not be read.")
      return(state())
    }
    size <- file.info(path[1])$size
    if (!is.na(size) && size > 2 * 1024 * 1024) {
      errors <- c(errors, "The experimental groups file is unexpectedly large (over 2 MB).")
      return(state())
    }

    lines <- tryCatch(readLines(path[1], warn = FALSE, n = 5000),
                      error = function(...) character(0))
    lines <- lines[nzchar(trimws(lines)) & !grepl("^\\s*#", lines)]
    if (!length(lines)) {
      errors <- c(errors, "The experimental groups file has no usable rows.")
      return(state())
    }

    # Sniffed from the first surviving row, never the physical first line: the
    # template ships with eight comment lines above the data, and a single
    # "[\t,;]|\\s+" pattern cannot tell a comma inside a label from a delimiter.
    first <- lines[[1]]
    delim <- if (grepl("\t", first, fixed = TRUE)) "\t"
             else if (grepl(",", first, fixed = TRUE)) ","
             else NA_character_
    # scan() rather than strsplit() for the two real delimiters, because
    # load_platemap() reads them with csv.reader and csv.reader honours quoting:
    # a label written "GF, M1" is one field there and was two here, so the
    # preview named a sample that would never appear in the results. quote is
    # the double quote alone -- scan's default also treats an apostrophe as one,
    # and 5' and 3' turn up in labels. The whitespace branch matches Python's
    # line.split(), which does no quoting either.
    split_row <- function(line) {
      if (is.na(delim)) return(trimws(strsplit(trimws(line), "[ \t]+")[[1]]))
      trimws(tryCatch(
        scan(text = line, what = "", sep = delim, quote = "\"", quiet = TRUE,
             strip.white = TRUE, blank.lines.skip = FALSE),
        error = function(...) strsplit(line, delim, fixed = TRUE)[[1]]))
    }
    fields <- lapply(lines, split_row)
    if (median(lengths(fields)) < 2) {
      errors <- c(errors, paste("The experimental groups file needs at least two columns:",
                                "the plate or well, then the label."))
      return(state())
    }

    # The same token set load_platemap() uses, plus plate_or_well, which the
    # downloaded template used to write and the parser used to reject.
    header_tokens <- c("well", "plate", "sample", "bio_id", "bio.id", "group",
                       "plate_well", "plate_id", "plate_or_well")
    # In the order load_platemap()'s get() tries the names, not in the order the
    # columns happen to sit. A file headed "plate,well,sample" is a well map to
    # the pipeline, which asks for "well" first; the form took the leftmost
    # matching column and previewed it as a whole-plate design.
    key_tokens <- c("well", "plate_well", "plate_id", "plate_or_well", "plate")
    head_fields <- tolower(trimws(fields[[1]]))
    has_header <- any(head_fields %in% header_tokens)
    key_col <- 1L
    sample_col <- 2L
    group_col <- 3L
    data_rows <- fields
    if (has_header) {
      first_named <- function(wanted) {
        for (name in wanted) {
          hit <- which(head_fields == name)
          if (length(hit)) return(hit[[1]])
        }
        NA_integer_
      }
      key_col <- first_named(key_tokens)
      if (is.na(key_col)) {
        errors <- c(errors, paste("The first row looks like a header but names no key column.",
                                  "Rename the first column to 'well'."))
        return(state())
      }
      sample_col <- first_named(c("sample", "bio_id", "bio.id"))
      group_col <- first_named(c("group", "strain"))
      data_rows <- fields[-1]
      shape <- "full"
    }
    n_rows <- length(data_rows)

    field_at <- function(row, i) {
      if (is.na(i) || length(row) < i) return("")
      trimws(as.character(row[[i]]))
    }
    keys_raw <- vapply(data_rows, field_at, character(1), key_col)
    sample_raw <- vapply(data_rows, field_at, character(1), sample_col)
    group_raw <- vapply(data_rows, field_at, character(1), group_col)
    # load_platemap() falls back to the group per ROW, and only in the header
    # branch: a design exported with the sample column left empty and the label
    # written under group is read correctly by the pipeline, and the form used
    # to block it on "rows name a plate or well with no label". The headerless
    # branch has no such fallback there, so it has none here either -- a blank
    # second column really is a row the pipeline would pool into a sample named
    # after the run.
    labels <- if (has_header) ifelse(nzchar(sample_raw), sample_raw, group_raw)
              else sample_raw

    # A row with a label and no key is dropped silently by load_platemap(), and
    # silently is the problem: it is a row the submitter meant to matter.
    keyless <- sum(!nzchar(keys_raw))
    if (keyless > 0L) {
      errors <- c(errors, sprintf("%d row(s) carry a label with no plate or well.", keyless))
    }
    keep <- nzchar(keys_raw)
    keys_raw <- keys_raw[keep]
    labels <- labels[keep]
    group_raw <- group_raw[keep]

    keys <- vapply(keys_raw, canonical_well_key, character(1), USE.NAMES = FALSE)
    is_plate <- grepl("^P[0-9]{2}$", keys)
    is_well <- grepl("^P[0-9]{2}[A-H][0-9]{2}$", keys) | grepl("^[A-H][0-9]{2}$", keys)
    unknown_keys <- unique(keys[!(is_plate | is_well)])
    if (length(unknown_keys)) {
      errors <- c(errors, sprintf(
        "These keys are not a plate or a well: %s. A plate is P01, a well is P01A01, A01P01 or 01P01A.",
        paste(head(unknown_keys, 5L), collapse = ", ")))
    }
    plate_keys <- sort(unique(keys[is_plate]))
    well_keys <- sort(unique(keys[is_well]))
    if (!has_header) {
      shape <- if (sum(is_well) > sum(is_plate)) "long"
               else if (sum(is_plate) > sum(is_well)) "short"
               else "unknown"
    }

    if (isTRUE(geometry$ok)) {
      # A well key is checked against the set's actual well list rather than
      # against its rows and columns separately, because a barcode set is not
      # always a full plate: LukaP defines four wells and every one of P01A01 to
      # P01H12 passed a row-and-column check.
      plateless <- well_keys[!grepl("^P", well_keys)]
      # A key with no plate is not a well of the only plate. The well the
      # demultiplexer names is P01A03, platemap_for_well() looks up P01A03 and
      # then the plate key P01, and A03 is neither -- so this passed validation
      # as a warning, every well then landed in the unmapped bucket, and the
      # summary still reported platemap_used yes.
      if (length(plateless) && any(grepl("^P", geometry$well_names))) {
        errors <- c(errors, sprintf(
          "These keys name no plate: %s. Barcode set %s names its wells like %s, so write the plate in.",
          paste(head(plateless, 5L), collapse = ", "), set_label,
          geometry$well_names[[1]]))
      }
      bad_wells <- setdiff(setdiff(well_keys, geometry$well_names), plateless)
      if (length(bad_wells)) {
        errors <- c(errors, sprintf(
          "Barcode set %s cannot produce %s. It defines %d well(s), from %s to %s.",
          set_label, paste(head(bad_wells, 5L), collapse = ", "), geometry$wells,
          geometry$well_names[[1]], geometry$well_names[[length(geometry$well_names)]]))
      }
      bad_plates <- setdiff(plate_keys, geometry$plates)
      if (length(bad_plates)) {
        errors <- c(errors, sprintf(
          "Barcode set %s has no %s. It spans %s.",
          set_label, paste(head(bad_plates, 5L), collapse = ", "),
          paste(geometry$plates, collapse = ", ")))
      }
    } else {
      # No geometry, no constraint. Inventing one here would block a design that
      # is correct for a barcode set this app simply could not read.
      warnings <- c(warnings, sprintf(
        "Barcode set %s could not be read, so the plate and well names were not checked against it.",
        set_label))
    }

    duplicates <- unique(keys[duplicated(keys)])
    same_label <- character(0)
    conflicts <- 0L
    for (key in duplicates) {
      seen <- unique(labels[keys == key])
      seen <- seen[nzchar(seen)]
      if (length(seen) > 1L) {
        conflicts <- conflicts + 1L
        if (conflicts <= 5L) {
          errors <- c(errors, sprintf("Key %s is listed twice with different labels (%s and %s).",
                                      key, seen[[1]], seen[[2]]))
        }
      } else {
        same_label <- c(same_label, key)
      }
    }
    if (conflicts > 5L) {
      errors <- c(errors, sprintf("%d further key(s) are listed with conflicting labels.",
                                  conflicts - 5L))
    }
    if (length(same_label)) {
      warnings <- c(warnings, sprintf("%d key(s) are listed more than once with the same label.",
                                      length(same_label)))
    }

    blank_labels <- sum(!nzchar(labels))
    if (blank_labels > 0L) {
      errors <- c(errors, sprintf("%d row(s) name a plate or well with no label.", blank_labels))
    }

    samples <- sort(unique(labels[nzchar(labels)]))
    effective_group <- ifelse(nzchar(group_raw), group_raw, labels)
    groups <- sort(unique(effective_group[nzchar(effective_group)]))

    # Counted on the keys rather than on wells_covered, which counts named
    # wells only and is therefore zero for a perfectly good short-form design
    # that lists nothing but plates.
    if (!length(plate_keys) && !length(well_keys)) {
      errors <- c(errors, "The experimental groups file names no wells.")
    }

    state()
  }

  # ok = FALSE with no errors is the no-file case, and that is the whole point:
  # platemap_attached in sample.json is isTRUE(platemap_state()$ok), and no file
  # means no design was attached rather than a design that failed.
  platemap_state <- reactive({
    upload <- input$platemap_file
    if (is.null(upload) || !length(upload$datapath) || !nzchar(upload$datapath[[1]])) {
      return(list(ok = FALSE, errors = character(0), warnings = character(0),
                  rows = 0L, shape = "unknown", plate_keys = character(0),
                  well_keys = character(0), wells_covered = 0L,
                  design_plates = character(0),
                  samples = character(0), groups = character(0),
                  duplicates = character(0), unknown_keys = character(0)))
    }
    validate_platemap(upload$datapath[[1]], input$primer)
  })

  plural <- function(n) if (identical(as.integer(n), 1L)) "" else "s"

  # One sentence for what a design covers, shared by the preview and by the
  # notification the staging step raises, so the two cannot drift apart. Plates
  # and named wells are counted separately because they are different claims: a
  # well row covers one well, a plate row covers whatever the run recovers on
  # that plate, which nothing knows until the library is demultiplexed.
  platemap_coverage <- function(state) {
    parts <- character(0)
    n_wells <- state$wells_covered
    n_plates <- length(state$plate_keys)
    if (n_wells > 0L) {
      parts <- c(parts, sprintf("%d named well%s", n_wells, plural(n_wells)))
    }
    if (n_plates > 0L) {
      parts <- c(parts, sprintf("every well recovered on %d plate%s",
                                n_plates, plural(n_plates)))
    }
    if (!length(parts)) return("no wells")
    touched <- state$design_plates
    paste0(paste(parts, collapse = " and "),
           if (length(touched))
             sprintf(" (plate%s %s)", plural(length(touched)),
                     paste(c(head(touched, 8L),
                             if (length(touched) > 8L)
                               sprintf("and %d more", length(touched) - 8L)),
                           collapse = ", "))
           else "")
  }

  # What the design will do, said before the run rather than after it. The
  # no-file wording is deliberately the longest of the four: a MIBS run with no
  # design still produces a full set of numbers, and they are the numbers for
  # one pooled sample.
  output$platemapPreview <- renderUI({
    upload <- input$platemap_file
    if (is.null(upload) || !length(upload$datapath)) {
      return(tags$div(
        class = "card-note", style = "margin-top:6px;",
        paste("No experimental groups file attached. Every well will be treated as one",
              "sample. Incidence is only meaningful within a sample, so a library holding",
              "more than one mouse will have its clones merged, and the per-sample and",
              "per-group figures will be missing from the report.")))
    }

    state <- platemap_state()
    if (!isTRUE(state$ok)) {
      # The counts are suppressed on purpose. Reporting "3 samples" for a file
      # that will not be used is the failure this preview exists to prevent.
      return(tags$div(
        class = "card-note", style = "margin-top:6px;color:#A31621;",
        tags$div("This file was not attached. Fix it and upload it again:"),
        tags$ul(lapply(state$errors, tags$li))))
    }

    geometry <- barcode_geometry(input$primer)
    summary_line <- sprintf(
      "%s form, %d row%s: %d sample%s in %d group%s, covering %s.",
      state$shape, state$rows, plural(state$rows),
      length(state$samples), plural(length(state$samples)),
      length(state$groups), plural(length(state$groups)),
      platemap_coverage(state))
    geometry_line <- if (isTRUE(geometry$ok)) {
      sprintf("Barcode set %s can produce %d well%s across %d plate%s (%d rows, %d columns).",
              input$primer, geometry$wells, plural(geometry$wells),
              length(geometry$plates), plural(length(geometry$plates)),
              length(geometry$rows), length(geometry$cols))
    } else NULL
    samples_line <- sprintf(
      "Samples: %s.",
      paste(c(head(sort(state$samples), 6L),
              if (length(state$samples) > 6L)
                sprintf("and %d more", length(state$samples) - 6L)),
            collapse = ", "))

    tagList(
      tags$div(class = "card-note", style = "margin-top:6px;",
               tags$div(summary_line),
               if (!is.null(geometry_line)) tags$div(geometry_line),
               tags$div(samples_line)),
      if (length(state$warnings))
        tags$div(class = "card-note", style = "color:#8a6d3b;",
                 tags$ul(lapply(state$warnings, tags$li)))
    )
  })

  # The run folder is rsynced to the cluster wholesale and becomes INPUTPATH, so
  # writing the file beside sample.json is all it takes to reach the pipeline.
  # Nothing here is silent: a file that is not copied is a file the submitter
  # believes is in the run, and the run would look normal without it.
  stage_platemap <- function(run_path) {
    upload <- input$platemap_file
    if (is.null(upload) || !length(upload$datapath) || !nzchar(upload$datapath[[1]])) {
      return(invisible(FALSE))
    }
    if (!identical(input$method, "tiago_sulkseq")) {
      showNotification(paste("The experimental groups file was not attached to this run:",
                             "only the MIBS pipeline reads it."),
                       type = "warning", duration = 12)
      return(invisible(FALSE))
    }
    state <- platemap_state()
    if (!isTRUE(state$ok)) {
      showNotification(
        paste0("The experimental groups file was not attached: ",
               if (length(state$errors)) state$errors[[1]] else "it could not be read."),
        type = "error", duration = 14)
      return(invisible(FALSE))
    }
    # file.copy() reports failure in its return value and nowhere else. It was
    # discarded, so a full disk or an unwritable run folder produced a success
    # notification, a sample.json saying the run has a design, and a cluster job
    # stopping on the guard in rules/tiago_sulkseq.smk with nobody expecting it.
    if (!isTRUE(file.copy(upload$datapath[[1]], file.path(run_path, "platemap.tsv"),
                          overwrite = TRUE))) {
      showNotification(paste("The experimental groups file could not be written into the",
                             "run folder, so it was not attached."),
                       type = "error", duration = 14)
      return(invisible(FALSE))
    }
    showNotification(sprintf("Experimental groups attached: %d sample%s covering %s.",
                             length(state$samples), plural(length(state$samples)),
                             platemap_coverage(state)),
                     type = "message", duration = 8)
    invisible(TRUE)
  }

  # ---- preflight ---------------------------------------------------------
  # scripts/preflight.bash emits one finding per line as
  # SEVERITY|CODE|MESSAGE|HINT and decides nothing itself. Splitting the parse
  # from the policy keeps the checks testable from a shell.
  parse_preflight <- function(lines) {
    out <- list(fail = list(), warn = list(), ok = list())
    for (line in lines) {
      parts <- strsplit(line, "|", fixed = TRUE)[[1]]
      if (length(parts) < 3) next
      finding <- list(code = parts[[2]], message = parts[[3]],
                      hint = if (length(parts) >= 4) parts[[4]] else "")
      key <- switch(toupper(parts[[1]]), FAIL = "fail", WARN = "warn", OK = "ok", NULL)
      if (!is.null(key)) out[[key]] <- c(out[[key]], list(finding))
    }
    out
  }

  preflight_html <- function(findings) {
    paste0("<ul style='text-align:left;margin:8px 0 0;padding-left:20px'>",
           paste0(vapply(findings, function(f) {
             paste0("<li style='margin-bottom:6px'>", htmltools::htmlEscape(f$message),
                    if (nzchar(f$hint))
                      paste0("<br><span style='opacity:.75'>",
                             htmltools::htmlEscape(f$hint), "</span>")
                    else "", "</li>")
           }, character(1)), collapse = ""),
           "</ul>")
  }

  run_preflight <- function() {
    script <- app_script("preflight.bash")
    if (!file.exists(script)) return(list(fail = list(), warn = list(), ok = list()))
    # Hard cap. Every check inside is itself bounded, but the app must not be
    # hostage to one that is not.
    result <- run_script(script, timeout_sec = 90, args = c(
      "-s", input$server_name %||% "",
      "-i", gsub("^/+|/+$", "", trimws(input$input_folder %||% "")),
      "-p", input$method %||% "",
      "-b", input$primer %||% "",
      "-t", if (isTRUE(as.logical(input$trim_primers))) (input$trim_sequence %||% "") else "",
      "-n", current_run_name() %||% "",
      "-u", input$user_id %||% ""
    ))
    findings <- parse_preflight(result$output)
    # A preflight that cannot run -- or that ran out of time -- must not block a
    # run. It is a courtesy check, and the pipeline's own guards are still in
    # place behind it.
    if (result$status != 0 && !length(findings$fail) && !length(findings$warn)) {
      return(list(fail = list(), warn = list(list(
        code = "preflight_unavailable",
        message = "The pre-submission check could not run, so this run was not checked in advance.",
        hint = "The run was submitted anyway; the pipeline's own checks still apply.")),
        ok = list()))
    }
    findings
  }

  # timeout_sec bounds the wall clock. Shiny is single threaded, so a script
  # that blocks does not just delay its own submission -- it freezes the app for
  # everyone connected. That happened: a preflight integrity check ran `gzip -t`
  # over multi-gigabyte FASTQs on the NAS and held the process for six minutes.
  # 0 keeps the old unbounded behaviour for the long-running dispatch scripts,
  # which are expected to take minutes and have nothing else to block.
  run_script <- function(script_path, args = character(), timeout_sec = 0) {
    command_name <- script_path
    command_args <- args

    if (grepl("\\.bash$", script_path, ignore.case = TRUE)) {
      command_name <- "bash"
      command_args <- c(script_path, args)
    }

    output <- tryCatch(
      suppressWarnings(system2(command_name, args = command_args, stdout = TRUE,
                               stderr = TRUE, timeout = timeout_sec)),
      error = function(e) structure(conditionMessage(e), status = 1)
    )

    status <- attr(output, "status")
    if (is.null(status)) {
      status <- 0
    }

    list(
      status = status,
      output = paste(output, collapse = "\n")
    )
  }

  parse_numeric_marker <- function(output_text, marker_name) {
    if (is.null(output_text) || !nzchar(output_text)) {
      return(NA_real_)
    }

    marker_match <- regmatches(
      output_text,
      regexpr(sprintf("%s=([0-9]+)", marker_name), output_text, perl = TRUE)
    )

    if (length(marker_match) == 0 || !nzchar(marker_match)) {
      return(NA_real_)
    }

    as.numeric(sub(sprintf("^%s=", marker_name), "", marker_match))
  }

  # Snakemake prints a long DAG/job-table preamble before doing any work, so the
  # first N characters of a failed run are almost never the reason it failed.
  # Keep the tail instead, and write the whole thing to a per-run log file.
  BACKEND_OUTPUT_KEEP_CHARS <- 4000

  tail_chars <- function(text, n = BACKEND_OUTPUT_KEEP_CHARS) {
    text <- as.character(text)
    if (length(text) < 1 || is.na(text[[1]])) {
      return("")
    }
    text <- text[[1]]
    if (nchar(text) <= n) {
      return(text)
    }
    paste0("[... truncated, full log in runs/<user>/<run>/backend.log ...]\n",
           substr(text, nchar(text) - n + 1, nchar(text)))
  }

  # What the run history and the emails say about a submission, read once when
  # it is submitted. The run now finishes long after the request that started
  # it, often after its session has gone, when input$ no longer describes it.
  submission_snapshot <- function() {
    list(
      user_id = input$user_id,
      email = input$email,
      run_name = current_run_name(),
      input_method = input$type_input,
      server_name = input$server_name,
      pipeline = input$method,
      # The run log shows chain and barcode set for every historical run; live
      # runs had neither recorded, so those columns were blank for anything
      # submitted through the app.
      primer_set = input$primer,
      igblast_panel = input$igblast_panel,
      # Two runs of the same library with and without a design give different
      # clone counts, and nothing in the log said which had one. Read exactly
      # as sample.json reads it, pipeline included: the other three never
      # receive the file, so recording true for them would be a log entry no
      # run folder backs.
      platemap_attached = identical(input$method, "tiago_sulkseq") &&
                            isTRUE(platemap_state()$ok),
      sequencing_kit = input$sequencing_kit %||% "",
      archive_format = input$archive_format
    )
  }

  append_run_history <- function(status_label, runtime_sec, input_bytes = NA_real_, backend_output = "",
                                 snap = submission_snapshot()) {
    plates_used <- parse_numeric_marker(backend_output, "TBRC_PLATES_USED")
    if (!length(plates_used) || !is.finite(plates_used)) plates_used <- NA_real_
    # Persist only a compact summary per run; the full workflow log is in
    # backend.log next to the run metadata.
    append_run_stat(
      data.frame(
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        user_id = snap$user_id,
        run_name = snap$run_name,
        input_method = snap$input_method,
        server_name = snap$server_name,
        pipeline = snap$pipeline,
        primer_set = snap$primer_set,
        igblast_panel = snap$igblast_panel,
        platemap_attached = snap$platemap_attached,
        sequencing_kit = snap$sequencing_kit,
        # Counted from the delivered wells by final.script.sh, not typed on the
        # form. A submitter does not reliably know this and should not have to.
        plates_used = as.numeric(plates_used),
        archive_format = snap$archive_format,
        input_bytes = as.numeric(input_bytes),
        runtime_sec = as.numeric(runtime_sec),
        status = as.character(status_label),
        backend_output = tail_chars(backend_output),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    )
  }

  estimated_input_bytes <- reactive({
    input_mode <- input$type_input
    if (is.null(input_mode) || !identical(input_mode, "Upload") || is.null(input$input_files$datapath)) {
      return(NA_real_)
    }

    upload_sizes <- suppressWarnings(file.info(input$input_files$datapath)$size)
    upload_sizes <- upload_sizes[is.finite(upload_sizes) & !is.na(upload_sizes)]
    if (!length(upload_sizes)) {
      return(NA_real_)
    }

    sum(upload_sizes)
  })

  run_metrics <- reactive({
    metrics <- compute_run_metrics(
      read_run_stats(),
      input_bytes = estimated_input_bytes(),
      manual_minutes_per_100mb = get_metric_value(servers_config, "manual_minutes_per_100mb", 29),
      manual_hours_per_run = get_metric_value(servers_config, "manual_hours_per_run", 0),
      run_index_offset = get_metric_value(servers_config, "run_index_offset", 0)
    )
    observed_count <- historical_run_count()

    # The structured runtime ledger is new, so fall back to the number of
    # historical run folders until enough completed runs populate run_stats.tsv.
    if (isTRUE(observed_count > metrics$dataset_count)) {
      metrics$dataset_count <- observed_count
      if (is.na(metrics$total_time_saved_sec) || metrics$total_time_saved_sec <= 0) {
        metrics$total_time_saved_sec <- observed_count * 6 * 60
      }
    }

    metrics
  })

  send_mail_message <- function(recipient, subject, body, context = "general") {
    if (is.null(recipient) || length(recipient) < 1 || is.na(recipient[[1]])) {
      recipient <- ""
    }
    recipient <- trimws(as.character(recipient[[1]]))
    if (!nzchar(recipient)) {
      append_notification_log(recipient, subject, "skipped", context, "No recipient provided")
      return(list(status = 0, message = "No email recipient provided."))
    }

    if (nzchar(Sys.which("mail"))) {
      output <- tryCatch(
        suppressWarnings(system2("mail", args = c("-s", subject, recipient), input = body, stdout = TRUE, stderr = TRUE)),
        error = function(e) structure(conditionMessage(e), status = 1)
      )
      status <- attr(output, "status")
      if (is.null(status)) {
        status <- 0
      }
      message_text <- paste(output, collapse = "\n")
      append_notification_log(recipient, subject, if (status == 0) "accepted" else "failed", context, message_text)
      return(list(status = status, message = message_text))
    }

    if (nzchar(Sys.which("sendmail"))) {
      message_text <- sprintf("To: %s\nSubject: %s\n\n%s\n", recipient, subject, body)
      output <- tryCatch(
        suppressWarnings(system2("sendmail", args = c("-t"), input = message_text, stdout = TRUE, stderr = TRUE)),
        error = function(e) structure(conditionMessage(e), status = 1)
      )
      status <- attr(output, "status")
      if (is.null(status)) {
        status <- 0
      }
      output_text <- paste(output, collapse = "\n")
      append_notification_log(recipient, subject, if (status == 0) "accepted" else "failed", context, output_text)
      return(list(status = status, message = output_text))
    }

    append_notification_log(recipient, subject, "failed", context, "Neither mail nor sendmail is available on the Shiny host.")
    list(status = 1, message = "Neither mail nor sendmail is available on the Shiny host.")
  }

  send_run_notifications <- function(status_label, sample_name, storage_name, archive_format, submitter_name, submitter_email) {
    admin_email <- get_cluster_value(servers_config, "notification_admin_email")
    user_subject <- sprintf("TBRC run %s %s", sample_name, status_label)
    user_body <- paste(
      sprintf("Your TBRC run %s is now %s.", sample_name, status_label),
      sprintf("Storage target: %s.", storage_name),
      sprintf("Archive format: %s.", archive_format),
      sep = "\n"
    )
    admin_subject <- sprintf("TBRC usage: %s ran %s (%s)", submitter_name, sample_name, status_label)
    admin_body <- paste(
      sprintf("TBRC was used by: %s", submitter_name),
      sprintf("User email: %s", submitter_email),
      sprintf("Run name: %s", sample_name),
      sprintf("Status: %s", status_label),
      sprintf("Storage target: %s", storage_name),
      sprintf("Archive format: %s", archive_format),
      sep = "\n"
    )

    list(
      user = send_mail_message(submitter_email, user_subject, user_body, context = "user_completion"),
      admin = send_mail_message(admin_email, admin_subject, admin_body, context = "admin_usage")
    )
  }

  # Starts goServer.bash or goUpload.bash in the background and returns at once,
  # so the app stays usable for everyone while the run goes. What used to follow
  # the blocking call -- run history, emails, the closing message -- is finish(),
  # called by the shared poller in global.R when the script exits.
  launch_submission <- function(script, cmd_args, path, started_at, input_bytes = NA_real_) {
    snap <- submission_snapshot()
    key <- paste(snap$user_id, snap$run_name, sep = "/")
    if (background_run_active(key)) {
      shinyalert(
        title = "Already running",
        text = sprintf(paste("Run %s from an earlier submission has not finished yet.",
                             "Wait for it, or submit under another run name."), snap$run_name),
        type = "error"
      )
      return(invisible(FALSE))
    }

    finish <- function(status, output) {
      bytes <- if (is.finite(input_bytes)) input_bytes else parse_numeric_marker(output, "TBRC_INPUT_BYTES")
      label <- if (identical(status, 0L)) "completed" else "failed"
      append_run_history(
        status_label = label,
        runtime_sec = as.numeric(difftime(Sys.time(), started_at, units = "secs")),
        input_bytes = bytes,
        backend_output = output,
        snap = snap
      )
      email_result <- send_run_notifications(
        status_label = label,
        sample_name = snap$run_name,
        storage_name = snap$server_name,
        archive_format = snap$archive_format,
        submitter_name = snap$user_id,
        submitter_email = snap$email
      )
      # The page, only if the submitter is still on it.
      if (session$isClosed()) return(invisible())
      withReactiveDomain(session, isolate({
        if (identical(label, "failed")) {
          failure_text <- summarize_backend_failure(
            output,
            "The server-side workflow failed. Please check the backend logs."
          )
          shinyalert(
            title = if (grepl("Folder not found on", output, fixed = TRUE)) "Folder not found" else "Submission failed",
            text = failure_text,
            type = "error"
          )
          set_progress_detail(sprintf("Run %s failed.", snap$run_name), failure_text)
        } else {
          set_progress_detail(
            sprintf("Run %s finished.", snap$run_name),
            "Look for pre-trim and post-trim histograms, plate QC, and the .imgt.ready.fasta file in the result bundle.",
            sprintf(
              "If the selected storage server is in pull mode, Synology should fetch the %s archive from the HPC results folder.",
              snap$archive_format
            )
          )
          pgPaneUpdate("thispg", "workflow", 100)
          pgPaneUpdate("thispg", "export", 100)
          if (email_result$user$status != 0 || email_result$admin$status != 0) {
            set_progress_detail(progress_detail(), "Shiny-side completion email could not be sent.")
          }
        }
      }))
    }

    start_background_script(key, script, cmd_args,
                            log_path = file.path(path, "backend.log"),
                            status_path = file.path(path, "backend.exit"),
                            on_done = finish)
    pgPaneUpdate("thispg", "dispatch", 100)
    pgPaneUpdate("thispg", "workflow", 20)
    set_progress_detail(
      sprintf("Run %s submitted.", snap$run_name),
      "It is running on the cluster now. This page stays usable, and you can close it.",
      if (nzchar(snap$email %||% "")) sprintf("An email goes to %s when it finishes.", snap$email)
      else "No email address was given, so watch the results folder or the run history for it."
    )
    invisible(TRUE)
  }

  write_submission_metadata <- function(path) {
    metadata_df <- submission_row()
    metadata_list <- as.list(metadata_df[1, , drop = TRUE])
    jsonlite::write_json(
      metadata_list,
      path = file.path(path, "sample.json"),
      auto_unbox = TRUE,
      pretty = TRUE,
      null = "null"
    )
  }

  summarize_backend_failure <- function(output_text, default_message) {
    if (is.null(output_text) || !nzchar(output_text)) {
      return(default_message)
    }

    if (grepl("Folder not found on", output_text, fixed = TRUE)) {
      return("The requested folder was not found on the selected server.")
    }

    if (grepl("Missing required arguments\\.", output_text) ||
        grepl("Error in rule", output_text, fixed = TRUE) ||
        grepl("Submitted batch job", output_text, fixed = TRUE) ||
        grepl("Missing output files:", output_text, fixed = TRUE) ||
        grepl("Complete log:", output_text, fixed = TRUE) ||
        grepl("server-side submission failed", output_text, fixed = TRUE)) {
      return(default_message)
    }

    output_text
  }

  output$run_local_control <- renderUI({
    if (isTRUE(LOCAL_BACKEND$available)) {
      app_switch_input("run_local", "Run on this server (not HPC)", FALSE, status = "success")
    } else {
      tags$div(
        class = "card-note",
        style = "margin-bottom:10px;",
        sprintf("Runs go to the HPC cluster. Running on this server is unavailable: %s.", LOCAL_BACKEND$reason)
      )
    }
  })

  # Follow the barcode set: selecting a T-cell set must not leave a B-cell handle
  # in the trim field. Six recorded T-cell runs trimmed with the Ig handle and
  # discarded their libraries, and the app gave no indication.
  observeEvent(input$primer, {
    default_seq <- barcode_default(servers_config, input$primer, "trim_seq")
    if (!is.na(default_seq)) {
      updateTextInput(session, "trim_sequence", value = default_seq)
    }
    # And the IgBLAST scope, but only across the B/T divide: a T set left on
    # an Ig panel is caught by the backend's receptor check, but the form
    # should not offer the mistake. A narrower choice within the right class
    # (igh, trb) is the submitter's to make and is left alone.
    chain <- barcode_default(servers_config, input$primer, "chain")
    panel <- isolate(input$igblast_panel) %||% ""
    tcr_panels <- c("tra", "trb", "tcr_all")
    if (identical(chain, "T") && !(panel %in% c(tcr_panels, "all_receptors"))) {
      app_update_picker_input(session, "igblast_panel", selected = "tcr_all")
    } else if (identical(chain, "B") && panel %in% tcr_panels) {
      app_update_picker_input(session, "igblast_panel", selected = "ig_all")
    }
  }, ignoreInit = FALSE)

  output$chainHint <- renderText({
    chain <- barcode_default(servers_config, input$primer, "chain")
    if (is.na(chain)) {
      "Unrecognised barcode set - check the trim sequence yourself."
    } else if (identical(chain, "T")) {
      "T cell set. The trim sequence has been set to the TCR handle and the IgBLAST scope to TCR."
    } else {
      "B cell set. The trim sequence has been set to the immunoglobulin handle."
    }
  })

  # Light-chain rescue used to force classic_luka here. IgBlast annotation
  # reversed that: on RUN10, tiago_vsearch is 79.0% productive against 61.1%, with
  # J called on 95.8% of sequences against 69.2% and a CDR3 on 93.5% against
  # 68.8% -- and the gain holds at every level of read support, not just the thin
  # wells. Its one real cost is 7 wells with 50+ reads that classic delivers and
  # it does not, so the choice is left open and the trade stated.
  observeEvent(input$luka_light, {
    if (identical(input$luka_light, "luka_light") && identical(input$method, "tiago_vsearch")) {
      showNotification(
        paste(
          "Light-chain rescue with tiago_vsearch: more productive sequences",
          "(79% vs 61%) but it drops a handful of well-supported wells that",
          "classic_luka keeps. Switch to classic_luka if well recovery matters more."
        ),
        type = "message", duration = 12
      )
    }
  }, ignoreInit = TRUE)

  runlogServer("runlog")
  changelogServer("changelog")

  # Update the choices for the folder dropdown
  observe({
    app_update_picker_input(session, "folder_dropdown", choices = folder_list())
  })

  observe({
    user_choices <- names(known_users())
    updateSelectizeInput(session, "user_id", choices = user_choices, selected = isolate(input$user_id), server = TRUE)
  })

  observeEvent(input$user_id, {
    selected_user <- ""
    if (!is.null(input$user_id)) {
      selected_user <- trimws(input$user_id)
    }
    remembered_email <- unname(known_users()[selected_user])

    if (length(remembered_email) > 0 && !is.na(remembered_email) && nzchar(remembered_email)) {
      updateTextInput(session, "email", value = remembered_email)
    }
  }, ignoreNULL = FALSE)

  observeEvent(input$run_clonality, {
    if (isTRUE(input$run_clonality) && !isTRUE(input$run_igblast)) {
      updateCheckboxInput(session, "run_igblast", value = TRUE)
    }
  }, ignoreInit = TRUE)
  
  # Reactive to read content of the selected barcode files
  barcode_content <- reactive({
    if (!is.null(input$folder_dropdown) && nzchar(input$folder_dropdown)) {
      folder_path <- file.path(BC_DIR, input$folder_dropdown)
      
      # Check if both files exist in the selected folder
      req(file.exists(file.path(folder_path, "BC1.txt")), file.exists(file.path(folder_path, "BC2.txt")))
      
      list(
        BC1 = readLines(file.path(folder_path, "BC1.txt"), warn = FALSE),
        BC2 = readLines(file.path(folder_path, "BC2.txt"), warn = FALSE)
      )
    } else {
      list(BC1 = NULL, BC2 = NULL)
    }
  })
  
  #Render barcode content files in screen
  output$barcode_display1 <- renderText({
    bc <- barcode_content()
    paste(c("Barcode 1:", bc$BC1), collapse = "\n")
  })
  
  output$barcode_display2 <- renderText({
    bc <- barcode_content()
    paste(c("Barcode 2:", bc$BC2), collapse = "\n")
  })
  
  
  output$log <- renderText({
    log()
  })
  
  output$progressDetail <- renderText({
    progress_detail()
  })
  
  observe({
    # If any fastq file is uploaded, set type_input to 'Upload' and set server_name to 'upload'
    if(!is.null(input$input_files$datapath)){
      app_update_radio_buttons(
        session,
        inputId = "type_input",
        choices = c("Upload", "Server"),
        selected = "Upload"
      )
      app_update_radio_buttons(
        session,
        inputId = "server_name",
        choices = c(server_names, "upload"),
        selected = "upload"
      )
    }
    
    # Depending on the value of type_input, update the text input and server_name input
    if(input$type_input == "Upload"){
      updateTextInput(
        session,
        inputId = "input_folder",
        value = "Not used when uploading..."
      )
      app_update_radio_buttons(
        session,
        inputId = "server_name",
        choices = c(server_names, "upload"),
        selected = "upload"
      )
    } else if(input$type_input == "Server") {
      
      updateTextInput(
        session,
        inputId = "input_folder",
        value = ""
      )
      app_update_radio_buttons(
        session,
        inputId = "server_name",
        choices = c(server_names, "upload")
      )
    }
  })
  
  
  # Observer to update the primer set selection if barcodes are uploaded
  observe({
    app_update_picker_input(
      session,
      inputId = "primer",
      choices = bc_names(),
      selected = "B2"
    )
  })
  
  output$datasetCountDisplay <- renderText({
    m <- run_metrics()
    idx <- m$run_index
    if (is.null(idx) || !is.finite(idx)) {
      sprintf("%d datasets tracked", m$dataset_count)
    } else {
      sprintf("Run #%s", format(idx, big.mark = ",", trim = TRUE))
    }
  })

  output$datasetCountNote <- renderText({
    m <- run_metrics()
    sprintf("%d completed runs on record", m$dataset_count)
  })

  # ---- library cost ------------------------------------------------------
  # Sequencing off the core facility's rate card, plus what the bench work cost
  # per barcode plate. Both are estimates; the card says so rather than
  # implying an invoice.
  observe({
    choices <- cost_sequencing_choices(servers_config)
    if (!length(choices)) return()
    # A long-read run was not sequenced on a MiSeq, so start it on the service
    # this lab actually uses rather than making the user notice and change it.
    selected <- if (identical(input$method, "tiago_nanopore") &&
                    "plasmidsaurus_custom" %in% choices) {
      "plasmidsaurus_custom"
    } else {
      as.character(cost_config(servers_config)$default_sequencing %||% choices[[1]])
    }
    app_update_picker_input(session, "sequencing_kit",
                            choices = as.list(choices), selected = selected)
  })

  # Plate counts as measured by the backend, keyed by library name.
  plates_by_run <- reactive({
    stats <- tryCatch(read_run_stats(), error = function(...) NULL)
    if (is.null(stats) || !nrow(stats) || is.null(stats$plates_used)) return(NULL)
    values <- suppressWarnings(as.numeric(stats$plates_used))
    keys <- trimws(tolower(as.character(stats$run_name)))
    keep <- nzchar(keys) & is.finite(values) & values > 0
    if (!any(keep)) return(NULL)
    # A library re-run keeps its largest measured plate count.
    tapply(values[keep], keys[keep], max)
  })

  library_cost <- reactive({
    log <- tryCatch(read_run_log(), error = function(...) NULL)
    cost_all_runs(servers_config, log, input$sequencing_kit, plates_by_run())
  })

  output$libraryCostDisplay <- renderText({
    summary <- library_cost()
    currency <- as.character(cost_config(servers_config)$currency %||% "USD")
    if (is.null(summary) || !is.finite(summary$total)) {
      return("Estimated cost: not configured")
    }
    sprintf("%s across %s libraries", format_money(summary$total, currency),
            format(summary$distinct_runs, big.mark = ","))
  })

  output$libraryCostNote <- renderText({
    summary <- library_cost()
    currency <- as.character(cost_config(servers_config)$currency %||% "USD")
    if (is.null(summary) || !is.finite(summary$total)) {
      return("No prices set \u2014 fill in costs in server_config.json")
    }
    sprintf("%s sequencing + %s bench, over %s barcode plates",
            format_money(summary$sequencing, currency),
            format_money(summary$prep, currency),
            format(summary$total_plates, big.mark = ""))
  })

  # The marginal number: what one more barcode plate costs once the flowcell is
  # shared. Prep is spent per plate whatever happens; the sequencing charge
  # divides across the plates pooled onto the run.
  output$libraryCostPerPlate <- renderText({
    summary <- library_cost()
    currency <- as.character(cost_config(servers_config)$currency %||% "USD")
    # Only meaningful when the sequencing charge is per flowcell. A per-project
    # service is not bought by the plate, and adding a per-library price to a
    # per-plate one would print a number in mixed units.
    if (is.null(summary) || !is.finite(summary$per_plate) ||
        !identical(summary$pricing, "flat")) return("")
    sanger <- cost_sanger_equivalent(servers_config)
    line <- sprintf("%s per plate at %g plates pooled per flowcell",
                    format_money(summary$per_plate, currency), summary$pooled)
    if (is.finite(sanger)) {
      line <- paste0(line, sprintf("; Sanger on the same 96 wells: %s",
                                   format_money(sanger, currency)))
    }
    line
  })

  output$libraryCostAllRuns <- renderText({
    summary <- library_cost()
    if (is.null(summary)) return("")
    if (summary$measured_runs >= summary$distinct_runs) {
      return("Plate counts measured from the delivered wells \u2014 estimate only")
    }
    sprintf("Plates counted for %s of %s libraries; the rest assumed to be one, so this is a floor",
            format(summary$measured_runs, big.mark = ","),
            format(summary$distinct_runs, big.mark = ","))
  })

  # Render form table
  # The Run button, rendered from readiness rather than always live.
  #
  # A button you can always press, that then tells you what is wrong, teaches
  # people to press it to find out. A button that states its own readiness --
  # and, when not ready, exactly what is missing -- means the answer is visible
  # before anyone clicks.
  output$run_control <- renderUI({
    checklist <- submission_checklist()
    missing <- Filter(function(x) !x$done, checklist)

    rows <- lapply(checklist, function(x) {
      tags$div(
        class = if (x$done) "rdy-item rdy-done" else "rdy-item rdy-todo",
        tags$span(class = "rdy-mark", if (x$done) "\u2713" else "\u25CB"),
        tags$span(class = "rdy-label", x$label),
        if (!x$done) tags$span(class = "rdy-todo-text", x$todo) else NULL
      )
    })

    tagList(
      tags$div(class = "rdy-box", rows),
      if (length(missing) == 0) {
        app_action_button("run", "Run", style = "material-flat",
                          color = "success", size = "md", block = TRUE)
      } else {
        tagList(
          tags$button(class = "btn btn-block rdy-blocked", disabled = NA,
                      sprintf("%d thing%s still needed",
                              length(missing), if (length(missing) == 1) "" else "s")),
          tags$div(class = "rdy-hint",
                   "The Run button becomes available once the list above is complete.")
        )
      }
    )
  })

  # Rendered as a static list rather than an interactive table. Sorting was
  # available and actively harmful here: these settings have a meaningful order,
  # and re-sorting them alphabetically by value tells the reader nothing while
  # destroying the grouping. Paging and search are noise on a list this short.
  output$submission_form <- DT::renderDataTable({
    DT::datatable(
      form(),
      rownames = FALSE,
      colnames = c("Setting", "Value"),
      selection = "none",
      class = "compact stripe",
      options = list(
        ordering = FALSE,   # the order is the point; do not let it be changed
        paging = FALSE,
        searching = FALSE,
        info = FALSE,
        scrollY = "420px",
        scrollCollapse = TRUE,
        columnDefs = list(
          list(targets = 0, width = "230px"),
          list(targets = "_all", className = "dt-left")
        )
      )
    )
  })
  
  #Barcode Uploader
  output$confirmation <- renderText({
    "Upload barcode files."
  })
  
  observeEvent(input$bcupload,{
    
    # Check if the two BC files are uploaded
    if (is.null(input$BC_files) | length(input$BC_files$datapath) != 2) {
      
      output$confirmation <- renderText({
        "Please upload both BC1 and BC2 files."
      })
      return()
    }
    
    if (!is_valid_path_segment(input$folder_name)) {
      output$confirmation <- renderText({
        invalid_segment_message("Folder name")
      })
      return()
    }

    if (!is_valid_path_segment(input$bc_name)) {
      output$confirmation <- renderText({
        invalid_segment_message("Barcode name")
      })
      return()
    }

    base_dir <- BC_DIR
    new_bc_dir <- file.path(base_dir, input$folder_name)

    if (dir.exists(new_bc_dir)) {
      output$confirmation <- renderText({
        paste0("A barcode folder named ", input$folder_name, " already exists. Choose another name.")
      })
      return()
    }
    
    # Save the uploaded BC files to the folder
    # We will just use names BC1.txt and BC2.txt for simplicity. Ensure users upload them in the correct order.
    # Validate the CONTENT before anything is written. These bytes are stored on
    # the cluster and rsynced to lab storage, so accepting an arbitrary file here
    # would make this a way to park anything anywhere in the lab's infrastructure.
    # Checking the name would not help: the copy is always written as BC1.txt or
    # BC2.txt, so the uploaded name is already discarded.
    bc_problems <- c(
      validate_barcode_file(input$BC_files$datapath[1], "The first barcode file (BC1)"),
      validate_barcode_file(input$BC_files$datapath[2], "The second barcode file (BC2)")
    )
    if (length(bc_problems)) {
      output$confirmation <- renderText({
        paste("Barcode files rejected.", paste(bc_problems, collapse = " "))
      })
      return()
    }

    # Only now, once the files are known to be good. Created before the check,
    # a rejected upload left an empty folder that the barcode menu then offered
    # and that blocked the same name from being uploaded again.
    dir.create(new_bc_dir)
    file.copy(input$BC_files$datapath[1], file.path(new_bc_dir, 'BC1.txt'))
    file.copy(input$BC_files$datapath[2], file.path(new_bc_dir, 'BC2.txt'))
    
    # A set whose barcodes sit within 2 positions of each other cannot tolerate a
    # sequencing error: the demultiplexer breaks the tie by file order, silently.
    warnings <- character()
    for (which_bc in c("BC1", "BC2")) {
      d <- barcode_min_distance(file.path(new_bc_dir, paste0(which_bc, ".txt")))
      if (!is.na(d) && d < 3) {
        warnings <- c(warnings, sprintf(
          "%s barcodes are only %d apart; reads carrying a sequencing error in them cannot be assigned safely, so they will be demultiplexed with exact matching.",
          which_bc, d
        ))
      }
    }

    # Append barcode name and folder name to bc_names.txt
    cat(paste('"', input$bc_name, '" "', input$folder_name, '"\n', sep = ""), file = file.path(base_dir, 'bc_names.txt'), append = TRUE)
    
      output$confirmation <- renderText({
        base <- paste("Uploaded", input$bc_name, "with folder name", input$folder_name)
        if (length(warnings)) paste(c(base, warnings), collapse = " ") else base
      })
    
    
    
    # run_script() reports a real exit status; the previous system2() call
    # returned output, so the status check below could never fire.
    sync_result <- run_script(app_script("updatebarcodes.bash"))

    if (sync_result$status != 0) {
      output$confirmation <- renderText({
        paste0(
          "Barcodes saved locally, but syncing them to the pipeline failed. ",
          "They will not be available to runs until this is resolved."
        )
      })
      return()
    }
    
  })
  
  # Observer to perform actions when the Run button is clicked
  observeEvent(input$run, {
    if (isTRUE(input$run_clonality) && !isTRUE(input$run_igblast)) {
      shinyalert(
        title = "IGBlast required",
        text = "Clonality analysis needs IGBlast output. Please enable Run IGBlast annotation.",
        type = "error"
      )
      return()
    }

    if (!is_valid_path_segment(input$user_id)) {
      shinyalert(
        title = "Check the user name",
        text = invalid_segment_message("User name"),
        type = "error"
      )
      return()
    }

    # The address is handed to the mail command as a recipient. Left unchecked,
    # anyone could send mail from this server to any address they liked.
    if (!is_valid_email(input$email)) {
      shinyalert(
        title = "Check the email address",
        text = paste("Enter a valid email address, or leave the field empty",
                     "to skip the notification."),
        type = "error"
      )
      return()
    }

    if (input$type_input == "Server" && !nzchar(trimws(input$input_folder))) {
      shinyalert(
        title = "Oops!",
        text = "Folder in server cannot be empty for Server submissions.",
        type = "error"
      )
      return()
    }

    if (input$type_input == "Server" && !is_valid_path_segment(input$input_folder)) {
      shinyalert(
        title = "Check the folder name",
        text = invalid_segment_message("Folder in server"),
        type = "error"
      )
      return()
    }

    # It becomes a folder on the storage server, inside a remote command.
    run_group <- trimws(input$run_group %||% "")
    if (nzchar(run_group) && !is_valid_path_segment(run_group)) {
      shinyalert(
        title = "Check the results subfolder",
        text = invalid_segment_message("Results subfolder"),
        type = "error"
      )
      return()
    }

    # Check if there are missing arguments in the form
    required_form_fields <- required_fields_for(input$method, input$type_input,
                                                names(submission_row()),
                                                input$trim_primers)

    required_values <- submission_row()[1, required_form_fields, drop = TRUE]

    if(any(is.na(required_values) | required_values == "")){
      # Named, not counted. "There are arguments missing" sent people back to a
      # form of twenty fields to find which one, and the MIBS numeric boxes made
      # that worse because an empty box looks filled in until you look at it.
      # The labels are the ones the preview uses, so a field is called the same
      # thing in both places.
      missing_fields <- required_form_fields[is.na(required_values) | required_values == ""]
      missing_labels <- ifelse(missing_fields %in% names(PREVIEW_FIELDS),
                               PREVIEW_FIELDS[missing_fields],
                               gsub("_", " ", missing_fields))
      shinyalert(
        title = "Oops!",
        text = paste0("There are arguments missing. Please complete the form: ",
                      paste(unname(missing_labels), collapse = ", "), "."),
        type = "error"
      )
      return()
    }

    # The check above only tests for NA and for "", which a typed 0 is neither
    # of. A MIBS number outside its range is discarded by the Snakefile and the
    # configured default is used in its place, with nothing in sample.json, the
    # log or summary.tsv recording the substitution -- so it has to be refused
    # here, where the submitter can still see what they typed.
    if (identical(input$method, "tiago_sulkseq")) {
      bad_mibs <- mibs_bad_numbers()
      if (length(bad_mibs)) {
        shinyalert(
          title = "Oops!",
          text = paste0("These MIBS options are empty or outside the range they allow: ",
                        paste(unname(PREVIEW_FIELDS[bad_mibs]), collapse = ", "),
                        ". The pipeline replaces a value it cannot use with its own ",
                        "default rather than refusing it, so fix them before running."),
          type = "error"
        )
        return()
      }
    }
    pgPaneUpdate("thispg", "validate", 100)
    pgPaneUpdate("thispg", "stage", 0)
    pgPaneUpdate("thispg", "dispatch", 0)
    pgPaneUpdate("thispg", "workflow", 0)
    pgPaneUpdate("thispg", "export", 0)
    
    # If Upload mode is on
    if(input$type_input == "Upload") {
      run_started_at <- Sys.time()
      set_progress_detail(
        "Validation complete.",
        "Preparing upload run folder.",
        "Copying uploaded FASTQ files into the local run directory.",
        "Submitting workflow launch command.",
        "Live rule-by-rule Snakemake status is not available in this UI yet."
      )
      pgPaneUpdate("thispg", "stage", 25)
      pgPaneUpdate("thispg", "stage", 65)
      pgPaneUpdate("thispg", "stage", 100)
      
      # Alert if no fastq files are uploaded
      if(is.null(input$input_files)){
        shinyalert(title = "Oops!", text = "No fastq files detected! Upload or change input method.", type = "error")
        return()
      }

      if (length(input$input_files$name) != 2) {
        shinyalert(title = "Oops!", text = "Upload mode requires exactly two FASTQ files: one R1 and one R2.", type = "error")
        return()
      }

      # The name of an uploaded file is chosen by whoever uploads it and used to
      # reach a command run over ssh on the login node. goUpload.bash re-checks
      # this; neither layer should be the only one standing there.
      bad_names <- input$input_files$name[
        !vapply(input$input_files$name, is_valid_fastq_filename, logical(1))
      ]
      if (length(bad_names)) {
        shinyalert(
          title = "Oops!",
          text = paste0(
            "These file names cannot be accepted: ",
            paste(bad_names, collapse = ", "),
            ". Use only letters, numbers, dot, dash and underscore, ending in ",
            ".fastq, .fq, .fastq.gz or .fq.gz."
          ),
          type = "error"
        )
        return()
      }

      upload_roles <- vapply(input$input_files$name, classify_uploaded_read, character(1))
      if (!all(c("R1", "R2") %in% upload_roles)) {
        shinyalert(title = "Oops!", text = "Uploaded files must include one R1 and one R2 FASTQ file.", type = "error")
        return()
      }

      r1_index <- which(upload_roles == "R1")[1]
      r2_index <- which(upload_roles == "R2")[1]
      
      # processing path and creating directory
      path <- file.path(RUNS_DIR, input$user_id, current_run_name())
      dir.create(path, showWarnings = T, recursive = T)
      
      # Writing the submission metadata on the run folder
      write_submission_metadata(path)
      # sample.json already says whether this run has a design, and the backend
      # stops a MIBS run that says yes and arrives without one. When the copy
      # itself failed there is nothing to gain from dispatching a run that is
      # going to stop at DAG build with a message about a missing file.
      # Staged first and tested afterwards, never inside the condition: && stops
      # at the first FALSE, and a run that is not MIBS has to reach
      # stage_platemap() to be told its file will not be used.
      staged <- stage_platemap(path)
      if (isTRUE(submission_row()$platemap_attached) && !isTRUE(staged)) {
        shinyalert(
          title = "Oops!",
          text = paste("The experimental groups file could not be written into the run",
                       "folder. This submission records that the run has one and the",
                       "pipeline stops when it does not arrive, so the run was not started."),
          type = "error"
        )
        return()
      }
      # Shiny deletes an upload's temporary files when its session ends, which
      # can now come before goUpload.bash has copied them. Put them where that
      # script would copy them to (runs/<user>/<R1 name minus its extension>/)
      # before it starts; it skips a copy onto the same file.
      upload_dir <- file.path(RUNS_DIR, input$user_id,
                              sub("\\.[^.]*$", "", input$input_files$name[r1_index]))
      dir.create(upload_dir, showWarnings = FALSE, recursive = TRUE)
      staged_reads <- file.path(upload_dir, input$input_files$name[c(r1_index, r2_index)])
      if (!all(file.copy(input$input_files$datapath[c(r1_index, r2_index)], staged_reads,
                         overwrite = TRUE))) {
        shinyalert(title = "Oops!", text = "The uploaded FASTQ files could not be staged for the run.",
                   type = "error")
        return()
      }
      cmd_args <- c(
        "-u", input$user_id,
        "-d", staged_reads[1],
        "-e", staged_reads[2],
        "-n", input$input_files$name[r1_index],
        "-m", input$input_files$name[r2_index]
      )

      launch_submission(app_script("goUpload.bash"), cmd_args, path, run_started_at,
                        input_bytes = estimated_input_bytes())
    }
    
    
    # If Server mode is on and all fields were filled, run the pipeline
    if(input$type_input == "Server" & !any(is.na(required_values) | required_values == "")) {
      run_started_at <- Sys.time()

      # Preflight, against the real storage server, before anything is
      # dispatched. Every guard this pipeline had until now fired AFTER
      # dispatch: the user waited, got a job failure, and had to read a cluster
      # log to learn the folder name was misspelled. A FAIL here stops the
      # submission and says what to fix; a WARN is shown and the run proceeds.
      set_progress_detail("Checking the run before submitting it.")
      preflight <- run_preflight()
      if (length(preflight$fail)) {
        shinyalert(
          title = "This run would not have worked",
          text = preflight_html(preflight$fail),
          type = "error",
          html = TRUE
        )
        set_progress_detail("Submission stopped by the preflight check.",
                            paste(vapply(preflight$fail, function(f) f$message,
                                         character(1)), collapse = " "))
        return()
      }
      if (length(preflight$warn)) {
        showNotification(
          HTML(paste0("<strong>Worth knowing before this runs</strong>",
                      preflight_html(preflight$warn))),
          type = "warning", duration = 20
        )
      }

      set_progress_detail(
        "Validation complete.",
        "Preparing remote submission.",
        "Writing sample metadata into the local run folder.",
        "The server job will pull FASTQ files, run Snakemake, then export results.",
        sprintf(
          "Final packaging will use %s.",
          submission_row()[["archive_format"]][[1]]
        ),
        "Live rule-by-rule Snakemake status is not available in this UI yet."
      )
      # Run progress updates and folder creation code here
      fname <- current_run_name() #folder on server name (Run name)
      path <- file.path(RUNS_DIR, input$user_id, fname) #full path to the local run
      dir.create(path, showWarnings = T, recursive = T) #create the folder
      pgPaneUpdate("thispg", "stage", 25)
      pgPaneUpdate("thispg", "stage", 65)
      pgPaneUpdate("thispg", "stage", 100)
      # Writing the submission metadata on the run folder
      write_submission_metadata(path)
      # Server mode is the default path, and the design was only ever staged on
      # the Upload branch -- so the file the submitter attached reached the
      # cluster for nobody. goServer.bash rsyncs this folder and makes it
      # INPUTPATH, so writing it here is the whole transfer.
      # sample.json already says whether this run has a design, and the backend
      # stops a MIBS run that says yes and arrives without one. When the copy
      # itself failed there is nothing to gain from dispatching a run that is
      # going to stop at DAG build with a message about a missing file.
      # Staged first and tested afterwards, never inside the condition: && stops
      # at the first FALSE, and a run that is not MIBS has to reach
      # stage_platemap() to be told its file will not be used.
      staged <- stage_platemap(path)
      if (isTRUE(submission_row()$platemap_attached) && !isTRUE(staged)) {
        shinyalert(
          title = "Oops!",
          text = paste("The experimental groups file could not be written into the run",
                       "folder. This submission records that the run has one and the",
                       "pipeline stops when it does not arrive, so the run was not started."),
          type = "error"
        )
        return()
      }

      # Running the pipeline script with arguments
      cmd_args <- c(
        "-u", input$user_id,
        "-s", input$server_name,
        # Trimmed, and stripped of slashes. The value goes straight into a remote
        # "test -d <base>/<folder>", so a trailing space or a pasted trailing slash
        # makes the check fail on a folder that plainly exists -- and the user is
        # told the folder was not found, which is worse than useless because it is
        # wrong. current_run_name() already trims; this path did not.
        "-f", gsub("^/+|/+$", "", trimws(input$input_folder %||% "")),
        "-n", fname
      )
      if (isTRUE(input$run_local) && isTRUE(LOCAL_BACKEND$available)) {
        cmd_args <- c(cmd_args, "-l")
      }
      
      
      launch_submission(app_script("goServer.bash"), cmd_args, path, run_started_at)
    }
    
  })
}
