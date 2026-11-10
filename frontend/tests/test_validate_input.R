# Tests for validate_input.R and run_script.R: the form refuses bad values, and
# nothing the form sends is handed to a shell as text. Also the run ledger and the
# completion steps in global.R.
#
# Run from frontend/ (testthat is not in the image, so these are base R):
#   docker run --rm -v "$WT:/repo" -w /repo/frontend \
#     -e TBRC_STATE_DIR=/tmp/s -e TBRC_RUNS_DIR=/tmp/r -e TBRC_BC_DIR=/tmp/b \
#     tbrc-frontend:latest Rscript tests/test_validate_input.R
#
# The shared good/bad values live in ../tests/validation_cases.json and are also
# read by the cluster-side tests, so the two validators cannot drift apart.

# The app creates its state folders on load. Point them at scratch whatever the
# caller's environment says, so no test can touch real runs or barcodes.
scratch <- tempfile("tbrc_test_")
dir.create(scratch)
Sys.setenv(
  TBRC_STATE_DIR = file.path(scratch, "state"),
  TBRC_RUNS_DIR  = file.path(scratch, "runs"),
  TBRC_BC_DIR    = file.path(scratch, "bc")
)

# The image prints "built under R version" warnings for its own packages.
suppressWarnings(suppressMessages(source("global.R")))
cases <- jsonlite::fromJSON("../tests/validation_cases.json", simplifyVector = TRUE)

# The app reads the list of barcode sets when a session starts.
writeLines(c('"bcname"', '"B cell2" "B2"', '"Real set" "real_set"'),
           file.path(BC_DIR, "bc_names.txt"))

SERVERS <- c("storage1", "storage2")
KITS <- c("miseq_nano_2x250", "plasmidsaurus_custom")
choices <- submission_choices(SERVERS, KITS)

check <- function(condition, what) {
  if (!isTRUE(condition)) stop("FAILED: ", what, call. = FALSE)
  invisible(TRUE)
}

# A complete submission as the form sends it by default.
good_fields <- function() {
  list(
    user_id = "tiago", email = "tiago@example.org", type_input = "Server",
    server_name = "storage1", input_folder = "RUN15", primer = "B2",
    sequencing_kit = "miseq_nano_2x250", keep_intermediate = "FALSE",
    archive_scope = "light", run_group = "", trim_primers = "TRUE",
    run_local = FALSE, run_igblast = TRUE, run_clonality = TRUE,
    archive_format = "zip", igblast_species = "mouse", igblast_panel = "ig_all",
    receptor_check = "auto", trim_sequence = "GGGAATTCGAGGTGCAGCTGCAGGAGTCTGG",
    n_read = 3, method = "tiago_vsearch", luka_light = "standard",
    mibs_junction_max_distance = 4, mibs_gene_grouping = "strict",
    mibs_clone_scope = "sample", mibs_min_reads_per_well = 6,
    mibs_min_wells_per_clone = 2, mibs_productive_only = FALSE,
    mibs_min_reads_per_cluster = 2
  )
}

# `[<-` with list() keeps NULL, a vector, or a list as the value of the field.
with_field <- function(fields, id, value) {
  fields[id] <- list(value)
  fields
}

check_field <- function(id, value, fields = good_fields()) {
  validate_submission(with_field(fields, id, value), choices)
}

names_the_field <- function(result, id) {
  any(grepl(SUBMISSION_FIELD_LABELS[[id]], result$errors, fixed = TRUE))
}

# Fields that carry text or a pick from a list: the ones a hostile string aims at.
TEXT_AND_CHOICE_FIELDS <- c(
  "user_id", "email", "type_input", "server_name", "input_folder", "primer",
  "sequencing_kit", "keep_intermediate", "archive_scope", "run_group",
  "trim_primers", "archive_format", "igblast_species", "igblast_panel",
  "receptor_check", "trim_sequence", "method", "luka_light",
  "mibs_gene_grouping", "mibs_clone_scope"
)

test_good_submission_is_accepted <- function() {
  result <- validate_submission(good_fields(), choices)
  check(isTRUE(result$ok), "a complete default submission is accepted")
  check(length(result$errors) == 0, "an accepted submission has no errors")
  check(identical(result$clean$method, "tiago_vsearch"), "clean carries the method")
}

test_each_field_accepts_good_and_refuses_bad <- function() {
  segment_fields <- c("user_id", "input_folder", "run_group", "primer")
  for (id in segment_fields) {
    for (value in cases$path_segment$good) {
      result <- check_field(id, value)
      check(result$ok, sprintf("%s accepts %s", id, encodeString(value)))
    }
    for (value in cases$path_segment$bad) {
      # run_group is the one segment that may be empty.
      if (id == "run_group" && value == "") next
      result <- check_field(id, value)
      check(!result$ok, sprintf("%s refuses %s", id, encodeString(value)))
      check(names_the_field(result, id),
            sprintf("the refusal of %s for %s names the field", encodeString(value), id))
    }
  }
  # The form trims user, folder and subfolder before it uses them, so stray
  # whitespace is cleaned, not refused. A barcode set is a pick and is not trimmed.
  for (id in c("user_id", "input_folder", "run_group")) {
    for (value in cases$path_segment$bad_untrimmed) {
      result <- check_field(id, value)
      check(result$ok, sprintf("%s cleans %s", id, encodeString(value)))
      check(identical(result$clean[[id]], "a"), sprintf("%s cleans to a", id))
    }
  }
  for (value in cases$path_segment$bad_untrimmed) {
    result <- check_field("primer", value)
    check(!result$ok, sprintf("primer refuses %s (not trimmed)", encodeString(value)))
    check(names_the_field(result, "primer"), "the primer refusal names the field")
  }

  for (value in cases$trim_sequence$good) {
    check(check_field("trim_sequence", value)$ok, paste("trim sequence accepts", value))
  }
  for (value in cases$trim_sequence$bad) {
    result <- check_field("trim_sequence", value)
    check(!result$ok, paste("trim sequence refuses", value))
    check(names_the_field(result, "trim_sequence"), "the trim refusal names the field")
  }

  for (value in cases$method$good) {
    check(check_field("method", value)$ok, paste("method accepts", value))
  }
  for (value in cases$method$bad) {
    result <- check_field("method", value)
    check(!result$ok, sprintf("method refuses %s", encodeString(value)))
    check(names_the_field(result, "method"), "the method refusal names the field")
  }
}

test_hostile_strings_refused_in_every_field <- function() {
  for (id in TEXT_AND_CHOICE_FIELDS) {
    for (value in cases$hostile) {
      result <- check_field(id, value)
      check(!result$ok, sprintf("%s refuses %s", id, encodeString(value)))
      check(names_the_field(result, id),
            sprintf("the refusal of %s for %s names the field", encodeString(value), id))
    }
  }
  # An empty trim sequence is only allowed with the trim off, and then a hostile
  # one is still refused, because it goes into sample.json either way.
  off <- with_field(good_fields(), "trim_primers", "FALSE")
  for (value in cases$hostile) {
    check(!check_field("trim_sequence", value, off)$ok,
          sprintf("trim sequence refuses %s with the trim off", encodeString(value)))
  }
  # A hostile value in a field must also not be echoed back into the message.
  result <- check_field("method", "$(touch canary)")
  check(!any(grepl("touch", result$errors, fixed = TRUE)), "refusals do not echo the value")
}

test_forged_shapes_refused <- function() {
  # Shiny lets a browser send any JSON for any input id.
  check(!check_field("primer", c("B2", "; touch x"))$ok, "a two-element primer is refused")
  check(!check_field("primer", list("B2"))$ok, "a list primer is refused")
  check(!check_field("primer", 42)$ok, "a numeric primer is refused")
  check(!check_field("primer", NULL)$ok, "a missing primer is refused")
  check(!check_field("primer", NA_character_)$ok, "an NA primer is refused")
  check(!check_field("primer", character(0))$ok, "an empty primer vector is refused")

  string_fields <- setdiff(TEXT_AND_CHOICE_FIELDS, c("run_group", "email"))
  for (id in union(string_fields, c("run_group", "email"))) {
    good <- good_fields()[[id]]
    for (forged in list(c(good, "; touch x"), list(good), 42, TRUE, NA, character(0))) {
      result <- check_field(id, forged)
      check(!result$ok, sprintf("%s refuses a forged %s", id, class(forged)[1]))
      check(names_the_field(result, id), sprintf("the forged %s refusal names %s", class(forged)[1], id))
    }
  }

  # Missing is refused for what the form always sends.
  required <- c("user_id", "type_input", "server_name", "input_folder", "primer",
                "keep_intermediate", "archive_scope", "trim_primers", "archive_format",
                "igblast_species", "igblast_panel", "trim_sequence", "n_read",
                "method", "luka_light", "run_igblast", "run_clonality")
  for (id in required) {
    result <- check_field(id, NULL)
    check(!result$ok, sprintf("%s refuses to be missing", id))
    check(names_the_field(result, id), sprintf("the missing-%s refusal names the field", id))
  }
  # Truth values must be real TRUE or FALSE, not a string a shell could read.
  for (id in c("run_igblast", "run_clonality", "run_local", "mibs_productive_only")) {
    for (forged in list("; touch x", "TRUE", 1, c(TRUE, FALSE), NA, list(TRUE))) {
      check(!check_field(id, forged)$ok, sprintf("%s refuses %s", id, deparse(forged)))
    }
    check(check_field(id, TRUE)$ok && check_field(id, FALSE)$ok, paste(id, "accepts TRUE and FALSE"))
  }
  check(!validate_submission("not a list", choices)$ok, "a non-list submission is refused")
  check(!validate_submission(list(), choices)$ok, "an empty submission is refused")
}

test_non_ascii_refused <- function() {
  for (id in c("user_id", "input_folder")) {
    for (value in c("café", "１２", "tiago​")) {
      result <- check_field(id, value)
      check(!result$ok, sprintf("%s refuses %s", id, encodeString(value)))
      check(names_the_field(result, id), "the non-ASCII refusal names the field")
    }
  }
  # toupper() maps U+017F (long s) to an ASCII S: it must be refused before cleaning
  # makes it look like a valid IUPAC letter.
  check(!check_field("trim_sequence", "ACGTACGTACſ")$ok, "a long s is not turned into an S")
  check(!check_field("email", "café@example.org")$ok, "a non-ASCII email is refused")
  # Bytes that are not text at all.
  check(!check_field("user_id", "\xff\xfe")$ok, "invalid UTF-8 is refused")
}

# The shared email cases (tests/validation_cases.json). A local part may not start with "-":
# an address that did would be an option if it were ever put on a mail program's command line.
test_the_shared_email_cases <- function() {
  for (value in cases$email$good) {
    check(isTRUE(is_valid_email(value)), sprintf("is_valid_email accepts %s", encodeString(value)))
    check(check_field("email", value)$ok, sprintf("the form accepts the email %s", encodeString(value)))
  }
  for (value in cases$email$bad) {
    check(!isTRUE(is_valid_email(value)), sprintf("is_valid_email refuses %s", encodeString(value)))
    result <- check_field("email", value)
    check(!result$ok, sprintf("the form refuses the email %s", encodeString(value)))
    check(names_the_field(result, "email"), "the email refusal names the field")
  }
}

# The form check trims an email before it tests it, so "a@b.org " passes the form. The row that
# is written to sample.json and the snapshot that feeds the run history and the emails must hold
# the trimmed address too: the cluster refuses a padded one, and that kills the run after the raw
# pull. submission_row() and submission_snapshot() live inside server() and need a Shiny session,
# so this test takes the right-hand side of each line of server.R that reads input$email, as
# written, and evaluates it against a fake input.
test_a_padded_email_reaches_the_row_trimmed <- function() {
  check(isTRUE(is_valid_email("a@b.org ")), "the form check accepts a padded address")
  result <- check_field("email", " a@b.org\t")
  check(result$ok, "the form accepts a padded address")
  check(identical(result$clean$email, "a@b.org"), "the form's own cleaned value is trimmed")

  lines <- readLines("server.R", encoding = "UTF-8")
  assigned <- grep("^\\s*[A-Za-z_]+\\s*=\\s*.*input\\$email", lines, value = TRUE)
  check(length(assigned) == 2,
        sprintf("server.R assigns input$email in two places (the row and the snapshot), found %d",
                length(assigned)))
  for (line in assigned) {
    rhs <- sub(",\\s*$", "", sub("^\\s*[A-Za-z_]+\\s*=\\s*", "", line))
    value_for <- function(email) {
      input <- if (is.null(email)) list() else list(email = email)
      eval(parse(text = rhs), list(input = input), environment())
    }
    check(identical(value_for("a@b.org "), "a@b.org"), sprintf("`%s` trims trailing spaces", rhs))
    check(identical(value_for(" \ta@b.org\n"), "a@b.org"), sprintf("`%s` trims all padding", rhs))
    check(identical(value_for("a@b.org"), "a@b.org"), sprintf("`%s` keeps a clean address", rhs))
    check(identical(value_for(""), ""), sprintf("`%s` keeps an empty address empty", rhs))
    check(identical(value_for(NULL), ""), sprintf("`%s` turns a missing address into empty", rhs))
  }
}

test_user_id_and_trim_are_cleaned <- function() {
  result <- check_field("user_id", " tiago\n")
  check(result$ok && identical(result$clean$user_id, "tiago"), "user id is trimmed")
  result <- check_field("trim_sequence", " acgtacgtacgt\n")
  check(result$ok && identical(result$clean$trim_sequence, "ACGTACGTACGT"),
        "trim sequence is trimmed and upper-cased")
  check(identical(clean_trim_sequence(" acgtacgtacgt\n"), "ACGTACGTACGT"),
        "clean_trim_sequence trims and upper-cases")
  # Case matters for a user: Tiago and tiago are different people.
  check(identical(check_field("user_id", "Tiago")$clean$user_id, "Tiago"), "user id keeps its case")
  result <- check_field("input_folder", "  RUN15\n")
  check(result$ok && identical(result$clean$input_folder, "RUN15"), "folder is trimmed")
  # A leading newline must not survive to a shell as a command separator.
  check(!check_field("user_id", "\ntouch canary")$ok, "a leading newline then a command is refused")
}

test_empty_trim_only_when_trim_off <- function() {
  on <- with_field(good_fields(), "trim_primers", "TRUE")
  off <- with_field(good_fields(), "trim_primers", "FALSE")
  for (empty in list("", "   ", "\n", NULL)) {
    result <- check_field("trim_sequence", empty, on)
    check(!result$ok, "an empty trim sequence is refused with the trim on")
    check(names_the_field(result, "trim_sequence"), "the empty-trim refusal names the field")
  }
  for (empty in list("", "   ", "\n", NULL)) {
    result <- check_field("trim_sequence", empty, off)
    check(result$ok, "an empty trim sequence is accepted with the trim off")
    check(identical(result$clean$trim_sequence, ""), "and it is cleaned to an empty string")
  }
  # With the trim off a sequence may still be there; it is checked like any other.
  check(check_field("trim_sequence", "ACGTACGTACGT", off)$ok, "a good sequence with the trim off")
  check(!check_field("trim_sequence", "ACGT", off)$ok, "a short sequence with the trim off")
}

test_numbers_in_range <- function() {
  for (bad in list(0, 501, "6;x", "6", 1.5, NA, NULL, Inf, -1, c(1, 2), TRUE)) {
    result <- check_field("n_read", bad)
    check(!result$ok, sprintf("n_read refuses %s", deparse(bad)))
    check(names_the_field(result, "n_read"), "the n_read refusal names the field")
  }
  for (good in list(1, 500, 3L, 250)) {
    check(check_field("n_read", good)$ok, sprintf("n_read accepts %s", deparse(good)))
  }

  ranges <- list(
    mibs_min_reads_per_well = c(1, 10000), mibs_min_wells_per_clone = c(1, 96),
    mibs_junction_max_distance = c(0, 30), mibs_min_reads_per_cluster = c(1, 1000)
  )
  for (id in names(ranges)) {
    lo <- ranges[[id]][1]; hi <- ranges[[id]][2]
    check(check_field(id, lo)$ok && check_field(id, hi)$ok, paste(id, "accepts its limits"))
    for (bad in list(lo - 1, hi + 1, "6;x", 2.5, c(lo, hi))) {
      result <- check_field(id, bad)
      check(!result$ok, sprintf("%s refuses %s", id, deparse(bad)))
      check(names_the_field(result, id), sprintf("the %s refusal names the field", id))
    }
    # An emptied box matters only to the pipeline that reads it.
    mibs <- with_field(good_fields(), "method", "tiago_sulkseq")
    check(!check_field(id, NA, mibs)$ok, paste(id, "refuses an empty box for sulkseq"))
    check(check_field(id, NA, good_fields())$ok, paste(id, "tolerates an empty box elsewhere"))
  }
}

test_primer_accepts_cluster_only_set <- function() {
  # The nanopore set LongMS_IgA exists only on the cluster and has real runs.
  check(check_field("primer", "LongMS_IgA")$ok, "LongMS_IgA is accepted")
}

test_server_name_and_upload_mode <- function() {
  check(check_field("server_name", "storage2")$ok, "a configured server is accepted")
  check(!check_field("server_name", "upload")$ok, "upload is not a server in Server mode")
  check(!check_field("server_name", "evil.example.org")$ok, "an unknown server is refused")
  # In Upload mode the form itself sets the server to "upload" and the folder box
  # to a note saying it is not used.
  upload <- with_field(good_fields(), "type_input", "Upload")
  upload$server_name <- "upload"
  upload$input_folder <- UPLOAD_FOLDER_PLACEHOLDER
  result <- validate_submission(upload, choices)
  check(result$ok, "the form's own Upload-mode values are accepted")
  check(identical(result$clean$input_folder, ""), "the placeholder cleans to no folder")
  upload$input_folder <- ""
  check(validate_submission(upload, choices)$ok, "an empty folder is accepted for Upload")
  upload$input_folder <- "; touch canary"
  check(!validate_submission(upload, choices)$ok, "a hostile folder is refused for Upload too")
  # In Server mode an empty folder is not allowed.
  check(!check_field("input_folder", "")$ok, "an empty folder is refused for Server")
}

test_sequencing_kit_and_choices <- function() {
  check(check_field("sequencing_kit", "")$ok, "an empty sequencing run is accepted")
  check(check_field("sequencing_kit", NULL)$ok, "a missing sequencing run is accepted")
  check(check_field("sequencing_kit", "plasmidsaurus_custom")$ok, "a configured kit is accepted")
  check(!check_field("sequencing_kit", "novaseq")$ok, "an unknown kit is refused")
  check(!check_field("igblast_panel", "ig_light")$ok, "a panel the widget no longer offers is refused")
  for (panel in unname(igblast_panel_choices)) {
    check(check_field("igblast_panel", panel)$ok, paste("panel", panel))
  }
  check(identical(unname(METHOD_CHOICES), c("tiago_vsearch", "classic_luka", "tiago_nanopore", "tiago_sulkseq")),
        "the method choices are the four pipelines")
  for (name in c("method", "archive_scope", "archive_format", "igblast_species", "igblast_panel",
                 "luka_light", "receptor_check", "mibs_gene_grouping", "mibs_clone_scope",
                 "keep_intermediate", "trim_primers", "type_input", "server_name", "sequencing_kit")) {
    check(is.character(choices[[name]]) && is.null(names(choices[[name]])),
          paste("choices for", name, "are a plain character vector"))
  }
}

test_path_segment_helper_is_strict <- function() {
  check(is_valid_path_segment("a") && is_valid_path_segment(" a\n"), "good segments, trimmed")
  check(!is_valid_path_segment(c("a", "b")), "a vector is not a segment")
  check(!is_valid_path_segment(list("a")), "a list is not a segment")
  check(!is_valid_path_segment(NULL) && !is_valid_path_segment(NA_character_), "NULL and NA")
  check(!is_valid_path_segment("a\nb") && !is_valid_path_segment("-x") && !is_valid_path_segment(".."),
        "bad segments")
  # PCRE's $ also matches before a final newline; the helper must not be fooled.
  check(!is_path_segment_text("a\n"), "an untrimmed trailing newline is not a segment")
  check(!is_path_segment_text("\xff"), "invalid UTF-8 is not a segment")
}

test_run_script_passes_hostile_args_verbatim <- function() {
  dir <- file.path(tempfile("run script "), "with space")
  dir.create(dir, recursive = TRUE)
  canary <- file.path(tempdir(), "canary")
  script <- file.path(dir, "echo args.bash")
  writeLines(c("#!/bin/bash", "printf '<%s>\\n' \"$@\""), script)

  args <- c(cases$hostile, sprintf("; touch %s", canary), sprintf("$(touch %s)", canary),
            sprintf("`touch %s`", canary), "", "a b", "it's", "*", "~", "$HOME", "\\")
  result <- run_script(script, args = args)
  expected <- paste0("<", args, ">", collapse = "\n")
  check(identical(result$status, 0), "the script ran")
  check(identical(result$output, expected),
        paste("the script received exactly the arguments; got:", encodeString(result$output)))
  check(!file.exists(canary), "no canary file was created")

  # With a time limit the arguments still arrive the same way.
  result <- run_script(script, args = c("-t", "$(touch x)"), timeout_sec = 20)
  check(identical(result$output, "<-t>\n<$(touch x)>"), "arguments are intact under a timeout")

  # A script that is not there reports a failure instead of raising.
  result <- run_script(file.path(dir, "missing.bash"), args = "x")
  check(result$status != 0, "a missing script is a failure, not an error")
}

# The app looks up run_script() and shinyalert() when a submission arrives, so
# the server tests replace both for the length of one testServer() call. They are
# put where global.R puts its own functions, and put back afterwards, because the
# other tests in this file use the real run_script().
install_stubs <- function(calls) {
  assign("run_script", function(script_path, args = character(), timeout_sec = 0) {
    calls$scripts[[length(calls$scripts) + 1]] <- list(script = script_path, args = args)
    list(status = 0, output = "OK|stub|stubbed|")
  }, envir = globalenv())
  assign("shinyalert", function(...) {
    calls$alerts[[length(calls$alerts) + 1]] <- list(...)
    invisible(NULL)
  }, envir = globalenv())
}

remove_stubs <- function() {
  if (exists("shinyalert", envir = globalenv(), inherits = FALSE)) {
    rm("shinyalert", envir = globalenv())
  }
  source("run_script.R", local = FALSE)
}

test_nothing_reaches_a_script_on_bad_input <- function() {
  canary <- file.path(tempdir(), "server_canary")
  calls <- new.env()
  calls$scripts <- list()
  calls$alerts <- list()
  on.exit(remove_stubs())

  shiny::testServer(shiny::shinyAppDir("."), {
    # After the app has loaded: loading global.R puts the real functions back.
    install_stubs(calls)

    base <- good_fields()
    base$run <- NULL
    run_number <- 0
    submit <- function(...) {
      inputs <- utils::modifyList(base, list(...), keep.null = TRUE)
      run_number <<- run_number + 1
      inputs$run <- run_number
      do.call(session$setInputs, inputs)
    }

    hostile <- c(cases$hostile, sprintf("$(touch %s)", canary), sprintf("; touch %s", canary))
    for (id in c("trim_sequence", "method", "primer", "user_id", "input_folder", "run_group",
                 "igblast_panel", "igblast_species", "luka_light", "archive_scope",
                 "receptor_check", "server_name")) {
      for (value in hostile) {
        before <- length(calls$scripts)
        alerts_before <- length(calls$alerts)
        args <- list(value)
        names(args) <- id
        do.call(submit, args)
        if (length(calls$scripts) != before) {
          stop("FAILED: a script was run for ", id, " = ", encodeString(value), call. = FALSE)
        }
        if (length(calls$alerts) != alerts_before + 1) {
          stop("FAILED: no alert for ", id, " = ", encodeString(value), call. = FALSE)
        }
        alert <- calls$alerts[[length(calls$alerts)]]
        if (!identical(alert$title, "Check the form") ||
            !grepl(SUBMISSION_FIELD_LABELS[[id]], alert$text, fixed = TRUE)) {
          stop("FAILED: the alert for ", id, " does not name it: ", alert$text, call. = FALSE)
        }
      }
    }
    # The same fields with the checks on the inputs the brief names.
    submit(trim_sequence = "$(touch canary)")
    submit(trim_sequence = "ACGT", trim_primers = "TRUE")
    submit(user_id = "tiago\n; touch canary")
    submit(method = "classic_luka; touch canary")
    if (length(calls$scripts) != 0) stop("FAILED: a script ran on bad input", call. = FALSE)
  })
  check(!file.exists(canary) && !file.exists(file.path(getwd(), "canary")), "no canary file exists")
  check(length(calls$scripts) == 0, "zero scripts were run on bad input")
}

test_a_good_submission_reaches_preflight_with_clean_values <- function() {
  calls <- new.env()
  calls$scripts <- list()
  calls$alerts <- list()
  on.exit(remove_stubs())

  shiny::testServer(shiny::shinyAppDir("."), {
    install_stubs(calls)

    inputs <- good_fields()
    inputs$user_id <- " tiago\n"
    inputs$trim_sequence <- " gggaattcgaggtgcagctgcaggagtctgg\n"
    inputs$run <- 1
    try(do.call(session$setInputs, inputs), silent = TRUE)
  })
  preflight <- Filter(function(call) grepl("preflight.bash$", call$script), calls$scripts)
  check(length(preflight) >= 1, "a good submission reaches preflight")
  args <- preflight[[1]]$args
  check(identical(args[which(args == "-u") + 1], "tiago"), "preflight gets the trimmed user name")
  check(identical(args[which(args == "-t") + 1], "GGGAATTCGAGGTGCAGCTGCAGGAGTCTGG"),
        "preflight gets the cleaned trim sequence")
  check(!any(vapply(calls$alerts, function(a) identical(a$title, "Check the form"), logical(1))),
        "a good submission is not stopped by the form check")
}

test_folder_dropdown_cannot_leave_the_barcode_folder <- function() {
  bc <- Sys.getenv("TBRC_BC_DIR")
  dir.create(file.path(bc, "real_set"), recursive = TRUE)
  writeLines("A01\tACGT", file.path(bc, "real_set", "BC1.txt"))
  writeLines("A01\tACGT", file.path(bc, "real_set", "BC2.txt"))
  outside <- file.path(dirname(bc), "outside_set")
  dir.create(outside, recursive = TRUE)
  writeLines("SECRET", file.path(outside, "BC1.txt"))
  writeLines("SECRET", file.path(outside, "BC2.txt"))

  results <- new.env()
  shiny::testServer(shiny::shinyAppDir("."), {
    session$setInputs(type_input = "Server", folder_dropdown = "real_set")
    results$real <- barcode_content()
    for (escape in c("../outside_set", "..", ".", "real_set/../../outside_set", "/etc")) {
      session$setInputs(folder_dropdown = escape)
      results[[escape]] <- barcode_content()
    }
  })
  check(identical(results$real$BC1, "A01\tACGT"), "a real barcode folder is read")
  for (escape in c("../outside_set", "..", ".", "real_set/../../outside_set", "/etc")) {
    check(is.null(results[[escape]]$BC1) && is.null(results[[escape]]$BC2),
          paste("folder_dropdown refuses", escape))
  }
}

test_ui_uses_the_choice_constants <- function() {
  page <- suppressWarnings(suppressMessages(source("ui.R", local = new.env())$value))
  html <- as.character(htmltools::renderTags(page)$html)
  for (value in c(unname(METHOD_CHOICES), unname(ARCHIVE_SCOPE_CHOICES), unname(ARCHIVE_FORMAT_CHOICES),
                  unname(IGBLAST_SPECIES_CHOICES), unname(RECEPTOR_CHECK_CHOICES),
                  unname(LUKA_LIGHT_CHOICES), unname(MIBS_GENE_GROUPING_CHOICES),
                  unname(MIBS_CLONE_SCOPE_CHOICES))) {
    check(grepl(sprintf("value=\"%s\"", value), html, fixed = TRUE),
          paste("the page offers", value))
  }
}

# A state folder the app cannot write to must stop the start with a message that
# says what to do, not let every later write fail silently (a root-owned volume
# left over from when the app ran as root made the barcode list and the run ledger
# unwritable with no error anywhere).
test_unwritable_state_dir_stops_with_a_clear_message <- function() {
  stopped_with <- function(path) {
    Sys.setenv(TBRC_TEST_STATE_DIR = path)
    tryCatch({ app_state_dir("TBRC_TEST_STATE_DIR", "unused"); NA_character_ },
             error = function(e) conditionMessage(e))
  }
  # A folder that cannot be created, whoever runs this: its parent is a file.
  blocker <- file.path(scratch, "a_file")
  writeLines("x", blocker)
  msg <- stopped_with(file.path(blocker, "state"))
  check(!is.na(msg), "app_state_dir stops when the folder cannot be created")
  check(grepl("TBRC_TEST_STATE_DIR", msg, fixed = TRUE), "the message names the variable")
  check(grepl("not writable", msg, fixed = TRUE), "the message says the folder is not writable")
  check(grepl("chown", msg, fixed = TRUE), "the message says how to fix it")
  # It has to serve any site: say what is needed (write access, or another folder), and show the Docker
  # chown as an example with a placeholder id and image, not this lab's user id or a named image.
  check(grepl("write access", msg, fixed = TRUE), "the message asks for write access to the folder")
  check(grepl("point TBRC_TEST_STATE_DIR at a writable one", msg, fixed = TRUE),
        "the message offers a writable folder through the variable")
  check(grepl("chown -R <uid>:<uid>", msg, fixed = TRUE), "the example uses <uid>")
  check(!grepl("501", msg, fixed = TRUE), "the example names no particular user id")
  check(!grepl("caddy|the Docker host", msg, ignore.case = TRUE), "the example names no image or site")

  # A folder that exists but is read-only. Root can write anywhere (access(2)
  # always says yes), so this half only means something for an ordinary user,
  # which is what the container runs as.
  if (!identical(Sys.info()[["effective_user"]], "root")) {
    read_only <- file.path(scratch, "read_only")
    dir.create(read_only)
    Sys.chmod(read_only, "555")
    on.exit(Sys.chmod(read_only, "755"), add = TRUE)
    check(!is.na(stopped_with(read_only)), "app_state_dir stops for a read-only folder")
  }

  # And a writable one is accepted and returned as before.
  fine <- file.path(scratch, "fine")
  Sys.setenv(TBRC_TEST_STATE_DIR = fine)
  check(identical(app_state_dir("TBRC_TEST_STATE_DIR", "unused"), normalizePath(fine, winslash = "/")),
        "a writable folder is created and returned")
}

test_barcode_preview_reads_illumina_and_nanopore_sets <- function() {
  store <- tempfile("bc_store")
  dir.create(store)
  on.exit(unlink(store, recursive = TRUE), add = TRUE)
  make <- function(name, files) {
    folder <- file.path(store, name)
    dir.create(folder)
    for (f in files) writeLines("x", file.path(folder, f))
    folder
  }

  illumina <- barcode_preview_files(make("illumina", c("BC1.txt", "BC2.txt")))
  check(identical(basename(illumina$paths), c("BC1.txt", "BC2.txt")), "an Illumina set previews BC1 and BC2")
  check(identical(illumina$labels, c("Barcode 1:", "Barcode 2:")), "an Illumina set keeps its labels")

  nanopore <- barcode_preview_files(make("nanopore", c("plate_barcodes.tsv", "BCFR1_anchors.tsv")))
  check(identical(basename(nanopore$paths), c("plate_barcodes.tsv", "BCFR1_anchors.tsv")),
        "a nanopore set previews its plate barcodes and anchors")
  check(all(nzchar(nanopore$labels)) && !any(grepl("^Barcode [12]:", nanopore$labels)),
        "a nanopore set says what its files are")

  check(is.null(barcode_preview_files(make("half", "BC1.txt"))), "half a set has no preview")
  check(is.null(barcode_preview_files(make("empty", character(0)))), "an empty folder has no preview")
}

# --------------------------------------------------------- run completion --
# The completion handler (finish() in server.R) records the run, emails, and then
# updates the page. A staging run showed that when the history step threw, the page
# never got its last update: Workflow and Export stayed unfinished on a run that
# had succeeded. complete_run() (global.R) is what finish() hands those three steps
# to, so each can be tested without a Shiny session.

# Runs complete_run() and returns what happened, in order, plus what the page got.
play_complete_run <- function(history = function() invisible(), emails = NULL) {
  events <- character(0)
  page <- NULL
  logged <- character(0)
  if (is.null(emails)) {
    emails <- function() list(user = list(status = 0), admin = list(status = 0))
  }
  withCallingHandlers(
    complete_run(
      record_history = function() { events <<- c(events, "history"); history() },
      send_notifications = function() { events <<- c(events, "email"); emails() },
      update_page = function(email_result) { events <<- c(events, "page"); page <<- email_result }
    ),
    message = function(m) {
      logged <<- c(logged, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(events = events, page = page, logged = logged)
}

test_a_run_is_recorded_emailed_and_then_the_page_is_updated <- function() {
  out <- play_complete_run()
  check(identical(out$events, c("history", "email", "page")), "history, then email, then the page")
  check(identical(out$page$user$status, 0) && identical(out$page$admin$status, 0),
        "the page gets the email result as the email step returned it")
  check(length(out$logged) == 0, "nothing is logged when nothing fails")
}

test_a_failing_history_step_does_not_stop_the_page_update <- function() {
  out <- play_complete_run(history = function() stop("replacement has 1 row, data has 0"))
  check(identical(out$events, c("history", "email", "page")),
        "a history error is followed by the email step and the page update")
  check(identical(out$page$user$status, 0), "the email result is untouched by a history error")
  check(any(grepl("history", out$logged, ignore.case = TRUE) &
              grepl("replacement has 1 row, data has 0", out$logged, fixed = TRUE)),
        "the history error is logged with message(), with its reason")
}

test_a_failing_email_step_does_not_stop_the_page_update <- function() {
  out <- play_complete_run(emails = function() stop("sendmail exploded"))
  check(identical(out$events, c("history", "email", "page")),
        "an email error is followed by the page update")
  check(is.list(out$page) && !is.null(out$page$user$status) && !is.null(out$page$admin$status),
        "the page still gets email_result$user$status and email_result$admin$status")
  check(out$page$user$status != 0 && out$page$admin$status != 0,
        "both email statuses are non-zero after an email error")
  check(any(grepl("email", out$logged, ignore.case = TRUE) &
              grepl("sendmail exploded", out$logged, fixed = TRUE)),
        "the email error is logged with message(), with its reason")
}

test_both_steps_failing_still_updates_the_page <- function() {
  out <- play_complete_run(history = function() stop("h"), emails = function() stop("e"))
  check(identical(out$events, c("history", "email", "page")), "the page update still happens")
  check(length(out$logged) == 2, "each failure is logged once")
  check(out$page$user$status != 0, "the page is told the email failed")
}

test_an_email_result_of_the_wrong_shape_is_treated_as_failed <- function() {
  for (odd in list(NULL, "sent", list(user = list(status = 0)))) {
    out <- play_complete_run(emails = function() odd)
    check(identical(out$events, c("history", "email", "page")), "the page update still happens")
    check(out$page$user$status != 0 && out$page$admin$status != 0,
          "a result the page cannot read is reported as a failed email")
  }
}

# finish() must hand the history and email steps to complete_run() as closures; a
# direct call there would put them back in front of the page update.
test_finish_in_server_r_runs_its_steps_through_complete_run <- function() {
  exprs <- parse(file = "server.R", keep.source = FALSE)
  find_function <- function(e, name) {
    if (is.call(e)) {
      if ((identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("="))) &&
          identical(e[[2]], as.name(name)) && is.call(e[[3]]) &&
          identical(e[[3]][[1]], as.name("function"))) return(e[[3]])
      for (i in seq_along(e)) {
        found <- find_function(e[[i]], name)
        if (!is.null(found)) return(found)
      }
    }
    NULL
  }
  finish <- NULL
  for (e in exprs) { finish <- find_function(e, "finish"); if (!is.null(finish)) break }
  check(!is.null(finish), "server.R defines finish()")

  # Calls to the two steps that sit directly in finish()'s body, not inside a
  # function(...) literal.
  direct <- character(0)
  walk <- function(e, inside_closure) {
    if (!is.call(e)) return(invisible())
    head <- e[[1]]
    if (is.name(head) && as.character(head) %in% c("append_run_history", "send_run_notifications") &&
        !inside_closure) direct <<- c(direct, as.character(head))
    now_inside <- inside_closure || identical(head, as.name("function"))
    for (i in seq_along(e)) walk(e[[i]], now_inside)
  }
  walk(finish[[3]], FALSE)
  check(length(direct) == 0,
        sprintf("finish() calls %s outside a closure handed to complete_run()", paste(direct, collapse = ", ")))
  uses <- FALSE
  find_call <- function(e) {
    if (!is.call(e)) return(invisible())
    if (identical(e[[1]], as.name("complete_run"))) uses <<- TRUE
    for (i in seq_along(e)) find_call(e[[i]])
  }
  find_call(finish[[3]])
  check(uses, "finish() calls complete_run()")
}

# ----------------------------------------------------- stage progress motion --
# A stage that has started and not finished moves (CSS in ui.R, keyed on the inline
# width spsComps writes: 0% is not started, 100% is done). What a browser does with
# the rule was checked in headless Chrome; these tests keep the rule's shape.
test_stage_motion_is_a_css_rule_that_respects_reduced_motion <- function() {
  ui_source <- paste(readLines("ui.R", warn = FALSE), collapse = "\n")
  check(grepl("@keyframes", ui_source, fixed = TRUE), "ui.R defines the keyframes")

  # The text of the @media block that guards the motion, by brace matching.
  start <- regexpr("@media (prefers-reduced-motion: no-preference)", ui_source, fixed = TRUE)
  check(start > 0, "the motion sits under prefers-reduced-motion: no-preference")
  chars <- strsplit(substring(ui_source, start), "")[[1]]
  depth <- 0
  end <- length(chars)
  for (i in seq_along(chars)) {
    if (chars[i] == "{") depth <- depth + 1
    if (chars[i] == "}") { depth <- depth - 1; if (depth == 0) { end <- i; break } }
  }
  guarded <- paste(chars[seq_len(end)], collapse = "")
  check(grepl("animation:", guarded, fixed = TRUE), "the animation is declared inside that block")
  check(!grepl("animation:", sub(guarded, "", ui_source, fixed = TRUE), fixed = TRUE),
        "no animation is declared outside it")
  check(grepl(".progress-bar.stage-running", guarded, fixed = TRUE), "the stripes move only on a bar marked as running")
  # Safari showed the stripes stuck: rules that read the inline width did not follow a script's change to it.
  check(!grepl("[style*=", ui_source, fixed = TRUE), "no rule reads the inline width")
  check(!grepl(":has(", ui_source, fixed = TRUE), "no rule depends on :has()")
  check(grepl("w > 0 && w < 100", ui_source, fixed = TRUE), "a stage counts as running only between 0% and 100%")
  check(grepl("MutationObserver", ui_source, fixed = TRUE), "the mark follows every change to a bar's width")
  check(grepl("[id$='-pg']", ui_source, fixed = TRUE), "the script watches the five stage bars, not the overall bar")
}

# After a failed run the Workflow bar was left at 20%. With motion on 20%, that
# would keep the stripes going under a "Run failed" message.
# Workflow waits at 20%, a strip 45 pixels wide in a panel at 70% opacity, so in a real session the
# moving stripes could not be seen. A running stage also sweeps its whole track and shows a spinner.
test_a_running_stage_shows_across_its_whole_track <- function() {
  ui_source <- paste(readLines("ui.R", warn = FALSE), collapse = "\n")
  check(grepl("@keyframes stage-sweep", ui_source, fixed = TRUE), "ui.R defines the sweep")
  check(grepl("@keyframes stage-spin", ui_source, fixed = TRUE), "ui.R defines the spinner's turn")
  motion <- substring(ui_source, regexpr("@media (prefers-reduced-motion: no-preference)", ui_source, fixed = TRUE)[1])
  check(grepl("::after", motion, fixed = TRUE) && grepl("animation: stage-sweep", motion, fixed = TRUE),
        "the sweep moves only for people who have not asked for less motion")
  check(grepl(".progress.stage-running::after", ui_source, fixed = TRUE), "the sweep is drawn on the track of a running stage")
  check(grepl("content: '\\\\f110'", ui_source, fixed = TRUE),
        "a running stage shows a spinner in place of the red cross")
}

test_a_failed_run_stops_the_workflow_bar <- function() {
  exprs <- parse(file = "server.R", keep.source = FALSE)
  branches <- list()
  walk <- function(e) {
    if (!is.call(e)) return(invisible())
    if (identical(e[[1]], as.name("if")) && length(e) >= 4 &&
        identical(deparse(e[[2]]), 'identical(label, "failed")')) branches[[length(branches) + 1]] <<- e
    for (i in seq_along(e)) walk(e[[i]])
  }
  for (e in exprs) walk(e)
  check(length(branches) == 1, "server.R has one place that tells a failed run from a finished one")
  updates <- function(branch) {
    found <- character(0)
    collect <- function(e) {
      if (!is.call(e)) return(invisible())
      if (identical(e[[1]], as.name("pgPaneUpdate")))
        found <<- c(found, paste(vapply(as.list(e)[-1], deparse, ""), collapse = ", "))
      for (i in seq_along(e)) collect(e[[i]])
    }
    collect(branch)
    found
  }
  check('"thispg", "workflow", 0' %in% updates(branches[[1]][[3]]), "a failed run sets Workflow back to 0")
  check('"thispg", "workflow", 100' %in% updates(branches[[1]][[4]]), "a finished run sets Workflow to 100")
}

# ------------------------------------------------------------- run ledger --
# append_run_stat() appends one row to run_stats.tsv, or rewrites the file when the
# row has columns the header lacks. It reads the global RUN_STATS_FILE, so each test
# points that at a ledger in its own scratch folder and puts it back afterwards.
in_scratch_folder <- function(code) {
  folder <- tempfile("ledger")
  dir.create(folder)
  original <- RUN_STATS_FILE
  assign("RUN_STATS_FILE", file.path(folder, "run_stats.tsv"), envir = globalenv())
  on.exit({
    assign("RUN_STATS_FILE", original, envir = globalenv())
    unlink(folder, recursive = TRUE)
  }, add = TRUE)
  force(code)
}

read_ledger <- function() {
  utils::read.table(RUN_STATS_FILE, sep = "\t", header = TRUE, quote = "\"",
                    stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
}

test_a_header_only_ledger_with_other_columns_takes_the_new_row <- function() {
  in_scratch_folder({
    writeLines('"old_a"\t"old_b"', RUN_STATS_FILE)
    row <- data.frame(timestamp = "2026-11-09 10:00:00", run_name = "r1", status = "completed",
                      stringsAsFactors = FALSE)
    append_run_stat(row)
    got <- read_ledger()
    check(nrow(got) == 1, "the ledger holds exactly the one new row")
    check(identical(names(got), names(row)), "the ledger has the new row's columns")
    check(identical(got$run_name, "r1") && identical(got$status, "completed"), "the row's values survive")
  })
}

test_a_header_only_ledger_with_the_same_columns_takes_the_new_row <- function() {
  in_scratch_folder({
    row <- data.frame(timestamp = "t", run_name = "r1", stringsAsFactors = FALSE)
    writeLines('"timestamp"\t"run_name"', RUN_STATS_FILE)
    append_run_stat(row)
    got <- read_ledger()
    check(nrow(got) == 1 && identical(got$run_name, "r1"), "one row, no repeated header")
  })
}

test_a_ledger_with_rows_still_appends_and_still_grows_its_columns <- function() {
  in_scratch_folder({
    first <- data.frame(timestamp = "t1", run_name = "r1", stringsAsFactors = FALSE)
    append_run_stat(first)                                   # no file yet: created
    append_run_stat(data.frame(timestamp = "t2", run_name = "r2", stringsAsFactors = FALSE))
    got <- read_ledger()
    check(nrow(got) == 2 && identical(got$run_name, c("r1", "r2")), "same columns: the row is appended")

    append_run_stat(data.frame(timestamp = "t3", run_name = "r3", extra = "x", stringsAsFactors = FALSE))
    got <- read_ledger()
    check(identical(names(got), c("timestamp", "run_name", "extra")), "a new column is added to the header")
    check(nrow(got) == 3 && identical(got$run_name, c("r1", "r2", "r3")), "old rows are kept, in order")
    check(identical(got$extra, c("", "", "x")), "old rows get an empty value in the new column")

    append_run_stat(data.frame(timestamp = "t4", stringsAsFactors = FALSE))
    got <- read_ledger()
    check(nrow(got) == 4 && identical(got$run_name, c("r1", "r2", "r3", "")),
          "a row lacking a column gets an empty value there")
  })
}

tests <- Filter(function(name) startsWith(name, "test_") && is.function(get(name)), ls())
failures <- character(0)
for (name in tests) {
  outcome <- tryCatch({ get(name)(); NULL }, error = function(e) conditionMessage(e))
  if (is.null(outcome)) {
    cat("PASS", name, "\n")
  } else {
    cat("FAIL", name, "-", outcome, "\n")
    failures <- c(failures, name)
  }
}
cat(sprintf("\n%d of %d tests passed\n", length(tests) - length(failures), length(tests)))
if (length(failures)) quit(status = 1)
