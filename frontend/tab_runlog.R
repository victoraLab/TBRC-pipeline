# ---------------------------------------------------------------------------
# tab_runlog.R -- self-contained Shiny module rendering docs/run-log.tsv
#
# The run log is the authoritative list of every TBRC run we can still evidence:
# the runs that carry pipeline metadata (sample.json / run_stats.tsv) plus the
# older ones reconstructed from result archives and raw-data folders on
# storage1 and storage2. It is a static file rebuilt by the storage audit,
# not app state, so this module only ever reads it.
#
# Wiring is the caller's job -- see the note at the bottom of this file.
# ---------------------------------------------------------------------------

RUN_LOG_COLUMNS <- c(
  "run_index", "run_name", "date", "user", "chain", "primer_set",
  "storage_server", "era", "status", "input_bytes", "runtime_sec",
  "evidence", "path", "notes"
)

RUN_LOG_DERIVED <- c("run_index_num", "date_parsed", "input_bytes_num", "runtime_sec_num")

# Labels standing in for an empty cell, so a filter can still select those rows.
RUN_LOG_BLANK <- list(
  user = "(not established)",
  chain = "unknown",
  primer_set = "(not recorded)",
  era = "(unknown)"
)

# The log lives in the repo beside the app rather than in the writable state
# directory, but a container may mount it elsewhere -- allow an override.
run_log_file <- function() {
  configured <- Sys.getenv("TBRC_RUN_LOG_FILE", unset = "")
  if (nzchar(configured)) {
    return(configured)
  }

  candidates <- c(
    if (exists("app_path", mode = "function")) app_path("..", "docs", "run-log.tsv"),
    if (exists("STATE_DIR")) file.path(STATE_DIR, "run-log.tsv"),
    file.path("..", "docs", "run-log.tsv"),
    file.path("docs", "run-log.tsv")
  )

  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = FALSE))
    }
  }
  if (length(candidates)) candidates[[1]] else "run-log.tsv"
}

run_log_empty <- function() {
  columns <- c(RUN_LOG_COLUMNS, RUN_LOG_DERIVED)
  frame <- as.data.frame(
    setNames(replicate(length(columns), character(0), simplify = FALSE), columns),
    stringsAsFactors = FALSE
  )
  frame$run_index_num <- integer(0)
  frame$date_parsed <- as.Date(character(0))
  frame$input_bytes_num <- numeric(0)
  frame$runtime_sec_num <- numeric(0)
  frame
}

# Runs submitted through the app are recorded in the writable state directory,
# not in docs/run-log.tsv -- which is a historical audit, is committed to the
# repository, and is mounted read-only in the container. The tab read only the
# audit, so nothing anyone submitted ever appeared in it. Read both.
runlog_live_file <- function() {
  configured <- Sys.getenv("TBRC_RUN_STATS_FILE", unset = "")
  if (nzchar(configured)) return(configured)
  if (exists("RUN_STATS_FILE")) return(get("RUN_STATS_FILE"))
  if (exists("STATE_DIR")) return(file.path(get("STATE_DIR"), "run_stats.tsv"))
  "run_stats.tsv"
}

# The two files were written for different purposes and do not share a schema.
# Map the live one onto the audit's columns rather than the other way round: the
# audit is the format the tab, its filters and its summary tiles already speak.
runlog_read_live <- function(path = runlog_live_file()) {
  if (!file.exists(path)) return(NULL)
  live <- tryCatch(
    utils::read.table(path, sep = "\t", header = TRUE, quote = "\"",
                      comment.char = "", stringsAsFactors = FALSE,
                      check.names = FALSE, colClasses = "character",
                      fill = TRUE, na.strings = character(0)),
    error = function(...) NULL
  )
  if (is.null(live) || !nrow(live)) return(NULL)

  pick <- function(column, default = "") {
    if (is.null(live[[column]])) rep(default, nrow(live))
    else {
      v <- as.character(live[[column]]); v[is.na(v)] <- default; trimws(v)
    }
  }
  stamp <- pick("timestamp")
  out <- data.frame(
    run_index      = rep("", nrow(live)),
    run_name       = pick("run_name"),
    date           = substr(stamp, 1, 10),
    user           = pick("user_id"),
    chain          = pick("igblast_panel"),
    primer_set     = pick("primer_set"),
    storage_server = pick("server_name"),
    era            = rep("submitted through the app", nrow(live)),
    status         = pick("status"),
    input_bytes    = pick("input_bytes"),
    runtime_sec    = pick("runtime_sec"),
    evidence       = pick("pipeline"),
    path           = rep("", nrow(live)),
    notes          = pick("backend_output"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  # Long backend output would swamp a table cell; the full log is written beside
  # the run metadata by write_backend_log().
  out$notes <- substr(gsub("[\r\n\t]+", " ", out$notes), 1, 300)
  out
}

read_run_log <- function(path = run_log_file()) {
  if (!file.exists(path)) {
    return(runlog_finalise(runlog_bind_live(run_log_empty())))
  }

  # quote = "" because the notes column carries apostrophes and quotation marks
  # that would otherwise swallow the rest of the file; comment.char = "" because
  # a path may contain a hash.
  log <- tryCatch(
    utils::read.table(
      path, sep = "\t", header = TRUE, quote = "", comment.char = "",
      stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character",
      na.strings = character(0)
    ),
    error = function(...) NULL
  )
  if (is.null(log) || nrow(log) == 0) {
    return(runlog_finalise(runlog_bind_live(run_log_empty())))
  }

  for (column in RUN_LOG_COLUMNS) {
    values <- if (is.null(log[[column]])) "" else as.character(log[[column]])
    values[is.na(values)] <- ""
    log[[column]] <- trimws(values)
  }

  log$run_index_num <- suppressWarnings(as.integer(log$run_index))
  log$date_parsed <- suppressWarnings(as.Date(log$date, format = "%Y-%m-%d"))
  log$input_bytes_num <- suppressWarnings(as.numeric(log$input_bytes))
  log$runtime_sec_num <- suppressWarnings(as.numeric(log$runtime_sec))

  runlog_finalise(runlog_bind_live(log))
}

runlog_finalise <- function(log) {
  if (is.null(log) || !nrow(log)) return(run_log_empty())
  log$run_index_num <- suppressWarnings(as.integer(log$run_index))
  log$date_parsed <- suppressWarnings(as.Date(log$date, format = "%Y-%m-%d"))
  log$input_bytes_num <- suppressWarnings(as.numeric(log$input_bytes))
  log$runtime_sec_num <- suppressWarnings(as.numeric(log$runtime_sec))
  log <- log[order(log$run_index_num, na.last = TRUE), , drop = FALSE]
  log[, c(RUN_LOG_COLUMNS, RUN_LOG_DERIVED), drop = FALSE]
}

# Append live runs after the audit, numbering them on from its last index so the
# tab's newest-first ordering keeps working. A run present in both -- the audit
# was rebuilt after the run happened -- is kept once, from the audit.
runlog_bind_live <- function(log) {
  live <- runlog_read_live()
  if (is.null(live) || !nrow(live)) return(log)

  for (column in RUN_LOG_COLUMNS) {
    values <- if (is.null(live[[column]])) "" else as.character(live[[column]])
    values[is.na(values)] <- ""
    live[[column]] <- trimws(values)
  }
  live <- live[, RUN_LOG_COLUMNS, drop = FALSE]

  seen <- paste(log$run_name, log$date, sep = "\r")
  live <- live[!(paste(live$run_name, live$date, sep = "\r") %in% seen), , drop = FALSE]
  if (!nrow(live)) return(log)

  last <- suppressWarnings(max(as.integer(log$run_index), na.rm = TRUE))
  if (!is.finite(last)) last <- 0L
  live$run_index <- as.character(seq.int(last + 1L, length.out = nrow(live)))

  rbind(log[, RUN_LOG_COLUMNS, drop = FALSE], live)
}

# --------------------------------------------------------------- formatting --

runlog_format_bytes <- function(bytes) {
  vapply(bytes, function(value) {
    if (!is.finite(value) || value <= 0) {
      return("")
    }
    units <- c("B", "kB", "MB", "GB", "TB")
    power <- min(floor(log(value, 1024)), length(units) - 1)
    scaled <- value / (1024^power)
    sprintf(if (power == 0 || scaled >= 100) "%.0f %s" else "%.1f %s", scaled, units[[power + 1]])
  }, character(1), USE.NAMES = FALSE)
}

runlog_format_runtime <- function(seconds) {
  vapply(seconds, function(value) {
    if (!is.finite(value) || value < 0) {
      return("")
    }
    if (value < 60) {
      return(sprintf("%.0f s", value))
    }
    sprintf("%d m %02d s", as.integer(value %/% 60), as.integer(round(value %% 60)))
  }, character(1), USE.NAMES = FALSE)
}

RUNLOG_DASH <- "\u2014"

# Most columns are handed to DT as plain text and escaped there, so the built-in
# search matches the real value and never the markup around it.
runlog_text_cell <- function(values) {
  ifelse(nzchar(values), values, RUNLOG_DASH)
}

# Only the three badge columns carry markup, and they are excluded from search --
# the pickers above the table are the way to select on them.
runlog_badge <- function(values, kind) {
  slug <- gsub("[^a-z0-9]+", "-", tolower(values))
  ifelse(
    nzchar(values),
    paste0("<span class=\"rl-badge rl-", kind, "-", slug, "\">",
           htmltools::htmlEscape(values), "</span>"),
    paste0("<span class=\"rl-muted\">", RUNLOG_DASH, "</span>")
  )
}

# Paths and notes are long enough to wreck the table. DataTables clips them at
# draw time and hangs the whole string off the title attribute; the cell value
# itself stays plain text, so search and the TSV download see the full string.
runlog_clip_js <- function(class, width) {
  DT::JS(sprintf(
    "function(data, type, row) {
       if (type !== 'display' || data === null) { return data; }
       var full = String(data);
       var short = full;
       if (full.length > %d) {
         short = full.substr(0, %d).replace(/&[^;]*$/, '') + '\u2026';
       }
       return '<span class=\"%s\" title=\"' + full.replace(/\"/g, '&quot;') +
              '\">' + short + '</span>';
     }", width, width - 1, class))
}

runlog_labelled <- function(values, blank_label) {
  ifelse(nzchar(values), values, blank_label)
}

# -------------------------------------------------------------------- style --

runlog_css <- function() {
  htmltools::tags$style(htmltools::HTML("
    .rl-wrap { --rl-sage-100: #edf4ea; --rl-sage-300: #cdddcf; --rl-sage-400: #8fad8d;
               --rl-sage-700: #35523c; --rl-soft: #5f7965; }
    .rl-metrics { display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
                  gap: 14px; }
    .rl-filters { display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr));
                  gap: 0 16px; align-items: end; }
    .rl-filters .form-group { margin-bottom: 10px; }
    .rl-actions { display: flex; gap: 10px; flex-wrap: wrap; align-items: center;
                  margin: 2px 0 12px 0; }
    .rl-summary { font-size: 13px; color: var(--rl-soft); line-height: 1.6; margin: 0; }
    .rl-summary b { color: var(--rl-sage-700); font-weight: 600; }
    .rl-userlist { margin-top: 10px; }
    .rl-chip { display: inline-block; background: var(--rl-sage-100);
      border: 1px solid var(--rl-sage-300); border-radius: 999px; padding: 1px 10px;
      margin: 0 6px 6px 0; color: var(--rl-sage-700); font-size: 12px; white-space: nowrap; }
    .rl-badge { display: inline-block; border-radius: 999px; padding: 1px 9px;
      font-size: 11px; font-weight: 600; letter-spacing: 0.02em;
      border: 1px solid transparent; white-space: nowrap; }
    .rl-muted { color: #9bab9d; }
    .rl-chain-b { background: #e4eef8; border-color: #c3d8ec; color: #2c4a68; }
    .rl-chain-t { background: #f6ecdf; border-color: #e6d3ba; color: #6b4c22; }
    .rl-chain-unknown { background: #f0f1ef; border-color: #dcdedb; color: #77807a; }
    .rl-era-pipeline { background: #e2efdf; border-color: #bcd6b8; color: #35523c; }
    .rl-era-pre-pipeline { background: #f2f1ec; border-color: #ddd9cc; color: #6f6a58; }
    .rl-status-completed { background: #dfeedd; border-color: #b3d2ae; color: #2f5c33; }
    .rl-status-failed { background: #fae6e3; border-color: #ecc3bc; color: #8a3b2e; }
    .rl-note-cell { font-size: 12px; color: var(--rl-soft); }
    .rl-path-cell { font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
                    font-size: 11.5px; color: var(--rl-soft); }
    .rl-footnote { font-size: 12px; color: var(--rl-soft); margin-top: 10px; }
    .rl-wrap table.dataTable { font-size: 13px; }
    .rl-wrap table.dataTable thead th { color: var(--rl-sage-700);
      border-bottom: 2px solid var(--rl-sage-300); }
    .rl-wrap .dataTables_wrapper .dataTables_filter input {
      border-radius: 12px; border: 1px solid var(--rl-sage-300); padding: 3px 10px; }
    @media (max-width: 767px) { .rl-filters { grid-template-columns: 1fr; } }
  "))
}

# ------------------------------------------------------------- input helpers --
# Use the app's picker wrappers when global.R is loaded, plain Shiny otherwise,
# so this module also works when sourced on its own.

runlog_has_widgets <- function() {
  isTRUE(get0("HAS_SHINYWIDGETS", envir = globalenv(), ifnotfound = FALSE)) &&
    exists("app_picker_input", mode = "function")
}

runlog_picker <- function(input_id, label, choices = NULL) {
  if (runlog_has_widgets()) {
    return(app_picker_input(
      input_id, label, choices = choices, selected = choices, multiple = TRUE,
      options = list(
        `live-search` = TRUE, `actions-box` = TRUE, `container` = "body",
        `selected-text-format` = "count > 2", `none-selected-text` = "All"
      )
    ))
  }
  shiny::selectizeInput(
    input_id, label, choices = choices, selected = choices, multiple = TRUE,
    options = list(placeholder = "All", plugins = list("remove_button"))
  )
}

runlog_update_picker <- function(session, input_id, choices, selected) {
  if (runlog_has_widgets()) {
    return(app_update_picker_input(session, input_id, choices = choices, selected = selected))
  }
  shiny::updateSelectizeInput(session, input_id, choices = choices, selected = selected)
}

# ----------------------------------------------------------------------- UI --

runlogUI <- function(id) {
  # Populate the filters here rather than only from the server. A pickerInput
  # built with choices = NULL reports "" back from the browser -- length 1, so it
  # slips past a NULL guard and filters every row away before the server has had
  # a chance to fill it in.
  initial <- local({
    log <- tryCatch(read_run_log(), error = function(...) run_log_empty())
    out <- lapply(names(RUN_LOG_BLANK), function(column) {
      if (is.null(log[[column]])) character(0)
      else sort(unique(runlog_labelled(log[[column]], RUN_LOG_BLANK[[column]])))
    })
    names(out) <- names(RUN_LOG_BLANK)
    out
  })

  ns <- shiny::NS(id)

  htmltools::tagList(
    runlog_css(),
    htmltools::tags$div(
      class = "rl-wrap",

      htmltools::tags$div(
        class = "rl-metrics",
        htmltools::tags$div(
          class = "metric-card",
          htmltools::tags$div(class = "metric-label", "Runs shown"),
          htmltools::tags$div(class = "metric-value",
                              shiny::textOutput(ns("metric_count"), inline = TRUE)),
          htmltools::tags$div(class = "metric-note",
                              shiny::textOutput(ns("metric_count_note"), inline = TRUE))
        ),
        htmltools::tags$div(
          class = "metric-card",
          htmltools::tags$div(class = "metric-label", "Date range"),
          htmltools::tags$div(class = "metric-value",
                              shiny::textOutput(ns("metric_range"), inline = TRUE)),
          htmltools::tags$div(class = "metric-note",
                              shiny::textOutput(ns("metric_range_note"), inline = TRUE))
        ),
        htmltools::tags$div(
          class = "metric-card",
          htmltools::tags$div(class = "metric-label", "User names"),
          htmltools::tags$div(class = "metric-value",
                              shiny::textOutput(ns("metric_users"), inline = TRUE)),
          htmltools::tags$div(class = "metric-note",
                              "Distinct names people submitted under")
        )
      ),

      htmltools::tags$div(
        class = "panel-card",
        htmltools::tags$div(class = "card-title", "Filter the log"),
        htmltools::tags$div(
          class = "card-note",
          paste(
            "Every run still evidenced on storage1 or storage2. Pipeline-era runs",
            "carry full metadata; older ones are reconstructed from result archives and",
            "raw-data folders, so their chain, primer set and outcome are unknown."
          )
        ),
        htmltools::tags$div(
          class = "rl-filters",
          runlog_picker(ns("f_user"), "User", initial$user),
          runlog_picker(ns("f_chain"), "Chain", initial$chain),
          runlog_picker(ns("f_primer"), "Primer set", initial$primer_set),
          runlog_picker(ns("f_era"), "Era", initial$era)
        ),
        htmltools::tags$div(
          class = "rl-actions",
          shiny::actionButton(ns("reset"), "Reset filters", class = "btn btn-default btn-sm"),
          shiny::downloadButton(ns("download"), "Download this view (TSV)",
                                class = "btn btn-default btn-sm")
        ),
        shiny::uiOutput(ns("summary")),
        htmltools::tags$div(class = "rl-userlist", shiny::uiOutput(ns("user_counts")))
      ),

      htmltools::tags$div(
        class = "panel-card",
        htmltools::tags$div(class = "card-title", "Run log"),
        htmltools::tags$div(
          class = "card-note",
          paste(
            "Search matches every column, the evidence path and the notes included.",
            "Click a header to sort; hover a clipped path or note to read it in full."
          )
        ),
        DT::dataTableOutput(ns("table")),
        htmltools::tags$div(class = "rl-footnote", shiny::uiOutput(ns("source_note")))
      )
    )
  )
}

# ------------------------------------------------------------------- SERVER --

runlogServer <- function(id, log_path = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    path <- if (is.null(log_path)) run_log_file() else log_path

    # Poll so a rebuilt audit appears without restarting the app; fall back to a
    # one-shot read when the file is not there yet.
    # Poll both files. reactiveFileReader watches a single path, and a new
    # submission only ever touches the live one -- so the tab sat on the audit's
    # mtime and never refreshed, which reads exactly like a log that has stopped
    # updating.
    run_log <- shiny::reactivePoll(
      5000, session,
      checkFunc = function() {
        stamp <- function(f) if (!is.null(f) && file.exists(f))
          paste(file.info(f)$mtime, file.info(f)$size) else ""
        paste(stamp(path), stamp(runlog_live_file()), sep = "|")
      },
      valueFunc = function() read_run_log(path)
    )

    filter_levels <- shiny::reactive({
      log <- run_log()
      lapply(names(RUN_LOG_BLANK), function(column) {
        sort(unique(runlog_labelled(log[[column]], RUN_LOG_BLANK[[column]])))
      }) -> levels
      names(levels) <- names(RUN_LOG_BLANK)
      levels
    })

    apply_all_levels <- function() {
      levels <- filter_levels()
      runlog_update_picker(session, "f_user", levels$user, levels$user)
      runlog_update_picker(session, "f_chain", levels$chain, levels$chain)
      runlog_update_picker(session, "f_primer", levels$primer_set, levels$primer_set)
      runlog_update_picker(session, "f_era", levels$era, levels$era)
    }

    shiny::observeEvent(filter_levels(), apply_all_levels())
    shiny::observeEvent(input$reset, apply_all_levels())

    filtered <- shiny::reactive({
      log <- run_log()
      if (nrow(log) == 0) {
        return(log)
      }

      # A NULL selection means the picker has not been populated yet -- show
      # everything rather than flashing an empty table on first paint.
      keep <- function(column, selected) {
        selected <- selected[nzchar(selected)]
        if (length(selected) == 0) {
          return(rep(TRUE, nrow(log)))
        }
        values <- runlog_labelled(log[[column]], RUN_LOG_BLANK[[column]])
        # No overlap at all means the control has not been populated yet, not
        # that the user deselected every option.
        if (!any(selected %in% values)) {
          return(rep(TRUE, nrow(log)))
        }
        values %in% selected
      }

      log[
        keep("user", input$f_user) &
          keep("chain", input$f_chain) &
          keep("primer_set", input$f_primer) &
          keep("era", input$f_era),
        , drop = FALSE
      ]
    })

    output$metric_count <- shiny::renderText(format(nrow(filtered()), big.mark = ","))

    output$metric_count_note <- shiny::renderText({
      total <- nrow(run_log())
      shown <- nrow(filtered())
      if (total == 0) "The run log file was not found"
      else if (shown == total) sprintf("the whole log (%s runs)", format(total, big.mark = ","))
      else sprintf("filtered from %s", format(total, big.mark = ","))
    })

    output$metric_range <- shiny::renderText({
      dates <- filtered()$date_parsed
      dates <- dates[!is.na(dates)]
      if (!length(dates)) {
        return("—")
      }
      sprintf("%s – %s", format(min(dates), "%b %Y"), format(max(dates), "%b %Y"))
    })

    output$metric_range_note <- shiny::renderText({
      undated <- sum(is.na(filtered()$date_parsed))
      if (nrow(filtered()) == 0) "Nothing to date"
      else if (undated == 0) "every run in this view is dated"
      else if (undated == 1) "1 run carries no date"
      else sprintf("%d runs carry no date", undated)
    })

    output$metric_users <- shiny::renderText({
      users <- filtered()$user
      format(length(unique(users[nzchar(users)])), big.mark = ",")
    })

    output$summary <- shiny::renderUI({
      log <- filtered()
      if (nrow(log) == 0) {
        message <- if (nrow(run_log()) == 0) {
          sprintf("No run log at %s — rebuild it with the storage audit.",
                  htmltools::htmlEscape(path))
        } else {
          "No run matches these filters."
        }
        return(htmltools::tags$p(class = "rl-summary", htmltools::HTML(message)))
      }

      dates <- log$date_parsed[!is.na(log$date_parsed)]
      span <- if (length(dates)) {
        sprintf("%s to %s", format(min(dates), "%d %b %Y"), format(max(dates), "%d %b %Y"))
      } else {
        "no dated run"
      }
      count <- function(n) format(n, big.mark = ",")
      pipeline <- sum(log$era == "pipeline")
      reconstructed <- nrow(log) - pipeline
      outcome <- sum(nzchar(log$status))

      sentences <- sprintf("<b>%s run%s</b> spanning <b>%s</b>.",
                           count(nrow(log)), if (nrow(log) == 1) "" else "s", span)

      # Only say a thing when there is something of that kind in the view.
      if (pipeline > 0 && reconstructed > 0) {
        sentences <- c(sentences, sprintf(
          paste0("<b>%s</b> carry pipeline metadata; <b>%s</b> are reconstructed ",
                 "from what survives on the storage servers."),
          count(pipeline), count(reconstructed)))
      } else if (pipeline > 0) {
        sentences <- c(sentences, "All of them carry pipeline metadata.")
      } else {
        sentences <- c(sentences, paste(
          "All of them are reconstructed from result archives and raw-data folders",
          "on the storage servers."))
      }

      # Two different reasons a run has no outcome, and the sentence should not
      # claim the wrong one: a ledger run may hold configuration and never have
      # been executed, while a reconstructed run simply predates the ledger.
      silent <- nrow(log) - outcome
      silent_pipeline <- sum(log$era == "pipeline" & !nzchar(log$status))
      why_silent <- if (silent_pipeline == silent) {
        sprintf("carry configuration but no execution record")
      } else if (silent_pipeline == 0) {
        sprintf("predate the run ledger, so only their evidence survives")
      } else {
        sprintf(paste("break down as %s with configuration but no execution record",
                      "and %s known only from storage"),
                count(silent_pipeline), count(silent - silent_pipeline))
      }

      if (outcome == 0) {
        sentences <- c(sentences, sprintf("None has a recorded outcome: they %s.", why_silent))
      } else {
        clause <- sprintf("An outcome is recorded for <b>%s</b> of them — %s completed, %s failed",
                          count(outcome), sum(log$status == "completed"),
                          sum(log$status == "failed"))
        clause <- if (silent > 0) {
          sprintf("%s; the other <b>%s</b> %s.", clause, count(silent), why_silent)
        } else {
          paste0(clause, ".")
        }
        sentences <- c(sentences, clause)
      }

      htmltools::tags$p(class = "rl-summary",
                        htmltools::HTML(paste(sentences, collapse = " ")))
    })

    output$user_counts <- shiny::renderUI({
      log <- filtered()
      users <- log$user[nzchar(log$user)]
      if (!length(users)) {
        return(htmltools::tags$span(class = "rl-muted",
                                    "No user name is established for these runs."))
      }

      counts <- sort(table(users), decreasing = TRUE)
      top <- utils::head(counts, 12)
      chips <- lapply(seq_along(top), function(i) {
        htmltools::tags$span(class = "rl-chip",
                             sprintf("%s · %d", names(top)[[i]], as.integer(top[[i]])))
      })
      remaining <- length(counts) - length(top)

      htmltools::tagList(
        chips,
        if (remaining > 0) {
          htmltools::tags$span(class = "rl-muted", sprintf(
            "and %d more name%s", remaining, if (remaining == 1) "" else "s"))
        }
      )
    })

    output$source_note <- shiny::renderUI({
      htmltools::HTML(sprintf(
        "Source: <code>%s</code> — %s rows, read-only. Rebuilt by the storage audit.",
        htmltools::htmlEscape(path), format(nrow(run_log()), big.mark = ",")
      ))
    })

    output$table <- DT::renderDataTable({
      log <- filtered()

      display <- data.frame(
        `#` = log$run_index,
        Run = runlog_text_cell(log$run_name),
        Date = runlog_text_cell(log$date),
        User = runlog_text_cell(log$user),
        Chain = runlog_badge(log$chain, "chain"),
        Primer = runlog_text_cell(log$primer_set),
        Server = runlog_text_cell(log$storage_server),
        Era = runlog_badge(log$era, "era"),
        Status = runlog_badge(log$status, "status"),
        Input = runlog_format_bytes(log$input_bytes_num),
        Runtime = runlog_format_runtime(log$runtime_sec_num),
        Evidence = runlog_text_cell(log$evidence),
        Path = runlog_text_cell(log$path),
        Notes = runlog_text_cell(log$notes),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )

      # Indices derived from the built table, never hardcoded. They were written
      # for a 15-column layout; dropping the email column left targets = 14 out
      # of range, DataTables failed to render, and the table reported "No
      # matching records found" over a full 308 rows.
      badge_cols <- match(c("Chain", "Era", "Status"), names(display))
      path_col <- match("Path", names(display))
      note_col <- match("Notes", names(display))
      # Most recent run first. Dates are ISO (YYYY-MM-DD) so they sort correctly
      # as text; the run index breaks ties, since many runs share a date.
      # Derived by name for the same reason the escape and columnDefs targets
      # are: a hardcoded position silently went out of range when a column was
      # removed, and the table reported "No matching records found" across all
      # 308 rows.
      date_col <- match("Date", names(display))
      index_col <- match("#", names(display))
      default_order <- if (is.na(date_col)) list(list(0, "desc")) else
        list(list(as.numeric(date_col - 1), "desc"))
      if (!is.na(index_col)) {
        default_order <- c(default_order,
                           list(list(as.numeric(index_col - 1), "desc")))
      }

      DT::datatable(
        display,
        # Escape every column except the three that carry badge markup.
        escape = -badge_cols,
        rownames = FALSE,
        selection = "none",
        class = "stripe hover row-border",
        options = list(
          pageLength = 25,
          lengthMenu = list(c(10, 25, 50, 100, -1), c("10", "25", "50", "100", "All")),
          order = default_order,
          scrollX = TRUE,
          autoWidth = FALSE,
          search = list(regex = FALSE, caseInsensitive = TRUE),
          columnDefs = list(
            list(targets = 0, width = "44px"),
            # Badge columns hold markup, so searching them would match class
            # names -- the pickers above the table cover them instead.
            # DataTables targets are 0-based; match() returns 1-based.
            list(targets = as.numeric(badge_cols - 1), searchable = FALSE),
            list(targets = as.numeric(path_col - 1), render = runlog_clip_js("rl-path-cell", 70)),
            list(targets = as.numeric(note_col - 1), render = runlog_clip_js("rl-note-cell", 110))
          ),
          dom = "<'row'<'col-sm-6'l><'col-sm-6'f>>tip",
          language = list(
            search = "Search:",
            emptyTable = "No run matches these filters.",
            info = "Runs _START_ to _END_ of _TOTAL_",
            infoEmpty = "No run to show",
            lengthMenu = "Show _MENU_ runs"
          )
        )
      )
    }, server = TRUE)

    output$download <- shiny::downloadHandler(
      filename = function() sprintf("tbrc-run-log-%s.tsv", format(Sys.Date(), "%Y-%m-%d")),
      content = function(file) {
        utils::write.table(
          filtered()[, RUN_LOG_COLUMNS, drop = FALSE],
          file, sep = "\t", quote = FALSE, row.names = FALSE, na = ""
        )
      }
    )

    # Handy for the caller: the currently filtered rows.
    shiny::reactive(filtered())
  })
}

# ---------------------------------------------------------------------------
# Integration (two lines, added by the caller -- not done here):
#
#   ui.R      tabPanel("Run log", { source(app_path("tab_runlog.R")); runlogUI("runlog") })
#   server.R  runlogServer("runlog")
#
# If you would rather keep the source() out of ui.R, put
#   source(app_path("tab_runlog.R"))
# at the end of global.R and the ui.R line becomes just
#   tabPanel("Run log", runlogUI("runlog"))
# ---------------------------------------------------------------------------
