# "What changed" tab -- the project's own history, read from docs/changelog.tsv.
#
# That file is generated from the git log by scripts/build_changelog.sh rather
# than maintained by hand: a hand-kept changelog drifts the moment someone
# forgets to update it, and the commit history is the only record that cannot.

CHANGELOG_COLUMNS <- c("date", "commit", "kind", "findings", "title", "detail")

changelog_file <- function() {
  configured <- Sys.getenv("TBRC_CHANGELOG_FILE", unset = "")
  if (nzchar(configured)) {
    return(configured)
  }
  candidates <- c(
    if (exists("app_path", mode = "function")) app_path("..", "docs", "changelog.tsv"),
    file.path("..", "docs", "changelog.tsv"),
    file.path("docs", "changelog.tsv")
  )
  for (candidate in candidates) {
    if (!is.null(candidate) && file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = FALSE))
    }
  }
  if (length(candidates)) candidates[[1]] else "changelog.tsv"
}

read_changelog <- function(path = changelog_file()) {
  empty <- as.data.frame(
    setNames(replicate(length(CHANGELOG_COLUMNS), character(0), simplify = FALSE), CHANGELOG_COLUMNS),
    stringsAsFactors = FALSE
  )
  if (!file.exists(path)) {
    return(empty)
  }
  out <- tryCatch(
    utils::read.table(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
                      stringsAsFactors = FALSE, colClasses = "character", fill = TRUE),
    error = function(...) NULL
  )
  if (is.null(out) || !all(CHANGELOG_COLUMNS %in% names(out))) {
    return(empty)
  }
  out[order(out$date, decreasing = TRUE), CHANGELOG_COLUMNS, drop = FALSE]
}

changelog_kind_badge <- function(kind) {
  cls <- switch(kind,
    fix = "cl-badge cl-fix",
    feature = "cl-badge cl-feature",
    housekeeping = "cl-badge cl-house",
    "cl-badge cl-change"
  )
  sprintf('<span class="%s">%s</span>', cls, htmltools::htmlEscape(kind))
}

changelogUI <- function(id) {
  ns <- shiny::NS(id)
  htmltools::tagList(
    htmltools::tags$style(htmltools::HTML("
      .cl-wrap .cl-entry { border-left: 3px solid var(--sage-300, #cdddcf);
        padding: 10px 0 10px 14px; margin-bottom: 4px; }
      .cl-wrap .cl-entry.is-fix { border-left-color: #2a78d6; }
      .cl-wrap .cl-entry.is-feature { border-left-color: #eb6834; }
      .cl-wrap .cl-head { display: flex; flex-wrap: wrap; align-items: baseline; gap: 10px; }
      .cl-wrap .cl-date { font-family: ui-monospace, Menlo, monospace; font-size: 12px;
        color: #5f7965; min-width: 84px; }
      .cl-wrap .cl-title { font-weight: 600; color: #294034; flex: 1 1 320px; }
      .cl-wrap .cl-detail { font-size: 13px; color: #4d6954; margin-top: 4px; line-height: 1.5; }
      .cl-badge { font-size: 10px; letter-spacing: .08em; text-transform: uppercase;
        padding: 2px 7px; border-radius: 3px; background: #edf4ea; color: #35523c; }
      .cl-badge.cl-fix { background: #e3edfa; color: #1c4f8f; }
      .cl-badge.cl-feature { background: #fbe9e1; color: #9c3d16; }
      .cl-badge.cl-house { background: #eef0f1; color: #55605c; }
      .cl-finding { font-family: ui-monospace, Menlo, monospace; font-size: 11px; color: #7a6a2f;
        background: #faf3dd; padding: 1px 6px; border-radius: 3px; }
      .cl-day { font-family: ui-monospace, Menlo, monospace; font-size: 11px; letter-spacing: .09em;
        text-transform: uppercase; color: #7d8f82; margin: 18px 0 6px; }
    ")),
    htmltools::tags$div(
      class = "cl-wrap",
      htmltools::tags$div(
        class = "panel-card",
        htmltools::tags$div(class = "card-title", "What changed"),
        htmltools::tags$div(
          class = "card-note",
          paste(
            "Generated from the commit history, newest first -- not maintained by hand,",
            "so it cannot drift. Entries tagged F1-F16 refer to findings from the",
            "August 2026 audit."
          )
        ),
        htmltools::tags$div(
          class = "rl-filters",
          style = "display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;margin-top:12px;",
          shiny::selectizeInput(ns("kind"), "Kind", choices = NULL, selected = NULL,
                                multiple = TRUE, options = list(placeholder = "All")),
          shiny::textInput(ns("q"), "Search", placeholder = "e.g. barcode, light chain, F10")
        ),
        shiny::uiOutput(ns("summary"))
      ),
      htmltools::tags$div(class = "panel-card", shiny::uiOutput(ns("entries")))
    )
  )
}

changelogServer <- function(id, path = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    log <- shiny::reactive({
      if (is.null(path)) read_changelog() else read_changelog(path)
    })

    shiny::observeEvent(log(), {
      kinds <- sort(unique(log()$kind))
      shiny::updateSelectizeInput(session, "kind", choices = kinds, selected = kinds)
    })

    filtered <- shiny::reactive({
      d <- log()
      if (nrow(d) == 0) return(d)
      sel <- input$kind
      sel <- sel[nzchar(sel)]
      if (length(sel) && any(sel %in% d$kind)) {
        d <- d[d$kind %in% sel, , drop = FALSE]
      }
      q <- trimws(input$q %||% "")
      if (nzchar(q)) {
        hay <- tolower(paste(d$title, d$detail, d$findings))
        d <- d[grepl(tolower(q), hay, fixed = TRUE), , drop = FALSE]
      }
      d
    })

    output$summary <- shiny::renderUI({
      d <- filtered(); all <- log()
      if (nrow(all) == 0) {
        return(htmltools::tags$div(class = "card-note", "No changelog found. Run scripts/build_changelog.sh."))
      }
      fixes <- sum(all$kind == "fix")
      htmltools::tags$div(
        class = "card-note", style = "margin-top:6px;",
        sprintf("Showing %d of %d entries. %d fixes since %s.",
                nrow(d), nrow(all), fixes, min(all$date[nzchar(all$date)]))
      )
    })

    output$entries <- shiny::renderUI({
      d <- filtered()
      if (nrow(d) == 0) {
        return(htmltools::tags$div(class = "card-note", "Nothing matches."))
      }
      blocks <- list(); last_date <- ""
      for (i in seq_len(nrow(d))) {
        row <- d[i, ]
        if (!identical(row$date, last_date)) {
          blocks[[length(blocks) + 1]] <- htmltools::tags$div(class = "cl-day", row$date)
          last_date <- row$date
        }
        findings <- trimws(row$findings)
        blocks[[length(blocks) + 1]] <- htmltools::tags$div(
          class = paste0("cl-entry is-", row$kind),
          htmltools::tags$div(
            class = "cl-head",
            htmltools::HTML(changelog_kind_badge(row$kind)),
            htmltools::tags$span(class = "cl-title", row$title),
            if (nzchar(findings)) htmltools::tags$span(class = "cl-finding", findings),
            htmltools::tags$span(class = "cl-date", row$commit)
          ),
          if (nzchar(trimws(row$detail)))
            htmltools::tags$div(class = "cl-detail", row$detail)
        )
      }
      htmltools::tagList(blocks)
    })
  })
}

`%||%` <- function(a, b) if (is.null(a)) b else a
