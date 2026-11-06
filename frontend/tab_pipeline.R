# "Pipelines" tab -- the Mermaid diagrams in docs/pipeline-diagrams.md.
#
# The diagrams live in that one Markdown file so GitHub renders the same source
# the app shows; drawing them here would give two copies to keep in step. The
# file is read per session, like the changelog, so a merged edit reaches the
# next visitor without a container restart.
#
# Mermaid sizes each box from its label measured in the font it was given, so
# that font must be the one the SVG inherits: inside a <pre> the labels came out
# in monospace, wider than measured, and were cut off. In a hidden tab every
# size is zero, so the diagrams are drawn only once their tab is visible: on the
# output's arrival (Shiny holds hidden outputs until shown) and on every tab
# change, each pass picking up only blocks not drawn yet.

MERMAID_JS <- "https://cdn.jsdelivr.net/npm/mermaid@11.4.1/dist/mermaid.min.js"

pipeline_diagrams_file <- function() {
  configured <- Sys.getenv("TBRC_DIAGRAMS_FILE", unset = "")
  if (nzchar(configured)) {
    return(configured)
  }
  candidates <- c(
    if (exists("app_path", mode = "function")) app_path("..", "docs", "pipeline-diagrams.md"),
    file.path("..", "docs", "pipeline-diagrams.md"),
    file.path("docs", "pipeline-diagrams.md")
  )
  for (candidate in candidates) {
    if (!is.null(candidate) && file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = FALSE))
    }
  }
  if (length(candidates)) candidates[[1]] else "pipeline-diagrams.md"
}

# Returns list(intro = <text>, sections = list(list(title, note, code), ...)).
# A section is a "## " heading, the prose under it, and its first ```mermaid
# block; a section without a block is dropped.
read_pipeline_diagrams <- function(path = pipeline_diagrams_file()) {
  out <- list(intro = "", sections = list())
  if (!file.exists(path)) {
    return(out)
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  intro <- character(0); sections <- list(); cur <- NULL; in_code <- FALSE
  for (line in lines) {
    if (in_code) {
      if (grepl("^```\\s*$", line)) {
        in_code <- FALSE
      } else {
        cur$code <- c(cur$code, line)
      }
    } else if (grepl("^## ", line)) {
      if (!is.null(cur)) sections[[length(sections) + 1]] <- cur
      cur <- list(title = trimws(sub("^## ", "", line)), note = character(0), code = character(0))
    } else if (is.null(cur)) {
      if (!grepl("^# ", line)) intro <- c(intro, line)
    } else if (grepl("^```mermaid\\s*$", line) && length(cur$code) == 0) {
      in_code <- TRUE
    } else if (length(cur$code) == 0) {
      cur$note <- c(cur$note, line)
    }
  }
  if (!is.null(cur)) sections[[length(sections) + 1]] <- cur
  sections <- Filter(function(s) length(s$code) > 0, sections)
  out$intro <- trimws(paste(intro, collapse = " "))
  out$sections <- lapply(sections, function(s) list(
    title = s$title,
    note = trimws(paste(s$note, collapse = " ")),
    code = paste(s$code, collapse = "\n")
  ))
  out
}

# The notes are Markdown prose; only `code` spans are worth rendering.
pipeline_note_html <- function(text) {
  htmltools::HTML(gsub("`([^`]+)`", "<code>\\1</code>", htmltools::htmlEscape(text)))
}

pipeline_anchor <- function(title) {
  paste0("pd-", gsub("[^a-z0-9]+", "-", tolower(title)))
}

pipelineUI <- function(id) {
  ns <- shiny::NS(id)
  htmltools::tagList(
    htmltools::tags$script(src = MERMAID_JS),
    htmltools::tags$style(htmltools::HTML("
      .pd-wrap .pd-jump { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 10px; }
      .pd-wrap .pd-jump a { font-size: 12px; padding: 3px 10px; border-radius: 999px;
        background: #edf4ea; color: #35523c; border: 1px solid #cdddcf; text-decoration: none; }
      .pd-wrap .pd-jump a:hover { background: #dde8db; }
      .pd-wrap .pd-figure { overflow-x: auto; padding: 8px 0 4px; }
      .pd-wrap .tbrc-mermaid { text-align: center; }
      .pd-wrap .tbrc-mermaid:not([data-processed]) { text-align: left; white-space: pre;
        font-family: ui-monospace, Menlo, monospace; font-size: 11px; color: #5f7965; }
      .pd-wrap details { margin-top: 8px; font-size: 12px; color: #5f7965; }
      .pd-wrap details pre { font-size: 11px; max-height: 320px; overflow: auto; }
    ")),
    htmltools::tags$script(htmltools::HTML(sprintf("
      (function () {
        var started = false;
        function render() {
          if (!window.mermaid) return;
          if (!started) {
            mermaid.initialize({
              startOnLoad: false, securityLevel: 'strict', theme: 'base',
              themeVariables: {
                fontFamily: getComputedStyle(document.body).fontFamily, fontSize: '14px',
                primaryColor: '#edf4ea', primaryBorderColor: '#8fad8d', primaryTextColor: '#294034',
                lineColor: '#5f7965', secondaryColor: '#f3f7f1', tertiaryColor: '#f7faf6',
                clusterBkg: '#f7faf6', clusterBorder: '#cdddcf', edgeLabelBackground: '#ffffff'
              },
              flowchart: { useMaxWidth: true, htmlLabels: true }
            });
            started = true;
          }
          var nodes = Array.prototype.filter.call(
            document.querySelectorAll('.tbrc-mermaid:not([data-processed])'),
            function (n) { return n.offsetParent !== null; });
          if (nodes.length) mermaid.run({ nodes: nodes });
        }
        $(document).on('shown.bs.tab', render);
        $(document).on('shiny:value', function (e) {
          if (e.name === '%s') setTimeout(render, 0);
        });
      })();
    ", ns("diagrams")))),
    htmltools::tags$div(class = "pd-wrap", shiny::uiOutput(ns("diagrams")))
  )
}

pipelineServer <- function(id, path = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    output$diagrams <- shiny::renderUI({
      d <- if (is.null(path)) read_pipeline_diagrams() else read_pipeline_diagrams(path)
      if (length(d$sections) == 0) {
        return(htmltools::tags$div(
          class = "panel-card",
          htmltools::tags$div(class = "card-title", "Pipelines"),
          htmltools::tags$div(class = "card-note", "No diagrams found: docs/pipeline-diagrams.md is missing.")
        ))
      }
      header <- htmltools::tags$div(
        class = "panel-card",
        htmltools::tags$div(class = "card-title", "Pipelines"),
        if (nzchar(d$intro)) htmltools::tags$div(class = "card-note", pipeline_note_html(d$intro)),
        htmltools::tags$div(
          class = "pd-jump",
          lapply(d$sections, function(s) htmltools::tags$a(href = paste0("#", pipeline_anchor(s$title)), s$title))
        )
      )
      cards <- lapply(d$sections, function(s) htmltools::tags$div(
        class = "panel-card", id = pipeline_anchor(s$title),
        htmltools::tags$div(class = "card-title", s$title),
        if (nzchar(s$note)) htmltools::tags$div(class = "card-note", pipeline_note_html(s$note)),
        htmltools::tags$div(class = "pd-figure", htmltools::tags$div(class = "tbrc-mermaid", s$code)),
        htmltools::tags$details(
          htmltools::tags$summary("Mermaid source"),
          htmltools::tags$pre(s$code)
        )
      ))
      htmltools::tagList(header, cards)
    })
  })
}
