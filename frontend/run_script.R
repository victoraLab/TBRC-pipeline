# Run one of the bundled scripts and wait for it. Used for the pre-submission
# check and for syncing barcodes; the long-running dispatch scripts go through
# start_background_script() in global.R instead.
#
# system2() hands its command line to /bin/sh and quotes only the command name,
# not the arguments. Passing a form value as an argument therefore let anything in
# it (a trim sequence, a pipeline name, a user name) run as a command, as root,
# inside the container. Every argument, and the script path, is now quoted
# with shQuote() before it gets there, so the shell sees each one as a single
# word and never as syntax.
#
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
    suppressWarnings(system2(command_name, args = shQuote(command_args), stdout = TRUE,
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
