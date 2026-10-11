args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript plate_depth_qc.R <input_fasta> <output_dir>", call. = FALSE)
}

suppressPackageStartupMessages(library(ggplot2))

# Victora-lab plot styling, vendored alongside this script.
local({
  theme_path <- file.path(
    dirname(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
    "lab_theme.R"
  )
  if (file.exists(theme_path)) source(theme_path, local = FALSE)
})
if (!exists("gvic_theme")) gvic_theme <- ggplot2::theme()
if (!exists("gvic_theme_plate")) gvic_theme_plate <- ggplot2::theme()

if (!requireNamespace("ggplate", quietly = TRUE)) {
  stop("plate_depth_qc.R requires ggplate.", call. = FALSE)
}

input_fasta <- args[1]
output_dir <- args[2]

if (!file.exists(input_fasta)) {
  stop(sprintf("FASTA file not found: %s", input_fasta), call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fasta_lines <- readLines(input_fasta, warn = FALSE)
headers <- sub("^>", "", fasta_lines[grepl("^>", fasta_lines)])

if (!length(headers)) {
  stop("No FASTA headers found for plate depth QC.", call. = FALSE)
}

header_df <- data.frame(cell_id = headers, stringsAsFactors = FALSE)

parse_header <- function(cell_id) {
  match <- regexec("^(.*?)([A-H][0-9]{2}P[0-9]{2}|P[0-9]{2}[A-H][0-9]{2})_([0-9]+)-([0-9]+)$", cell_id)
  parts <- regmatches(cell_id, match)[[1]]

  if (length(parts) != 5) {
    stop(sprintf("Failed to parse plate metadata from FASTA header: %s", cell_id), call. = FALSE)
  }

  plate_well <- parts[3]
  plate <- sub(".*(P[0-9]{2}).*", "\\1", plate_well)
  well <- sub(".*([A-H][0-9]{2}).*", "\\1", plate_well)

  data.frame(
    cell_id = cell_id,
    plate = plate,
    well = well,
    ggplate_well = sub("^([A-H])0?([0-9]+)$", "\\1\\2", well),
    contig_number = as.integer(parts[4]),
    read_depth = as.integer(parts[5]),
    stringsAsFactors = FALSE
  )
}

metadata <- do.call(rbind, lapply(headers, parse_header))
metadata$plate_numeric <- as.integer(sub("^P", "", metadata$plate))
metadata$plate_label <- sprintf("Plate %02d", metadata$plate_numeric)

depth_by_well <- aggregate(
  read_depth ~ plate + plate_numeric + plate_label + ggplate_well,
  data = metadata,
  FUN = sum
)

depth_by_well <- depth_by_well[order(depth_by_well$plate_numeric, depth_by_well$ggplate_well), ]

write.table(
  depth_by_well,
  file = file.path(output_dir, "plate_read_depth.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plate_palette <- c("#255f6e", "#78aebb", "#f3efe4", "#d39a6a", "#8a4b2d")
plates <- unique(depth_by_well$plate)

make_plate_plot <- function(plate_id) {
  plate_df <- depth_by_well[depth_by_well$plate == plate_id, c("ggplate_well", "read_depth")]
  colnames(plate_df) <- c("position", "value")
  plate_df$label <- plate_df$value

  ggplate::plate_plot(
    data = plate_df,
    position = position,
    value = value,
    label = label,
    plate_size = 96,
    plate_type = "round",
    colour = plate_palette,
    title = sprintf("Read depth: %s", plate_id),
    show_legend = TRUE,
    silent = TRUE
  ) +
    theme(plot.title = element_text(face = "bold"))
}

plot_list <- lapply(plates, make_plate_plot)
names(plot_list) <- plates

panel_png <- file.path(output_dir, "plate_read_depth.panel.png")
panel_pdf <- file.path(output_dir, "plate_read_depth.panel.pdf")

draw_panel <- function(device_fn, path) {
  n_plates <- length(plot_list)
  n_cols <- min(2, n_plates)
  n_rows <- ceiling(n_plates / n_cols)

  if (identical(device_fn, png)) {
    device_fn(path, width = 8.5 * n_cols, height = 6.5 * n_rows, units = "in", res = 180)
  } else {
    device_fn(path, width = 8.5 * n_cols, height = 6.5 * n_rows)
  }

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n_rows, n_cols)))

  for (idx in seq_along(plot_list)) {
    row_idx <- ceiling(idx / n_cols)
    col_idx <- ((idx - 1) %% n_cols) + 1
    print(
      plot_list[[idx]],
      vp = grid::viewport(layout.pos.row = row_idx, layout.pos.col = col_idx)
    )
  }

  invisible(grDevices::dev.off())
}

draw_panel(png, panel_png)
draw_panel(pdf, panel_pdf)
