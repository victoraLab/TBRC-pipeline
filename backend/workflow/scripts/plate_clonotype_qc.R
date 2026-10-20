args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript plate_clonotype_qc.R <clonality_or_igblast_tsv> <output_dir>", call. = FALSE)
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
  stop("plate_clonotype_qc.R requires ggplate.", call. = FALSE)
}

input_tsv <- args[1]
output_dir <- args[2]

if (!file.exists(input_tsv)) {
  stop(sprintf("Input TSV not found: %s", input_tsv), call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_df <- read.table(
  input_tsv,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!nrow(input_df)) {
  stop("Input TSV has no rows for plate clonotype QC.", call. = FALSE)
}

pick_col <- function(candidates) {
  match_idx <- match(tolower(candidates), tolower(names(input_df)))
  found <- names(input_df)[match_idx[!is.na(match_idx)]]
  if (length(found)) found[[1]] else NA_character_
}

sequence_id_col <- pick_col(c("sequence_id", "sequence_id_1"))
if (is.na(sequence_id_col)) {
  stop("Input TSV is missing sequence_id, so plate clonotype QC cannot parse wells.", call. = FALSE)
}

parse_header <- function(cell_id) {
  match <- regexec("^(.*?)([A-H][0-9]{2}P[0-9]{2}|P[0-9]{2}[A-H][0-9]{2})_([0-9]+)-([0-9]+)$", cell_id)
  parts <- regmatches(cell_id, match)[[1]]

  if (length(parts) != 5) {
    stop(sprintf("Failed to parse plate metadata from sequence_id: %s", cell_id), call. = FALSE)
  }

  plate_well <- parts[3]
  plate <- sub(".*(P[0-9]{2}).*", "\\1", plate_well)
  well <- sub(".*([A-H][0-9]{2}).*", "\\1", plate_well)

  data.frame(
    sequence_id = cell_id,
    plate = plate,
    well = well,
    ggplate_well = sub("^([A-H])0?([0-9]+)$", "\\1\\2", well),
    contig_number = as.integer(parts[4]),
    read_depth = as.numeric(parts[5]),
    stringsAsFactors = FALSE
  )
}

header_meta <- do.call(rbind, lapply(input_df[[sequence_id_col]], parse_header))

clonality_col <- pick_col(c("clonality"))
v_call_col <- pick_col(c("v_call", "best_v_hit", "V-GENE and allele", "v_gene"))
j_call_col <- pick_col(c("j_call", "best_j_hit", "J-GENE and allele", "j_gene"))
junction_col <- pick_col(c("junction", "junction_aa", "cdr3", "CDR3-IMGT"))

if (!is.na(clonality_col)) {
  clonality_value <- trimws(as.character(input_df[[clonality_col]]))
  clonality_primary <- sub("\\..*$", "", clonality_value)
  clonality_primary[!nzchar(clonality_primary)] <- "unassigned"
  clonotype_label <- clonality_value
  clonotype_family <- clonality_primary
  source_mode <- "clonality"
} else {
  if (anyNA(c(v_call_col, j_call_col, junction_col))) {
    stop("Input TSV is missing both clonality labels and the IGBlast columns needed to derive clonotypes.", call. = FALSE)
  }

  # Gene of the first call, as the clonality package keys clones: the raw
  # strings carry tie lists and allele numbers, which split one clone whenever
  # a tie resolves differently in two wells (test results/igkj1-2026-10-19/).
  gene_only <- function(x) sub("\\*.*$", "", sub(",.*$", "", x))
  v_val <- gene_only(trimws(as.character(input_df[[v_call_col]])))
  j_val <- gene_only(trimws(as.character(input_df[[j_call_col]])))
  junction_val <- gsub("[^A-Za-z]", "", trimws(as.character(input_df[[junction_col]])))

  clonotype_family <- paste(v_val, j_val, sep = "|")
  clonotype_family[!nzchar(v_val) | !nzchar(j_val)] <- "unassigned"
  clonotype_label <- paste(clonotype_family, junction_val, sep = "|")
  clonotype_label[!nzchar(junction_val)] <- clonotype_family[!nzchar(junction_val)]
  source_mode <- "igblast"
}

qc_df <- data.frame(
  sequence_id = as.character(input_df[[sequence_id_col]]),
  clonotype_label = clonotype_label,
  clonotype_family = clonotype_family,
  stringsAsFactors = FALSE
)
qc_df <- merge(qc_df, header_meta, by = "sequence_id", all.x = TRUE, sort = FALSE)
qc_df <- qc_df[order(qc_df$plate, qc_df$well, -qc_df$read_depth, qc_df$contig_number), , drop = FALSE]

score_well <- function(well_df) {
  well_df <- well_df[order(-well_df$read_depth, well_df$contig_number), , drop = FALSE]
  total_depth <- sum(well_df$read_depth, na.rm = TRUE)
  if (!is.finite(total_depth) || total_depth <= 0) {
    total_depth <- nrow(well_df)
    well_df$read_depth <- 1
  }

  dominant <- well_df[1, , drop = FALSE]
  dominant_family <- dominant$clonotype_family[[1]]
  dominant_label <- dominant$clonotype_label[[1]]
  dominant_fraction <- dominant$read_depth[[1]] / total_depth

  other_family_mask <- well_df$clonotype_family != dominant_family
  contaminant_depth <- sum(well_df$read_depth[other_family_mask], na.rm = TRUE)
  contaminant_fraction <- contaminant_depth / total_depth
  alt_same_family_depth <- sum(well_df$read_depth[!other_family_mask & well_df$clonotype_label != dominant_label], na.rm = TRUE)
  alt_same_family_fraction <- alt_same_family_depth / total_depth
  unique_family_count <- length(unique(well_df$clonotype_family))
  top_two_fraction <- if (nrow(well_df) >= 2) well_df$read_depth[[2]] / total_depth else 0

  status <- "green"
  rationale <- "single dominant clonotype"

  if (unique_family_count > 1 && contaminant_fraction >= 0.35) {
    status <- "red"
    rationale <- sprintf("multiple clonotype families at comparable depth (%.0f%% secondary)", contaminant_fraction * 100)
  } else if (unique_family_count > 1 && contaminant_fraction >= 0.10) {
    status <- "yellow"
    rationale <- sprintf("secondary clonotype family detected (%.0f%% secondary)", contaminant_fraction * 100)
  } else if (alt_same_family_fraction >= 0.35 || (unique_family_count == 1 && top_two_fraction >= 0.35)) {
    status <- "yellow"
    rationale <- sprintf("same-family isoform mixture above alert threshold (%.0f%% secondary)", max(alt_same_family_fraction, top_two_fraction) * 100)
  }

  risk_score <- c(green = 1, yellow = 2, red = 3)[[status]]
  plot_label <- c(green = "OK", yellow = "!", red = "!!")[[status]]

  data.frame(
    plate = well_df$plate[[1]],
    well = well_df$well[[1]],
    ggplate_well = well_df$ggplate_well[[1]],
    source_mode = source_mode,
    dominant_clonotype = dominant_label,
    dominant_family = dominant_family,
    contig_count = nrow(well_df),
    total_read_depth = total_depth,
    dominant_fraction = dominant_fraction,
    contaminant_fraction = contaminant_fraction,
    alt_same_family_fraction = alt_same_family_fraction,
    unique_family_count = unique_family_count,
    risk_status = status,
    risk_score = risk_score,
    plot_label = plot_label,
    rationale = rationale,
    stringsAsFactors = FALSE
  )
}

well_keys <- unique(qc_df[c("plate", "well")])
well_scores <- do.call(
  rbind,
  lapply(seq_len(nrow(well_keys)), function(i) {
    key <- well_keys[i, , drop = FALSE]
    score_well(qc_df[qc_df$plate == key$plate & qc_df$well == key$well, , drop = FALSE])
  })
)

well_scores$plate_numeric <- as.integer(sub("^P", "", well_scores$plate))
well_scores$plate_label <- sprintf("Plate %02d", well_scores$plate_numeric)
well_scores <- well_scores[order(well_scores$plate_numeric, well_scores$ggplate_well), , drop = FALSE]

write.table(
  well_scores,
  file = file.path(output_dir, "plate_clonotype_contamination.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

make_plate_plot <- function(plate_id) {
  plate_df <- well_scores[well_scores$plate == plate_id, c("ggplate_well", "risk_score", "plot_label")]
  colnames(plate_df) <- c("position", "value", "label")

  ggplate::plate_plot(
    data = plate_df,
    position = position,
    value = value,
    label = label,
    plate_size = 96,
    plate_type = "round",
    colour = c("#4f9d69", "#e3bf4f", "#c75a57"),
    title = sprintf("Clonotype contamination risk: %s", plate_id),
    show_legend = TRUE,
    silent = TRUE
  ) +
    scale_fill_gradientn(
      colours = c("#4f9d69", "#e3bf4f", "#c75a57"),
      limits = c(1, 3),
      breaks = c(1, 2, 3),
      labels = c("Green", "Yellow", "Red")
    ) +
    theme(plot.title = element_text(face = "bold"))
}

plot_list <- lapply(unique(well_scores$plate), make_plate_plot)
names(plot_list) <- unique(well_scores$plate)

panel_png <- file.path(output_dir, "plate_clonotype_contamination.panel.png")
panel_pdf <- file.path(output_dir, "plate_clonotype_contamination.panel.pdf")

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
