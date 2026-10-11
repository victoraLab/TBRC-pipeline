args <- commandArgs(trailingOnly = TRUE)

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

if (length(args) < 3) {
  stop("Usage: Rscript fasta_qc.R <input_fasta> <output_dir> <label>", call. = FALSE)
}

input_fasta <- args[1]
output_dir <- args[2]
label <- args[3]

# Read sequence lengths directly from FASTA records so the QC step stays
# independent from any Bioconductor/FASTA parser installation.
read_fasta_lengths <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("FASTA file not found: %s", path), call. = FALSE)
  }

  lines <- readLines(path, warn = FALSE)
  if (!length(lines)) {
    return(integer())
  }

  header_idx <- grep("^>", lines)
  if (!length(header_idx)) {
    return(integer())
  }

  sequence_lengths <- integer(length(header_idx))
  for (i in seq_along(header_idx)) {
    start_idx <- header_idx[i] + 1
    end_idx <- if (i < length(header_idx)) header_idx[i + 1] - 1 else length(lines)
    sequence_lines <- lines[start_idx:end_idx]
    sequence_lines <- sequence_lines[nzchar(sequence_lines)]
    sequence_lengths[i] <- nchar(paste(sequence_lines, collapse = ""))
  }

  sequence_lengths
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

lengths <- read_fasta_lengths(input_fasta)
safe_label <- gsub("[^A-Za-z0-9._-]+", "_", label)
tsv_path <- file.path(output_dir, paste0(safe_label, ".contig_lengths.tsv"))
png_path <- file.path(output_dir, paste0(safe_label, ".contig_lengths.png"))
pdf_path <- file.path(output_dir, paste0(safe_label, ".contig_lengths.pdf"))

summary_df <- data.frame(
  label = label,
  sequence_count = length(lengths),
  min_length = if (length(lengths)) min(lengths) else NA_integer_,
  median_length = if (length(lengths)) median(lengths) else NA_real_,
  mean_length = if (length(lengths)) mean(lengths) else NA_real_,
  max_length = if (length(lengths)) max(lengths) else NA_integer_
)

write.table(summary_df, file = tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)

plot_lengths <- function(device_fn, path) {
  if (identical(device_fn, png)) {
    device_fn(path, width = 9, height = 5, units = "in", res = 180)
  } else {
    device_fn(path, width = 9, height = 5)
  }

  if (length(lengths)) {
    lengths_df <- data.frame(length = lengths)
    x_range <- range(lengths)
    # Keep the default plot shape stable across runs by targeting ~30 bins.
    binwidth <- max(1, ceiling((x_range[2] - x_range[1]) / 30))

    plot_obj <- ggplot(lengths_df, aes(x = length)) +
      geom_histogram(
        binwidth = binwidth,
        fill = "#9ab49a",
        color = "#294034",
        linewidth = 0.45,
        alpha = 0.92,
        boundary = 0
      ) +
      # Scale the density curve onto the histogram count axis so both shapes
      # can be read together in a single panel.
      geom_density(
        aes(y = after_stat(count * binwidth)),
        color = "#163028",
        linewidth = 1.05,
        adjust = 1.05
      ) +
      labs(
        title = paste("Contig length distribution:", label),
        subtitle = sprintf("n=%s | median=%.1f | mean=%.1f", length(lengths), median(lengths), mean(lengths)),
        x = "Contig length (nt)",
        y = "Count"
      ) +
      theme_bw(base_size = 12) +
      gvic_theme +
      theme(
        plot.title = element_text(face = "bold", color = "#294034"),
        plot.subtitle = element_text(color = "#4d6954"),
        axis.title = element_text(color = "#294034"),
        axis.text = element_text(color = "#35523c"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
      )

    print(plot_obj)
  } else {
    empty_plot <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No sequences found in FASTA", size = 5, color = "#35523c") +
      labs(title = paste("Contig length distribution:", label)) +
      theme_void(base_size = 12) +
      gvic_theme_plate +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#294034"))

    print(empty_plot)
  }

  invisible(dev.off())
}

plot_lengths(png, png_path)
plot_lengths(pdf, pdf_path)
