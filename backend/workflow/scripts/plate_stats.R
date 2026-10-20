args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript plate_stats.R <input_prefix> <output_dir> [sample_name]", call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

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

has_ggplate <- requireNamespace("ggplate", quietly = TRUE)
has_platetools <- requireNamespace("platetools", quietly = TRUE)

input_prefix <- args[1]
output_dir <- args[2]
sample_name <- if (length(args) >= 3) args[3] else ""
counts_file <- paste0(input_prefix, "BC.stats.txt")
quality_file <- paste0(input_prefix, "collapsed.stats.txt")

if (!file.exists(counts_file) || !file.exists(quality_file)) {
  stop("Plate QC inputs not found. Expected BC.stats.txt and collapsed.stats.txt.", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Which plates exist is a property of the run, not a constant. The layout used to
# be hardcoded to 12, so a run using plates 13-16 -- as RUN13 does -- had them
# silently dropped from every plate figure and table.
build_plate_layout <- function(plate_ids) {
  plate_numbers <- sort(unique(suppressWarnings(as.integer(sub("^plate_", "", plate_ids)))))
  plate_numbers <- plate_numbers[is.finite(plate_numbers)]

  if (length(plate_numbers) == 0) {
    plate_numbers <- 1:12
  }

  levels_out <- sprintf("plate_%d", plate_numbers)
  n <- length(levels_out)

  wells <- if (has_platetools) {
    rep(platetools::num_to_well(1:96), n)
  } else {
    rep(paste0(rep(LETTERS[1:8], each = 12), rep(1:12, times = 8)), n)
  }

  list(
    levels = levels_out,
    layout = data.frame(
      plate_id = rep(levels_out, each = 96),
      well = wells,
      plate_well_id = paste0(rep(levels_out, each = 96), "_", wells),
      stringsAsFactors = FALSE
    )
  )
}

plate_layout <- data.frame(
  plate_id = character(0),
  well = character(0),
  plate_well_id = character(0),
  stringsAsFactors = FALSE
)

# Split files are named <sample><BC1><BC2>.fasta, and the barcode sets disagree
# on which of the two is the plate: B2 is well-then-plate (A01P01) while CBC is
# plate-then-well (P01A01). The previous regex matched only plate-first, so the
# default barcode set silently produced empty plates.
#
# Strip the sample name first -- it can itself contain a plate-like token, e.g.
# RUN10 -- then take the plate from whichever P<n> remains and the well from the
# row letter and column digits beside it.
split_token <- function(token, sample_name = "") {
  base <- sub("\\.fasta$", "", basename(token))

  # Plain prefix strip: the sample name is arbitrary user text and must not be
  # interpreted as a regular expression.
  if (nzchar(sample_name) && startsWith(base, sample_name)) {
    base <- substring(base, nchar(sample_name) + 1L)
  }

  plate_match <- regmatches(base, regexpr("P[0-9]{1,2}", base))
  plate_id <- if (length(plate_match) == 1L && nzchar(plate_match)) {
    sprintf("plate_%d", as.integer(sub("^P", "", plate_match)))
  } else {
    NA_character_
  }

  remainder <- sub("P[0-9]{1,2}", "", base)
  row_letter <- regmatches(remainder, regexpr("[A-Ha-h]", remainder))
  col_digits <- regmatches(remainder, regexpr("[0-9]{1,2}", remainder))

  well <- if (length(row_letter) == 1L && length(col_digits) == 1L) {
    sprintf("%s%d", toupper(row_letter), as.integer(col_digits))
  } else {
    NA_character_
  }

  list(plate_id = plate_id, well = well)
}

parse_count_table <- function(path) {
  df_counts <- read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(df_counts) <- c("barcode", "count", "wells")

  parsed <- lapply(df_counts$wells, split_token, sample_name = sample_name)
  df_counts$well_token <- basename(df_counts$wells)
  df_counts$plate_id <- vapply(parsed, function(x) x$plate_id, character(1))
  df_counts$well <- vapply(parsed, function(x) x$well, character(1))
  df_counts$plate_well_id <- paste0(df_counts$plate_id, "_", df_counts$well)

  df_counts
}

parse_quality_table <- function(path) {
  df_quality <- read.table(path, fill = TRUE, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  parsed_cols <- lapply(colnames(df_quality), split_token, sample_name = sample_name)
  colnames(df_quality) <- vapply(parsed_cols, function(x) {
    if (is.na(x$plate_id) || is.na(x$well)) NA_character_ else paste0(x$plate_id, "_", x$well)
  }, character(1))

  df_quality_pct <- apply(df_quality, MARGIN = 2, FUN = function(x) {
    sum(x[1:3], na.rm = TRUE) / sum(x, na.rm = TRUE)
  })

  keep <- !is.na(names(df_quality_pct))
  df_quality_pct <- df_quality_pct[keep]

  data.frame(
    score = as.numeric(df_quality_pct),
    plate_id = sub("_[A-H][0-9]+$", "", names(df_quality_pct)),
    well = sub("^.*_", "", names(df_quality_pct)),
    plate_well_id = names(df_quality_pct),
    stringsAsFactors = FALSE
  )
}

demux_counts <- parse_count_table(counts_file)
quality_scores <- parse_quality_table(quality_file)

observed <- build_plate_layout(c(demux_counts$plate_id, quality_scores$plate_id))
plate_levels <- observed$levels
plate_layout <- observed$layout
message("Plate QC covering ", length(plate_levels), " plates: ", paste(plate_levels, collapse = ", "))

demux_unmatched <- demux_counts[grep("unmatched", demux_counts$well, ignore.case = TRUE), c("barcode", "count", "well_token")]
demux_plates <- left_join(plate_layout, demux_counts[, c("count", "plate_well_id")], by = "plate_well_id")
quality_plates <- left_join(plate_layout, quality_scores[, c("score", "plate_well_id")], by = "plate_well_id")

demux_plates$plate_id <- factor(demux_plates$plate_id, levels = plate_levels)
quality_plates$plate_id <- factor(quality_plates$plate_id, levels = plate_levels)

write.table(demux_plates, file = file.path(output_dir, "plate_demultiplexed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(quality_plates, file = file.path(output_dir, "plate_quality.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(demux_unmatched, file = file.path(output_dir, "plate_unmatched.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

plot_plate_series <- function(data, value_col, label_col = NULL, title_prefix, file_prefix, colour_scale) {
  # The layout spans a fixed 12 plates, but a run rarely uses them all. Rendering
  # every empty plate produced ~24 blank figures per run and buried the real ones.
  populated <- plate_levels[vapply(plate_levels, function(p) {
    any(!is.na(data[data$plate_id == p, value_col, drop = TRUE]))
  }, logical(1))]

  if (length(populated) == 0) {
    message("No populated plates for ", file_prefix, "; skipping plate figures.")
    return(invisible(NULL))
  }

  pdf(file.path(output_dir, paste0(file_prefix, ".pdf")), width = 8.5, height = 6.5)

  for (current_plate in populated) {
    plate_data <- data[data$plate_id == current_plate, , drop = FALSE]

    # A well is a small circle. A raw score like 0.9704225 overflows it and
    # collides with its neighbour, so shorten what is drawn inside the well.
    if (!is.null(label_col)) {
      plate_data$plate_label <- if (is.numeric(plate_data[[label_col]])) {
        ifelse(
          is.na(plate_data[[label_col]]), NA_character_,
          format(round(plate_data[[label_col]], 2), trim = TRUE, nsmall = 0)
        )
      } else {
        as.character(plate_data[[label_col]])
      }
    }

    if (has_ggplate) {
      if (is.null(label_col)) {
        plate_plot <- ggplate::plate_plot(
          data = plate_data,
          position = well,
          value = !!rlang::sym(value_col),
          plate_size = 96,
          plate_type = "round",
          colour = colour_scale,
          title = paste(title_prefix, current_plate),
          show_legend = TRUE,
          silent = TRUE
        )
      } else {
        plate_plot <- ggplate::plate_plot(
          data = plate_data,
          position = well,
          value = !!rlang::sym(value_col),
          label = plate_label,
          plate_size = 96,
          plate_type = "round",
          colour = colour_scale,
          title = paste(title_prefix, current_plate),
          show_legend = TRUE,
          silent = TRUE
        )
      }
      plate_plot <- plate_plot + gvic_theme_plate
      print(plate_plot)
      ggsave(
        filename = file.path(output_dir, paste0(file_prefix, ".", current_plate, ".png")),
        plot = plate_plot,
        width = 8.5,
        height = 6.5,
        dpi = 180
      )
    } else if (has_platetools) {
      fallback_plot <- platetools::raw_grid(
        data = plate_data[[value_col]],
        well = plate_data$well,
        plate = 96
      ) +
        ggtitle(paste(title_prefix, current_plate)) +
        gvic_theme_plate
      print(fallback_plot)
      ggsave(
        filename = file.path(output_dir, paste0(file_prefix, ".", current_plate, ".png")),
        plot = fallback_plot,
        width = 8.5,
        height = 6.5,
        dpi = 180
      )
    }
  }

  dev.off()
}

if (has_ggplate || has_platetools) {
  demux_palette <- c("#edf4ea", "#b9cfb5", "#8fad8d", "#5f7965", "#35523c")
  quality_palette <- c("#b94b55", "#e6c75b", "#8fad8d", "#35523c")

  plot_plate_series(
    data = demux_plates,
    value_col = "count",
    label_col = "count",
    title_prefix = "Demultiplexed sequences",
    file_prefix = "plate_demultiplexed",
    colour_scale = demux_palette
  )

  plot_plate_series(
    data = quality_plates,
    value_col = "score",
    label_col = "score",
    title_prefix = "Collapsed quality score",
    file_prefix = "plate_quality",
    colour_scale = quality_palette
  )
}
