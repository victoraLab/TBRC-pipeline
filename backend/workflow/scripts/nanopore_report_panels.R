# Nanopore-specific panels for the run report.
#
# The shared report answers "did this run work" for any platform: wells
# recovered, plate maps, length distributions, trim retention. Nanopore needs
# more, because its failure modes are different in kind rather than degree.
#
# Four things matter here that have no Illumina equivalent:
#
#   read quality      Illumina reads are near-uniform in quality; nanopore reads
#                     are not, and a library's mean-Q distribution decides how
#                     much consensus depth is needed per well.
#   assignment        A read can fail to reach a well in several distinct ways --
#                     no anchor, ambiguous anchor, or a missing plate cassette --
#                     and they mean different things. A library sequenced one
#                     plate per file has NO plate cassette by design, and that
#                     must read as a library mode rather than as loss.
#   dominant fraction What share of a well's reads belong to the winning clone.
#                     This is the doublet and barcode-misassignment signal, and
#                     it exists only because clone identity is now decided by
#                     clustering. Nothing in the Illumina reports corresponds.
#   orientation       Nanopore reads arrive in both orientations. A strong skew
#                     is a library-prep signal worth seeing.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr)
})

np_read_tsv <- function(path) {
  if (!file.exists(path)) return(NULL)
  out <- tryCatch(suppressWarnings(readr::read_tsv(path, show_col_types = FALSE,
                                                   progress = FALSE)),
                  error = function(e) NULL)
  if (is.null(out) || !nrow(out)) NULL else out
}

# Reads may be assigned, unassigned for a stated reason, or unassigned because
# the library has no plate cassette at all. The last is not a loss and must not
# be drawn as one.
np_assignment <- function(summary_tbl) {
  if (is.null(summary_tbl) || !all(c("metric", "value") %in% names(summary_tbl))) return(NULL)
  g <- summary_tbl %>% filter(.data$metric %in% c(
    "raw_reads", "well_assigned", "no_well_anchor", "ambiguous_well",
    "ambiguous_plate", "well_only_missing_plate_cassette", "partial")) %>%
    mutate(value = suppressWarnings(as.numeric(.data$value))) %>%
    group_by(.data$metric) %>% summarise(value = sum(.data$value, na.rm = TRUE), .groups = "drop")
  if (!nrow(g)) return(NULL)
  get1 <- function(k) { v <- g$value[g$metric == k]; if (length(v)) v[1] else 0 }
  raw <- get1("raw_reads")
  well_only <- get1("well_only_missing_plate_cassette")
  list(
    raw = raw,
    assigned = get1("well_assigned"),
    no_anchor = get1("no_well_anchor"),
    ambiguous = get1("ambiguous_well") + get1("ambiguous_plate"),
    well_only = well_only,
    # A library is well-only when every assigned read lacks a plate cassette.
    mode = if (well_only > 0 && well_only >= get1("well_assigned") * 0.99)
             "one plate per file (no plate cassette)" else "plate-multiplexed"
  )
}

np_plot_quality <- function(read_qc, theme_extra = NULL) {
  if (is.null(read_qc) || !"mean_q" %in% names(read_qc)) return(NULL)
  d <- read_qc %>% mutate(mean_q = suppressWarnings(as.numeric(.data$mean_q))) %>%
    filter(!is.na(.data$mean_q))
  if (!nrow(d)) return(NULL)
  med <- stats::median(d$mean_q)
  p <- ggplot(d, aes(x = .data$mean_q)) +
    geom_histogram(bins = 60, fill = "#2F5FA8", colour = NA) +
    geom_vline(xintercept = med, linetype = "dashed", colour = "#A31621") +
    labs(x = "Mean quality per read (Phred)", y = "Reads",
         title = "Read quality",
         subtitle = sprintf("median Q%.1f; consensus depth is what converts this into an accurate sequence", med)) +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p
}

np_plot_lengths <- function(read_qc, theme_extra = NULL) {
  if (is.null(read_qc)) return(NULL)
  col <- intersect(c("insert_length", "read_length"), names(read_qc))
  if (!length(col)) return(NULL)
  d <- read_qc %>% mutate(len = suppressWarnings(as.numeric(.data[[col[1]]]))) %>%
    filter(!is.na(.data$len), .data$len > 0)
  if (!nrow(d)) return(NULL)
  p <- ggplot(d, aes(x = .data$len)) +
    geom_histogram(bins = 60, fill = "#1B6B4A", colour = NA) +
    labs(x = sprintf("%s (nt)", gsub("_", " ", col[1])), y = "Reads",
         title = "Read length",
         subtitle = "A full-length amplicon spans the whole VDJ plus the constant-region tail") +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p
}

# The panel with no Illumina counterpart. A well whose dominant clone accounts
# for most of its reads is a clean single cell; one where it accounts for half
# is a doublet or a barcode collision, and its consensus should be distrusted.
np_plot_dominance <- function(cons_tbl, theme_extra = NULL) {
  if (is.null(cons_tbl) || !"dominant_fraction" %in% names(cons_tbl)) return(NULL)
  d <- cons_tbl %>%
    mutate(dominant_fraction = suppressWarnings(as.numeric(.data$dominant_fraction))) %>%
    filter(!is.na(.data$dominant_fraction), .data$consensus_length > 0)
  if (!nrow(d)) return(NULL)
  flagged <- sum(d$dominant_fraction < 0.7)
  p <- ggplot(d, aes(x = .data$dominant_fraction)) +
    geom_histogram(bins = 40, fill = "#7A5AA8", colour = NA) +
    geom_vline(xintercept = 0.7, linetype = "dashed", colour = "#A31621") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Share of the well's reads in the dominant clone", y = "Wells",
         title = "Clone dominance per well",
         subtitle = sprintf("%d well%s below 70%%: likely doublets or barcode collisions, not sequencing error",
                            flagged, ifelse(flagged == 1, "", "s"))) +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p
}

np_plot_depth <- function(cons_tbl, theme_extra = NULL) {
  if (is.null(cons_tbl) || !"reads_used" %in% names(cons_tbl)) return(NULL)
  d <- cons_tbl %>% mutate(reads_used = suppressWarnings(as.numeric(.data$reads_used))) %>%
    filter(!is.na(.data$reads_used), .data$reads_used > 0)
  if (!nrow(d)) return(NULL)
  p <- ggplot(d, aes(x = .data$reads_used)) +
    geom_histogram(bins = 50, fill = "#A85B10", colour = NA) +
    scale_x_log10(labels = scales::comma) +
    labs(x = "Reads entering the consensus (log scale)", y = "Wells",
         title = "Consensus depth per well",
         subtitle = "Nanopore accuracy comes from depth: error averages out, real mutations do not") +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p
}

np_plot_orientation <- function(read_qc, theme_extra = NULL) {
  if (is.null(read_qc) || !"orientation" %in% names(read_qc)) return(NULL)
  d <- read_qc %>% filter(!is.na(.data$orientation), .data$orientation != "") %>%
    count(.data$orientation, name = "reads")
  if (nrow(d) < 2) return(NULL)
  names(d)[1] <- "orientation"
  p <- ggplot(d, aes(x = .data$orientation, y = .data$reads)) +
    geom_col(fill = "#42566B", width = 0.6) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Orientation as sequenced", y = "Reads",
         title = "Read orientation",
         subtitle = "Both strands are expected; a strong skew points at library prep") +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p
}
