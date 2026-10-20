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

# Library structure (nanopore_library_qc.py). A run with few wells can be a run
# with too few reads, or a run whose PCR made mostly primer dimers -- they look
# the same in every other panel.
NP_LIB_LABELS <- c(
  full_length = "Full length", short_product = "Short product", primer_dimer = "Primer dimer",
  concatemer = "Concatemer", mixed_strands = "Mixed strands", no_3p_cassette = "No 3' cassette",
  no_5p_handle = "No 5' handle", unrecognised = "Unrecognised")
NP_LIB_COLOURS <- c(
  full_length = "#1B6B4A", short_product = "#E0A030", primer_dimer = "#A31621",
  concatemer = "#7A5AA8", mixed_strands = "#B07AA1", no_3p_cassette = "#8DA0B3",
  no_5p_handle = "#B8C4D0", unrecognised = "#DADFE4")

# Headline numbers for the report text, or NULL when the table is absent.
np_library_summary <- function(lib_tbl) {
  if (is.null(lib_tbl) || !all(c("scope", "class", "reads") %in% names(lib_tbl))) return(NULL)
  a <- lib_tbl %>% filter(.data$scope == "all") %>%
    mutate(reads = suppressWarnings(as.numeric(.data$reads)))
  get1 <- function(k) { v <- a$reads[a$class == k]; if (length(v)) v[1] else 0 }
  products <- get1("full_length") + get1("short_product") + get1("primer_dimer")
  list(reads = sum(a$reads, na.rm = TRUE), products = products,
       full = get1("full_length"), short = get1("short_product"), dimer = get1("primer_dimer"),
       concatemer = get1("concatemer"),
       dimer_share = if (products > 0) get1("primer_dimer") / products else NA_real_,
       full_share = if (products > 0) get1("full_length") / products else NA_real_)
}

np_plot_library <- function(lib_tbl, theme_extra = NULL) {
  if (is.null(lib_tbl) || !all(c("scope", "class", "reads") %in% names(lib_tbl))) return(NULL)
  d <- lib_tbl %>% mutate(reads = suppressWarnings(as.numeric(.data$reads))) %>%
    filter(!is.na(.data$reads), .data$class %in% names(NP_LIB_LABELS)) %>%
    group_by(.data$scope) %>% mutate(share = .data$reads / sum(.data$reads)) %>% ungroup() %>%
    mutate(class = factor(.data$class, levels = rev(names(NP_LIB_LABELS))),
           scope = factor(ifelse(.data$scope == "all", "All reads", .data$scope),
                          levels = c("All reads", setdiff(unique(.data$scope), "all"))))
  if (!nrow(d)) return(NULL)
  s <- np_library_summary(lib_tbl)
  p <- ggplot(d, aes(x = .data$share, y = .data$scope, fill = .data$class)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = NP_LIB_COLOURS, labels = NP_LIB_LABELS, breaks = names(NP_LIB_LABELS),
                      name = NULL) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    labs(x = "Share of reads", y = NULL, title = "Library structure",
         subtitle = if (!is.null(s) && is.finite(s$dimer_share))
           sprintf("Of %s reads with both primers: %.0f%% full length, %.0f%% primer dimer",
                   format(s$products, big.mark = ","), 100 * s$full_share, 100 * s$dimer_share)
         else "Reads classified by the primers they carry") +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p + theme(legend.position = "bottom")
}

np_plot_inserts <- function(ins_tbl, theme_extra = NULL, dimer_max = 100, min_insert = 250) {
  if (is.null(ins_tbl) || !all(c("insert_bin_start", "reads") %in% names(ins_tbl))) return(NULL)
  # The run's own thresholds, written beside the histogram, win over the defaults.
  if ("dimer_max" %in% names(ins_tbl)) dimer_max <- suppressWarnings(as.numeric(ins_tbl$dimer_max[1]))
  if ("min_insert" %in% names(ins_tbl)) min_insert <- suppressWarnings(as.numeric(ins_tbl$min_insert[1]))
  d <- ins_tbl %>% mutate(x = suppressWarnings(as.numeric(.data$insert_bin_start)),
                          reads = suppressWarnings(as.numeric(.data$reads))) %>%
    filter(!is.na(.data$x), !is.na(.data$reads)) %>%
    mutate(class = ifelse(.data$x < dimer_max, "primer_dimer",
                   ifelse(.data$x < min_insert, "short_product", "full_length")))
  if (!nrow(d)) return(NULL)
  p <- ggplot(d, aes(x = .data$x, y = .data$reads, fill = .data$class)) +
    geom_col(width = 10, just = 0) +
    geom_vline(xintercept = c(dimer_max, min_insert), linetype = "dashed", colour = "grey30") +
    scale_fill_manual(values = NP_LIB_COLOURS, labels = NP_LIB_LABELS, name = NULL) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Insert between the primers (nt; last bin is 1000+)", y = "Reads",
         title = "Insert length",
         subtitle = sprintf("Under %d nt: primer dimer; %d-%d nt: short product; the pipeline uses %d nt and over",
                            dimer_max, dimer_max, min_insert - 1, min_insert)) +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p + theme(legend.position = "bottom")
}
