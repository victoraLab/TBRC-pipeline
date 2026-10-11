# MIBS / sulkseq panels: incidence, not depth.
#
# The shared report answers "did this run work" and the repertoire panels answer
# "what is in it". Neither answers the question a MIBS experiment is run to
# settle, which is whether any clone is expanded enough to be sampled from
# several independent wells. That is a distribution over wells per clone, and
# nothing else in the report shows it.
#
# One panel is deliberately absent. Every other pipeline reports well dominance
# -- the share of a well's reads in its top sequence -- as a doublet check. A
# MIBS well holds a hundred cells and is supposed to be polyclonal, so low
# dominance there is the assay working. Showing it would invite exactly the
# wrong reading.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr)
})

sulk_read <- function(path) {
  if (!file.exists(path)) return(NULL)
  out <- tryCatch(suppressWarnings(readr::read_tsv(path, show_col_types = FALSE,
                                                   progress = FALSE,
                                                   guess_max = 10000)),
                  error = function(e) NULL)
  if (is.null(out) || !nrow(out)) NULL else out
}

sulk_load <- function(qc_dir, sample) {
  f <- function(suffix) file.path(qc_dir, paste0(sample, suffix))
  list(
    clones  = sulk_read(f(".clones.tsv")),
    by_well = sulk_read(f(".clone_by_well.tsv")),
    public  = sulk_read(f(".public_clonotypes.tsv")),
    summary = sulk_read(f(".summary.tsv")),
    sample_stats = sulk_read(f(".sample_stats.tsv"))
  )
}

# metric -> value, for the summary table the report renders as tiles.
sulk_metric <- function(summary, key, default = NA) {
  if (is.null(summary) || !all(c("metric", "value") %in% names(summary))) return(default)
  hit <- summary$value[match(key, summary$metric)]
  if (length(hit) < 1 || is.na(hit)) default else hit
}

# The central figure. A clone present in one well is indistinguishable from a
# PCR artefact, which is why the method discards it; the shape to the right of
# the cut is the expansion the experiment is looking for.
sulk_plot_incidence <- function(clones, min_wells = 2, theme_extra = NULL) {
  if (is.null(clones) || !"wells" %in% names(clones)) return(NULL)
  d <- clones %>%
    mutate(wells = suppressWarnings(as.integer(.data$wells))) %>%
    filter(!is.na(.data$wells)) %>%
    count(.data$wells, name = "clones")
  if (!nrow(d)) return(NULL)
  p <- ggplot(d, aes(x = .data$wells, y = .data$clones)) +
    geom_col(fill = "#2F5FA8", width = 0.8) +
    geom_vline(xintercept = min_wells - 0.5, linetype = "dashed",
               colour = "#A31621", linewidth = 0.4) +
    scale_x_continuous(breaks = function(x) unique(floor(pretty(x)))) +
    scale_y_continuous(trans = "log1p",
                       breaks = c(0, 1, 10, 100, 1000, 10000)) +
    labs(x = "Wells the clone appears in", y = "Clones (log scale)",
         title = "Clone incidence",
         subtitle = sprintf(paste("Everything left of the dashed line is discarded:",
                                  "a clone seen in fewer than %d wells cannot be",
                                  "told from a PCR artefact"), min_wells)) +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p
}

# Ranked clone sizes: how much of the repertoire sits in the few largest clones.
sulk_plot_ranked <- function(clones, top_n = 60, theme_extra = NULL) {
  if (is.null(clones) || !all(c("wells", "verdict") %in% names(clones))) return(NULL)
  d <- clones %>%
    filter(.data$verdict == "kept") %>%
    mutate(wells = suppressWarnings(as.integer(.data$wells))) %>%
    filter(!is.na(.data$wells)) %>%
    arrange(desc(.data$wells))
  if (!nrow(d)) return(NULL)
  d <- head(d, top_n) %>% mutate(rank = seq_len(dplyr::n()))
  has_sample <- "sample" %in% names(d) && length(unique(d$sample)) > 1
  p <- ggplot(d, aes(x = .data$rank, y = .data$wells)) +
    geom_col(fill = "#1B6B4A", width = 0.85) +
    labs(x = "Clone, ranked by incidence", y = "Wells",
         title = sprintf("The %d most widely sampled clones", nrow(d)),
         subtitle = if (has_sample)
           "Across all samples in the run; clones are assigned within a sample"
         else "Within this sample") +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p
}

# Wells recovered per sample. A sample well short of its 16-32 is a
# demultiplexing problem, and incidence is scaled by the number of wells a clone
# could have appeared in -- so this has to be read before the panel above.
sulk_plot_wells_per_sample <- function(by_well, theme_extra = NULL) {
  if (is.null(by_well) || !all(c("sample", "well") %in% names(by_well))) return(NULL)
  d <- by_well %>% distinct(.data$sample, .data$well) %>%
    count(.data$sample, name = "wells") %>% arrange(desc(.data$wells))
  if (nrow(d) < 2) return(NULL)
  d$sample <- factor(d$sample, levels = rev(d$sample))
  p <- ggplot(d, aes(x = .data$sample, y = .data$wells)) +
    geom_col(fill = "#4A525E") + coord_flip() +
    labs(x = NULL, y = "Wells recovered",
         title = "Wells per sample",
         subtitle = "Incidence is only comparable between samples with similar well counts") +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p + theme(axis.text.y = element_text(size = 8))
}

# Where the read filter falls. The published threshold is 6 reads for a clone in
# a well; showing the distribution says whether that cut lands in a gap or
# through the middle of the data.
sulk_plot_reads <- function(by_well, min_reads = 6, theme_extra = NULL) {
  if (is.null(by_well) || !"reads" %in% names(by_well)) return(NULL)
  d <- by_well %>% mutate(reads = suppressWarnings(as.numeric(.data$reads))) %>%
    filter(!is.na(.data$reads), .data$reads > 0)
  if (!nrow(d)) return(NULL)
  p <- ggplot(d, aes(x = .data$reads)) +
    geom_histogram(bins = 60, fill = "#2F5FA8") +
    geom_vline(xintercept = min_reads, linetype = "dashed",
               colour = "#A31621", linewidth = 0.4) +
    scale_x_continuous(trans = "log10") +
    labs(x = "Reads for a clone in one well (log scale)", y = "Clone-well observations",
         title = "Read support per clone and well",
         subtitle = sprintf("Dashed line: the %d-read threshold a well must clear to count", min_reads)) +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p
}

# The private view: how each sample's own repertoire is distributed across
# clone sizes. This is the paper's main private-repertoire figure, and it is not
# the complement of the public one -- public asks which junctions are shared
# between mice, private asks how much of one mouse sits in expanded clones.
#
# Bins are theirs exactly: [1-2) [2-3) [3-4) [4-7) [7-13) [13-33), left-closed,
# on the number of wells a clone occupies. Each sample is normalised to 1 so
# samples with different well counts can be compared, which is the whole point
# of binning rather than plotting raw counts.
SULK_INCIDENCE_BREAKS <- c(1, 2, 3, 4, 7, 13, 33)
SULK_INCIDENCE_TAGS <- c("1", "2", "3", "4-6", "7-12", "13+")

sulk_private_composition <- function(clones) {
  if (is.null(clones) || !all(c("wells", "sample") %in% names(clones))) return(NULL)
  d <- clones %>%
    mutate(wells = suppressWarnings(as.integer(.data$wells))) %>%
    filter(!is.na(.data$wells), .data$wells >= 1)
  if (!nrow(d)) return(NULL)
  d$bin <- cut(d$wells, breaks = SULK_INCIDENCE_BREAKS, include.lowest = TRUE,
               right = FALSE, labels = SULK_INCIDENCE_TAGS)
  d <- d[!is.na(d$bin), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  keep <- c("sample", "bin",
            intersect(c("group", "tissue", "mouse"), names(d)))
  out <- d[, keep, drop = FALSE] %>%
    count(dplyr::across(dplyr::all_of(setdiff(keep, character(0)))), name = "clones") %>%
    group_by(.data$sample) %>%
    mutate(fraction = .data$clones / sum(.data$clones)) %>%
    ungroup()
  out
}

sulk_plot_private <- function(clones, theme_extra = NULL) {
  d <- sulk_private_composition(clones)
  if (is.null(d) || !nrow(d)) return(NULL)
  n_samples <- length(unique(d$sample))
  if (n_samples < 2) return(NULL)

  # Order samples by how much of the repertoire is in the largest bins, so the
  # expanded ones group together rather than scattering alphabetically.
  expanded <- d %>% filter(.data$bin %in% c("4-6", "7-12", "13+")) %>%
    group_by(.data$sample) %>% summarise(x = sum(.data$fraction), .groups = "drop")
  ord <- expanded$sample[order(expanded$x)]
  ord <- c(setdiff(unique(d$sample), ord), ord)
  d$sample <- factor(d$sample, levels = ord)

  p <- ggplot(d, aes(x = .data$sample, y = 100 * .data$fraction, fill = .data$bin)) +
    geom_col(width = 0.85) + coord_flip() +
    scale_fill_brewer(palette = "YlGnBu", direction = 1, name = "Wells") +
    labs(x = NULL, y = "Percent of the sample's clones",
         title = "Clone size within each sample",
         subtitle = paste("Each sample normalised to 100%, so samples with",
                          "different well counts are comparable")) +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p + theme(axis.text.y = element_text(size = 7.5))
}

# The same composition collapsed to the experimental groups, which is the
# comparison the experiment exists to make.
sulk_plot_private_groups <- function(clones, theme_extra = NULL) {
  d <- sulk_private_composition(clones)
  if (is.null(d) || !nrow(d) || !"group" %in% names(d)) return(NULL)
  d <- d[nzchar(as.character(d$group)), , drop = FALSE]
  if (!nrow(d) || length(unique(d$group)) < 2) return(NULL)
  g <- d %>% group_by(.data$group, .data$bin) %>%
    summarise(fraction = mean(.data$fraction), .groups = "drop") %>%
    group_by(.data$group) %>% mutate(fraction = .data$fraction / sum(.data$fraction)) %>%
    ungroup()
  p <- ggplot(g, aes(x = .data$group, y = 100 * .data$fraction, fill = .data$bin)) +
    geom_col(width = 0.7) +
    scale_fill_brewer(palette = "YlGnBu", direction = 1, name = "Wells") +
    labs(x = NULL, y = "Percent of clones (mean across samples)",
         title = "Clone size by group",
         subtitle = "Averaged over the samples in each group, not pooled, so a sample with many clones does not dominate") +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p
}

# The per-sample table. Every figure above answers a question about the run;
# this answers it about each mouse, which is the comparison the plate map was
# attached to make. Ordered by wells recovered rather than by clone count,
# because incidence is scaled by the number of wells a clone could have turned
# up in and a sample short of its 16-32 has to be read first.
#
# It suppresses itself below two rows, the same rule the private panels use. A
# run with no plate map has exactly one sample and a single row would only
# restate the tiles the report already shows above it, under a heading implying
# a comparison that is not there.
#
# mouse and tissue stay in the TSV and out of here: the table has to remain
# readable at forty samples.
sulk_sample_table <- function(sample_stats) {
  if (is.null(sample_stats)) return(NULL)
  # clones_below_min_wells, not clones_single_well: the column counts clones
  # that failed either filter, and the wells filter is a number the submission
  # form now sets. A results directory written before the rename has the old
  # name, and this returns NULL for it, which the report renders as a missing
  # table rather than a wrong one.
  need <- c("sample", "group", "wells", "sequences", "reads", "clones_called",
            "clones_kept", "clones_below_min_wells", "max_wells_per_clone",
            "mean_wells_per_clone", "public_clonotypes")
  if (!all(need %in% names(sample_stats))) return(NULL)
  if (nrow(sample_stats) < 2) return(NULL)
  d <- sample_stats[order(-suppressWarnings(as.numeric(sample_stats$wells)),
                          sample_stats$sample), need, drop = FALSE]
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  rownames(d) <- NULL
  d
}

# Public clonotypes: the same junction in more than one individual. Only
# meaningful when the run holds several mice.
sulk_plot_public <- function(public, theme_extra = NULL) {
  if (is.null(public) || !"individuals" %in% names(public)) return(NULL)
  d <- public %>% mutate(individuals = suppressWarnings(as.integer(.data$individuals))) %>%
    filter(!is.na(.data$individuals)) %>% count(.data$individuals, name = "junctions")
  if (!nrow(d)) return(NULL)
  p <- ggplot(d, aes(x = factor(.data$individuals), y = .data$junctions)) +
    geom_col(fill = "#A85B10") +
    labs(x = "Individuals sharing the junction", y = "Junctions",
         title = "Public clonotypes",
         subtitle = "Exact junction match across mice, among clones that passed both filters") +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p
}
