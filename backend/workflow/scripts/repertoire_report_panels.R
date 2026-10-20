# Repertoire panels: what the run actually found, not just whether it ran.
#
# The rest of the report answers "did this work" -- wells recovered, plate maps,
# read depth. These answer "what is in it": which V genes were used, how mutated
# the sequences are, how long the CDR3s are, and what fraction is productive.
#
# One convention is enforced here rather than left to whoever reads the table
# afterwards: only the first allele of a call is kept. IgBLAST returns
# comma-separated alleles by descending score, so counting the raw string
# treats "IGHV1-9*01" and "IGHV1-9*01,IGHV1-9*02" as different genes.
#
# The MiniMouse IGHV3-6CHI1-64 -> IGHV3-6 collapse is deliberately NOT applied.
# It exists because the 10x cellranger reference and bulk annotation name that
# gene differently; TBRC annotates with IgBLAST against one IMGT set, where it
# has one name, so the rename never fired (0 rows on all four RUN15 tables) while
# the plot and its caption said it had. It belongs in the 10x analysis code.
#
# Somatic hypermutation is taken as 100 - v_identity, which is the V-region
# divergence from germline that IgBLAST reports directly. On nanopore that number
# carries sequencing error as well as biology, and the report says so rather than
# presenting it as a clean measurement.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr)
})

rep_read_airr <- function(path) {
  if (!file.exists(path)) return(NULL)
  out <- tryCatch(suppressWarnings(readr::read_tsv(path, show_col_types = FALSE,
                                                   progress = FALSE,
                                                   guess_max = 5000)),
                  error = function(e) NULL)
  if (is.null(out) || !nrow(out)) NULL else out
}

# One gene per call: first allele, allele suffix dropped.
rep_gene <- function(x) {
  g <- sub(",.*$", "", ifelse(is.na(x), "", as.character(x)))
  sub("\\*.*$", "", g)
}

rep_truthy <- function(x) toupper(substr(as.character(x), 1, 1)) == "T"

rep_summary <- function(airr) {
  if (is.null(airr)) return(NULL)
  n <- nrow(airr)
  pick <- function(col) if (col %in% names(airr)) airr[[col]] else rep(NA, n)
  shm <- suppressWarnings(as.numeric(pick("v_identity")))
  # IgBLAST reports v_identity as a percentage in outfmt 19, but some versions
  # emit a fraction. Normalise rather than silently reporting 99% mutation.
  if (length(shm[!is.na(shm)]) && stats::median(shm, na.rm = TRUE) <= 1) shm <- shm * 100
  list(
    n = n,
    productive = sum(rep_truthy(pick("productive")), na.rm = TRUE),
    in_frame   = sum(rep_truthy(pick("vj_in_frame")), na.rm = TRUE),
    stop_codon = sum(rep_truthy(pick("stop_codon")), na.rm = TRUE),
    cdr3       = sum(nzchar(ifelse(is.na(pick("cdr3_aa")), "", pick("cdr3_aa")))),
    v_genes    = length(unique(rep_gene(pick("v_call"))[nzchar(rep_gene(pick("v_call")))])),
    shm_median = if (any(!is.na(shm))) stats::median(100 - shm, na.rm = TRUE) else NA_real_
  )
}

rep_plot_vgene <- function(airr, top_n = 20, theme_extra = NULL) {
  if (is.null(airr) || !"v_call" %in% names(airr)) return(NULL)
  d <- data.frame(gene = rep_gene(airr$v_call), stringsAsFactors = FALSE) %>%
    filter(nzchar(.data$gene)) %>% count(.data$gene, name = "n") %>%
    arrange(desc(.data$n))
  if (!nrow(d)) return(NULL)
  names(d)[1] <- "gene"
  total <- sum(d$n)
  d <- head(d, top_n) %>% mutate(pct = 100 * .data$n / total,
                                 gene = factor(.data$gene, levels = rev(.data$gene)))
  p <- ggplot(d, aes(x = .data$gene, y = .data$pct)) +
    geom_col(fill = "#2F5FA8") + coord_flip() +
    labs(x = NULL, y = "Percent of sequences",
         title = sprintf("V gene usage (top %d of %d genes)", nrow(d), total_genes <- length(unique(rep_gene(airr$v_call)))),
         subtitle = "First allele of each call") +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p + theme(axis.text.y = element_text(size = 9))
}

rep_plot_jgene <- function(airr, theme_extra = NULL) {
  if (is.null(airr) || !"j_call" %in% names(airr)) return(NULL)
  d <- data.frame(gene = rep_gene(airr$j_call), stringsAsFactors = FALSE) %>%
    filter(nzchar(.data$gene)) %>% count(.data$gene, name = "n")
  if (!nrow(d)) return(NULL)
  names(d)[1] <- "gene"
  d$pct <- 100 * d$n / sum(d$n)
  p <- ggplot(d, aes(x = stats::reorder(.data$gene, -.data$pct), y = .data$pct)) +
    geom_col(fill = "#1B6B4A") +
    labs(x = NULL, y = "Percent of sequences", title = "J gene usage") +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

rep_plot_shm <- function(airr, nanopore = FALSE, theme_extra = NULL) {
  if (is.null(airr) || !"v_identity" %in% names(airr)) return(NULL)
  v <- suppressWarnings(as.numeric(airr$v_identity))
  if (all(is.na(v))) return(NULL)
  if (stats::median(v, na.rm = TRUE) <= 1) v <- v * 100
  d <- data.frame(shm = 100 - v) %>% filter(!is.na(.data$shm), .data$shm >= 0, .data$shm < 40)
  if (!nrow(d)) return(NULL)
  med <- stats::median(d$shm)
  sub <- if (nanopore)
    sprintf("median %.1f%%. On nanopore this includes residual sequencing error, not biology alone", med)
  else
    sprintf("median %.1f%% divergence from germline V", med)
  p <- ggplot(d, aes(x = .data$shm)) +
    geom_histogram(bins = 50, fill = "#7A5AA8", colour = NA) +
    geom_vline(xintercept = med, linetype = "dashed", colour = "#A31621") +
    labs(x = "V-region divergence from germline (%)", y = "Sequences",
         title = "Somatic hypermutation", subtitle = sub) +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p
}

rep_plot_cdr3 <- function(airr, theme_extra = NULL) {
  if (is.null(airr) || !"cdr3_aa" %in% names(airr)) return(NULL)
  d <- data.frame(len = nchar(ifelse(is.na(airr$cdr3_aa), "", airr$cdr3_aa))) %>%
    filter(.data$len > 0, .data$len < 60)
  if (!nrow(d)) return(NULL)
  p <- ggplot(d, aes(x = .data$len)) +
    geom_histogram(binwidth = 1, fill = "#A85B10", colour = NA) +
    labs(x = "CDR3 length (amino acids)", y = "Sequences",
         title = "CDR3 length distribution",
         subtitle = sprintf("median %d aa", stats::median(d$len))) +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p
}

# Clone sizes say whether the run caught expansion or a flat repertoire.
rep_plot_clones <- function(final_tbl, col = "Clonality.90", theme_extra = NULL) {
  if (is.null(final_tbl) || !col %in% names(final_tbl)) return(NULL)
  d <- final_tbl %>% filter(!is.na(.data[[col]])) %>% count(.data[[col]], name = "size")
  if (!nrow(d)) return(NULL)
  tab <- as.data.frame(table(d$size), stringsAsFactors = FALSE)
  names(tab) <- c("size", "clones")
  tab$size <- as.integer(tab$size)
  expanded <- sum(d$size > 1)
  p <- ggplot(tab, aes(x = .data$size, y = .data$clones)) +
    geom_col(fill = "#0E7C7B") +
    scale_x_continuous(breaks = function(x) unique(floor(pretty(x)))) +
    labs(x = "Cells per clone", y = "Clones",
         title = "Clone size distribution (90% identity)",
         subtitle = sprintf("%d of %d clones contain more than one cell", expanded, nrow(d))) +
    theme_bw()
  if (!is.null(theme_extra)) p <- p + theme_extra
  p
}
