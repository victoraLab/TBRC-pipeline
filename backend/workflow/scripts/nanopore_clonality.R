#!/usr/bin/env Rscript
# Clonal grouping across wells, plus the final per-well table.
#
# Two sequences belong to the same clone when they share a V gene, a J gene and
# a CDR3 length, and their CDR3s fall within a distance threshold. That grouping
# rule and the three thresholds (75/80/90) are taken from the assay's own
# implementation rather than chosen here, so results stay comparable with
# everything analysed before.
#
# Reported twice, and the distinction matters biologically: within a sample, and
# globally across the run. The same clone appearing twice in one mouse is clonal
# expansion; appearing in two mice is either a public clone or a contamination
# signal, and conflating them would hide both.
#
# Wells that IgBLAST could not resolve are carried through with NA clonality
# rather than dropped, so the final table has one row per well that produced a
# consensus and the count reconciles with the consensus table.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringdist)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: nanopore_clonality.R <airr.tsv> <isotype.tsv> <consensus.tsv> <out.tsv>")
}
airr_path <- args[1]; iso_path <- args[2]; cons_path <- args[3]; out_path <- args[4]

read_safe <- function(p) {
  if (!file.exists(p)) return(NULL)
  out <- tryCatch(suppressWarnings(readr::read_tsv(p, show_col_types = FALSE,
                                                   progress = FALSE)),
                  error = function(e) NULL)
  if (is.null(out) || !nrow(out)) NULL else out
}

airr <- read_safe(airr_path)
iso  <- read_safe(iso_path)
cons <- read_safe(cons_path)

if (is.null(airr)) {
  # No annotation is a legitimate outcome for a failed run. Emit the header so
  # downstream steps see an empty table rather than a missing file.
  writeLines(paste("platewell", "v_gene", "j_gene", "cdr3_nt", "cdr3_aa",
                   "productive", "tail_call", "iga_locus",
                   "Clonality.75", "Clonality.80", "Clonality.90",
                   sep = "\t"), out_path)
  cat("no AIRR annotation; wrote empty table\n"); quit(status = 0)
}

well_of <- function(x) sub("\\|.*$", "", x)

# One gene per call: IgBLAST returns comma-separated alleles by descending score.
# The allele goes too, as in the report and the clonality package: tied alleles
# are listed alphabetically, so the same gene can come back as *01 in one well
# and *02 in another and split one clone (test results/igkj1-2026-10-19/).
first_gene <- function(x) sub("\\*.*$", "", sub(",.*$", "", ifelse(is.na(x), "", x)))

d <- airr %>%
  mutate(platewell = well_of(.data$sequence_id),
         v_gene = first_gene(.data$v_call),
         j_gene = first_gene(.data$j_call),
         cdr3_nt = ifelse(is.na(.data$cdr3), "", .data$cdr3),
         cdr3_aa = ifelse(is.na(.data$cdr3_aa), "", .data$cdr3_aa),
         productive = .data$productive) %>%
  select(any_of(c("platewell", "v_gene", "j_gene", "cdr3_nt", "cdr3_aa", "productive")))

assign_clonality <- function(df, threshold) {
  ans <- rep(NA_character_, nrow(df))
  valid <- !is.na(df$v_gene) & df$v_gene != "" &
           !is.na(df$j_gene) & df$j_gene != "" &
           !is.na(df$cdr3_nt) & df$cdr3_nt != ""
  if (!any(valid)) return(ans)
  composite <- paste(df$v_gene, df$j_gene, nchar(df$cdr3_nt), sep = "|")
  levels_here <- sort(unique(composite[valid]))
  for (gi in seq_along(levels_here)) {
    ix <- which(valid & composite == levels_here[gi])
    s <- df$cdr3_nt[ix]
    if (length(ix) == 1) {
      cl <- 1L
    } else {
      # Hamming is valid only because the group is fixed at one CDR3 length.
      dm <- stringdist::stringdistmatrix(s, method = "hamming")
      L <- nchar(s[1])
      hc <- stats::hclust(stats::as.dist(dm / L), method = "average")
      hc$height <- cummax(hc$height)
      raw <- stats::cutree(hc, h = 1 - threshold)
      # Renumber by the lexicographically smallest CDR3 so cluster ids are
      # stable between runs rather than depending on input order.
      reps <- tapply(s, raw, min)
      stable <- setNames(seq_along(order(reps)), names(sort(reps)))
      cl <- unname(stable[as.character(raw)])
    }
    ans[ix] <- paste0(gi, ".", cl)
  }
  ans
}

d$Clonality.75 <- assign_clonality(d, 0.75)
d$Clonality.80 <- assign_clonality(d, 0.80)
d$Clonality.90 <- assign_clonality(d, 0.90)

if (!is.null(iso)) d <- dplyr::left_join(d, iso, by = "platewell")
if (!is.null(cons)) {
  keep <- intersect(c("platewell", "reads_used", "unique_used",
                      "consensus_length", "status", "dominant_fraction"),
                    names(cons))
  d <- dplyr::left_join(d, cons[, keep, drop = FALSE], by = "platewell")
}

readr::write_tsv(d, out_path)

# CSV alongside the TSV, deliberately. The reference pipeline's final artefact is
# an Excel workbook; a binary that needs Excel to read is a poor archive format
# for results that should still be readable in ten years. CSV opens in Excel for
# anyone who wants it, and in everything else too.
csv_path <- sub("\\.tsv$", ".csv", out_path)
readr::write_csv(d, csv_path)

n_clones <- length(unique(stats::na.omit(d$Clonality.90)))
cat(sprintf("wells: %d   distinct clones at 90%%: %d   productive: %d\n",
            nrow(d), n_clones,
            sum(toupper(substr(as.character(d$productive), 1, 1)) == "T", na.rm = TRUE)))
if ("tail_call" %in% names(d)) {
  tb <- table(d$tail_call, useNA = "no")
  cat("tails:", paste(names(tb), tb, sep = "=", collapse = "  "), "\n")
}
if ("iga_locus" %in% names(d)) {
  tb <- table(d$iga_locus, useNA = "no")
  cat("loci :", paste(names(tb), tb, sep = "=", collapse = "  "), "\n")
}
