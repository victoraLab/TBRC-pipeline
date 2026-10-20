#!/usr/bin/env Rscript
#
# One report per run, shipped inside the archive.
#
# Until now a run delivered about twenty-five loose files -- PDFs, PNGs and TSVs
# scattered across qc/ and qc/plates/ -- and nothing you could open to answer
# "did this run work". The numbers that would answer it sat in TSVs nobody
# opens. This builds a single self-contained HTML plus a multi-page PDF from
# exactly those files, so the archive carries its own verdict.
#
# Deliberately dependency-free beyond what the cluster already has: ggplot2,
# dplyr, tidyr, readr, ggplate, scales. rmarkdown, knitr and base64enc are NOT
# installed there, so the HTML is assembled by hand and figures are inlined as
# SVG from the cairo device rather than base64 PNGs. Inline SVG is smaller,
# stays crisp at any zoom, and needs no encoder.
#
# Every input is optional. Older runs predate some of the QC tables, and a run
# that fails early may write none of them; each section renders what it has and
# says plainly when something is missing rather than aborting the report.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = "") {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1L]]
}

SAMPLE     <- get_arg("-s", "run")
RESULT_DIR <- get_arg("-d", ".")
FINAL_FASTA<- get_arg("-f", "")
OUT_PREFIX <- get_arg("-o", file.path(RESULT_DIR, SAMPLE))
PIPELINE   <- get_arg("-p", "unknown")
CHAIN      <- get_arg("-c", "")
PRIMER     <- get_arg("-b", "")
USER_ID    <- get_arg("-u", "")
SCRIPT_DIR <- get_arg("-S", ".")
NANOPORE   <- tolower(get_arg("-N", "")) %in% c("1", "true", "yes", "nanopore")

theme_file <- file.path(SCRIPT_DIR, "lab_theme.R")
if (file.exists(theme_file)) source(theme_file) else {
  gvic_theme <- theme(); gvic_theme_plate <- theme(); vicLabTheme <- theme()
}

QC     <- file.path(RESULT_DIR, "qc")
PLATES <- file.path(QC, "plates")
has_ggplate <- requireNamespace("ggplate", quietly = TRUE)

# ---------------------------------------------------------------- input helpers

read_tsv_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  out <- tryCatch(suppressWarnings(readr::read_tsv(path, show_col_types = FALSE,
                                                   progress = FALSE)),
                  error = function(e) NULL)
  if (is.null(out) || !nrow(out)) NULL else out
}

# Base R rather than Biostrings, which is not installed on the cluster.
fasta_lengths <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(integer(0))
  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
  if (!length(lines)) return(integer(0))
  is_hdr <- startsWith(lines, ">")
  if (!any(is_hdr)) return(integer(0))
  grp <- cumsum(is_hdr)
  nchar_seq <- nchar(lines)
  nchar_seq[is_hdr] <- 0L
  as.integer(tapply(nchar_seq, grp, sum))
}

fmt_int <- function(x) if (length(x) && is.finite(x)) formatC(x, big.mark = ",", format = "d") else "--"
fmt_pct <- function(x) if (length(x) && is.finite(x)) sprintf("%.1f%%", x) else "--"
esc <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

# ---------------------------------------------------------------- load the run

depth  <- read_tsv_safe(file.path(PLATES, "plate_read_depth.tsv"))
demux  <- read_tsv_safe(file.path(PLATES, "plate_demultiplexed.tsv"))
qual   <- read_tsv_safe(file.path(PLATES, "plate_quality.tsv"))
unmat  <- read_tsv_safe(file.path(PLATES, "plate_unmatched.tsv"))
len_pre  <- read_tsv_safe(file.path(QC, "pre_trim.contig_lengths.tsv"))
len_post <- read_tsv_safe(file.path(QC, "post_trim.contig_lengths.tsv"))
if (is.null(len_post)) len_post <- read_tsv_safe(file.path(QC, "final_output.contig_lengths.tsv"))

seq_len_vec <- fasta_lengths(FINAL_FASTA)

# Nanopore runs carry QC that has no Illumina counterpart, so the panels live in
# their own file and are folded in only when the run was nanopore. The shared
# sections above are unchanged: a nanopore run still gets plate maps, depth and
# length distributions drawn the same way as any other.
# Repertoire annotation, when the run produced any. Both pipelines can: the
# nanopore one always annotates its consensus, and an Illumina run does when
# IgBLAST was requested. The section is drawn from whichever AIRR table exists
# rather than from the pipeline name, so it appears whenever there is something
# to say and stays silent when there is not.
REP <- list(airr = NULL, final = NULL, summary = NULL)
rep_panels <- file.path(SCRIPT_DIR, "repertoire_report_panels.R")
if (file.exists(rep_panels)) {
  source(rep_panels)
  # Both pipelines that annotate internally write somewhere of their own, and a
  # sulkseq run was rendering with no Repertoire section at all despite carrying
  # a perfectly good AIRR table -- the report only knew about igblast/.
  airr_candidates <- c(
    file.path(QC, "nanopore", paste0(SAMPLE, ".long_consensus.airr.tsv")),
    file.path(QC, "sulkseq", paste0(SAMPLE, ".airr.tsv")),
    file.path(RESULT_DIR, "igblast", paste0(SAMPLE, ".igblast.tsv")),
    file.path(RESULT_DIR, "igblast", paste0(SAMPLE, ".airr.tsv"))
  )
  for (cand in airr_candidates) {
    if (is.null(REP$airr) && file.exists(cand)) REP$airr <- rep_read_airr(cand)
  }
  REP$final <- rep_read_airr(file.path(QC, "nanopore", paste0(SAMPLE, ".per_well_final.tsv")))
  REP$summary <- rep_summary(REP$airr)
}

# What happened to annotation, including nothing. A run with no igblast/ folder
# used to render exactly like one that had been annotated: plate maps, yield and
# read depth, and no hint that anything was missing. There was no way to tell
# whether it was skipped by choice, skipped for a missing database, or broken.
ANN <- NULL
local({
  tb <- read_tsv_safe(file.path(QC, "annotation", "status.tsv"))
  if (!is.null(tb) && all(c("metric", "value") %in% names(tb))) {
    v <- as.character(tb$value); names(v) <- as.character(tb$metric)
    ANN <<- v
  }
})

# MIBS / sulkseq. Detected by the tables existing rather than by the pipeline
# name, the same way the repertoire section is, so a rerun of an older result
# directory still renders what it has.
SK <- list(clones = NULL, by_well = NULL, public = NULL, summary = NULL)
local({
  panels <- file.path(SCRIPT_DIR, "sulkseq_report_panels.R")
  qc_dir <- file.path(QC, "sulkseq")
  if (file.exists(panels) && dir.exists(qc_dir)) {
    source(panels)
    SK <<- sulk_load(qc_dir, SAMPLE)
  }
})
IS_SULK <- !is.null(SK$clones) || !is.null(SK$summary)

# Receptor-type check. The panel a run was annotated with is a choice made on a
# form; this is the only place the report can say whether the sequences agreed
# with it. A run that predates the check writes nothing here and the section
# stays silent.
RC <- list(fields = NULL, comp = NULL)
local({
  chk <- read_tsv_safe(file.path(QC, "receptor", "locus_check.tsv"))
  if (!is.null(chk) && all(c("metric", "value") %in% names(chk))) {
    v <- as.character(chk$value); names(v) <- as.character(chk$metric)
    RC$fields <<- v
  }
  RC$comp <<- read_tsv_safe(file.path(QC, "receptor", "locus_composition.tsv"))
})

NP <- list(read_qc = NULL, summary = NULL, cons = NULL, assign = NULL,
           lib = NULL, lib_inserts = NULL, lib_sum = NULL)
if (NANOPORE) {
  panels <- file.path(SCRIPT_DIR, "nanopore_report_panels.R")
  if (file.exists(panels)) {
    source(panels)
    NP$read_qc <- np_read_tsv(file.path(RESULT_DIR, "nanopore",
                                        paste0(SAMPLE, ".platewell_read_qc.tsv")))
    NP$summary <- np_read_tsv(file.path(QC, "nanopore",
                                        paste0(SAMPLE, ".demux_summary.tsv")))
    NP$cons    <- np_read_tsv(file.path(QC, "nanopore",
                                        paste0(SAMPLE, ".consensus_per_well.tsv")))
    NP$assign  <- np_assignment(NP$summary)
    NP$lib     <- np_read_tsv(file.path(QC, "nanopore", paste0(SAMPLE, ".library_qc.tsv")))
    NP$lib_inserts <- np_read_tsv(file.path(QC, "nanopore",
                                            paste0(SAMPLE, ".library_qc_inserts.tsv")))
    if (exists("np_library_summary")) NP$lib_sum <- np_library_summary(NP$lib)
  }
}

MIN_DEPTH_FLAG <- 10L   # matches the "flag, do not drop" rule in pipelines.yaml

# ---------------------------------------------------------------- statistics

stats <- list(
  records      = length(seq_len_vec),
  len_median   = if (length(seq_len_vec)) stats::median(seq_len_vec) else NA_real_,
  len_min      = if (length(seq_len_vec)) min(seq_len_vec) else NA_real_,
  len_max      = if (length(seq_len_vec)) max(seq_len_vec) else NA_real_,
  wells        = NA_integer_, plates = NA_integer_,
  depth_median = NA_real_, depth_total = NA_real_,
  low_wells    = NA_integer_, empty_wells = NA_integer_
)

if (!is.null(depth) && all(c("plate", "ggplate_well", "read_depth") %in% names(depth))) {
  depth <- depth %>% mutate(read_depth = suppressWarnings(as.numeric(read_depth)))
  occupied <- depth %>% filter(!is.na(read_depth), read_depth > 0)
  stats$wells        <- nrow(occupied)
  stats$plates       <- dplyr::n_distinct(depth$plate)
  stats$depth_median <- if (nrow(occupied)) stats::median(occupied$read_depth) else NA_real_
  stats$depth_total  <- sum(depth$read_depth, na.rm = TRUE)
  stats$low_wells    <- sum(occupied$read_depth < MIN_DEPTH_FLAG, na.rm = TRUE)
  stats$empty_wells  <- nrow(depth) - nrow(occupied)
}

# Retention through trimming, the number that says whether the trim sequence was
# right. A collapse here is the signature of a wrong handle -- the failure mode
# that made luka_light discard 88% of a library without ever reporting it.
retention <- NA_real_
if (!is.null(len_pre) && !is.null(len_post) &&
    "sequence_count" %in% names(len_pre) && "sequence_count" %in% names(len_post)) {
  a <- suppressWarnings(as.numeric(len_pre$sequence_count[1]))
  b <- suppressWarnings(as.numeric(len_post$sequence_count[1]))
  if (is.finite(a) && a > 0 && is.finite(b)) retention <- 100 * b / a
}

# The count column has been spelled differently by different QC scripts over
# time -- plate_stats.R writes "count", earlier tables used "read_count". Take
# whichever is present rather than assuming one.
count_col <- function(tb) {
  if (is.null(tb)) return(NULL)
  hit <- intersect(c("read_count", "count", "reads", "n"), names(tb))
  if (length(hit)) hit[1] else NULL
}

# A match rate needs BOTH sides. plate_unmatched.tsv is written by the classic
# split_fasta.sh route; the vsearch route demultiplexes in demux_wells.sh and
# does not record unmatched barcodes at all, so that table arrives header-only.
# Treating an absent denominator as zero would turn "not measured" into a
# flattering "100% matched" -- the report must say it does not know instead.
matched_pct <- NA_real_
matched_known <- FALSE
dc <- count_col(demux); uc <- count_col(unmat)
if (!is.null(dc) && !is.null(uc) && nrow(unmat) > 0) {
  m <- sum(suppressWarnings(as.numeric(demux[[dc]])), na.rm = TRUE)
  u <- sum(suppressWarnings(as.numeric(unmat[[uc]])), na.rm = TRUE)
  if (is.finite(m) && is.finite(u) && (m + u) > 0) {
    matched_pct <- 100 * m / (m + u); matched_known <- TRUE
  }
}
# Reads that did land in a well is a fact we always have, and is worth showing
# even when the match rate is not measurable.
demux_reads <- if (!is.null(dc)) sum(suppressWarnings(as.numeric(demux[[dc]])), na.rm = TRUE) else NA_real_

# ---------------------------------------------------------------- figures

PLOTS <- list()   # name -> list(plot=, w=, h=, caption=)
add_plot <- function(key, p, w, h, caption) {
  if (!is.null(p)) PLOTS[[key]] <<- list(plot = p, w = w, h = h, caption = caption)
}

# Plate maps: the graphic people actually look at. ggplate draws real 96-well
# geometry, so a bad column or a dead corner is visible as a shape rather than
# something you have to infer from a table of well ids.
plate_plots <- list()
if (!is.null(depth) && has_ggplate && all(c("plate", "ggplate_well", "read_depth") %in% names(depth))) {
  for (pl in sort(unique(depth$plate))) {
    d <- depth %>% filter(plate == !!pl) %>%
         mutate(read_depth = suppressWarnings(as.numeric(read_depth)))
    if (!nrow(d)) next
    ttl <- if ("plate_label" %in% names(d) && !all(is.na(d$plate_label)))
             as.character(d$plate_label[1]) else paste("Plate", pl)
    p <- tryCatch(
      ggplate::plate_plot(data = d, position = ggplate_well, value = read_depth,
                          plate_size = 96, plate_type = "round",
                          colour = c("#F2F4F7", "#9EC5E8", "#2F5FA8", "#12305C"),
                          title = ttl, show_legend = TRUE, silent = TRUE) + gvic_theme_plate,
      error = function(e) NULL)
    if (!is.null(p)) plate_plots[[as.character(pl)]] <- p
  }
}

if (!is.null(depth) && "read_depth" %in% names(depth)) {
  d <- depth %>% mutate(read_depth = suppressWarnings(as.numeric(read_depth))) %>%
       filter(!is.na(read_depth), read_depth > 0)
  if (nrow(d)) {
    p <- ggplot(d, aes(x = read_depth)) +
      geom_histogram(bins = 50, fill = "#2F5FA8", colour = NA) +
      geom_vline(xintercept = MIN_DEPTH_FLAG, linetype = "dashed", colour = "#A31621") +
      scale_x_log10(labels = scales::comma) +
      labs(x = "Reads per well (log scale)", y = "Wells",
           title = "Read depth per well") +
      theme_bw() + gvic_theme
    add_plot("depth_hist", p, 9, 4.6,
      sprintf("Dashed line marks the %d-read flag threshold. Wells below it are reported, never dropped.", MIN_DEPTH_FLAG))
  }
}

if (length(seq_len_vec)) {
  p <- ggplot(data.frame(len = seq_len_vec), aes(x = len)) +
    geom_histogram(bins = 60, fill = "#1B6B4A", colour = NA) +
    labs(x = "Delivered sequence length (nt)", y = "Sequences",
         title = "Length distribution of delivered sequences") +
    theme_bw() + gvic_theme
  add_plot("len_hist", p, 9, 4.6,
    "Read off the expected amplicon size. A second mode well below it usually means adapter or primer-dimer carry-through.")
}

# Per-plate yield: which plates carried the run and which contributed nothing.
plate_summary <- NULL
if (!is.null(depth) && all(c("plate", "read_depth") %in% names(depth))) {
  plate_summary <- depth %>%
    mutate(read_depth = suppressWarnings(as.numeric(read_depth))) %>%
    group_by(plate) %>%
    summarise(
      wells_used   = sum(!is.na(read_depth) & read_depth > 0),
      wells_low    = sum(!is.na(read_depth) & read_depth > 0 & read_depth < MIN_DEPTH_FLAG),
      median_depth = suppressWarnings(stats::median(read_depth[!is.na(read_depth) & read_depth > 0])),
      total_reads  = sum(read_depth, na.rm = TRUE),
      .groups = "drop"
    ) %>% arrange(plate)

  if (nrow(plate_summary) > 1) {
    p <- ggplot(plate_summary, aes(x = factor(plate), y = wells_used)) +
      geom_col(fill = "#2F5FA8") +
      geom_col(aes(y = wells_low), fill = "#A85B10") +
      labs(x = "Plate", y = "Wells recovered",
           title = "Wells recovered per plate",
           subtitle = "Amber shows the portion below the flag threshold") +
      theme_bw() + gvic_theme +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10))
    add_plot("plate_bar", p, 9, 4.6,
      "A plate contributing far fewer wells than its neighbours usually points at a pipetting or barcode problem confined to that plate.")
  }
}

if (NANOPORE && exists("np_plot_quality")) {
  add_plot("np_quality", np_plot_quality(NP$read_qc, gvic_theme), 9, 4.4,
    "Nanopore read quality is broad where Illumina is narrow. What makes a base call trustworthy here is depth, not any single read.")
  add_plot("np_length", np_plot_lengths(NP$read_qc, gvic_theme), 9, 4.4,
    "A truncated mode well below the amplicon size means reads that will not span the region of interest.")
  add_plot("np_dominance", np_plot_dominance(NP$cons, gvic_theme), 9, 4.4,
    "The doublet signal. A well whose dominant clone holds only half its reads probably received two cells, and its consensus is a blend of two antibodies rather than a noisy copy of one.")
  add_plot("np_depth", np_plot_depth(NP$cons, gvic_theme), 9, 4.4,
    "How many reads each consensus was built from. Sequencing error averages away with depth; real somatic hypermutation does not.")
  add_plot("np_orient", np_plot_orientation(NP$read_qc, gvic_theme), 7, 3.8,
    "Reads arrive on both strands and are oriented during demultiplexing.")
  if (exists("np_plot_library")) {
    add_plot("np_library", np_plot_library(NP$lib, gvic_theme), 9, 3.8,
      "Every read, classed by the primers it carries. Grey classes are reads too broken to place; the coloured ones are what the PCR made.")
    add_plot("np_inserts", np_plot_inserts(NP$lib_inserts, gvic_theme), 9, 4.4,
      "Bases between the 5' handle and the reverse primer. A peak under the first line is primer dimer: no template between the primers.")
  }
}

if (!is.null(REP$airr) && exists("rep_plot_vgene")) {
  add_plot("rep_v", rep_plot_vgene(REP$airr, 20, gvic_theme), 9, 6.2,
    "Which V genes the run used, counting the first allele of each call as its gene.")
  add_plot("rep_j", rep_plot_jgene(REP$airr, gvic_theme), 8, 4.2,
    "J usage is narrower than V and mostly reflects the locus rather than selection.")
  add_plot("rep_shm", rep_plot_shm(REP$airr, NANOPORE, gvic_theme), 9, 4.4,
    "Divergence from the germline V. This is the number to treat most carefully on nanopore, because sequencing error and somatic hypermutation land in the same column.")
  add_plot("rep_cdr3", rep_plot_cdr3(REP$airr, gvic_theme), 9, 4.2,
    "CDR3 length is set at recombination and is a useful sanity check: a shifted distribution usually means mis-assembled junctions rather than unusual biology.")
  add_plot("rep_clone", rep_plot_clones(REP$final, "Clonality.90", gvic_theme), 9, 4.2,
    "Clone sizes say whether this run caught expansion or a flat repertoire.")
}

if (IS_SULK && exists("sulk_plot_incidence")) {
  min_wells <- suppressWarnings(as.integer(sulk_metric(SK$summary, "min_wells_per_clone", 2)))
  min_reads <- suppressWarnings(as.integer(sulk_metric(SK$summary, "min_reads_per_well", 6)))
  if (!is.finite(min_wells)) min_wells <- 2L
  if (!is.finite(min_reads)) min_reads <- 6L
  add_plot("sk_incidence", sulk_plot_incidence(SK$clones, min_wells, gvic_theme), 9, 4.4,
           paste("Wells per clone. This is the measurement: a clone is large because it",
                 "was sampled from many independent hundred-cell wells, not because it",
                 "has many reads."))
  add_plot("sk_ranked", sulk_plot_ranked(SK$clones, 60, gvic_theme), 9, 4.2,
           "The clones that survived both filters, ordered by how many wells they reached.")
  add_plot("sk_wells", sulk_plot_wells_per_sample(SK$by_well, gvic_theme), 9, 5.0,
           paste("Wells recovered per sample. Read this before the incidence panel:",
                 "a clone cannot appear in more wells than the sample has."))
  add_plot("sk_reads", sulk_plot_reads(SK$by_well, min_reads, gvic_theme), 9, 4.2,
           "Read support behind each clone-in-well observation, and where the threshold falls.")
  add_plot("sk_public", sulk_plot_public(SK$public, gvic_theme), 8, 4.0,
           "Junctions found in more than one individual.")
  # The private view. Both suppress themselves when the run carries no
  # experimental design, rather than drawing one bar and calling it a
  # comparison.
  add_plot("sk_private", sulk_plot_private(SK$clones, gvic_theme), 9, 5.4,
           paste("How much of each sample sits in expanded clones. Normalised per",
                 "sample, so samples with different well counts are comparable."))
  add_plot("sk_private_groups", sulk_plot_private_groups(SK$clones, gvic_theme), 8, 4.2,
           paste("The same composition by experimental group, averaged over samples",
                 "rather than pooled so one large sample cannot dominate."))
}

# ---------------------------------------------------------------- rendering

# Cairo SVG, inlined. Strips the XML prolog and DOCTYPE so the fragment can sit
# directly inside the HTML body; the <svg> element itself carries the geometry.
# Cairo does not emit <text>. It defines each glyph once as
# <symbol id="glyph0-1"> and draws it with <use xlink:href="#glyph0-1">, and it
# restarts that numbering at zero for every figure. Inlining several SVGs into
# one document therefore collides every id, and each figure after the first
# renders its letters from the FIRST figure's glyph table -- which looks like a
# broken font rather than like the id clash it actually is. Clip paths and
# embedded images collide the same way.
#
# So every id is namespaced per figure before inlining, along with every
# reference to it. The trailing quote and the closing paren in the patterns
# below matter: without them "glyph0-1" would also rewrite "glyph0-10".
.fig_counter <- 0L

namespace_svg_ids <- function(s, tag) {
  ids <- regmatches(s, gregexpr('id="[^"]+"', s))[[1]]
  ids <- unique(sub('"$', "", sub('^id="', "", ids)))
  for (id in ids) {
    new_id <- paste0(tag, "-", id)
    s <- gsub(paste0('id="', id, '"'),    paste0('id="', new_id, '"'),    s, fixed = TRUE)
    s <- gsub(paste0('#', id, '"'),       paste0('#', new_id, '"'),       s, fixed = TRUE)
    s <- gsub(paste0('url(#', id, ')'),   paste0('url(#', new_id, ')'),   s, fixed = TRUE)
  }
  s
}

svg_inline <- function(p, w, h) {
  f <- tempfile(fileext = ".svg")
  ok <- tryCatch({ svg(f, width = w, height = h, bg = "white"); print(p); dev.off(); TRUE },
                 error = function(e) { try(dev.off(), silent = TRUE); FALSE })
  if (!ok || !file.exists(f)) return("<p class='missing'>Figure could not be rendered.</p>")
  s <- paste(readLines(f, warn = FALSE), collapse = "\n")
  unlink(f)
  s <- sub("^<\\?xml[^>]*\\?>", "", s)
  s <- gsub("<!DOCTYPE[^>]*>", "", s)
  .fig_counter <<- .fig_counter + 1L
  s <- namespace_svg_ids(s, paste0("f", .fig_counter))
  paste0("<div class='fig'>", s, "</div>")
}

# One measured value, as a row in the results table.
#
# This used to be a grid of large-number cards. That reads as a dashboard, and a
# dashboard is the wrong genre: nobody watches a sequencing run in real time.
# What people do is read a result sheet -- the thing a core facility hands back,
# where the numbers line up in a column so the eye can run down them and the
# labels sit quietly to the left. Right-aligned tabular figures make that
# comparison possible; forty-point numerals in boxes do not.
tile <- function(value, label, tone = "") {
  mark <- switch(tone,
                 ok   = "<span class='mk mk-ok'>&#9679;</span>",
                 warn = "<span class='mk mk-warn'>&#9679;</span>",
                 crit = "<span class='mk mk-crit'>&#9679;</span>",
                 "<span class='mk'></span>")
  sprintf("<tr><th scope='row'>%s</th><td class='v'>%s</td><td class='s'>%s</td></tr>",
          label, value, mark)
}

section <- function(title, body, note = NULL) {
  paste0("<section><h2>", esc(title), "</h2>",
         if (!is.null(note)) paste0("<p class='note'>", note, "</p>") else "",
         body, "</section>")
}

miss <- function(what)
  paste0("<p class='missing'>", esc(what),
         " was not produced for this run, so this section is empty. ",
         "Older runs predate some QC tables.</p>")

# ---------------------------------------------------------------- HTML

css <- "
:root{--ink:#16191F;--ink2:#4A525E;--ink3:#79818E;--paper:#fff;--paper2:#F4F6F8;
--rule:#DDE2E8;--accent:#2F5FA8;--ok:#1B6B4A;--warn:#A85B10;--crit:#A31621;}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){
--ink:#E9EDF2;--ink2:#AEB8C4;--ink3:#7C8794;--paper:#12151A;--paper2:#181C23;
--rule:#2A313B;--accent:#7FA9E8;--ok:#6FC49C;--warn:#E8A45C;--crit:#F3868F;}}
*{box-sizing:border-box}
body{margin:0;background:var(--paper);color:var(--ink);
font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif;
font-size:15px;line-height:1.6;-webkit-font-smoothing:antialiased}
.wrap{max-width:1000px;margin:0 auto;padding:0 22px 80px}
header{border-bottom:2px solid var(--ink);padding:40px 0 18px;margin-bottom:28px}
h1{font-size:2rem;margin:0 0 8px;letter-spacing:-.02em}
.meta{display:flex;flex-wrap:wrap;gap:7px;margin-top:14px}
.meta span{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px;
background:var(--paper2);border:1px solid var(--rule);padding:3px 8px;border-radius:5px;color:var(--ink2)}
/* The results table. Hairline rules and alignment do the work that boxes and
   drop shadows were doing; the numbers sit in one column so they can be read
   down rather than hunted for, which is how a facility result sheet works. */
table.summary{width:100%;border-collapse:collapse;margin:4px 0 30px;
font-variant-numeric:tabular-nums}
table.summary th[scope=row]{text-align:left;font-weight:400;color:var(--ink2);
padding:7px 12px 7px 0;border-bottom:1px solid var(--rule);font-size:14.5px;width:62%}
table.summary td.v{text-align:right;padding:7px 10px 7px 0;
border-bottom:1px solid var(--rule);font-size:15px;font-weight:600;
font-family:ui-monospace,SFMono-Regular,Menlo,monospace;white-space:nowrap}
table.summary td.s{width:16px;padding:7px 0;border-bottom:1px solid var(--rule)}
table.summary tr:last-child th,table.summary tr:last-child td{border-bottom:none}
.mk{display:inline-block;width:9px;font-size:10px;line-height:1}
.mk-ok{color:var(--ok)}.mk-warn{color:var(--warn)}.mk-crit{color:var(--crit)}
h2{font-size:1.12rem;margin:40px 0 4px;padding-bottom:6px;border-bottom:2px solid var(--ink);letter-spacing:-.005em}
.note{color:var(--ink2);margin:12px 0 18px;max-width:70ch}
/* One callout style, used only when something needs to be said before the
   numbers are read. Left rule rather than a filled card: it has to interrupt
   without looking like a product notification. */
.alert{border-left:3px solid var(--rule);background:var(--paper2);
padding:14px 16px;margin:0 0 26px;font-size:14.5px}
.alert b{display:block;font-size:11px;text-transform:uppercase;letter-spacing:.07em;
margin-bottom:5px;font-weight:600}
.alert-crit{border-left-color:var(--crit)}.alert-crit b{color:var(--crit)}
.alert-warn{border-left-color:var(--warn)}.alert-warn b{color:var(--warn)}
.alert-ok{border-left-color:var(--ok)}.alert-ok b{color:var(--ok)}
.missing{color:var(--ink3);font-style:italic;background:var(--paper2);
border:1px dashed var(--rule);border-radius:6px;padding:12px 14px}
.fig{margin:14px 0 2px;overflow-x:auto}
.fig svg{max-width:100%;height:auto;display:block}
.cap{color:var(--ink3);font-size:13px;margin:-4px 0 20px;max-width:72ch}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:14px}
table{border-collapse:collapse;width:100%;font-size:13.5px;margin:14px 0}
th,td{border-bottom:1px solid var(--rule);padding:7px 9px;text-align:right;
font-variant-numeric:tabular-nums}
th:first-child,td:first-child{text-align:left}
th{font-size:11px;text-transform:uppercase;letter-spacing:.07em;color:var(--ink3);font-weight:600}
tr.flag td{color:var(--warn)}
footer{margin-top:52px;padding-top:20px;border-top:2px solid var(--ink);
color:var(--ink3);font-size:13px}
@media print{body{font-size:11px}.wrap{max-width:none}header{padding-top:0}
h2{page-break-after:avoid}.fig{page-break-inside:avoid}}
"

verdict_tone <- function(v, good, warn) {
  if (!is.finite(v)) "" else if (v >= good) "ok" else if (v >= warn) "warn" else "crit"
}

tiles <- paste0(
  tile(fmt_int(stats$records), "Sequences delivered", if (stats$records > 0) "ok" else "crit"),
  tile(fmt_int(stats$wells),   "Wells recovered"),
  tile(fmt_int(stats$plates),  "Plates"),
  tile(if (is.finite(stats$depth_median)) fmt_int(stats$depth_median) else "--", "Median reads / well"),
  tile(fmt_pct(retention), "Kept through trim", verdict_tone(retention, 60, 30)),
  if (matched_known)
    tile(fmt_pct(matched_pct), "Barcode matched", verdict_tone(matched_pct, 70, 40))
  else
    tile(fmt_int(demux_reads), "Reads assigned to wells"),
  tile(fmt_int(stats$low_wells), "Low-depth wells flagged",
       if (is.finite(stats$low_wells) && stats$low_wells > 0) "warn" else "")
)

# The receptor verdict, rendered ahead of everything else when it is not a
# plain agreement. A corrected panel changes what every downstream number in
# this report means, so it cannot sit at the bottom.
receptor_banner <- ""
receptor_section <- ""
if (!is.null(RC$fields)) {
  f <- RC$fields
  g <- function(k, d = "") {
    if (!(k %in% names(f))) return(d)
    v <- f[[k]]
    if (!length(v) || is.na(v[1]) || !nzchar(v[1])) d else v[1]
  }
  verdict <- g("verdict", "inconclusive")
  detail  <- g("detail", "")
  tones <- c(corrected = "crit", contradicted = "crit", widened = "crit",
             narrow_panel = "crit", mixed = "warn",
             inconclusive = "warn", failed = "warn", agree = "ok", skipped = "")
  titles <- c(corrected    = "Receptor type corrected",
              contradicted = "Receptor type contradicts the submission",
              widened      = "Panel widened: the plate carries more loci than submitted",
              narrow_panel = "The submitted panel misses loci that are on the plate",
              mixed        = "Both receptor classes present",
              inconclusive = "Receptor type could not be determined",
              failed       = "Receptor check did not run",
              agree        = "Receptor type confirmed",
              skipped      = "Receptor check disabled")
  tone <- if (verdict %in% names(tones)) tones[[verdict]] else "warn"
  ttl  <- if (verdict %in% names(titles)) titles[[verdict]] else "Receptor type"

  if (verdict != "agree" && verdict != "skipped" && nzchar(ttl)) {
    receptor_banner <- paste0(
      "<div class='alert alert-", tone, "'><b>", esc(ttl), "</b>",
      esc(detail), "</div>")
  }

  # The species, judged by the same probe. A wrong species is the one failure
  # worse than a wrong panel: every record is "annotated", at 75% identity,
  # to genes the animal does not have. A plate whose loci need different sets
  # -- a human V knock-in on mouse D/J beside mouse kappa -- is annotated per
  # locus and said out loud; one that could not be planned is kept as declared.
  sp_verdict <- g("species_verdict", "")
  sp_tones  <- c(corrected = "crit", contradicted = "crit", mixed = "warn",
                 per_locus = "warn",
                 unchecked = "", agree = "", skipped = "", failed = "warn")
  sp_titles <- c(corrected    = "Species corrected",
                 contradicted = "Species contradicts the submission",
                 mixed        = "More than one species on the plate",
                 per_locus    = "Loci annotated against different germline sets",
                 failed       = "Species check did not run")
  if (sp_verdict %in% names(sp_titles)) {
    receptor_banner <- paste0(receptor_banner,
      "<div class='alert alert-", sp_tones[[sp_verdict]], "'><b>",
      esc(sp_titles[[sp_verdict]]), "</b>", esc(g("species_detail", "")), "</div>")
  }

  comp_table <- ""
  if (!is.null(RC$comp) && nrow(RC$comp)) {
    has_species <- "species" %in% names(RC$comp)
    rows <- paste0(apply(RC$comp, 1, function(r) {
      paste0("<tr><td>", esc(r[["locus"]]), "</td>",
             if (has_species) paste0("<td>", esc(r[["species"]]), "</td>") else "",
             "<td>", fmt_int(as.integer(r[["sequences"]])), "</td><td>",
             fmt_pct(100 * as.numeric(r[["fraction"]])), "</td></tr>")
    }), collapse = "")
    comp_table <- paste0(
      "<table><thead><tr><th>Locus</th>",
      if (has_species) "<th>Species</th>" else "",
      "<th>Sequences</th><th>Share</th>",
      "</tr></thead><tbody>", rows, "</tbody></table>")
  }

  receptor_section <- section("Receptor type",
    paste0("<table class='summary'>",
           tile(esc(g("declared_panel", "--")), "Panel submitted"),
           tile(esc(g("applied_panel", "--")), "Panel used",
                if (identical(g("applied_panel"), g("declared_panel"))) "" else "warn"),
           tile(esc(g("observed_class", "--")), "Receptor class in the data"),
           tile(esc(g("species_declared", "--")), "Species submitted"),
           tile(esc(g("applied_species", g("species_declared", "--"))), "Species used",
                if (identical(g("applied_species"), g("species_declared")) ||
                    !nzchar(g("applied_species"))) "" else "warn"),
           tile(esc(if (nzchar(sp_verdict)) sp_verdict else "--"), "Species verdict",
                if (sp_verdict %in% names(sp_tones)) sp_tones[[sp_verdict]] else ""),
           tile(fmt_int(suppressWarnings(as.integer(g("sequences_tested", "0")))),
                "Sequences probed"),
           tile(fmt_int(suppressWarnings(as.integer(g("aligned", "0")))),
                "Aligned to a germline"),
           tile(esc(verdict), "Verdict", tone),
           "</table>", comp_table),
    paste0("A deterministic subsample of the delivered sequences was aligned against ",
           "the immunoglobulin and the T-cell receptor germline databases of every ",
           "species built on the server (", esc(g("species_probed", "--")), "), and each ",
           "sequence counted towards whichever scored highest. The point ",
           "is not to second-guess the submission but to catch the one failure that looks ",
           "like success: annotating a T-cell library against immunoglobulin germlines ",
           "produces an empty table and a clonality run that finds nothing, both exiting ",
           "cleanly."))
}

fig_block <- function(key) {
  it <- PLOTS[[key]]
  if (is.null(it)) return("")
  paste0(svg_inline(it$plot, it$w, it$h), "<p class='cap'>", esc(it$caption), "</p>")
}

# Said plainly, above the figures, because it changes what the rest of the
# report can be expected to contain.
annotation_note <- ""
if (!is.null(ANN)) {
  a <- function(k, d = "") {
    if (!(k %in% names(ANN))) return(d)
    v <- ANN[[k]]
    if (!length(v) || is.na(v[1]) || !nzchar(v[1])) d else v[1]
  }
  ig <- a("igblast_state", "not_requested")
  cl <- a("clonality_state", "not_requested")
  n_ann <- suppressWarnings(as.integer(a("annotated_sequences", "0")))

  if (!(ig %in% c("completed", "completed_by_pipeline"))) {
    tone <- switch(ig, failed = "crit", skipped = "warn", "")
    title <- switch(ig,
      not_requested = "This run was not annotated",
      skipped       = "Annotation was requested but did not run",
      failed        = "Annotation failed",
      "Annotation did not run")
    body <- if (ig == "not_requested")
      paste0("IgBLAST and clonality were not selected on the submission form, so this ",
             "archive has no <code>igblast/</code> or <code>clonality/</code> folder and ",
             "no repertoire panels below. Everything above is assembly quality only. ",
             "Tick <em>Run IgBLAST</em> when submitting to get V and J usage, CDR3 ",
             "lengths, mutation load and clone calls.")
    else paste0(esc(a("igblast_reason", "No reason was recorded.")),
                " Nothing else in this report is affected; the assembly ran normally.")
    annotation_note <- paste0("<div class='alert", if (nzchar(tone)) paste0(" alert-", tone) else "",
                              "'><b>", esc(title), "</b>", body, "</div>")
  } else if (cl != "completed" && cl != "not_requested") {
    annotation_note <- paste0(
      "<div class='alert alert-warn'><b>Clones were not called</b>",
      esc(a("clonality_reason", "")),
      " The IgBLAST annotation itself completed: ", fmt_int(n_ann), " sequences.</div>")
  }

  # Annotation that completed against the wrong germline set. Every record has
  # a V call, so nothing upstream notices; the median identity does.
  if (identical(a("igblast_identity_verdict", ""), "low")) {
    annotation_note <- paste0(annotation_note,
      "<div class='alert alert-crit'><b>The gene calls in this archive are not trustworthy</b>",
      "Median V identity is ", esc(a("igblast_median_v_identity", "?")), "%, below the ",
      esc(a("igblast_identity_floor", "85")), "% floor. The germline set used (",
      esc(a("igblast_species", "?")), " / ", esc(a("igblast_panel", "?")),
      ") does not match these sequences: wrong species, wrong strain, or a locus the ",
      "panel does not carry. V and J usage, mutation load and clone calls below all ",
      "inherit this. See the receptor and species check above.</div>")
  }

  # Heavy chains through the light-chain protocol. Its fixed 23-nt cut after
  # the handle is the kappa primer, but on a heavy read it is template (FR1 nt
  # 24-46), and nothing else in the run says so. A floor, so a handful of
  # forced calls on a light-only plate does not raise it.
  n_heavy <- suppressWarnings(as.integer(a("luka_light_heavy_records", "")))
  if (isTRUE(n_heavy > 0 && (n_heavy >= 20 || n_heavy >= 0.05 * n_ann))) {
    annotation_note <- paste0(annotation_note,
      "<div class='alert alert-warn'><b>Heavy chains on a light-chain (luka_light) run</b>",
      fmt_int(n_heavy), " of ", fmt_int(n_ann), " annotated records are heavy chains. ",
      "The light-chain protocol removes 23 nt after the 5' handle: the kappa primer on a ",
      "kappa read, but real sequence on a heavy read, so these records lack FR1 positions ",
      "24-46. Their FR1 and V identity are short by those bases. For a plate that carries ",
      "both chains, submit with the standard protocol and the ig_all panel.</div>")
  }
}

# MIBS summary and section. Built after fig_block, which it uses to embed
# the incidence figures.
sulk_tiles <- ""
sulk_section <- ""
if (IS_SULK) {
  m <- function(k, d = "--") {
    v <- sulk_metric(SK$summary, k, NA)
    if (is.na(v) || !nzchar(as.character(v))) d else as.character(v)
  }
  int <- function(k) {
    v <- suppressWarnings(as.integer(sulk_metric(SK$summary, k, NA)))
    if (is.finite(v)) fmt_int(v) else "--"
  }
  # The first of these metrics that a results directory actually has. A summary
  # written before a metric existed is not a summary with a missing value in it,
  # and an older report has to keep showing the number it showed before.
  int_any <- function(...) {
    for (k in c(...)) {
      v <- suppressWarnings(as.integer(sulk_metric(SK$summary, k, NA)))
      if (is.finite(v)) return(fmt_int(v))
    }
    "--"
  }
  unmapped_wells <- suppressWarnings(as.integer(
    sulk_metric(SK$summary, "wells_not_in_platemap", NA)))
  kept <- suppressWarnings(as.integer(sulk_metric(SK$summary, "clones_kept", NA)))
  called <- suppressWarnings(as.integer(sulk_metric(SK$summary, "clones_called", NA)))
  kept_pct <- if (is.finite(kept) && is.finite(called) && called > 0)
    100 * kept / called else NA_real_

  sulk_tiles <- paste0(
    "<table class='summary'>",
    # samples_mapped, not samples: '(unmapped wells)' is a diagnostic bucket
    # rather than a mouse, and counting it here over-reports the run by one
    # whenever a single well falls outside the design.
    tile(int_any("samples_mapped", "samples"), "Samples in the run"),
    tile(int("wells"), "Wells recovered"),
    # A well the design does not name is a well whose sequences are attributed
    # to nobody. It was reported nowhere, and the run looked complete.
    if (is.finite(unmapped_wells))
      tile(fmt_int(unmapped_wells), "Wells not in the plate map",
           if (unmapped_wells > 0) "crit" else "")
    else "",
    tile(m("wells_per_sample"), "Wells per sample"),
    tile(int("clones_called"), "Clones called"),
    tile(int("clones_kept"), "Clones passing both filters",
         if (is.finite(kept) && kept > 0) "ok" else "crit"),
    tile(fmt_pct(kept_pct), "Kept of called"),
    # Labelled for the filter rather than for its default: the submission form
    # sets the wells threshold, and at 3 this counts two-well clones too.
    tile(int("clones_single_well"), "Dropped below the wells filter"),
    tile(int("max_wells_per_clone"), "Widest clone (wells)"),
    tile(m("mean_wells_per_clone"), "Mean wells per kept clone"),
    tile(int("public_clonotypes"), "Public clonotypes"),
    "</table>")

  platemap_note <- if (identical(m("platemap_used"), "no"))
    paste0("<div class='alert alert-warn'><b>No plate map</b>",
           "Every well was treated as belonging to one sample, because no ",
           "platemap.tsv was found beside the FASTQ files. Incidence only means ",
           "something within a sample: if this library holds more than one mouse, ",
           "the clone counts above merge them. Add a platemap.tsv with columns ",
           "<code>well</code> and <code>sample</code> and rerun.</div>")
  else ""

  # Say what is missing rather than leaving a gap where the comparison should be.
  design_note <- if (identical(m("platemap_used"), "no")) paste0(
    "<p class='note'>The per-sample and per-group figures are not shown because ",
    "this run carries no experimental design. TBRC is given a folder of FASTQ ",
    "files and a barcode set; it can see which wells exist but not which mouse, ",
    "tissue or group each belongs to. A template listing every well found is in ",
    "<code>qc/sulkseq/", esc(SAMPLE), ".platemap_template.tsv</code> &mdash; fill ",
    "in sample, mouse, tissue and group, save it as <code>platemap.tsv</code> ",
    "beside the FASTQ files and rerun. Annotation, the slow step, is not ",
    "repeated.</p>") else ""

  # The per-sample table. Everything above it is a number for the run, and a
  # MIBS run is several mice: a run-level clone count averages over designs that
  # need not match, since one mouse with sixteen wells and another with
  # thirty-two do not contribute to incidence comparably. This is the comparison
  # the plate map is attached in order to get.
  #
  # sulk_sample_table() returns NULL for a single-sample run, where the row
  # would only restate the tiles above it, and for a results directory written
  # before the file existed -- which is the only case worth announcing.
  sample_table <- {
    d <- sulk_sample_table(SK$sample_stats)
    if (!is.null(d)) {
      rows <- paste0(apply(d, 1, function(r) {
        sprintf(paste0("<tr class='%s'><td>%s</td><td>%s</td><td>%s</td><td>%s</td>",
                       "<td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td>",
                       "<td>%s</td><td>%s</td></tr>"),
                if (identical(r[["sample"]], "(unmapped wells)")) "flag" else "",
                esc(r[["sample"]]), esc(r[["group"]]),
                fmt_int(as.numeric(r[["wells"]])),
                fmt_int(as.numeric(r[["sequences"]])),
                fmt_int(as.numeric(r[["reads"]])),
                fmt_int(as.numeric(r[["clones_called"]])),
                fmt_int(as.numeric(r[["clones_kept"]])),
                fmt_int(as.numeric(r[["clones_below_min_wells"]])),
                fmt_int(as.numeric(r[["max_wells_per_clone"]])),
                esc(r[["mean_wells_per_clone"]]),
                fmt_int(as.numeric(r[["public_clonotypes"]])))
      }), collapse = "")
      paste0("<table><thead><tr><th>Sample</th><th>Group</th><th>Wells</th>",
             "<th>Sequences</th><th>Reads</th><th>Clones called</th>",
             "<th>Clones kept</th><th>Below the wells filter</th>",
             "<th>Widest clone</th><th>Mean wells/clone</th><th>Public</th>",
             "</tr></thead><tbody>", rows, "</tbody></table>",
             "<p class='cap'>Per sample, ordered by wells recovered, because ",
             "incidence is scaled by the number of wells a clone could have ",
             "turned up in and a sample short of its 16 to 32 has to be read ",
             "first. A row flagged red is the bucket of wells the plate map does ",
             "not name; it is not a mouse.</p>")
    } else if (is.null(SK$sample_stats)) miss("Per-sample statistics") else ""
  }

  sulk_section <- paste0(platemap_note, section("Clones by incidence",
    paste0(sulk_tiles, sample_table,
           fig_block("sk_incidence"), fig_block("sk_ranked"),
           fig_block("sk_private"), fig_block("sk_private_groups"), design_note,
           fig_block("sk_wells"), fig_block("sk_reads"), fig_block("sk_public")),
    paste0("This library was processed as MIBS, where clone size is the number of ",
           "wells a clone turns up in rather than its read depth. Each well holds ",
           "roughly a hundred cells, so a clone reaching several wells was sampled ",
           "several times from an expanded population. A clone confined to one well ",
           "cannot be distinguished from a PCR artefact and is dropped -- that is the ",
           "method's error control, standing in for the depth thresholds used ",
           "elsewhere. Clones are assigned across wells within a sample: same V gene, ",
           "same J gene, junction within ", m("junction_max_distance"),
           " nucleotides.")))
}


plate_block <- if (length(plate_plots)) {
  paste0("<div class='grid'>",
         paste0(vapply(plate_plots, function(p) svg_inline(p, 5.6, 3.9), character(1)),
                collapse = ""),
         "</div>")
} else if (!has_ggplate) {
  miss("ggplate is not installed on this node, so plate maps")
} else miss("plate_read_depth.tsv")

# Wells judged per locus (well_verdicts_by_locus.py): on a mixed plate a well
# with one cell's heavy and kappa chains is two monoclonal loci, not a doublet,
# which is what the per-well verdict calls it. Rows are loci, columns verdicts.
locus_verdicts <- read_tsv_safe(file.path(QC, "assembly", "well_verdicts_by_locus.tsv"))
locus_verdict_table <- if (!is.null(locus_verdicts) && nrow(locus_verdicts)) {
  verdict_order <- c("monoclonal", "multiclonal", "low_depth", "unannotated")
  tab <- table(factor(locus_verdicts$locus),
               factor(locus_verdicts$verdict, levels = verdict_order))
  rows <- vapply(rownames(tab), function(locus) {
    # fmt_int formats one number. Handed the row's four counts, its `&&` only
    # warns on R 4.2 (the cluster) but stops the whole report from R 4.3 on.
    paste0("<tr><td>", esc(locus), "</td>",
           paste0("<td>", vapply(as.integer(tab[locus, ]), fmt_int, character(1)), "</td>",
                  collapse = ""), "</tr>")
  }, character(1))
  paste0("<table><thead><tr><th>Locus</th><th>Monoclonal</th><th>Multiclonal</th>",
         "<th>Low depth</th><th>Unannotated</th></tr></thead><tbody>",
         paste0(rows, collapse = ""), "</tbody></table>")
} else ""

plate_table <- if (!is.null(plate_summary)) {
  rows <- apply(plate_summary, 1, function(r) {
    flag <- suppressWarnings(as.numeric(r[["wells_low"]])) > 0
    sprintf("<tr class='%s'><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
            if (isTRUE(flag)) "flag" else "",
            esc(r[["plate"]]), fmt_int(as.numeric(r[["wells_used"]])),
            fmt_int(as.numeric(r[["wells_low"]])),
            fmt_int(as.numeric(r[["median_depth"]])),
            fmt_int(as.numeric(r[["total_reads"]])))
  })
  paste0("<table><thead><tr><th>Plate</th><th>Wells recovered</th><th>Flagged low</th>",
         "<th>Median depth</th><th>Total reads</th></tr></thead><tbody>",
         paste0(rows, collapse = ""), "</tbody></table>")
} else miss("Per-plate summary")

len_table <- {
  mk <- function(lbl, tb) {
    if (is.null(tb)) return("")
    sprintf("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
            esc(lbl), fmt_int(as.numeric(tb$sequence_count[1])),
            fmt_int(as.numeric(tb$min_length[1])),
            fmt_int(as.numeric(tb$median_length[1])),
            fmt_int(as.numeric(tb$max_length[1])))
  }
  body <- paste0(mk("Before trimming", len_pre), mk("After trimming", len_post))
  if (nzchar(body))
    paste0("<table><thead><tr><th>Stage</th><th>Sequences</th><th>Min</th>",
           "<th>Median</th><th>Max</th></tr></thead><tbody>", body, "</tbody></table>")
  else miss("Contig length tables")
}

html <- paste0(
  "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width,initial-scale=1'>",
  "<title>", esc(SAMPLE), " &mdash; TBRC run report</title><style>", css, "</style></head><body><div class='wrap'>",
  "<header><h1>", esc(SAMPLE), "</h1>",
  "<div style='color:var(--ink2)'>TBRC assembly run report</div>",
  "<div class='meta'>",
  paste0("<span>", esc(format(Sys.time(), "%Y-%m-%d %H:%M")), "</span>"),
  if (nzchar(USER_ID))  paste0("<span>user: ", esc(USER_ID), "</span>") else "",
  paste0("<span>pipeline: ", esc(PIPELINE), "</span>"),
  if (nzchar(CHAIN))  paste0("<span>chain: ", esc(CHAIN), "</span>") else "",
  if (nzchar(PRIMER)) paste0("<span>barcodes: ", esc(PRIMER), "</span>") else "",
  "</div></header>",

  receptor_banner,
  "<table class='summary'>", tiles, "</table>",

  if (NANOPORE && !is.null(NP$assign)) paste0(
    "<table class='summary'>",
    tile(fmt_int(NP$assign$raw), "Raw reads"),
    tile(fmt_pct(if (NP$assign$raw > 0) 100 * NP$assign$assigned / NP$assign$raw else NA),
         "Assigned to a well",
         verdict_tone(if (NP$assign$raw > 0) 100 * NP$assign$assigned / NP$assign$raw else NA, 50, 30)),
    tile(fmt_int(NP$assign$no_anchor), "No well anchor"),
    tile(fmt_int(NP$assign$ambiguous), "Ambiguous barcode"),
    "</table>",
    "<p class='note'>Library mode: <strong>", esc(NP$assign$mode), "</strong>. ",
    if (grepl("one plate", NP$assign$mode))
      paste0("Reads in this library carry a well anchor and no plate cassette. That is the design, ",
             "not a loss: identity is the well alone, and no read was discarded for it.")
    else
      "Reads carry both a plate cassette and a well anchor, so identity is plate plus well.",
    "</p>") else "",

  if (!is.null(REP$summary)) paste0(
    "<table class='summary'>",
    tile(fmt_pct(100 * REP$summary$productive / max(1, REP$summary$n)), "Productive",
         verdict_tone(100 * REP$summary$productive / max(1, REP$summary$n), 80, 60)),
    tile(fmt_pct(100 * REP$summary$in_frame / max(1, REP$summary$n)), "In frame"),
    tile(fmt_pct(100 * REP$summary$stop_codon / max(1, REP$summary$n)), "Stop codon",
         if (100 * REP$summary$stop_codon / max(1, REP$summary$n) > 10) "warn" else ""),
    tile(fmt_pct(100 * REP$summary$cdr3 / max(1, REP$summary$n)), "CDR3 recovered"),
    tile(fmt_int(REP$summary$v_genes), "Distinct V genes"),
    tile(if (is.finite(REP$summary$shm_median)) sprintf("%.1f%%", REP$summary$shm_median) else "--",
         "Median SHM"),
    if (NANOPORE && !is.null(NP$cons) && "dominant_fraction" %in% names(NP$cons)) {
      dd <- suppressWarnings(as.numeric(NP$cons$dominant_fraction))
      dd <- dd[!is.na(dd) & NP$cons$consensus_length > 0]
      pct <- if (length(dd)) 100 * sum(dd < 0.70) / length(dd) else NA_real_
      tile(fmt_pct(pct), "Probable doublets",
           if (is.finite(pct) && pct > 20) "warn" else "")
    } else "",
    if (NANOPORE && !is.null(NP$lib_sum) && is.finite(NP$lib_sum$dimer_share)) {
      tile(fmt_pct(100 * NP$lib_sum$dimer_share), "Primer dimers (of PCR products)",
           if (NP$lib_sum$dimer_share > 0.2) "warn" else "")
    } else "",
    "</table>") else "",

  annotation_note,
  sulk_section,
  receptor_section,

  if (!is.null(REP$airr)) section("Repertoire",
          paste0(fig_block("rep_v"), fig_block("rep_j"),
                 fig_block("rep_cdr3"), fig_block("rep_shm"),
                 fig_block("rep_clone")),
          paste0("What the run found, rather than whether it ran. Gene usage, junction ",
                 "length, mutation load and clonal structure, read from the same IgBLAST ",
                 "annotation that produced the productive rate above.")) else "",

  if (nzchar(locus_verdict_table)) section("Wells by locus", locus_verdict_table,
          paste0("Each well's delivered records, grouped by locus and judged within it: monoclonal ",
                 "when the second record of a locus holds at most a fifth of the first's reads. A well ",
                 "with one heavy and one kappa chain is two monoclonal loci, one cell. Only delivered ",
                 "records are counted, so a weak second clone below the per-well cut is not seen here.")) else "",

  if (NANOPORE && !is.null(NP$lib)) section("Library structure",
          paste0(fig_block("np_library"), fig_block("np_inserts")),
          paste0("What the PCR made, before any well is called. A run with few wells is either short ",
                 "of reads or short of product, and only this tells them apart: primer dimers and short ",
                 "products carry both primers and a well barcode, so the demultiplexer finds them and ",
                 "then drops them as too short. The per-well counts are in qc/nanopore/",
                 SAMPLE, ".library_qc_wells.tsv.")) else "",

  if (NANOPORE) section("Nanopore read quality",
          paste0(fig_block("np_quality"), fig_block("np_length"), fig_block("np_orient")),
          paste0("These have no Illumina counterpart. Illumina reads are near-uniform in quality; ",
                 "nanopore reads span a wide range, and it is consensus depth rather than any ",
                 "individual read that makes a base call trustworthy.")) else "",

  if (NANOPORE) section("Wells and clones",
          paste0(fig_block("np_dominance"), fig_block("np_depth")),
          paste0("Clone identity here comes from clustering each well's own reads rather than from ",
                 "annotating every one of them. That yields a number the previous pipeline could not ",
                 "report: the share of a well's reads belonging to its dominant clone, which is a ",
                 "direct doublet and barcode-collision signal.")) else "",

  section("Plate maps",
          plate_block,
          "Read depth per well, in real 96-well geometry. A dead column or a weak corner shows up as a shape here long before it shows up in a table."),

  section("Per-plate yield", paste0(fig_block("plate_bar"), plate_table)),

  section("Read depth", paste0(fig_block("depth_hist")),
          "Wells below the flag threshold are reported, never silently discarded."),

  section("Sequence lengths", paste0(fig_block("len_hist"), len_table),
          "Retention through trimming is the number to check first: a collapse here means the trim sequence did not match the library."),

  "<footer>Generated by TBRC <code>run_report.R</code> on ",
  esc(format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  ". Built from the QC tables in <code>qc/</code>; every figure is inline SVG, so this file is self-contained and needs no network. ",
  "Sources: <code>plate_read_depth.tsv</code>, <code>plate_demultiplexed.tsv</code>, <code>*.contig_lengths.tsv</code> and the delivered FASTA.",
  "</footer></div></body></html>"
)

html_path <- paste0(OUT_PREFIX, ".report.html")
writeLines(html, html_path)
cat("Wrote", html_path, "\n")

# ---------------------------------------------------------------- PDF

pdf_path <- paste0(OUT_PREFIX, ".report.pdf")
ok <- tryCatch({
  pdf(pdf_path, width = 11, height = 8.5, onefile = TRUE)
  # Cover page, drawn with base graphics so it needs no extra package.
  plot.new()
  text(0.5, 0.92, SAMPLE, cex = 2.2, font = 2)
  text(0.5, 0.85, "TBRC assembly run report", cex = 1.1, col = "grey30")
  lines <- c(
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M")),
    paste("Pipeline:", PIPELINE),
    if (nzchar(USER_ID)) paste("User:", USER_ID) else NULL,
    if (nzchar(CHAIN))   paste("Chain:", CHAIN) else NULL,
    if (nzchar(PRIMER))  paste("Barcode set:", PRIMER) else NULL,
    "",
    paste("Sequences delivered:", fmt_int(stats$records)),
    paste("Wells recovered:",     fmt_int(stats$wells)),
    paste("Plates:",              fmt_int(stats$plates)),
    paste("Median reads/well:",   if (is.finite(stats$depth_median)) fmt_int(stats$depth_median) else "--"),
    paste("Kept through trim:",   fmt_pct(retention)),
    paste("Barcode matched:",     fmt_pct(matched_pct)),
    paste("Low-depth wells flagged:", fmt_int(stats$low_wells))
  )
  for (i in seq_along(lines))
    text(0.5, 0.72 - (i - 1) * 0.042, lines[i], cex = 1.05,
         font = if (grepl("^Sequences delivered", lines[i])) 2 else 1)
  for (k in names(PLOTS)) print(PLOTS[[k]]$plot)
  for (p in plate_plots) print(p)
  dev.off(); TRUE
}, error = function(e) { try(dev.off(), silent = TRUE); message("PDF failed: ", conditionMessage(e)); FALSE })

if (ok) cat("Wrote", pdf_path, "\n")
