args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript run_clonality.R <igblast_tsv> <output_dir> <sample_name> [igblast_panel]", call. = FALSE)
}

Sys.setenv(OMP_NUM_THREADS = "1", OMP_THREAD_LIMIT = "1")
options(sd_num_thread = 1)

if (!requireNamespace("clonality", quietly = TRUE)) {
  stop("run_clonality.R requires the clonality R package.", call. = FALSE)
}

igblast_tsv <- args[1]
output_dir <- args[2]
sample_name <- args[3]
igblast_panel <- if (length(args) >= 4) args[4] else "ig_all"

cell_type_for_panel <- function(panel_name) {
  if (panel_name %in% c("igh", "ig_light", "ig_all")) {
    return("B")
  }

  if (panel_name %in% c("tra", "trb", "tcr_all")) {
    return("T")
  }

  NA_character_
}

cell_type <- cell_type_for_panel(igblast_panel)
if (is.na(cell_type)) {
  message(sprintf("Clonality is not configured for mixed receptor scope: %s", igblast_panel))
  quit(save = "no", status = 2)
}

if (!file.exists(igblast_tsv)) {
  stop(sprintf("IGBlast TSV not found: %s", igblast_tsv), call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

igblast_df <- read.table(
  igblast_tsv,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pick_col <- function(candidates) {
  match_idx <- match(tolower(candidates), tolower(names(igblast_df)))
  found <- names(igblast_df)[match_idx[!is.na(match_idx)]]
  if (length(found)) {
    found[[1]]
  } else {
    NA_character_
  }
}

sequence_id_col <- pick_col(c("sequence_id", "sequence_id_1", "Sequence ID"))
v_call_col <- pick_col(c("v_call", "best_v_hit", "V-GENE and allele", "v_gene"))
j_call_col <- pick_col(c("j_call", "best_j_hit", "J-GENE and allele", "j_gene"))
junction_col <- pick_col(c("junction", "junction_aa", "cdr3", "CDR3-IMGT"))
productive_col <- pick_col(c("productive", "functionality"))

required_map <- c(
  sequence_id = sequence_id_col,
  v_call = v_call_col,
  j_call = j_call_col,
  junction = junction_col
)

if (anyNA(required_map)) {
  missing_keys <- names(required_map)[is.na(required_map)]
  stop(
    sprintf(
      "IGBlast output is missing required clonality columns: %s",
      paste(missing_keys, collapse = ", ")
    ),
    call. = FALSE
  )
}

clonality_input <- data.frame(
  sequence_id = trimws(as.character(igblast_df[[sequence_id_col]])),
  v_call = trimws(as.character(igblast_df[[v_call_col]])),
  j_call = trimws(as.character(igblast_df[[j_call_col]])),
  junction = gsub("[^A-Za-z]", "", trimws(as.character(igblast_df[[junction_col]]))),
  productive = if (!is.na(productive_col)) trimws(as.character(igblast_df[[productive_col]])) else "",
  stringsAsFactors = FALSE
)

row_keep <- nzchar(clonality_input$sequence_id) &
  nzchar(clonality_input$v_call) &
  nzchar(clonality_input$j_call) &
  nzchar(clonality_input$junction)
clonality_input <- clonality_input[row_keep, , drop = FALSE]

productive_flags <- tolower(clonality_input$productive)
productive_keep <- productive_flags %in% c("", "true", "t", "productive", "yes")
clonality_input <- clonality_input[productive_keep, , drop = FALSE]

if (!nrow(clonality_input)) {
  message("No productive IGBlast records were available for clonality analysis.")
  quit(save = "no", status = 2)
}

clonality_output <- clonality::clonality(
  data = clonality_input,
  ident_col = "sequence_id",
  vgene_col = "v_call",
  jgene_col = "j_call",
  cdr3_col = "junction",
  cell = cell_type,
  output_original = TRUE,
  quiet = TRUE,
  project = sample_name
)

# clonality 0.10 wrote 19 columns; 0.15 joins the clone label onto the five input
# columns only and writes 9. Nine of the ten dropped columns repeat another one
# (CDR3 = junction, *_original = v_call/j_call, *_all = v_genes/j_genes) or are
# multi-annotation bookkeeping that IgBLAST's single call leaves NA/FALSE. CDR3L,
# the CDR3 length, is the one a reader might filter on, so keep it and keep the
# surviving columns in the order delivered archives have always had. A 0.10
# result already has CDR3L and is written exactly as before. Checked on five
# finished runs: CDR3L == nchar(junction) on every row.
if (!"CDR3L" %in% names(clonality_output) && "junction" %in% names(clonality_output)) {
  clonality_output$CDR3L <- nchar(clonality_output$junction)
  legacy_order <- c("sequence_id", "v_call", "j_call", "junction", "productive",
                    "v_genes", "j_genes", "CDR3L", "clonality", "project")
  clonality_output <- clonality_output[, c(intersect(legacy_order, names(clonality_output)),
                                           setdiff(names(clonality_output), legacy_order)),
                                       drop = FALSE]
}

write.table(
  clonality_output,
  file = file.path(output_dir, sprintf("%s.clonality.tsv", sample_name)),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_df <- as.data.frame(sort(table(clonality_output$clonality), decreasing = TRUE), stringsAsFactors = FALSE)
colnames(summary_df) <- c("clonality", "sequence_count")

write.table(
  summary_df,
  file = file.path(output_dir, sprintf("%s.clonality.summary.tsv", sample_name)),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
