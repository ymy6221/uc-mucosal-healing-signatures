#!/usr/bin/env Rscript
# 用法：
# Rscript tools/score_public_signature_panel.R expression.csv output_directory
# expression.csv：第一列为 HGNC gene symbol，后续每列为一个样本的表达量。
# 输出不是临床概率，只是可移植的签名分数、覆盖度和公共参考百分位。

all_arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", all_arguments, value = TRUE)
script_path <- if (length(file_argument)) sub("^--file=", "", file_argument[1L]) else NA_character_
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop("用法：Rscript tools/score_public_signature_panel.R expression.csv output_directory")
}
input_path <- normalizePath(arguments[1L], mustWork = TRUE)
output_directory <- arguments[2L]
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("data.table", quietly = TRUE)) stop("缺少 data.table 包。")
registry <- data.table::fread(file.path(project_root, "signatures",
                                        "published_signature_registry.csv"))
registry <- registry[primary_set == TRUE & definition_status == "exact_text"]
robust_ids <- c("IFX_ARJS_5", "OSM_OSMR", "IL22_NEUTROPHIL_6", "VDZ_NR_3")

input <- data.table::fread(input_path)
if (ncol(input) < 2L) stop("输入至少需要一列基因符号和一列样本表达。")
data.table::setnames(input, 1L, "gene")
input[, gene := toupper(trimws(as.character(gene)))]
input <- input[!is.na(gene) & gene != ""]
sample_columns <- setdiff(names(input), "gene")
for (column in sample_columns) {
  input[, (column) := suppressWarnings(as.numeric(get(column)))]
}
input <- input[, lapply(.SD, stats::median, na.rm = TRUE), by = gene,
               .SDcols = sample_columns]
expression <- as.matrix(input[, ..sample_columns])
rownames(expression) <- input$gene
if (any(!is.finite(expression))) {
  stop("表达矩阵含缺失或非数值；请先清理 NA/Inf。")
}

ranks <- apply(expression, 2L, rank, ties.method = "average") /
  (nrow(expression) + 1)
if (is.null(dim(ranks))) ranks <- matrix(ranks, ncol = 1L)
rownames(ranks) <- rownames(expression)
colnames(ranks) <- colnames(expression)

signature_ids <- unique(registry$signature_id)
scores <- data.table::data.table(sample = colnames(ranks))
coverage_rows <- list()
for (i in seq_along(signature_ids)) {
  sig_current <- signature_ids[i]
  definition <- registry[registry$signature_id == sig_current]
  present <- definition[gene %in% rownames(ranks)]
  required <- if (nrow(definition) <= 3L) nrow(definition) else ceiling(0.60 * nrow(definition))
  values <- rep(NA_real_, ncol(ranks))
  if (nrow(present) >= required) {
    positive <- present[component_weight > 0, gene]
    negative <- present[component_weight < 0, gene]
    values <- colMeans(ranks[positive, , drop = FALSE])
    if (length(negative) > 0L) {
      values <- values - colMeans(ranks[negative, , drop = FALSE])
    }
    values <- values * unique(present$risk_orientation)
  }
  scores[, (sig_current) := as.numeric(values)]
  coverage_rows[[i]] <- data.table::data.table(
    signature_id = sig_current,
    genes_total = nrow(definition),
    genes_present = nrow(present),
    required = required,
    coverage = nrow(present) / nrow(definition),
    passed = nrow(present) >= required,
    missing_genes = paste(setdiff(definition$gene, present$gene), collapse = ";")
  )
}
coverage <- data.table::rbindlist(coverage_rows)
if (all(robust_ids %in% names(scores))) {
  scores[, robust_consensus_rank := rowMeans(.SD, na.rm = TRUE), .SDcols = robust_ids]
}

reference_path <- file.path(project_root, "results", "cbm_upgrade",
                            "signature_reference_distribution.csv")
percentiles <- data.table::data.table(sample = scores$sample)
if (file.exists(reference_path)) {
  reference <- data.table::fread(reference_path)
  for (score_current in intersect(setdiff(names(scores), "sample"), unique(reference$score_id))) {
    reference_values <- reference[reference$score_id == score_current & is.finite(score), score]
    percentile_function <- stats::ecdf(reference_values)
    percentiles[, (paste0(score_current, "_reference_percentile")) :=
                  as.numeric(percentile_function(scores[[score_current]]))]
  }
}

data.table::fwrite(scores, file.path(output_directory, "signature_scores.csv"))
data.table::fwrite(coverage, file.path(output_directory, "signature_coverage.csv"))
data.table::fwrite(percentiles, file.path(output_directory, "reference_percentiles.csv"))
note <- c(
  "This tool does not estimate treatment response probability.",
  "Higher risk-oriented scores indicate greater expression of a published failure-associated module.",
  "Reference percentiles are descriptive positions in the pooled public baseline compendium.",
  "Do not use for clinical decisions without prospective analytical and clinical validation."
)
writeLines(note, file.path(output_directory, "IMPORTANT_README.txt"))
message("评分完成：", normalizePath(output_directory, mustWork = TRUE))
