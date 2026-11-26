# 自动验收升级阶段。技术复现失败会立即停止；生物学门槛失败不会改阈值，
# 而会把状态明确写成 DO_NOT_SUBMIT，防止把阴性结果包装成阳性结果。

required_results <- c(
  "benchmark_cohort_effects.csv",
  "benchmark_sample_counts.csv",
  "benchmark_signature_coverage.csv",
  "benchmark_active_meta.csv",
  "benchmark_placebo_meta.csv",
  "benchmark_leave_source_out_meta.csv",
  "matched_random_empirical_tests.csv",
  "benchmark_treatment_interactions.csv",
  "benchmark_interaction_meta.csv",
  "benchmark_severity_correlations.csv",
  "benchmark_severity_residual_effects.csv",
  "GSE72819_TNFnaive_sensitivity.csv",
  "signature_convergence_loadings.csv",
  "signature_convergence_variance.csv",
  "published_signature_evidence_matrix.csv",
  "high_impact_claim_gate.csv"
)
required_figure_stems <- c(
  "MainFigure1_upgraded_study_design",
  "MainFigure2_published_signature_benchmark",
  "MainFigure3_placebo_and_specificity",
  "MainFigure4_random_controls_and_convergence"
)
required_figures <- as.vector(outer(required_figure_stems, c(".pdf", ".tiff"), paste0))
result_paths <- file.path(DIR_RESULTS, required_results)
figure_paths <- file.path(DIR_FIGURES, required_figures)

sample_counts_check <- data.table::fread(file.path(DIR_RESULTS, "benchmark_sample_counts.csv"))
coverage_check <- data.table::fread(file.path(DIR_RESULTS, "benchmark_signature_coverage.csv"))
interactions_check <- data.table::fread(file.path(DIR_RESULTS, "benchmark_treatment_interactions.csv"))
leave_source_check <- data.table::fread(file.path(DIR_RESULTS, "benchmark_leave_source_out_meta.csv"))
empirical_check <- data.table::fread(file.path(DIR_RESULTS, "matched_random_empirical_tests.csv"))
claim_check <- data.table::fread(file.path(DIR_RESULTS, "high_impact_claim_gate.csv"))
registry_check <- read_signature_registry()
signature_status <- unique(registry_check[, .(signature_id, primary_set, definition_status)])
coverage_check <- merge(coverage_check, signature_status, by = "signature_id", all.x = TRUE)

treated_names_check <- c("ustekinumab", "golimumab", "infliximab", "vedolizumab", "etrolizumab")
active_rows_check <- sample_counts_check[treatment %in% treated_names_check]
placebo_rows_check <- sample_counts_check[treatment == "placebo"]
primary_exact_check <- signature_status[primary_set == TRUE]

random_files <- file.path(
  DIR_INTERMEDIATE,
  paste0("benchmark_", c("GSE206285", "GSE92415", "GSE16879", "GSE12251",
                         "GSE73661_IFX", "GSE73661_VDZ", "GSE23597", "GSE72819"),
         "_random_null.rds")
)
random_iterations_check <- if (all(file.exists(random_files))) {
  vapply(random_files, function(path) nrow(readRDS(path)), integer(1))
} else {
  integer()
}

download_paths <- file.path(DIR_RAW, BENCHMARK_FILES$filename)
download_sizes <- ifelse(file.exists(download_paths), file.info(download_paths)$size, 0)

technical_validation <- data.table::data.table(
  check = c(
    "all required result tables exist and are non-empty",
    "all final Main Figure 1-4 PDF and 600-dpi TIFF files exist",
    "three benchmark downloads exceed locked minimum sizes",
    "eight active cohorts contain exactly 629 samples",
    "three placebo arms contain exactly 225 samples",
    "seven primary signatures are exact-text definitions only",
    "all primary signatures have at least 75 percent gene coverage",
    "all eight cohorts contain exactly the requested random iterations",
    "three randomized trials have interaction models for every registered signature",
    "source-cohort-excluded Meta-analysis covers all seven primary signatures",
    "empirical tests cover all seven primary signatures",
    "biological claim gate contains six locked tests"
  ),
  observed = c(
    paste0(sum(file.exists(result_paths) & file.info(result_paths)$size > 100L), "/", length(result_paths)),
    paste0(sum(file.exists(figure_paths) & file.info(figure_paths)$size > 100L), "/", length(figure_paths)),
    paste(download_sizes, collapse = ";"),
    paste0(nrow(active_rows_check), " cohorts, n=", sum(active_rows_check$n)),
    paste0(nrow(placebo_rows_check), " arms, n=", sum(placebo_rows_check$n)),
    paste0(nrow(primary_exact_check), " exact primary signatures"),
    sprintf("minimum %.2f", min(coverage_check[primary_set == TRUE, coverage])),
    if (length(random_iterations_check)) paste(random_iterations_check, collapse = ";") else "missing",
    paste(interactions_check[, .N, by = cohort][, paste0(cohort, "=", N)], collapse = ";"),
    paste0(sum(primary_exact_check$signature_id %in% leave_source_check$signature_id), "/7"),
    paste0(sum(primary_exact_check$signature_id %in% empirical_check$signature_id), "/7"),
    paste0(nrow(claim_check), " gates")
  ),
  passed = c(
    all(file.exists(result_paths) & file.info(result_paths)$size > 100L),
    all(file.exists(figure_paths) & file.info(figure_paths)$size > 100L),
    all(download_sizes >= BENCHMARK_FILES$min_bytes),
    nrow(active_rows_check) == 8L && sum(active_rows_check$n) == 629L,
    nrow(placebo_rows_check) == 3L && sum(placebo_rows_check$n) == 225L,
    nrow(primary_exact_check) == 7L && all(primary_exact_check$definition_status == "exact_text"),
    min(coverage_check[primary_set == TRUE, coverage]) >= 0.75,
    length(random_iterations_check) == 8L && all(random_iterations_check == N_MATCHED_RANDOM),
    data.table::uniqueN(interactions_check$cohort) == 3L &&
      all(interactions_check[, .N, by = cohort]$N == data.table::uniqueN(registry_check$signature_id)),
    all(primary_exact_check$signature_id %in% leave_source_check$signature_id),
    all(primary_exact_check$signature_id %in% empirical_check$signature_id),
    nrow(claim_check) == 6L
  )
)
data.table::fwrite(
  technical_validation,
  file.path(DIR_RESULTS, "benchmark_technical_validation.csv")
)
print(technical_validation)

if (!all(technical_validation$passed)) {
  stop("升级阶段技术验收失败；请检查 results/benchmark_technical_validation.csv")
}

readiness <- data.table::data.table(
  technical_reproducibility = "PASS",
  fixed_target_claim_gate = if (all(claim_check$passed)) "PASS" else "FAIL",
  decision = if (all(claim_check$passed)) {
    "READY_TO_DRAFT_FOR_FIXED_HIGH_IMPACT_TARGET"
  } else {
    "DO_NOT_SUBMIT_FIXED_TARGET; RETAIN_LOCKED_NEGATIVE_RESULT"
  },
  random_iterations = N_MATCHED_RANDOM,
  active_samples = sum(active_rows_check$n),
  placebo_samples = sum(placebo_rows_check$n)
)
data.table::fwrite(readiness, file.path(DIR_RESULTS, "submission_readiness.csv"))
capture.output(sessionInfo(), file = file.path(DIR_RESULTS, "sessionInfo_stage5.txt"))
version_packages <- c("R", "data.table", "ggplot2", "pROC", "patchwork", "Matrix",
                      "scales", "digest", "GEOquery", "Biobase", "AnnotationDbi",
                      "org.Hs.eg.db")
package_versions <- data.table::data.table(
  package = version_packages,
  version = c(
    as.character(getRversion()),
    vapply(version_packages[-1L], function(pkg) as.character(utils::packageVersion(pkg)),
           character(1))
  )
)
data.table::fwrite(package_versions, file.path(DIR_RESULTS, "package_versions_stage5.csv"))
print(readiness)
message("升级阶段技术验收全部通过。投稿决策请以 results/submission_readiness.csv 为准。")
