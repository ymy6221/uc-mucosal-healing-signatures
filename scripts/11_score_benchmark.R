# 对 8 个治疗队列和 3 个安慰剂臂计算全部注册签名，
# 同时生成 10000 次大小、平均表达和变异度匹配的随机基因集零分布。

check_benchmark_packages()
registry <- read_signature_registry()
source_ids <- c("GSE206285", "GSE92415", "GSE16879", "GSE12251",
                "GSE73661", "GSE23597", "GSE72819")

all_effects <- list()
all_coverage <- list()
all_real_z <- list()
all_summaries <- list()
effect_index <- 0L
coverage_index <- 0L
z_index <- 0L
summary_index <- 0L

for (source_number in seq_along(source_ids)) {
  source_id <- source_ids[source_number]
  message("读取并评分：", source_id)
  cohort_objects <- prepare_source_cohorts(source_id)

  for (cohort_number in seq_along(cohort_objects)) {
    object <- cohort_objects[[cohort_number]]
    scored <- score_registered_signatures(object$expression, registry)
    phenotype <- merge(object$phenotype, scored$scores,
                       by = "geo_accession", all.x = TRUE, sort = FALSE)
    phenotype <- phenotype[match(colnames(object$expression), geo_accession)]
    saveRDS(
      phenotype,
      file.path(DIR_INTERMEDIATE, paste0("benchmark_", object$cohort, "_scored.rds")),
      compress = "xz"
    )

    coverage_index <- coverage_index + 1L
    one_coverage <- data.table::copy(scored$coverage)
    one_coverage[, cohort := object$cohort]
    all_coverage[[coverage_index]] <- one_coverage

    treatment_names <- unique(phenotype$treatment)
    for (treatment_name in treatment_names) {
      effect_index <- effect_index + 1L
      all_effects[[effect_index]] <- fit_all_signature_effects(
        phenotype, object$cohort, treatment_name, registry
      )
    }

    active_rows <- phenotype$treatment == object$active_treatment & !is.na(phenotype$failure)
    active_accessions <- phenotype$geo_accession[active_rows]
    active_ranks <- scored$ranks[, active_accessions, drop = FALSE]
    active_outcome <- phenotype$failure[active_rows]
    random_null <- matched_random_null(
      active_ranks,
      active_outcome,
      registry,
      n_random = N_MATCHED_RANDOM,
      seed_offset = source_number * 100L + cohort_number
    )
    saveRDS(
      random_null,
      file.path(DIR_INTERMEDIATE, paste0("benchmark_", object$cohort, "_random_null.rds")),
      compress = "xz"
    )

    signature_ids <- unique(registry$signature_id)
    real_z <- lapply(signature_ids, function(sig_id) {
      values <- phenotype[active_rows, get(sig_id)]
      if (all(is.na(values)) || stats::sd(values, na.rm = TRUE) == 0) return(NULL)
      data.table::data.table(
        cohort = object$cohort,
        treatment = object$active_treatment,
        signature_id = sig_id,
        fisher_z = score_matrix_to_fisher_z(matrix(values, nrow = 1L), active_outcome),
        n = length(active_outcome),
        failures = sum(active_outcome == 1L)
      )
    })
    z_index <- z_index + 1L
    all_real_z[[z_index]] <- data.table::rbindlist(real_z, fill = TRUE)

    summary_index <- summary_index + 1L
    one_summary <- phenotype[, .(
      n = .N,
      failures = sum(failure == 1L),
      successes = sum(failure == 0L)
    ), by = treatment]
    one_summary[, cohort := object$cohort]
    data.table::setcolorder(one_summary, c("cohort", "treatment", "n", "failures", "successes"))
    all_summaries[[summary_index]] <- one_summary

    rm(scored, phenotype, active_ranks, random_null)
    invisible(gc())
    message("完成：", object$cohort)
  }
  rm(cohort_objects)
  invisible(gc())
}

effects <- data.table::rbindlist(all_effects, fill = TRUE)
coverage <- data.table::rbindlist(all_coverage, fill = TRUE)
real_z <- data.table::rbindlist(all_real_z, fill = TRUE)
summaries <- data.table::rbindlist(all_summaries, fill = TRUE)
data.table::fwrite(effects, file.path(DIR_RESULTS, "benchmark_cohort_effects.csv"))
data.table::fwrite(coverage, file.path(DIR_RESULTS, "benchmark_signature_coverage.csv"))
data.table::fwrite(real_z, file.path(DIR_RESULTS, "benchmark_real_fisher_z.csv"))
data.table::fwrite(summaries, file.path(DIR_RESULTS, "benchmark_sample_counts.csv"))
print(summaries)
