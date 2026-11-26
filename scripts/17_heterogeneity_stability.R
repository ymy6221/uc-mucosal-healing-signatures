# CBM 升级 2：在仅有 8 个队列的条件下，用预测区间、逐队列删除、
# 平台/试验敏感性和影响度拆解异质性；不做低功效的 Meta 回归。

DIR_CBM_RESULTS <- file.path(DIR_RESULTS, "cbm_upgrade")
dir.create(DIR_CBM_RESULTS, recursive = TRUE, showWarnings = FALSE)

registry <- read_signature_registry()
primary_ids <- unique(registry[primary_set == TRUE, signature_id])
effects <- data.table::fread(file.path(DIR_RESULTS, "benchmark_cohort_effects.csv"))
treated_names <- c("ustekinumab", "golimumab", "infliximab", "vedolizumab", "etrolizumab")
treated <- effects[treatment %in% treated_names & signature_id %in% primary_ids &
                     is.finite(log_or) & is.finite(se) & se > 0]
treated[, platform := ifelse(cohort == "GSE72819", "RNA-seq", "microarray")]
treated[, randomized_trial := cohort %in% c("GSE206285", "GSE92415", "GSE23597")]

meta_with_prediction <- function(one, label) {
  if (nrow(one) < 2L) return(NULL)
  meta <- random_effects_meta(one$log_or, one$se)
  k <- nrow(one)
  multiplier <- if (k >= 3L) stats::qt(0.975, df = k - 2L) else NA_real_
  prediction_se <- sqrt(meta$tau2 + meta$se^2)
  data.table::data.table(
    analysis = label,
    cohorts = k,
    total_n = sum(one$n),
    log_or = meta$log_or,
    se = meta$se,
    or = meta$or,
    ci_low = meta$ci_low,
    ci_high = meta$ci_high,
    p_value = meta$p_value,
    tau2 = meta$tau2,
    q = meta$q,
    i2 = meta$i2,
    prediction_low = if (k >= 3L) exp(meta$log_or - multiplier * prediction_se) else NA_real_,
    prediction_high = if (k >= 3L) exp(meta$log_or + multiplier * prediction_se) else NA_real_
  )
}

stability <- list()
stability_index <- 0L
leave_one_out <- list()
leave_index <- 0L
influence <- list()
influence_index <- 0L

for (sig_current in primary_ids) {
  # data.table 中避免变量名与列名冲突。
  one <- treated[treated$signature_id == sig_current]
  subsets <- list(
    all_active = one,
    microarray_only = one[platform == "microarray"],
    randomized_active_arms = one[randomized_trial == TRUE],
    excluding_largest_cohort = one[cohort != "GSE206285"]
  )
  for (label in names(subsets)) {
    estimate <- meta_with_prediction(subsets[[label]], label)
    if (is.null(estimate)) next
    stability_index <- stability_index + 1L
    stability[[stability_index]] <- cbind(
      data.table::data.table(signature_id = sig_current), estimate
    )
  }

  all_meta <- random_effects_meta(one$log_or, one$se)
  fixed_weight <- 1 / one$se^2
  fixed_mean <- sum(fixed_weight * one$log_or) / sum(fixed_weight)
  q_contribution <- fixed_weight * (one$log_or - fixed_mean)^2
  for (i in seq_len(nrow(one))) {
    reduced <- one[-i]
    reduced_meta <- random_effects_meta(reduced$log_or, reduced$se)
    leave_index <- leave_index + 1L
    leave_one_out[[leave_index]] <- data.table::data.table(
      signature_id = sig_current,
      omitted_cohort = one$cohort[i],
      cohorts = nrow(reduced),
      total_n = sum(reduced$n),
      or = reduced_meta$or,
      ci_low = reduced_meta$ci_low,
      ci_high = reduced_meta$ci_high,
      p_value = reduced_meta$p_value,
      tau2 = reduced_meta$tau2,
      i2 = reduced_meta$i2
    )
    influence_index <- influence_index + 1L
    influence[[influence_index]] <- data.table::data.table(
      signature_id = sig_current,
      cohort = one$cohort[i],
      q_contribution = q_contribution[i],
      delta_pooled_log_or = all_meta$log_or - reduced_meta$log_or,
      abs_delta_pooled_log_or = abs(all_meta$log_or - reduced_meta$log_or),
      i2_all = all_meta$i2,
      i2_without = reduced_meta$i2,
      delta_i2 = all_meta$i2 - reduced_meta$i2
    )
  }
}

stability <- data.table::rbindlist(stability, fill = TRUE)
leave_one_out <- data.table::rbindlist(leave_one_out, fill = TRUE)
influence <- data.table::rbindlist(influence, fill = TRUE)
stability[, fdr := stats::p.adjust(p_value, method = "BH"), by = analysis]

data.table::fwrite(stability,
                   file.path(DIR_CBM_RESULTS, "meta_prediction_and_sensitivity.csv"))
data.table::fwrite(leave_one_out,
                   file.path(DIR_CBM_RESULTS, "meta_leave_one_cohort_out.csv"))
data.table::fwrite(influence,
                   file.path(DIR_CBM_RESULTS, "meta_influence_baujat.csv"))

message("异质性稳健性分析完成。")
print(stability[signature_id %in% c("IFX_ARJS_5", "OSM_OSMR",
                                   "IL22_NEUTROPHIL_6", "VDZ_NR_3")])
