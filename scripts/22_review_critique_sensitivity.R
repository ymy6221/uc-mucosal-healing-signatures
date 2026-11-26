# 审稿意见定向升级：交互效应精度、终点异质性和单细胞逐患者剔除。
#
# 本脚本只读取已有派生结果 CSV，不读写或修改任何公共原始数据。
# 0.67–1.50 交互 OR 区间是审稿回应中的事后、探索性精度界限，
# 不是经临床验证的最小重要差异，也不将未显著结果解释为等效。

DIR_CBM_RESULTS <- file.path(DIR_RESULTS, "cbm_upgrade")
dir.create(DIR_CBM_RESULTS, recursive = TRUE, showWarnings = FALSE)

registry <- read_signature_registry()
primary_ids <- unique(registry[primary_set == TRUE, signature_id])
robust_ids <- c("IFX_ARJS_5", "OSM_OSMR", "IL22_NEUTROPHIL_6", "VDZ_NR_3")

# -----------------------------------------------------------------------------
# 1. 随机试验 treatment×signature 交互的精度/等效敏感性审计
# -----------------------------------------------------------------------------
interaction <- data.table::fread(file.path(DIR_RESULTS,
                                            "benchmark_interaction_meta.csv"))
interaction <- interaction[signature_id %in% primary_ids]
sample_counts <- data.table::fread(file.path(DIR_RESULTS,
                                              "benchmark_sample_counts.csv"))
randomized_cohorts <- c("GSE206285", "GSE92415", "GSE23597")
randomized_counts <- sample_counts[cohort %in% randomized_cohorts]

total_randomized_n <- sum(randomized_counts$n)
active_randomized_n <- sum(randomized_counts[treatment != "placebo", n])
placebo_randomized_n <- sum(randomized_counts[treatment == "placebo", n])

equivalence_low <- 2 / 3
equivalence_high <- 1.5
log_margin <- log(equivalence_high)

interaction_precision <- data.table::copy(interaction)
interaction_precision[, `:=`(
  total_randomized_n = total_randomized_n,
  active_randomized_n = active_randomized_n,
  placebo_randomized_n = placebo_randomized_n,
  illustrative_margin_low = equivalence_low,
  illustrative_margin_high = equivalence_high,
  ci_contains_null = ci_low <= 1 & ci_high >= 1,
  ci_within_illustrative_margin = ci_low >= equivalence_low &
    ci_high <= equivalence_high,
  ci_allows_positive_important_interaction = ci_high > equivalence_high,
  ci_allows_negative_important_interaction = ci_low < equivalence_low,
  normal_approx_mde_or_80 = exp((stats::qnorm(0.975) + stats::qnorm(0.80)) * se)
)]

# TOST 仅作事后精度审计：H0 为 beta<=-log(1.5) 或 beta>=log(1.5)。
interaction_precision[, `:=`(
  tost_lower_p = 1 - stats::pnorm((log_or + log_margin) / se),
  tost_upper_p = 1 - stats::pnorm((log_margin - log_or) / se)
)]
interaction_precision[, tost_p := pmax(tost_lower_p, tost_upper_p)]
interaction_precision[, equivalence_supported_05 := tost_p < 0.05]
interaction_precision[, precision_conclusion := data.table::fcase(
  equivalence_supported_05, "within exploratory 0.67-1.50 margin",
  ci_contains_null & (ci_allows_positive_important_interaction |
                        ci_allows_negative_important_interaction),
  "inconclusive: null and important interaction both compatible",
  ci_contains_null, "inconclusive: compatible with null",
  default = "interaction estimate excludes null"
)]

interaction_precision <- merge(
  interaction_precision,
  unique(registry[, .(signature_id, display_name)]),
  by = "signature_id", all.x = TRUE, sort = FALSE
)
data.table::setcolorder(interaction_precision, c(
  "signature_id", "display_name", "trials", "total_randomized_n",
  "active_randomized_n", "placebo_randomized_n", "log_or", "se", "or",
  "ci_low", "ci_high", "p_value", "fdr", "i2",
  "illustrative_margin_low", "illustrative_margin_high", "tost_p",
  "equivalence_supported_05", "ci_contains_null",
  "ci_within_illustrative_margin", "ci_allows_positive_important_interaction",
  "ci_allows_negative_important_interaction", "normal_approx_mde_or_80",
  "precision_conclusion"
))
data.table::fwrite(interaction_precision,
                   file.path(DIR_CBM_RESULTS, "interaction_precision_audit.csv"))

# -----------------------------------------------------------------------------
# 2. 临床终点调和登记表与分层 Meta 分析
# -----------------------------------------------------------------------------
# 分类完全基于本项目已锁定的队列结局定义。GSE12251 同时包含
# 应答与黏膜愈合，因此明确标记为“客观或混合黏膜终点”，不伪装成纯客观终点。
endpoint_registry <- data.table::data.table(
  cohort = c(
    "GSE206285", "GSE92415", "GSE16879", "GSE12251",
    "GSE73661_IFX", "GSE73661_VDZ", "GSE23597", "GSE72819"
  ),
  geo_accession = c(
    "GSE206285", "GSE92415", "GSE16879", "GSE12251",
    "GSE73661", "GSE73661", "GSE23597", "GSE72819"
  ),
  endpoint_definition = c(
    "week-8 mucosal healing", "week-6 response",
    "endoscopic and histologic response",
    "week-8 response and mucosal healing",
    "endoscopic mucosal healing", "endoscopic mucosal healing",
    "week-8 clinical response", "week-10 clinical remission"
  ),
  endpoint_family = c(
    "objective_or_composite_mucosal", "clinical_only",
    "objective_or_composite_mucosal", "objective_or_composite_mucosal",
    "objective_or_composite_mucosal", "objective_or_composite_mucosal",
    "clinical_only", "clinical_only"
  ),
  assessment_window = c(
    "week_6_8", "week_6_8", "other_post_induction",
    "week_6_8", "other_post_induction", "other_post_induction",
    "week_6_8", "other_post_induction"
  ),
  source = "locked project cohort phenotype harmonisation"
)
data.table::fwrite(endpoint_registry,
                   file.path(DIR_CBM_RESULTS, "endpoint_harmonization_registry.csv"))

effects <- data.table::fread(file.path(DIR_RESULTS, "benchmark_cohort_effects.csv"))
treated_names <- c("ustekinumab", "golimumab", "infliximab", "vedolizumab",
                   "etrolizumab")
active <- effects[treatment %in% treated_names & signature_id %in% primary_ids &
                    is.finite(log_or) & is.finite(se) & se > 0]
active <- merge(active, endpoint_registry, by = "cohort", all.x = TRUE, sort = FALSE)
if (anyNA(active$endpoint_family)) {
  stop("存在未登记的活性治疗队列终点。")
}

meta_one_group <- function(one, dimension, stratum) {
  if (nrow(one) < 2L) return(NULL)
  fit <- random_effects_meta(one$log_or, one$se)
  cbind(
    data.table::data.table(
      dimension = dimension,
      stratum = stratum,
      cohorts = nrow(one),
      total_n = sum(one$n),
      direction_positive_cohorts = sum(one$log_or > 0),
      direction_consistency = mean(one$log_or > 0)
    ),
    fit
  )
}

endpoint_meta_list <- list()
meta_index <- 0L
for (sig_current in primary_ids) {
  one_signature <- active[signature_id == sig_current]
  strata <- list(
    all_active = one_signature,
    objective_or_composite_mucosal =
      one_signature[endpoint_family == "objective_or_composite_mucosal"],
    clinical_only = one_signature[endpoint_family == "clinical_only"],
    week_6_8 = one_signature[assessment_window == "week_6_8"],
    other_post_induction =
      one_signature[assessment_window == "other_post_induction"]
  )
  dimensions <- c("overall", "endpoint_family", "endpoint_family",
                  "assessment_window", "assessment_window")
  for (i in seq_along(strata)) {
    estimate <- meta_one_group(strata[[i]], dimensions[i], names(strata)[i])
    if (is.null(estimate)) next
    meta_index <- meta_index + 1L
    endpoint_meta_list[[meta_index]] <- cbind(
      data.table::data.table(signature_id = sig_current), estimate
    )
  }
}
endpoint_meta <- data.table::rbindlist(endpoint_meta_list, fill = TRUE)
endpoint_meta[, fdr_within_stratum := stats::p.adjust(p_value, method = "BH"),
              by = .(dimension, stratum)]
endpoint_meta <- merge(
  endpoint_meta,
  unique(registry[, .(signature_id, display_name)]),
  by = "signature_id", all.x = TRUE, sort = FALSE
)
data.table::fwrite(endpoint_meta,
                   file.path(DIR_CBM_RESULTS, "endpoint_stratified_meta.csv"))

# 两个互不重叠终点家族的探索性组间差异检验。队列数仅 5 vs 3，
# 因此这是单因素敏感性检查，不称为多变量 Meta 回归。
endpoint_difference <- data.table::rbindlist(lapply(primary_ids, function(sig_current) {
  objective <- endpoint_meta[
    signature_id == sig_current & stratum == "objective_or_composite_mucosal"
  ]
  clinical <- endpoint_meta[
    signature_id == sig_current & stratum == "clinical_only"
  ]
  if (nrow(objective) != 1L || nrow(clinical) != 1L) return(NULL)
  difference <- objective$log_or - clinical$log_or
  difference_se <- sqrt(objective$se^2 + clinical$se^2)
  data.table::data.table(
    signature_id = sig_current,
    objective_log_or = objective$log_or,
    clinical_log_or = clinical$log_or,
    difference_log_or = difference,
    difference_se = difference_se,
    ratio_of_ors = exp(difference),
    ci_low = exp(difference - 1.96 * difference_se),
    ci_high = exp(difference + 1.96 * difference_se),
    p_value = 2 * stats::pnorm(-abs(difference / difference_se))
  )
}), fill = TRUE)
endpoint_difference[, fdr := stats::p.adjust(p_value, method = "BH")]
endpoint_difference <- merge(
  endpoint_difference,
  unique(registry[, .(signature_id, display_name)]),
  by = "signature_id", all.x = TRUE, sort = FALSE
)
data.table::fwrite(endpoint_difference,
                   file.path(DIR_CBM_RESULTS, "endpoint_subgroup_difference.csv"))

# -----------------------------------------------------------------------------
# 3. 单细胞定位的逐患者剔除（leave-one-patient-out, LOPO）
# -----------------------------------------------------------------------------
localization_scores <- data.table::fread(file.path(
  DIR_CBM_RESULTS, "singlecell_fine_localization_scores.csv"
))
if ("fine_annotations_MK_V4" %in% names(localization_scores) &&
    !"fine_cell_type" %in% names(localization_scores)) {
  data.table::setnames(localization_scores, "fine_annotations_MK_V4",
                       "fine_cell_type")
}
localization_tests <- data.table::fread(file.path(
  DIR_CBM_RESULTS, "singlecell_fine_localization_tests.csv"
))
patients <- sort(unique(localization_scores$patient_short))

retest_localization <- function(data) {
  answer <- data[stats::complete.cases(centered_score), {
    current_median <- stats::median(centered_score)
    if (.N < 5L || stats::sd(centered_score) == 0) {
      list(n_patients = .N, median_centered_score = current_median,
           localization_p = NA_real_)
    } else {
      test <- suppressWarnings(stats::wilcox.test(
        centered_score, mu = 0, alternative = "greater", exact = FALSE
      ))
      list(n_patients = .N, median_centered_score = current_median,
           localization_p = test$p.value)
    }
  }, by = .(signature_id, fine_cell_type, coarse_cell_type)]
  answer[, localization_fdr := stats::p.adjust(localization_p, method = "BH")]
  answer
}

lopo_full <- data.table::rbindlist(lapply(patients, function(omitted) {
  one <- localization_scores[patient_short != omitted, .(
    signature_id,
    fine_cell_type,
    coarse_cell_type,
    centered_score
  )]
  retested <- retest_localization(one)
  retested[, omitted_patient := omitted]
  retested
}), fill = TRUE)
data.table::setcolorder(lopo_full, c(
  "omitted_patient", "signature_id", "fine_cell_type", "coarse_cell_type",
  "n_patients", "median_centered_score", "localization_p", "localization_fdr"
))
data.table::fwrite(lopo_full, file.path(
  DIR_CBM_RESULTS, "singlecell_leave_one_patient_out.csv"
))

original_significant <- localization_tests[
  is.finite(localization_fdr) & localization_fdr < 0.05
]
lopo_significant <- merge(
  lopo_full,
  original_significant[, .(
    signature_id, fine_cell_type, coarse_cell_type,
    original_n_patients = n_patients,
    original_median_centered_score = median_centered_score,
    original_localization_p = localization_p,
    original_localization_fdr = localization_fdr
  )],
  by = c("signature_id", "fine_cell_type", "coarse_cell_type"),
  all = FALSE, sort = FALSE
)

lopo_summary <- lopo_significant[, .(
  omissions = .N,
  valid_test_fraction = mean(is.finite(localization_p)),
  positive_median_fraction = mean(is.finite(median_centered_score) &
                                    median_centered_score > 0),
  nominal_significance_fraction = mean(localization_p < 0.05, na.rm = TRUE),
  fdr_significance_fraction = mean(localization_fdr < 0.05, na.rm = TRUE),
  min_median_centered_score = min(median_centered_score, na.rm = TRUE),
  max_median_centered_score = max(median_centered_score, na.rm = TRUE),
  max_localization_p = max(localization_p, na.rm = TRUE),
  max_localization_fdr = max(localization_fdr, na.rm = TRUE)
), by = .(
  signature_id, fine_cell_type, coarse_cell_type,
  original_n_patients, original_median_centered_score,
  original_localization_p, original_localization_fdr
)]
lopo_summary[, `:=`(
  sign_stable_all_omissions = positive_median_fraction == 1,
  nominal_stable_all_omissions = nominal_significance_fraction == 1,
  fdr_stable_all_omissions = fdr_significance_fraction == 1
)]
lopo_summary <- merge(
  lopo_summary,
  unique(registry[, .(signature_id, display_name)]),
  by = "signature_id", all.x = TRUE, sort = FALSE
)
data.table::fwrite(lopo_summary, file.path(
  DIR_CBM_RESULTS, "singlecell_leave_one_patient_out_summary.csv"
))

message("审稿意见定向敏感性分析完成；未修改任何原始数据。")
print(interaction_precision[, .(
  signature_id, or, ci_low, ci_high, tost_p, equivalence_supported_05,
  precision_conclusion
)])
print(endpoint_meta[signature_id %in% robust_ids &
                      stratum %in% c("objective_or_composite_mucosal", "clinical_only"),
                    .(signature_id, stratum, cohorts, total_n, or, ci_low, ci_high,
                      p_value, fdr_within_stratum)])
print(lopo_summary[, .(
  originally_significant_tests = .N,
  sign_stable = sum(sign_stable_all_omissions),
  fdr_stable_all = sum(fdr_stable_all_omissions)
)])
