# 解析每个试验的基线结局，并完成主效应、临床校正和治疗交互分析。

check_packages()

# ---------- UNIFI：Ustekinumab 与安慰剂 ----------
unifi <- data.table::as.data.table(readRDS(file.path(DIR_INTERMEDIATE, "GSE206285_scored.rds")))
unifi[, treatment := data.table::fcase(
  grepl("Ustekinumab", get("treatment:ch1")), "ustekinumab",
  grepl("Placebo", get("treatment:ch1")), "placebo",
  default = NA_character_
)]
unifi[, failure := data.table::fcase(
  get("mucosal healing at week 8:ch1") == "N", 1L,
  get("mucosal healing at week 8:ch1") == "Y", 0L,
  default = NA_integer_
)]
unifi_clinical <- data.table::fread(file.path(DIR_RAW, "GSE206285_additional_sample_metadata.tsv.gz"))
data.table::setnames(
  unifi_clinical,
  c("Accession", "Field[CRP]", "Field[Calprotectin]", "Field[Mayo_score]",
    "Field[Disease_duration]", "Field[Anti_TNF_history]", "Field[Neutrophils_LP]",
    "Field[Neutrophils_Epithelium]"),
  c("geo_accession", "crp", "calprotectin", "mayo", "disease_duration",
    "anti_tnf_history", "neutrophils_lp", "neutrophils_epithelium")
)
unifi <- merge(unifi, unifi_clinical, by = "geo_accession", all.x = TRUE, sort = FALSE)
for (column in c("crp", "calprotectin", "mayo", "disease_duration",
                 "neutrophils_lp", "neutrophils_epithelium")) {
  data.table::set(unifi, j = column, value = as_numeric_safely(unifi[[column]]))
}
unifi[, eri_z := standardize_score(eri)]

# UST 和安慰剂分别报告，防止把预后效应误叫作药物预测效应。
result_unifi_active <- fit_simple_cohort(unifi[treatment == "ustekinumab"], "GSE206285", "ustekinumab")
result_unifi_placebo <- fit_simple_cohort(unifi[treatment == "placebo"], "GSE206285", "placebo")

# 主校正模型只在 UST 组内运行；连续炎症指标先 log1p 降低极端值影响。
unifi_active_adjusted <- unifi[
  treatment == "ustekinumab" &
    stats::complete.cases(failure, eri_z, mayo, crp, calprotectin, anti_tnf_history)
]
unifi_active_adjusted[, anti_tnf_history := factor(anti_tnf_history)]
model_unifi_adjusted <- stats::glm(
  failure ~ eri_z + mayo + log1p(crp) + log1p(calprotectin) + anti_tnf_history,
  family = stats::binomial(),
  data = unifi_active_adjusted
)
result_unifi_adjusted <- extract_glm_term(
  model_unifi_adjusted, "eri_z", "GSE206285", "ustekinumab", "clinical_adjusted"
)

# 更严格的敏感性模型再校正固有层和上皮内中性粒细胞评分。
unifi_histology_adjusted <- unifi[
  treatment == "ustekinumab" &
    stats::complete.cases(failure, eri_z, mayo, crp, calprotectin, anti_tnf_history,
                          neutrophils_lp, neutrophils_epithelium)
]
unifi_histology_adjusted[, anti_tnf_history := factor(anti_tnf_history)]
model_unifi_histology <- stats::glm(
  failure ~ eri_z + mayo + log1p(crp) + log1p(calprotectin) + anti_tnf_history +
    neutrophils_lp + neutrophils_epithelium,
  family = stats::binomial(),
  data = unifi_histology_adjusted
)
result_unifi_histology <- extract_glm_term(
  model_unifi_histology, "eri_z", "GSE206285", "ustekinumab", "histology_adjusted"
)

# 比较“仅临床模型”和“临床+ERI 模型”的 AUC、Brier 分数和似然比检验。
model_unifi_clinical_only <- stats::glm(
  failure ~ mayo + log1p(crp) + log1p(calprotectin) + anti_tnf_history,
  family = stats::binomial(),
  data = unifi_active_adjusted
)
prediction_clinical <- stats::predict(model_unifi_clinical_only, type = "response")
prediction_augmented <- stats::predict(model_unifi_adjusted, type = "response")
roc_clinical <- pROC::roc(unifi_active_adjusted$failure, prediction_clinical, quiet = TRUE)
roc_augmented <- pROC::roc(unifi_active_adjusted$failure, prediction_augmented, quiet = TRUE)
roc_comparison <- pROC::roc.test(roc_clinical, roc_augmented, paired = TRUE, method = "delong")
lrt_comparison <- stats::anova(model_unifi_clinical_only, model_unifi_adjusted, test = "LRT")
incremental_performance <- data.table::data.table(
  n = nrow(unifi_active_adjusted),
  auc_clinical = as.numeric(pROC::auc(roc_clinical)),
  auc_clinical_plus_eri = as.numeric(pROC::auc(roc_augmented)),
  delta_auc = as.numeric(pROC::auc(roc_augmented) - pROC::auc(roc_clinical)),
  delong_p = roc_comparison$p.value,
  brier_clinical = mean((unifi_active_adjusted$failure - prediction_clinical)^2),
  brier_clinical_plus_eri = mean((unifi_active_adjusted$failure - prediction_augmented)^2),
  likelihood_ratio_p = lrt_comparison$`Pr(>Chi)`[2L]
)

# 治疗×ERI 是“预测”与“仅预后”的正式区分检验。
unifi_interaction_data <- unifi[
  treatment %in% c("placebo", "ustekinumab") &
    stats::complete.cases(failure, eri_z, mayo, crp, calprotectin, anti_tnf_history)
]
unifi_interaction_data[, treatment := stats::relevel(factor(treatment), ref = "placebo")]
unifi_interaction_data[, anti_tnf_history := factor(anti_tnf_history)]
model_unifi_interaction <- stats::glm(
  failure ~ treatment * eri_z + mayo + log1p(crp) + log1p(calprotectin) + anti_tnf_history,
  family = stats::binomial(),
  data = unifi_interaction_data
)
result_unifi_interaction <- extract_glm_term(
  model_unifi_interaction,
  "treatmentustekinumab:eri_z",
  "GSE206285",
  "ustekinumab_vs_placebo",
  "treatment_interaction"
)

# ---------- PURSUIT-SC：Golimumab 与安慰剂 ----------
glm_trial <- data.table::as.data.table(readRDS(file.path(DIR_INTERMEDIATE, "GSE92415_scored.rds")))
glm_trial <- glm_trial[get("visit:ch1") == "Week 0" & get("treatment:ch1") %in% c("golimumab", "Placebo")]
glm_trial[, treatment := tolower(get("treatment:ch1"))]
glm_trial[, failure := data.table::fcase(
  get("wk6response:ch1") == "No", 1L,
  get("wk6response:ch1") == "Yes", 0L,
  default = NA_integer_
)]
glm_trial[, mayo := as_numeric_safely(get("mayo score:ch1"))]
glm_trial[, age := as_numeric_safely(get("age:ch1"))]
glm_trial[, eri_z := standardize_score(eri)]
result_glm_active <- fit_simple_cohort(glm_trial[treatment == "golimumab"], "GSE92415", "golimumab")
result_glm_placebo <- fit_simple_cohort(glm_trial[treatment == "placebo"], "GSE92415", "placebo")
glm_interaction_data <- glm_trial[stats::complete.cases(failure, eri_z, mayo, age)]
glm_interaction_data[, treatment := stats::relevel(factor(treatment), ref = "placebo")]
model_glm_interaction <- stats::glm(
  failure ~ treatment * eri_z + mayo + age,
  family = stats::binomial(),
  data = glm_interaction_data
)
result_glm_interaction <- extract_glm_term(
  model_glm_interaction,
  "treatmentgolimumab:eri_z",
  "GSE92415",
  "golimumab_vs_placebo",
  "treatment_interaction"
)

# ---------- GSE16879：Infliximab ----------
ifx_16879 <- data.table::as.data.table(readRDS(file.path(DIR_INTERMEDIATE, "GSE16879_scored.rds")))
ifx_16879 <- ifx_16879[
  get("disease:ch1") == "UC" & grepl("Before", get("before or after first infliximab treatment:ch1"))
]
ifx_16879[, treatment := "infliximab"]
ifx_16879[, failure := data.table::fcase(
  get("response to infliximab:ch1") == "No", 1L,
  get("response to infliximab:ch1") == "Yes", 0L,
  default = NA_integer_
)]
result_ifx_16879 <- fit_simple_cohort(ifx_16879, "GSE16879", "infliximab")

# ---------- GSE12251：Infliximab，P13 技术重复只保留第一个 ----------
ifx_12251 <- data.table::as.data.table(readRDS(file.path(DIR_INTERMEDIATE, "GSE12251_scored.rds")))
ifx_12251[, patient_id := sub("/.*$", "", title)]
ifx_12251 <- ifx_12251[!duplicated(patient_id)]
ifx_12251[, treatment := "infliximab"]
ifx_12251[, failure := data.table::fcase(
  get("WK8RSPHM:ch1") == "No", 1L,
  get("WK8RSPHM:ch1") == "Yes", 0L,
  default = NA_integer_
)]
result_ifx_12251 <- fit_simple_cohort(ifx_12251, "GSE12251", "infliximab")

# ---------- GSE73661：Infliximab 与 Vedolizumab ----------
mixed_trial <- data.table::as.data.table(readRDS(file.path(DIR_INTERMEDIATE, "GSE73661_scored.rds")))
mixed_trial[, patient_id := trimws(as.character(get("study individual number:ch1")))]
mixed_trial[, regimen := as.character(get("induction therapy_maintenance therapy:ch1"))]
mixed_trial[, week := as.character(get("week (w):ch1"))]
mixed_trial[, response_marker := data.table::fcase(
  grepl("_UC NR ", title), "NR",
  grepl("_UC R ", title), "R",
  default = NA_character_
)]
mixed_trial[, week_number := as_numeric_safely(gsub("[^0-9]", "", week))]
mixed_trial[, key := paste(regimen, patient_id, sep = "__")]
response_map <- mixed_trial[!is.na(response_marker)][order(week_number), .SD[1L], by = key]
response_map <- response_map[, .(key, assigned_response = response_marker)]
mixed_trial <- merge(mixed_trial, response_map, by = "key", all.x = TRUE, sort = FALSE)
mixed_trial[, treatment := data.table::fcase(
  regimen == "IFX", "infliximab",
  grepl("vdz", regimen) & regimen != "plac_plac", "vedolizumab",
  default = NA_character_
)]
mixed_trial[, failure := data.table::fcase(
  assigned_response == "NR", 1L,
  assigned_response == "R", 0L,
  default = NA_integer_
)]
mixed_baseline <- mixed_trial[week == "W0" & !is.na(treatment) & !is.na(failure)]
result_ifx_73661 <- fit_simple_cohort(mixed_baseline[treatment == "infliximab"], "GSE73661_IFX", "infliximab")
result_vdz_73661 <- fit_simple_cohort(mixed_baseline[treatment == "vedolizumab"], "GSE73661_VDZ", "vedolizumab")

# 合并所有队列级结果并写出 CSV。
cohort_results <- data.table::rbindlist(list(
  result_unifi_active,
  result_unifi_placebo,
  result_glm_active,
  result_glm_placebo,
  result_ifx_16879,
  result_ifx_12251,
  result_ifx_73661,
  result_vdz_73661
), fill = TRUE)
interaction_results <- data.table::rbindlist(list(
  result_unifi_interaction,
  result_glm_interaction
), fill = TRUE)
adjusted_results <- data.table::rbindlist(list(result_unifi_adjusted, result_unifi_histology), fill = TRUE)
data.table::fwrite(cohort_results, file.path(DIR_RESULTS, "cohort_effects.csv"))
data.table::fwrite(interaction_results, file.path(DIR_RESULTS, "treatment_interactions.csv"))
data.table::fwrite(adjusted_results, file.path(DIR_RESULTS, "clinical_adjusted_effects.csv"))
data.table::fwrite(incremental_performance, file.path(DIR_RESULTS, "incremental_performance.csv"))

# 保存最小化的个体级作图表；不保存整张表达矩阵，文件很小。
baseline_samples <- data.table::rbindlist(list(
  unifi[treatment %in% c("ustekinumab", "placebo") & !is.na(failure),
        .(cohort = "GSE206285", treatment, failure, eri)],
  glm_trial[!is.na(failure),
            .(cohort = "GSE92415", treatment, failure, eri)],
  ifx_16879[!is.na(failure),
            .(cohort = "GSE16879", treatment, failure, eri)],
  ifx_12251[!is.na(failure),
            .(cohort = "GSE12251", treatment, failure, eri)],
  mixed_baseline[!is.na(failure),
                 .(cohort = ifelse(treatment == "infliximab", "GSE73661_IFX", "GSE73661_VDZ"),
                   treatment, failure, eri)]
), fill = TRUE)
baseline_samples[, eri_z := standardize_score(eri), by = .(cohort, treatment)]
data.table::fwrite(baseline_samples, file.path(DIR_RESULTS, "baseline_analysis_samples.csv"))

saveRDS(
  list(
    unifi_adjusted = model_unifi_adjusted,
    unifi_histology_adjusted = model_unifi_histology,
    unifi_clinical_only = model_unifi_clinical_only,
    unifi_interaction = model_unifi_interaction,
    golimumab_interaction = model_glm_interaction
  ),
  file.path(DIR_INTERMEDIATE, "fitted_models.rds")
)
print(cohort_results)
print(interaction_results)
print(adjusted_results)
print(incremental_performance)
