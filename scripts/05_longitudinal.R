# 验证 ERI 是否随黏膜改善而下降。所有检验均按患者配对，不把重复活检当独立样本。

check_packages()

# 一个通用函数：按 response 分层做配对 Wilcoxon，并比较两组 delta。
summarize_pairs <- function(pair_table, cohort_name, treatment_name) {
  pair_table <- pair_table[stats::complete.cases(eri_baseline, eri_followup, response)]
  pair_table[, delta := eri_followup - eri_baseline]
  output <- pair_table[, {
    test <- stats::wilcox.test(eri_followup, eri_baseline, paired = TRUE, exact = FALSE)
    list(
      n = .N,
      median_baseline = stats::median(eri_baseline),
      median_followup = stats::median(eri_followup),
      median_delta = stats::median(delta),
      paired_p = test$p.value
    )
  }, by = response]
  between_p <- if (length(unique(pair_table$response)) == 2L) {
    stats::wilcox.test(delta ~ response, data = pair_table, exact = FALSE)$p.value
  } else {
    NA_real_
  }
  output[, `:=`(cohort = cohort_name, treatment = treatment_name, delta_between_response_p = between_p)]
  list(summary = output, pairs = pair_table)
}

# ---------- GSE16879：IFX 治疗前后 ----------
ifx <- data.table::as.data.table(readRDS(file.path(DIR_INTERMEDIATE, "GSE16879_scored.rds")))
ifx <- ifx[get("disease:ch1") == "UC"]
ifx[, patient_id := sub("_(beforeT|afterT)$", "", title)]
ifx[, timepoint := ifelse(grepl("beforeT$", title), "baseline", "followup")]
ifx[, response := ifelse(get("response to infliximab:ch1") == "Yes", "R", "NR")]
ifx_pairs <- data.table::dcast(
  ifx[, .(patient_id, response, timepoint, eri)],
  patient_id + response ~ timepoint,
  value.var = "eri"
)
data.table::setnames(ifx_pairs, c("baseline", "followup"), c("eri_baseline", "eri_followup"))
ifx_summary <- summarize_pairs(ifx_pairs, "GSE16879", "infliximab")

# ---------- GSE92415：Golimumab/安慰剂 W0 到 W6 ----------
glm_trial <- data.table::as.data.table(readRDS(file.path(DIR_INTERMEDIATE, "GSE92415_scored.rds")))
glm_trial <- glm_trial[get("disease:ch1") == "Ulcerative Colitis (UC)" &
                         get("visit:ch1") %in% c("Week 0", "Week 6")]
glm_trial[, patient_id := as.character(get("subject:ch1"))]
glm_trial[, treatment := tolower(get("treatment:ch1"))]
glm_trial[, response := ifelse(get("wk6response:ch1") == "Yes", "R", "NR")]
glm_trial[, timepoint := ifelse(get("visit:ch1") == "Week 0", "baseline", "followup")]
glm_pairs <- data.table::dcast(
  glm_trial[, .(patient_id, treatment, response, timepoint, eri)],
  patient_id + treatment + response ~ timepoint,
  value.var = "eri"
)
data.table::setnames(glm_pairs, c("baseline", "followup"), c("eri_baseline", "eri_followup"))
glm_active_summary <- summarize_pairs(glm_pairs[treatment == "golimumab"], "GSE92415", "golimumab")
glm_placebo_summary <- summarize_pairs(glm_pairs[treatment == "placebo"], "GSE92415", "placebo")

# ---------- GSE73661：IFX 或 VDZ 的基线到首次有结局标记的活检 ----------
mixed <- data.table::as.data.table(readRDS(file.path(DIR_INTERMEDIATE, "GSE73661_scored.rds")))
mixed[, patient_id := trimws(as.character(get("study individual number:ch1")))]
mixed[, regimen := as.character(get("induction therapy_maintenance therapy:ch1"))]
mixed[, week := as.character(get("week (w):ch1"))]
mixed[, week_number := as_numeric_safely(gsub("[^0-9]", "", week))]
mixed[, response := data.table::fcase(
  grepl("_UC NR ", title), "NR",
  grepl("_UC R ", title), "R",
  default = NA_character_
)]
mixed[, treatment := data.table::fcase(
  regimen == "IFX", "infliximab",
  grepl("vdz", regimen) & regimen != "plac_plac", "vedolizumab",
  default = NA_character_
)]
mixed[, key := paste(regimen, patient_id, sep = "__")]
mixed_baseline <- mixed[week == "W0" & !is.na(treatment), .(key, patient_id, treatment, eri_baseline = eri)]
mixed_followup <- mixed[!is.na(response) & !is.na(treatment)][order(week_number), .SD[1L], by = key]
mixed_followup <- mixed_followup[, .(key, response, followup_week = week_number, eri_followup = eri)]
mixed_pairs <- merge(mixed_baseline, mixed_followup, by = "key", all = FALSE)
mixed_ifx_summary <- summarize_pairs(mixed_pairs[treatment == "infliximab"], "GSE73661_IFX", "infliximab")
mixed_vdz_summary <- summarize_pairs(mixed_pairs[treatment == "vedolizumab"], "GSE73661_VDZ", "vedolizumab")

# 合并统计摘要和配对个体数据。
longitudinal_summary <- data.table::rbindlist(list(
  ifx_summary$summary,
  glm_active_summary$summary,
  glm_placebo_summary$summary,
  mixed_ifx_summary$summary,
  mixed_vdz_summary$summary
), fill = TRUE)
longitudinal_pairs <- data.table::rbindlist(list(
  data.table::copy(ifx_summary$pairs)[, `:=`(cohort = "GSE16879", treatment = "infliximab")],
  data.table::copy(glm_active_summary$pairs)[, cohort := "GSE92415"],
  data.table::copy(glm_placebo_summary$pairs)[, cohort := "GSE92415"],
  data.table::copy(mixed_ifx_summary$pairs)[, cohort := "GSE73661_IFX"],
  data.table::copy(mixed_vdz_summary$pairs)[, cohort := "GSE73661_VDZ"]
), fill = TRUE, use.names = TRUE)
data.table::fwrite(longitudinal_summary, file.path(DIR_RESULTS, "longitudinal_summary.csv"))
data.table::fwrite(longitudinal_pairs, file.path(DIR_RESULTS, "longitudinal_pairs.csv"))
print(longitudinal_summary)

