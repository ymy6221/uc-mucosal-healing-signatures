# 系统基准分析：随机效应 Meta、发现集排除、安慰剂效应、治疗交互、
# 严重度相关/残差化、匹配随机基因集和签名收敛轴。

check_benchmark_packages()
registry <- read_signature_registry()
signature_meta <- unique(registry[, .(
  signature_id, display_name, primary_set, source_doi, source_pmid,
  source_cohorts, definition_status, note
)])
treated_names <- c("ustekinumab", "golimumab", "infliximab", "vedolizumab", "etrolizumab")

effects <- data.table::fread(file.path(DIR_RESULTS, "benchmark_cohort_effects.csv"))
treated <- effects[treatment %in% treated_names]
placebo <- effects[treatment == "placebo"]

# ---------- 治疗队列和安慰剂队列的随机效应 Meta ----------
meta_by_signature <- function(effect_table, analysis_name) {
  output <- lapply(unique(effect_table$signature_id), function(sig_id) {
    one <- effect_table[signature_id == sig_id & is.finite(log_or) & is.finite(se) & se > 0]
    if (nrow(one) < 2L) return(NULL)
    pooled <- random_effects_meta(one$log_or, one$se)
    pooled[, `:=`(
      signature_id = sig_id,
      analysis = analysis_name,
      cohorts = nrow(one),
      total_n = sum(one$n),
      positive_cohorts = sum(one$log_or > 0)
    )]
    pooled
  })
  data.table::rbindlist(output, fill = TRUE)
}

active_meta <- meta_by_signature(treated, "all_treated")
placebo_meta <- meta_by_signature(placebo, "all_placebo")
active_meta[, fdr := stats::p.adjust(p_value, method = "BH")]
placebo_meta[, fdr := stats::p.adjust(p_value, method = "BH")]

# 只排除原文中用于发现或公开验证该签名的队列，防止把再分析当成新验证。
leave_source_out <- lapply(unique(treated$signature_id), function(sig_id) {
  source_string <- signature_meta[signature_id == sig_id, source_cohorts][1L]
  one <- treated[signature_id == sig_id & is.finite(log_or) & is.finite(se) & se > 0]
  base_cohort <- sub("_(IFX|VDZ)$", "", one$cohort)
  used_before <- vapply(base_cohort, function(x) grepl(x, source_string, fixed = TRUE), logical(1))
  external <- one[!used_before]
  if (nrow(external) < 2L) return(NULL)
  pooled <- random_effects_meta(external$log_or, external$se)
  pooled[, `:=`(
    signature_id = sig_id,
    analysis = "leave_source_cohorts_out",
    cohorts = nrow(external),
    total_n = sum(external$n),
    positive_cohorts = sum(external$log_or > 0),
    omitted_source_cohorts = source_string,
    included_cohorts = paste(external$cohort, collapse = ";")
  )]
  pooled
})
leave_source_meta <- data.table::rbindlist(leave_source_out, fill = TRUE)
leave_source_meta[, fdr := stats::p.adjust(p_value, method = "BH")]

data.table::fwrite(active_meta, file.path(DIR_RESULTS, "benchmark_active_meta.csv"))
data.table::fwrite(placebo_meta, file.path(DIR_RESULTS, "benchmark_placebo_meta.csv"))
data.table::fwrite(leave_source_meta, file.path(DIR_RESULTS, "benchmark_leave_source_out_meta.csv"))

# ---------- 10000 次匹配随机基因集经验 P 值 ----------
real_z <- data.table::fread(file.path(DIR_RESULTS, "benchmark_real_fisher_z.csv"))
treated_cohorts <- unique(real_z$cohort)
random_objects <- setNames(lapply(treated_cohorts, function(cohort_name) {
  readRDS(file.path(DIR_INTERMEDIATE, paste0("benchmark_", cohort_name, "_random_null.rds")))
}), treated_cohorts)

empirical_results <- lapply(unique(real_z$signature_id), function(sig_id) {
  available <- real_z[signature_id == sig_id & is.finite(fisher_z)]
  available <- available[vapply(cohort, function(x) {
    sig_id %in% colnames(random_objects[[x]]) &&
      any(is.finite(random_objects[[x]][, sig_id]))
  }, logical(1))]
  if (nrow(available) < 2L) return(NULL)
  real_stouffer <- sum(available$fisher_z) / sqrt(nrow(available))
  null_matrix <- do.call(cbind, lapply(available$cohort, function(x) random_objects[[x]][, sig_id]))
  null_stouffer <- rowSums(null_matrix) / sqrt(ncol(null_matrix))
  p_one_sided <- (1 + sum(null_stouffer >= real_stouffer)) / (length(null_stouffer) + 1)
  p_two_sided <- (1 + sum(abs(null_stouffer) >= abs(real_stouffer))) / (length(null_stouffer) + 1)

  source_string <- signature_meta[signature_id == sig_id, source_cohorts][1L]
  external_rows <- !vapply(sub("_(IFX|VDZ)$", "", available$cohort), function(x) {
    grepl(x, source_string, fixed = TRUE)
  }, logical(1))
  external_real <- NA_real_
  external_p <- NA_real_
  external_n <- sum(external_rows)
  if (external_n >= 2L) {
    external_real <- sum(available$fisher_z[external_rows]) / sqrt(external_n)
    external_null <- rowSums(null_matrix[, external_rows, drop = FALSE]) / sqrt(external_n)
    external_p <- (1 + sum(external_null >= external_real)) / (length(external_null) + 1)
  }
  data.table::data.table(
    signature_id = sig_id,
    cohorts = nrow(available),
    real_stouffer_z = real_stouffer,
    random_iterations = length(null_stouffer),
    empirical_p_one_sided = p_one_sided,
    empirical_p_two_sided = p_two_sided,
    external_cohorts = external_n,
    external_stouffer_z = external_real,
    external_empirical_p_one_sided = external_p
  )
})
empirical <- data.table::rbindlist(empirical_results, fill = TRUE)
empirical[, empirical_fdr := stats::p.adjust(empirical_p_one_sided, method = "BH")]
empirical[, external_empirical_fdr := stats::p.adjust(external_empirical_p_one_sided, method = "BH")]
data.table::fwrite(empirical, file.path(DIR_RESULTS, "matched_random_empirical_tests.csv"))
rm(random_objects)
invisible(gc())

# ---------- 三个 RCT 的 treatment×signature 交互 ----------
fit_trial_interactions <- function(cohort_name, active_name, adjustment) {
  d0 <- data.table::as.data.table(readRDS(
    file.path(DIR_INTERMEDIATE, paste0("benchmark_", cohort_name, "_scored.rds"))
  ))
  if ("crp" %in% names(d0)) d0[, log_crp := log1p(crp)]
  if ("calprotectin" %in% names(d0)) d0[, log_calprotectin := log1p(calprotectin)]
  signature_ids <- unique(registry$signature_id)
  results <- lapply(signature_ids, function(sig_id) {
    needed <- c("failure", "treatment", sig_id, adjustment)
    d <- d0[stats::complete.cases(d0[, ..needed])]
    if (nrow(d) < 20L || stats::sd(d[[sig_id]]) == 0) return(NULL)
    d[, score_z := standardize_score(get(sig_id))]
    d[, treatment := stats::relevel(factor(treatment), ref = "placebo")]
    if ("anti_tnf_history" %in% names(d)) d[, anti_tnf_history := factor(anti_tnf_history)]
    right_side <- paste(c("treatment * score_z", adjustment), collapse = " + ")
    model <- stats::glm(
      stats::as.formula(paste("failure ~", right_side)),
      family = stats::binomial(), data = d
    )
    term <- paste0("treatment", active_name, ":score_z")
    fit_interaction_term(model, term, sig_id, cohort_name)
  })
  data.table::rbindlist(results, fill = TRUE)
}

interaction_unifi <- fit_trial_interactions(
  "GSE206285", "ustekinumab",
  c("mayo", "log_crp", "log_calprotectin", "anti_tnf_history")
)
interaction_glm <- fit_trial_interactions("GSE92415", "golimumab", c("mayo", "age"))
interaction_ifx <- fit_trial_interactions("GSE23597", "infliximab", character())
interactions <- data.table::rbindlist(
  list(interaction_unifi, interaction_glm, interaction_ifx), fill = TRUE
)
interactions[, fdr_within_trial := stats::p.adjust(p_value, method = "BH"), by = cohort]
data.table::fwrite(interactions, file.path(DIR_RESULTS, "benchmark_treatment_interactions.csv"))

interaction_meta <- lapply(unique(interactions$signature_id), function(sig_id) {
  one <- interactions[signature_id == sig_id & is.finite(log_or) & is.finite(se) & se > 0]
  if (nrow(one) < 2L) return(NULL)
  pooled <- random_effects_meta(one$log_or, one$se)
  pooled[, `:=`(signature_id = sig_id, trials = nrow(one))]
  pooled
})
interaction_meta <- data.table::rbindlist(interaction_meta, fill = TRUE)
interaction_meta[, fdr := stats::p.adjust(p_value, method = "BH")]
data.table::fwrite(interaction_meta, file.path(DIR_RESULTS, "benchmark_interaction_meta.csv"))

# ---------- UNIFI 严重度相关和残差化 ----------
unifi <- data.table::as.data.table(readRDS(
  file.path(DIR_INTERMEDIATE, "benchmark_GSE206285_scored.rds")
))
unifi[, `:=`(log_crp = log1p(crp), log_calprotectin = log1p(calprotectin))]
severity_markers <- c("mayo", "log_crp", "log_calprotectin",
                      "neutrophils_lp", "neutrophils_epithelium")
severity_results <- list()
severity_index <- 0L
for (sig_id in unique(registry$signature_id)) {
  for (marker in severity_markers) {
    d <- unifi[stats::complete.cases(get(sig_id), get(marker))]
    test <- suppressWarnings(stats::cor.test(d[[sig_id]], d[[marker]],
                                             method = "spearman", exact = FALSE))
    severity_index <- severity_index + 1L
    severity_results[[severity_index]] <- data.table::data.table(
      signature_id = sig_id, marker = marker, n = nrow(d),
      spearman_rho = unname(test$estimate), p_value = test$p.value
    )
  }
}
severity <- data.table::rbindlist(severity_results)
severity[, fdr_global := stats::p.adjust(p_value, method = "BH")]
severity[, fdr_within_signature := stats::p.adjust(p_value, method = "BH"), by = signature_id]
data.table::fwrite(severity, file.path(DIR_RESULTS, "benchmark_severity_correlations.csv"))

residual_results <- list()
residual_index <- 0L
for (sig_id in unique(registry$signature_id)) {
  for (model_name in c("clinical_residual", "histology_residual")) {
    adjustment <- c("mayo", "log_crp", "log_calprotectin", "anti_tnf_history")
    if (model_name == "histology_residual") {
      adjustment <- c(adjustment, "neutrophils_lp", "neutrophils_epithelium")
    }
    needed <- c("failure", sig_id, adjustment)
    d <- unifi[treatment == "ustekinumab" & stats::complete.cases(unifi[, ..needed])]
    d[, score_z := standardize_score(get(sig_id))]
    d[, anti_tnf_history := factor(anti_tnf_history)]
    score_model <- stats::lm(
      stats::as.formula(paste("score_z ~", paste(adjustment, collapse = " + "))),
      data = d
    )
    d[, residual_z := standardize_score(stats::residuals(score_model))]
    outcome_model <- stats::glm(failure ~ residual_z, family = stats::binomial(), data = d)
    effect <- extract_glm_term(
      outcome_model, "residual_z", "GSE206285", "ustekinumab", model_name
    )
    residual_index <- residual_index + 1L
    residual_results[[residual_index]] <- cbind(
      data.table::data.table(signature_id = sig_id), effect
    )
  }
}
residual_effects <- data.table::rbindlist(residual_results, fill = TRUE)
residual_effects[, fdr := stats::p.adjust(p_value, method = "BH"), by = model]
data.table::fwrite(residual_effects, file.path(DIR_RESULTS, "benchmark_severity_residual_effects.csv"))

# GSE72819 抗 TNF 初治亚组是原文生物标志物分析中的关键人群，因此预先作为敏感性分析。
etro <- data.table::as.data.table(readRDS(
  file.path(DIR_INTERMEDIATE, "benchmark_GSE72819_scored.rds")
))
etro_naive <- etro[prior_anti_tnf == "Naive"]
etro_naive_effects <- fit_all_signature_effects(
  etro_naive, "GSE72819_TNFnaive", "etrolizumab", registry
)
data.table::fwrite(etro_naive_effects, file.path(DIR_RESULTS, "GSE72819_TNFnaive_sensitivity.csv"))

# ---------- 无结局监督的签名收敛轴 ----------
scored_files <- file.path(
  DIR_INTERMEDIATE,
  paste0("benchmark_", c("GSE206285", "GSE92415", "GSE16879", "GSE12251",
                         "GSE73661_IFX", "GSE73661_VDZ", "GSE23597", "GSE72819"),
         "_scored.rds")
)
primary_ids <- signature_meta[primary_set == TRUE, signature_id]
stacked <- data.table::rbindlist(lapply(scored_files, function(path) {
  d <- data.table::as.data.table(readRDS(path))
  cohort_name <- sub("^benchmark_|_scored\\.rds$", "", basename(path))
  d[, cohort := cohort_name]
  d[, c("cohort", "geo_accession", primary_ids), with = FALSE]
}), fill = TRUE)

for (sig_id in primary_ids) {
  stacked[, (sig_id) := {
    x <- get(sig_id)
    if (sum(is.finite(x)) < 2L || stats::sd(x, na.rm = TRUE) == 0) rep(NA_real_, .N)
    else as.numeric(scale(x))
  }, by = cohort]
}
pca_matrix <- as.matrix(stacked[, ..primary_ids])
pca_matrix[!is.finite(pca_matrix)] <- 0
pca <- stats::prcomp(pca_matrix, center = TRUE, scale. = TRUE)
if (sum(pca$rotation[, 1L]) < 0) {
  pca$rotation[, 1L] <- -pca$rotation[, 1L]
  pca$x[, 1L] <- -pca$x[, 1L]
}
variance <- pca$sdev^2 / sum(pca$sdev^2)
pca_loadings <- data.table::data.table(
  signature_id = rownames(pca$rotation),
  pc1_loading = pca$rotation[, 1L],
  pc2_loading = pca$rotation[, 2L]
)
pca_variance <- data.table::data.table(
  component = paste0("PC", seq_along(variance)), variance_explained = variance
)
pca_scores <- data.table::data.table(
  cohort = stacked$cohort,
  geo_accession = stacked$geo_accession,
  pc1 = pca$x[, 1L], pc2 = pca$x[, 2L]
)
correlation <- stats::cor(pca_matrix, method = "spearman", use = "pairwise.complete.obs")
correlation_long <- data.table::as.data.table(as.table(correlation))
data.table::setnames(correlation_long, c("signature_1", "signature_2", "spearman_rho"))
data.table::fwrite(pca_loadings, file.path(DIR_RESULTS, "signature_convergence_loadings.csv"))
data.table::fwrite(pca_variance, file.path(DIR_RESULTS, "signature_convergence_variance.csv"))
data.table::fwrite(pca_scores, file.path(DIR_RESULTS, "signature_convergence_sample_scores.csv"))
data.table::fwrite(correlation_long, file.path(DIR_RESULTS, "signature_score_correlations.csv"))

# ---------- 证据级别表和事先定义的高影响力主张门槛 ----------
target_trial <- data.table::data.table(
  signature_id = c("IFX_ARJS_5", "GLM_MPS_13", "OSM_OSMR",
                   "ETRO_2", "IL22_NEUTROPHIL_6", "VDZ_4", "VDZ_NR_3", "VDZ_R_5"),
  target_cohort = c("GSE23597", "GSE92415", "GSE23597",
                    NA, "GSE206285", NA, NA, NA)
)
evidence <- merge(signature_meta[primary_set == TRUE], active_meta, by = "signature_id", all.x = TRUE)
evidence <- merge(evidence, leave_source_meta[, .(
  signature_id, external_or = or, external_ci_low = ci_low,
  external_ci_high = ci_high, external_p = p_value, external_fdr = fdr,
  external_cohorts = cohorts
)], by = "signature_id", all.x = TRUE)
evidence <- merge(evidence, placebo_meta[, .(
  signature_id, placebo_or = or, placebo_ci_low = ci_low,
  placebo_ci_high = ci_high, placebo_p = p_value, placebo_fdr = fdr
)], by = "signature_id", all.x = TRUE)
evidence <- merge(evidence, empirical[, .(
  signature_id, empirical_p_one_sided, empirical_fdr,
  external_empirical_p_one_sided, external_empirical_fdr
)], by = "signature_id", all.x = TRUE)
evidence <- merge(evidence, target_trial, by = "signature_id", all.x = TRUE)

target_interactions <- interactions[target_trial, on = .(signature_id, cohort = target_cohort), nomatch = 0]
target_interactions <- target_interactions[, .(
  signature_id, target_interaction_or = or,
  target_interaction_p = p_value,
  target_interaction_fdr = fdr_within_trial
)]
evidence <- merge(evidence, target_interactions, by = "signature_id", all.x = TRUE)
severity_min <- severity[, .(
  severity_min_fdr = min(fdr_within_signature, na.rm = TRUE),
  severity_max_abs_rho = max(abs(spearman_rho), na.rm = TRUE)
), by = signature_id]
evidence <- merge(evidence, severity_min, by = "signature_id", all.x = TRUE)

evidence[, evidence_class := data.table::fcase(
  !is.na(target_interaction_fdr) & target_interaction_fdr < 0.05,
  "drug_specific_interaction_supported",
  fdr < 0.05 & empirical_fdr < 0.05 &
    ((placebo_fdr < 0.05 & placebo_or > 1) | severity_min_fdr < 0.05),
  "general_healing_propensity_supported",
  fdr < 0.05 & empirical_fdr < 0.05,
  "treated_association_without_specificity",
  default = "insufficient_or_inconsistent"
)]
data.table::fwrite(evidence, file.path(DIR_RESULTS, "published_signature_evidence_matrix.csv"))

sample_counts <- data.table::fread(file.path(DIR_RESULTS, "benchmark_sample_counts.csv"))
primary_n <- nrow(evidence)
active_supported_n <- evidence[fdr < 0.05 & empirical_fdr < 0.05, .N]
general_n <- evidence[evidence_class == "general_healing_propensity_supported", .N]
drug_specific_n <- evidence[evidence_class == "drug_specific_interaction_supported", .N]
pc1_variance <- pca_variance[component == "PC1", variance_explained]
claim_gate <- data.table::data.table(
  gate = c(
    "At least 7 exact primary published modules",
    "At least 4 primary modules beat matched random controls in treated cohorts",
    "At least 3 primary modules show general healing-propensity evidence",
    "No more than 1 target-trial interaction survives FDR",
    "Unsupervised primary-signature PC1 explains at least 35 percent",
    "Eight active cohorts and three placebo arms are present"
  ),
  observed = c(
    primary_n, active_supported_n, general_n, drug_specific_n,
    pc1_variance,
    paste0(sum(sample_counts$treatment %in% treated_names), "/",
           sum(sample_counts$treatment == "placebo"))
  ),
  passed = c(
    primary_n >= 7L,
    active_supported_n >= 4L,
    general_n >= 3L,
    drug_specific_n <= 1L,
    pc1_variance >= 0.35,
    sum(sample_counts$treatment %in% treated_names) == 8L &
      sum(sample_counts$treatment == "placebo") == 3L
  )
)
claim_gate[, overall_status := ifelse(all(passed),
                                      "READY_FOR_FIXED_HIGH_IMPACT_TARGET",
                                      "DO_NOT_SUBMIT_FIXED_TARGET_NOT_YET_MET")]
data.table::fwrite(claim_gate, file.path(DIR_RESULTS, "high_impact_claim_gate.csv"))
print(claim_gate)
