# Run from the project root after stage 5. The unit is the treated UNIFI participant.
# The end point is failure to achieve mucosal healing at week 8 (1 = failure).
suppressPackageStartupMessages(library(data.table))
out <- "results/review_clinical_auc"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
scores <- as.data.table(readRDS("data/intermediate/benchmark_GSE206285_scored.rds"))
patients <- scores[treatment == "ustekinumab"][order(geo_accession)]
stopifnot(nrow(patients) == 358L, !anyDuplicated(patients[["donor id:ch1"]]))
patients[, `:=`(log_crp = log1p(crp), log_calprotectin = log1p(calprotectin),
                anti_tnf_history = factor(anti_tnf_history))]
ids <- fread("signatures/published_signature_registry.csv")[primary_set == TRUE, unique(signature_id)]
stopifnot(length(ids) == 7L)
clinical <- c("mayo", "log_crp", "log_calprotectin", "anti_tnf_history")
histology <- c(clinical, "neutrophils_lp", "neutrophils_epithelium")
adjustments <- list(clinical_residual = clinical, histology_residual = histology)
complete_set <- function(vars) {
  patients[complete.cases(patients[, c("failure", ids, vars), with = FALSE])]
}

# Refit both stages of the original residual analysis using each complete-case set.
residuals_out <- rbindlist(lapply(names(adjustments), function(model_name) {
  covariates <- adjustments[[model_name]]
  d <- complete_set(covariates)
  rbindlist(lapply(ids, function(id) {
    d[, signature_z := as.numeric(scale(get(id)))]
    severity <- lm(reformulate(covariates, "signature_z"), data = d)
    d[, residual_z := as.numeric(scale(resid(severity)))]
    model <- glm(failure ~ residual_z, family = binomial(), data = d)
    term <- coef(summary(model))["residual_z", ]
    interval <- exp(term[1] + c(-1, 1) * qnorm(.975) * term[2])
    data.table(signature_id = id, model = model_name, n = nobs(model),
               failures = sum(d$failure), or = exp(term[1]), ci_low = interval[1],
               ci_high = interval[2], p_value = term[4])
  }))
}))
# Keep the original eleven-module FDR for provenance and state the family used in the revision.
historical <- fread("results/benchmark_severity_residual_effects.csv")
stopifnot(all(historical[, uniqueN(signature_id), by = model]$V1 == 11L))
residuals_out <- merge(residuals_out, historical[signature_id %in% ids,
      .(signature_id, model, historical_or = or, historical_p = p_value, historical_fdr_11 = fdr)],
      by = c("signature_id", "model"))
stopifnot(nrow(residuals_out) == 14L, max(abs(residuals_out$or - residuals_out$historical_or)) < 1e-8,
          max(abs(residuals_out$p_value - residuals_out$historical_p)) < 1e-8)
residuals_out[, fdr := p.adjust(p_value, "BH"), by = model]
residuals_out[, c("historical_or", "historical_p") := NULL]
fwrite(residuals_out, file.path(out, "residual_effects.csv"))

# Compare independently fitted nested models. Treat in-sample AUC as apparent performance.
# Paired DeLong inference is performed only with fixed coefficients on held-out participants,
# not on the data used to estimate the models (Demler et al., doi:10.1002/sim.5328).
dat <- complete_set(clinical)
stopifnot(nrow(dat) == 330L, sum(dat$failure) == 281L)
set.seed(20260826L)
val_idx <- unlist(lapply(split(seq_len(nrow(dat)), dat$failure), function(ix)
  sample(ix, size = round(length(ix) / 3), replace = FALSE)), use.names = FALSE)
train <- dat[-val_idx]
validation <- dat[val_idx]
stopifnot(all(table(validation$failure) > 0), length(intersect(train$geo_accession, validation$geo_accession)) == 0L)
make_roc <- function(fit, newdata) {
  pROC::roc(response = newdata$failure, predictor = as.numeric(predict(fit, newdata = newdata, type = "link")),
            levels = c(0, 1), direction = "<", quiet = TRUE)
}
base_formula <- reformulate(clinical, "failure")
base_full <- glm(base_formula, family = binomial(), data = dat)
base_train <- glm(base_formula, family = binomial(), data = train)
base_apparent_roc <- make_roc(base_full, dat)
base_test_roc <- make_roc(base_train, validation)
auc_out <- rbindlist(lapply(ids, function(id) {
  expanded_formula <- reformulate(c(clinical, id), "failure")
  expanded_full <- glm(expanded_formula, family = binomial(), data = dat)
  expanded_train <- glm(expanded_formula, family = binomial(), data = train)
  stopifnot(base_full$converged, base_train$converged, expanded_full$converged, expanded_train$converged,
            nobs(base_full) == nobs(expanded_full), nobs(base_train) == nobs(expanded_train))
  lrt <- anova(base_full, expanded_full, test = "LRT")[["Pr(>Chi)"]][2]
  expanded_apparent_roc <- make_roc(expanded_full, dat)
  expanded_test_roc <- make_roc(expanded_train, validation)
  test <- pROC::roc.test(expanded_test_roc, base_test_roc, paired = TRUE, method = "delong")
  data.table(signature_id = id, n_complete = nrow(dat), complete_failures = sum(dat$failure),
             auc_clinical_apparent = as.numeric(base_apparent_roc$auc),
             auc_plus_apparent = as.numeric(expanded_apparent_roc$auc),
             delta_auc_apparent = as.numeric(expanded_apparent_roc$auc) - as.numeric(base_apparent_roc$auc),
             likelihood_ratio_p = lrt, n_training = nrow(train),
             training_failures = sum(train$failure), n_holdout = nrow(validation),
             holdout_failures = sum(validation$failure), auc_clinical_holdout = as.numeric(base_test_roc$auc),
             auc_plus_holdout = as.numeric(expanded_test_roc$auc),
             delta_auc_holdout = as.numeric(expanded_test_roc$auc) - as.numeric(base_test_roc$auc),
             holdout_ci_low = as.numeric(test$conf.int[1]), holdout_ci_high = as.numeric(test$conf.int[2]),
             paired_delong_p = test$p.value)
}))
auc_out[, `:=`(paired_delong_fdr = p.adjust(paired_delong_p, "BH"),
               likelihood_ratio_fdr = p.adjust(likelihood_ratio_p, "BH"))]
stopifnot(all(auc_out$holdout_ci_low <= auc_out$delta_auc_holdout),
          all(auc_out$holdout_ci_high >= auc_out$delta_auc_holdout))
fwrite(auc_out, file.path(out, "incremental_auc.csv"))
