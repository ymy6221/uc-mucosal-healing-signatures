# Run after scripts/29_unifi_residual_auc.R from the project root.
# Internally validate eight fixed models on identical folds. The training fold
# supplies both the most variable probe per gene and all logistic coefficients.
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
source("scripts/helpers.R")
DIR_RAW <- "data/raw"
out <- "results/review_repeated_cv"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

registry <- fread("signatures/published_signature_registry.csv")[primary_set == TRUE]
ids <- unique(registry$signature_id)
clinical <- c("mayo", "log_crp", "log_calprotectin", "anti_tnf_history")
dat <- as.data.table(readRDS("data/intermediate/benchmark_GSE206285_scored.rds"))
dat <- dat[treatment == "ustekinumab"][order(geo_accession)]
dat[, `:=`(log_crp = log1p(crp), log_calprotectin = log1p(calprotectin),
           anti_tnf_history = factor(anti_tnf_history))]
dat <- dat[complete.cases(dat[, c("failure", clinical, ids), with = FALSE])]
stopifnot(nrow(dat) == 330L, sum(dat$failure) == 281L,
          !anyDuplicated(dat[["donor id:ch1"]]), length(ids) == 7L)
forms <- lapply(c(list(clinical), lapply(ids, function(id) c(clinical, id))),
                reformulate, response = "failure")
names(forms) <- c("clinical", ids)
fit_model <- function(formula, data) {
  fit <- glm(formula, data = data, family = binomial())
  stopifnot(fit$converged, nobs(fit) == nrow(data), all(is.finite(coef(fit))))
  fit
}

# Refit the full-data stage-29 models on their original, locked signature scores.
full <- lapply(forms, fit_model, data = dat)
lrt <- rbindlist(lapply(ids, function(id) {
  test <- anova(full$clinical, full[[id]], test = "LRT")
  data.table(signature_id = id, n = nrow(dat), failures = sum(dat$failure),
             df = test$Df[2L], likelihood_ratio = test$Deviance[2L],
             p_value = test[["Pr(>Chi)"]][2L])
}))
lrt[, fdr := p.adjust(p_value, "BH")]
prior <- fread("results/review_clinical_auc/incremental_auc.csv")[match(ids, signature_id)]
stopifnot(identical(prior$signature_id, lrt$signature_id), all(lrt$df == 1L),
          max(abs(prior$likelihood_ratio_p - lrt$p_value)) < 1e-12,
          max(abs(prior$likelihood_ratio_fdr - lrt$fdr)) < 1e-12)

# Read only the 330 participants' columns from the public processed GEO matrix.
annotation <- read_geo_annotation("GPL13158")
raw <- fread(file.path(DIR_RAW, "GSE206285_series_matrix.txt.gz"), skip = '"ID_REF"',
             select = c("ID_REF", dat$geo_accession), showProgress = FALSE)
annotation <- annotation[probe %in% raw$ID_REF]
expression <- as.matrix(raw[match(annotation$probe, ID_REF),
                            dat$geo_accession, with = FALSE])
rm(raw)
invisible(gc())
stopifnot(!anyDuplicated(annotation$probe), !anyNA(expression),
          identical(colnames(expression), dat$geo_accession))
symbols <- unique(annotation$symbol)
symbol_index <- match(annotation$symbol, symbols)
multis <- duplicated(symbol_index) | duplicated(symbol_index, fromLast = TRUE)
multi_rows <- which(multis)
registry_by_id <- lapply(ids, function(id) {
  definition <- registry[signature_id == id]
  definition <- definition[gene %in% symbols]
  stopifnot(nrow(definition) >= if (nrow(registry[signature_id == id]) <= 3L)
    nrow(registry[signature_id == id]) else ceiling(.6 * nrow(registry[signature_id == id])))
  list(rows = match(definition$gene, symbols), weight = definition$component_weight,
       direction = unique(definition$risk_orientation))
})

# Type-7 quartiles, identical to stats::IQR for complete expression values.
# Restrict order-statistic work to duplicated genes, without copying all probes.
train_quartile <- function(train_rows, prob) {
  h <- 1 + (length(train_rows) - 1) * prob
  low <- matrixStats::rowOrderStats(expression, rows = multi_rows, cols = train_rows,
                                    which = as.integer(floor(h)))
  if (h == floor(h)) return(low)
  high <- matrixStats::rowOrderStats(expression, rows = multi_rows, cols = train_rows,
                                     which = as.integer(ceiling(h)))
  fraction <- h - floor(h)
  (1 - fraction) * low + fraction * high
}
score_fold <- function(train_rows) {
  variability <- numeric(nrow(annotation))
  variability[multi_rows] <- train_quartile(train_rows, .75) - train_quartile(train_rows, .25)
  candidates <- order(symbol_index, -variability, seq_along(variability), method = "radix")
  selected <- candidates[!duplicated(symbol_index[candidates])]
  stopifnot(identical(symbol_index[selected], seq_along(symbols)))
  ranks <- matrixStats::colRanks(expression, rows = selected, ties.method = "average",
                                preserveShape = TRUE, useNames = FALSE)
  vapply(registry_by_id, function(spec) {
    ranks_subset <- ranks[spec$rows, , drop = FALSE] / (length(symbols) + 1)
    if (any(spec$weight < 0)) {
      score <- colMeans(ranks_subset[spec$weight > 0, , drop = FALSE]) -
        colMeans(ranks_subset[spec$weight < 0, , drop = FALSE])
    } else {
      score <- colMeans(ranks_subset)
    }
    score * spec$direction
  }, numeric(nrow(dat)))
}
auc <- function(y, predictions) {
  cases <- sum(y == 1L)
  controls <- sum(y == 0L)
  (sum(rank(predictions, ties.method = "average")[y == 1L]) -
    cases * (cases + 1) / 2) / (cases * controls)
}

# In each repeat, every patient has exactly one out-of-fold prediction per
# model. Compute one paired AUC difference per repeat across its 330 patients.
repetitions <- 50L
folds <- 5L
set.seed(20261126L)
rows <- seq_len(nrow(dat))
strata <- split(rows, dat$failure)
repeat_auc <- vector("list", repetitions)
for (rep_id in seq_len(repetitions)) {
  shuffled <- unlist(lapply(strata, sample), use.names = FALSE)
  fold_id <- integer(nrow(dat))
  fold_id[shuffled] <- rep_len(sample.int(folds), nrow(dat))
  fold_counts <- table(fold_id, dat$failure)
  stopifnot(all(rowSums(fold_counts) == 66L), all(fold_counts[, "0"] %in% c(9L, 10L)))
  predictions <- matrix(NA_real_, nrow(dat), length(forms), dimnames = list(NULL, names(forms)))
  heldout_scores <- matrix(NA_real_, nrow(dat), length(ids), dimnames = list(NULL, ids))
  for (fold in seq_len(folds)) {
    test_rows <- which(fold_id == fold)
    train_rows <- which(fold_id != fold)
    fold_scores <- score_fold(train_rows)
    heldout_scores[test_rows, ] <- fold_scores[test_rows, , drop = FALSE]
    fold_data <- copy(dat[, c("failure", clinical), with = FALSE])
    for (j in seq_along(ids)) set(fold_data, j = ids[j], value = fold_scores[, j])
    train <- fold_data[train_rows]
    test <- fold_data[test_rows]
    for (model_name in names(forms)) {
      model <- fit_model(forms[[model_name]], train)
      predictions[test_rows, model_name] <- predict(model, newdata = test, type = "link")
    }
  }
  stopifnot(all(is.finite(predictions)), all(is.finite(heldout_scores)))
  aucs <- apply(predictions, 2L, function(x) auc(dat$failure, x))
  repeat_auc[[rep_id]] <- data.table(
    repeat_id = rep_id, signature_id = ids, auc_clinical = unname(aucs["clinical"]),
    auc_plus = unname(aucs[ids]), delta_auc = unname(aucs[ids] - aucs["clinical"]),
    score_correlation_with_locked = vapply(ids, function(id)
      cor(heldout_scores[, id], dat[[id]]), numeric(1)),
    mean_abs_score_diff_from_locked = vapply(ids, function(id)
      mean(abs(heldout_scores[, id] - dat[[id]])), numeric(1)))
}
details <- rbindlist(repeat_auc)
# Percentiles describe sensitivity to the random folds in this one cohort.
# Overlapping repeats are not independent studies and these are not confidence limits.
summary <- details[, .(
  n = nrow(dat), failures = sum(dat$failure), healings = sum(dat$failure == 0L),
  repeats = .N, folds = folds, mean_auc_clinical = mean(auc_clinical),
  mean_auc_plus = mean(auc_plus), mean_delta_auc = mean(delta_auc),
  delta_split_p025 = quantile(delta_auc, .025, names = FALSE),
  delta_split_p975 = quantile(delta_auc, .975, names = FALSE),
  mean_score_correlation_with_locked = mean(score_correlation_with_locked),
  mean_abs_score_diff_from_locked = mean(mean_abs_score_diff_from_locked)
), by = signature_id]
stopifnot(nrow(details) == 350L, all(summary$repeats == repetitions),
          all(abs(summary$mean_auc_plus - summary$mean_auc_clinical - summary$mean_delta_auc) < 1e-12))
fwrite(lrt, file.path(out, "full_lrt.csv"))
fwrite(details, file.path(out, "auc_by_repeat.csv"))
fwrite(summary, file.path(out, "repeated_cv_summary.csv"))
