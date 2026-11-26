# Run from the project root after the overlap analysis (stage 26).
# The exploratory score averages the within-sample percentile ranks of the
# available Telesco panel genes. The published qPCR classifier uses a different
# measurement scale and algorithm.
library(data.table)
source("config.R")
source("scripts/helpers.R")
source("scripts/benchmark_helpers.R")
setDTthreads(1L)

out <- file.path(DIR_RESULTS, "review_telesco")
dir.create(out, showWarnings = FALSE)
id <- "GLM_MPS_13"
panel <- read_signature_registry()[signature_id == id]
stopifnot(nrow(panel) == 13L, !anyDuplicated(panel$gene),
          all(panel$component_weight == 1), all(panel$risk_orientation == 1))
selection <- fread(file.path(DIR_RESULTS, "overlap_reanalysis", "arm_selection.csv"))
selection <- selection[scenario == "primary" & retained == TRUE]
source_map <- fread(file.path(DIR_RESULTS, "overlap_reanalysis", "source_exclusion_map.csv"))
development <- strsplit(source_map[signature_id == id, excluded_arms], ";", fixed = TRUE)[[1]]
coverage <- fread(file.path(DIR_RESULTS, "benchmark_signature_coverage.csv"))[signature_id == id]
endpoints <- fread(file.path(DIR_RESULTS, "cbm_upgrade", "endpoint_harmonization_registry.csv"))

fits <- rbindlist(lapply(selection$cohort, function(arm) {
  samples <- as.data.table(readRDS(file.path(DIR_INTERMEDIATE, paste0("benchmark_", arm, "_scored.rds"))))
  if (arm == "GSE23597") {
    setorder(samples, geo_accession)
    stopifnot(identical(samples[duplicated(patient_id), geo_accession], "GSM578741"))
    samples <- samples[!duplicated(patient_id)]
  }
  stopifnot(all(complete.cases(samples[, c("failure", "treatment", id), with = FALSE])))
  rbindlist(lapply(sort(unique(samples$treatment)), function(therapy) {
    fit_all_signature_effects(samples, arm, therapy, panel)
  }))
}))
fits <- merge(fits, coverage[, .(cohort, genes_present, genes_total, genes_missing)], by = "cohort")
fits <- merge(fits, endpoints[, .(cohort, endpoint_definition)], by = "cohort")
fits[, development_trial := cohort %in% development]
active <- fits[treatment != "placebo"]
placebo <- fits[treatment == "placebo"]
stopifnot(nrow(active) == 6L, nrow(placebo) == 3L, sum(active$n) == sum(selection$n),
          all(active$n == selection$n[match(active$cohort, selection$cohort)]))

pools <- list(active_primary_arms = active,
              active_excluding_ACT1_and_PURSUIT = active[development_trial == FALSE],
              placebo = placebo)
pooled <- rbindlist(lapply(names(pools), function(name) {
  d <- pools[[name]]
  cbind(data.table(population = name, arms = nrow(d), n = sum(d$n), failures = sum(d$events),
                   cohorts = paste(d$cohort, collapse = ";")),
        random_effects_meta(d$log_or, d$se))
}))
stopifnot(all(pooled$failures > 0L & pooled$failures < pooled$n))

# Estimate the treatment contrast in the trial in which the panel was refined.
# The deposited endpoint is clinical response; Telesco et al. used mucosal healing.
pursuit <- as.data.table(readRDS(file.path(DIR_INTERMEDIATE, "benchmark_GSE92415_scored.rds")))
pursuit <- pursuit[complete.cases(failure, treatment, GLM_MPS_13, mayo, age)]
pursuit[, score_z := standardize_score(GLM_MPS_13)]
pursuit[, treatment := relevel(factor(treatment), ref = "placebo")]
model <- glm(failure ~ treatment * score_z + mayo + age, data = pursuit, family = binomial())
interaction <- fit_interaction_term(model, "treatmentgolimumab:score_z", id, "GSE92415")

# Compare against the historical benchmark everywhere the analysis population
# is comparable. The revised ACT-1 active arm has one fewer baseline array.
archived <- fread(file.path(DIR_RESULTS, "benchmark_cohort_effects.csv"))[signature_id == id]
audit <- merge(fits[!(cohort == "GSE23597" & treatment == "infliximab")], archived,
               by = c("cohort", "treatment"), suffixes = c("", "_archive"))
stopifnot(nrow(audit) == nrow(fits) - 1L, all(audit$n == audit$n_archive),
          max(abs(audit$log_or - audit$log_or_archive)) < 1e-8)
archived_int <- fread(file.path(DIR_RESULTS, "benchmark_treatment_interactions.csv"))[
  signature_id == id & cohort == "GSE92415"]
stopifnot(nrow(archived_int) == 1L, nrow(interaction) == 1L,
          abs(interaction$log_or - archived_int$log_or) < 1e-8)
fwrite(fits, file.path(out, "arm_effects.csv"))
fwrite(pooled, file.path(out, "pooled_effects.csv"))
fwrite(interaction, file.path(out, "pursuit_interaction.csv"))
