# Run from the R project root after stages 5–6. Input files are read only.
library(data.table)
source("config.R")
source("scripts/helpers.R")
source("scripts/benchmark_helpers.R")
setDTthreads(1)
out <- file.path(DIR_RESULTS, "overlap_reanalysis")
dir.create(out, showWarnings = FALSE)
save_result <- function(x, file) fwrite(x, file.path(out, paste0(file, ".csv")))
load_result <- function(file) fread(file.path(DIR_RESULTS, paste0(file, ".csv")))

genes <- read_signature_registry()
registry <- unique(genes[,
    .(signature_id, display_name, primary_set, source_doi, source_cohorts)])
primary_ids <- registry[primary_set == TRUE, signature_id]
effects <- load_result("benchmark_cohort_effects")[signature_id %in% primary_ids]
# Two baseline arrays for ACT-1 participant P13. Keep the smaller accession,
# a deterministic choice independent of score and response.
act_object <- prepare_source_cohorts("GSE23597")[[1L]]
act_pheno <- copy(act_object$phenotype)[order(geo_accession)]
stopifnot(identical(act_pheno[duplicated(patient_id), geo_accession], "GSM578741"))
act_pheno <- act_pheno[!duplicated(patient_id)]
act_scored <- score_registered_signatures(act_object$expression[,act_pheno$geo_accession,drop=FALSE], genes)
act_pheno <- merge(act_pheno, act_scored$scores, by="geo_accession", sort=FALSE)
act_active <- act_pheno[treatment=="infliximab"]
act_null <- matched_random_null(act_scored$ranks[,act_active$geo_accession,drop=FALSE],
                                act_active$failure, genes, n_random=10000L, seed_offset=601L)
act_z <- rbindlist(lapply(primary_ids, function(id)
  data.table(signature_id=id, cohort="GSE23597", treatment="infliximab",
    n=nrow(act_active),failures=sum(act_active$failure),
    fisher_z=score_matrix_to_fisher_z(matrix(act_active[[id]], nrow=1L),act_active$failure))))
effects <- rbindlist(list(effects[cohort!="GSE23597"],
  fit_all_signature_effects(act_pheno, "GSE23597", "infliximab", genes)[signature_id %in% primary_ids],
  fit_all_signature_effects(act_pheno, "GSE23597", "placebo", genes)[signature_id %in% primary_ids]))
active <- effects[treatment != "placebo"]
stopifnot(anyDuplicated(active, by = c("signature_id", "cohort")) == 0L,
          all(is.finite(active$log_or)), all(active$se > 0))
all_arms <- unique(active$cohort)
act <- c("GSE23597", "GSE12251")
leuven <- c("GSE16879", "GSE73661_IFX")
other <- setdiff(all_arms, c(act, leuven))
scenarios <- list(primary = c(other, act[1], leuven[1]),
                  switch_ACT = c(other, act[2], leuven[1]),
                  switch_Leuven = c(other, act[1], leuven[2]),
                  switch_both = c(other, act[2], leuven[2]))
counts <- load_result("benchmark_sample_counts")[treatment != "placebo"]
counts[cohort=="GSE23597", `:=`(n=nrow(act_active),failures=sum(act_active$failure),
                                 successes=sum(act_active$failure==0))]
stopifnot(counts[cohort == act[1], n] > counts[cohort == act[2], n],
          counts[cohort == leuven[1], n] > counts[cohort == leuven[2], n])
plan <- rbindlist(lapply(names(scenarios), function(s)
  counts[, .(scenario = s, cohort, treatment, n, retained = cohort %in% scenarios[[s]])]))
save_result(plan, "arm_selection")

# Conservative leave-publication-source-out: when a paper used either arm of
# GSE73661, exclude that entire deposit and the linked Leuven IFX records.
# GSE14580 shares UC baseline GSM accessions with GSE16879 (e.g. GSM364633).
# West et al. (10.1038/nm.4307) also used GSE16879 and NCT00487539 (PURSUIT).
# Mennillo et al. (10.1038/s41467-024-45665-6, Fig. 8) used GSE73661.
registry[, corrected_source := source_cohorts]
registry[signature_id == "OSM_OSMR",
         corrected_source := "GSE12251;GSE23597;GSE16879;GSE92415"]
registry[signature_id %in% c("VDZ_R_5","VDZ_NR_3"),
         corrected_source := paste0(source_cohorts, ";GSE73661")]
expand_sources <- function(s) {
  ids <- strsplit(s, ";", fixed = TRUE)[[1]]
  if (any(ids %in% act)) ids <- union(ids, act)
  if ("GSE73661" %in% ids) ids <- union(ids, c(leuven, "GSE73661_VDZ"))
  if (any(ids %in% c("GSE14580", leuven))) ids <- union(ids, leuven)
  intersect(all_arms, ids)
}
registry[, excluded_arms := vapply(corrected_source, function(x)
  paste(expand_sources(x), collapse = ";"), character(1))]
save_result(registry, "source_exclusion_map")

# The revised benchmark reports the seven primary signatures as its BH
# family. Non-estimable external tests count in that family's size.
pool <- function(d) {
  # The unpenalised GSE12251 OSM slope has no finite MLE (complete
  # separation). Its archive contains a numerical iterate of 1863.8.
  # Exclude that unestimable coefficient when computing pooled ORs.
  if (all(c("signature_id","cohort") %in% names(d)))
    d <- d[!(signature_id=="OSM_OSMR" & cohort=="GSE12251")]
  k <- nrow(d)
  fit <- if (k >= 2) random_effects_meta(d$log_or, d$se) else
    data.table(log_or=NA_real_, se=NA_real_, or=NA_real_, ci_low=NA_real_,
               ci_high=NA_real_, p_value=NA_real_, tau2=NA_real_, i2=NA_real_)
  if (k >= 3) {
    delta <- qt(0.975, k - 2) * sqrt(fit$tau2 + fit$se^2)
    fit[, `:=`(prediction_low = exp(log_or-delta), prediction_high = exp(log_or+delta))]
  } else fit[, `:=`(prediction_low = NA_real_, prediction_high = NA_real_)]
  cbind(data.table(cohorts = k, total_n = sum(d$n), positive_cohorts = sum(d$log_or > 0),
                   included_cohorts = paste(d$cohort, collapse = ";")), fit)
}
fisher <- load_result("benchmark_real_fisher_z")
fisher <- rbindlist(list(fisher[cohort!="GSE23597"], act_z),use.names=TRUE)
null <- setNames(lapply(all_arms, function(x)
  readRDS(file.path(DIR_INTERMEDIATE, paste0("benchmark_", x, "_random_null.rds")))), all_arms)
null[["GSE23597"]] <- act_null
empirical <- function(id, arms) {
  z <- fisher[signature_id == id & cohort %in% arms]
  stopifnot(nrow(z) == length(arms), anyDuplicated(z$cohort) == 0)
  if (nrow(z) < 2) return(data.table(random_iterations=NA_integer_, real_stouffer_z=NA_real_,
                                      empirical_p_one_sided=NA_real_, empirical_p_two_sided=NA_real_))
  m <- do.call(cbind, lapply(z$cohort, function(x) null[[x]][,id]))
  stopifnot(nrow(m)==10000L, all(is.finite(m)))
  observed <- sum(z$fisher_z) / sqrt(nrow(z))
  simulated <- rowSums(m) / sqrt(ncol(m))
  data.table(random_iterations = length(simulated), real_stouffer_z = observed,
             empirical_p_one_sided = (1 + sum(simulated >= observed)) / (1 + length(simulated)),
             empirical_p_two_sided = (1 + sum(abs(simulated) >= abs(observed))) / (1 + length(simulated)))
}
ids <- unique(active$signature_id)
combined <- list()
null_out <- list()
for (s in names(scenarios)) {
  for (id in ids) {
    blocked <- registry[signature_id == id, excluded_arms]
    exclude <- if (!length(blocked)) character() else strsplit(blocked, ";", fixed=TRUE)[[1]]
    subsets <- list(all = scenarios[[s]], source_excluded = setdiff(scenarios[[s]], exclude))
    for (analysis in names(subsets)) {
      arms <- subsets[[analysis]]
      key <- data.table(scenario=s, signature_id=id, analysis=analysis)
      d <- active[signature_id == id & cohort %in% arms]
      combined[[length(combined)+1L]] <- cbind(key, pool(d))
      null_out[[length(null_out)+1L]] <- cbind(key, data.table(cohorts=nrow(d)), empirical(id, arms))
    }
  }
}
combined <- rbindlist(combined, fill=TRUE)
combined[, fdr := p.adjust(p_value, "BH", n = length(ids)), by=.(scenario, analysis)]
save_result(combined[analysis=="all"], "active_meta")
save_result(combined[analysis=="source_excluded"], "source_excluded_meta")
null_out <- rbindlist(null_out)
null_out[, empirical_fdr := p.adjust(empirical_p_one_sided, "BH", n=length(ids)),
         by=.(scenario, analysis)]
save_result(null_out, "matched_null")

endpoint <- load_result("cbm_upgrade/endpoint_harmonization_registry")
retained <- merge(active[signature_id %in% primary_ids & cohort %in% scenarios$primary],
                  endpoint, by="cohort", sort=FALSE)
pool_by_signature <- function(d) d[, pool(.SD), by=signature_id]
strata <- function(subsets) rbindlist(lapply(names(subsets), function(s)
  pool_by_signature(subsets[[s]])[, analysis := s]), fill=TRUE)
stability <- strata(list(all_active=retained,
                         microarray_only=retained[cohort!="GSE72819"],
                         randomized_active_arms=retained[cohort %in% c("GSE206285","GSE23597","GSE92415")],
                         excluding_largest_cohort=retained[cohort!="GSE206285"],
                         omit_ACT_family=retained[!cohort %in% act],
                         omit_Leuven_IFX_family=retained[!cohort %in% leuven],
                         omit_both_IFX_families=retained[!cohort %in% c(act, leuven)],
                         omit_all_Leuven_arms=retained[!cohort %in% c(leuven, "GSE73661_VDZ")]))
stability[, fdr:=p.adjust(p_value, "BH"), by=analysis]
save_result(stability, "stability")
loo <- rbindlist(lapply(scenarios$primary, function(omit)
  pool_by_signature(retained[cohort != omit])[,omitted_cohort := omit]))
save_result(loo, "leave_one_out")
ep <- strata(list(all_active=retained,
                  objective_or_composite_mucosal=retained[endpoint_family=="objective_or_composite_mucosal"],
                  clinical_only=retained[endpoint_family=="clinical_only"],
                  week_6_8=retained[assessment_window=="week_6_8"],
                  other_post_induction=retained[assessment_window=="other_post_induction"]))
setnames(ep, "analysis","stratum")
ep[, fdr_within_stratum:=p.adjust(p_value, "BH"), by=stratum]
contrast <- merge(ep[stratum=="objective_or_composite_mucosal"], ep[stratum=="clinical_only"],
                  by="signature_id", suffixes=c("_mucosal","_clinical"))
contrast[, difference := log_or_mucosal - log_or_clinical]
contrast[, se_difference := sqrt(se_mucosal^2+se_clinical^2)]
contrast[, `:=`(ratio_of_ors=exp(difference), ci_low=exp(difference-1.96*se_difference),
                ci_high=exp(difference+1.96*se_difference),
                p_value=2*pnorm(-abs(difference/se_difference)))]
contrast[, fdr:=p.adjust(p_value,"BH")]
save_result(ep, "endpoint_stratified")
save_result(contrast[, .(signature_id, ratio_of_ors, ci_low, ci_high, p_value, fdr)], "endpoint_contrasts")

save_result(effects, "cohort_effects")
placebo <- pool_by_signature(effects[treatment=="placebo"])
placebo[, fdr:=p.adjust(p_value,"BH")]
save_result(placebo,"placebo_meta")
trial_treatments <- c(GSE206285="ustekinumab", GSE92415="golimumab", GSE23597="infliximab")
interaction <- rbindlist(lapply(names(trial_treatments), function(arm) {
  d0 <- if(arm=="GSE23597") copy(act_pheno) else
    as.data.table(readRDS(file.path(DIR_INTERMEDIATE, paste0("benchmark_",arm,"_scored.rds"))))
  adjust <- switch(arm, GSE206285=c("mayo","log_crp","log_calprotectin","anti_tnf_history"),
                        GSE92415=c("mayo","age"), GSE23597=character())
  if (arm=="GSE206285") d0[, `:=`(log_crp=log1p(crp),log_calprotectin=log1p(calprotectin),
                                   anti_tnf_history=as.factor(anti_tnf_history))]
  rbindlist(lapply(primary_ids, function(id) {
    vars <- c("failure","treatment",id,adjust)
    d <- d0[complete.cases(d0[,..vars])]
    d[, `:=`(score_z=standardize_score(get(id)),treatment=relevel(factor(treatment),"placebo"))]
    model <- glm(reformulate(c("treatment * score_z",adjust),"failure"), family=binomial(),data=d)
    fit <- fit_interaction_term(model,paste0("treatment",trial_treatments[arm],":score_z"),id,arm)
    fit[, `:=`(eligible_n=nrow(d0), active_n=sum(d$treatment!="placebo"),placebo_n=sum(d$treatment=="placebo"))]
    fit
  }))
}))
interaction[,fdr_within_trial:=p.adjust(p_value,"BH"),by=cohort]
save_result(interaction,"interaction_trial")
precision <- interaction[,cbind(pool(.SD),data.table(total_randomized_n=sum(eligible_n),
                      active_complete_n=sum(active_n),placebo_complete_n=sum(placebo_n))),by=signature_id]
precision[, `:=`(trials=cohorts,complete_case_n=total_n, fdr=p.adjust(p_value,"BH"),
                ci_within_illustrative_margin=ci_low>=2/3 & ci_high<=1.5,
                normal_approx_mde_or_80=exp((qnorm(.975)+qnorm(.80))*se),
                tost_p=pmax(pnorm((log_or+log(1.5))/se,lower.tail=FALSE),
                            pnorm((log(1.5)-log_or)/se,lower.tail=FALSE)))]
precision[, equivalence_supported_05:=tost_p<.05]
save_result(precision,"interaction_precision")

# Sanity-check archived cohort-level models by independently refitting the
# within-arm Bernoulli logistic slope from patient-level scores.
refit <- rbindlist(lapply(all_arms, function(arm) {
  scored <- as.data.table(readRDS(file.path(DIR_INTERMEDIATE, paste0("benchmark_",arm,"_scored.rds"))))
  d <- scored[treatment != "placebo"]
  if (arm=="GSE23597") d <- d[geo_accession != "GSM578741"]
  rbindlist(lapply(ids, function(id) {
    m <- glm(failure ~ scale(score), family=binomial(),
             data=data.frame(failure=d$failure,score=d[[id]]))
    data.table(cohort=arm, signature_id=id, n_refit=nobs(m),
               slope_refit=unname(coef(m)[2]), se_refit=sqrt(vcov(m)[2,2]))
  }))
}))
audit <- merge(refit, active[,.(cohort,signature_id,n,log_or,se)], by=c("cohort","signature_id"))
stopifnot(nrow(audit)==nrow(active), all(audit$n_refit==audit$n),
          max(abs(audit$slope_refit-audit$log_or))<1e-6,
          max(abs(audit$se_refit-audit$se))<1e-6)
save_result(audit, "cohort_refit_audit")
print(combined[scenario=="primary" & analysis=="all" & signature_id %in% primary_ids,
               .(signature_id, cohorts, total_n, or, ci_low, ci_high, i2, fdr)])
