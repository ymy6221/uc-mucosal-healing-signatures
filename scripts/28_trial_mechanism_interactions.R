# Run from the project root after 26_overlap_reanalysis.R.
# Clinical outcomes and treatment remain as coded in each source trial.
library(data.table)
library(ggplot2)
source("config.R")
source("scripts/helpers.R")
setDTthreads(1L)

previous <- file.path(DIR_RESULTS, "overlap_reanalysis")
out <- file.path(DIR_RESULTS, "review_clinical_trials")
dir.create(out, showWarnings = FALSE)
write_table <- function(x, name) fwrite(x, file.path(out, paste0(name, ".csv")))
prior <- fread(file.path(previous, "interaction_trial.csv"))
ids <- unique(prior$signature_id)
stopifnot(length(ids) == 7L)
names_by_trial <- c("UNIFI" = "GSE206285", "PURSUIT-SC" = "GSE92415", "ACT-1" = "GSE23597")
covariates <- list("UNIFI" = c("mayo", "log_crp", "log_calprotectin", "anti_tnf_history"),
                   "PURSUIT-SC" = c("mayo", "age"), "ACT-1" = character())
cohort <- list()
flow <- list()
missing_patterns <- list()
missing_variables <- list()
for (tr in names(names_by_trial)) {
  acc <- names_by_trial[[tr]]
  a <- as.data.table(readRDS(file.path(DIR_INTERMEDIATE, paste0("benchmark_", acc, "_scored.rds"))))
  if (tr == "ACT-1") {
    setorder(a, geo_accession)
    stopifnot(identical(a[duplicated(patient_id), geo_accession], "GSM578741"))
    a <- a[!duplicated(patient_id)]
  }
  if (tr == "UNIFI") {
    a[, log_crp := log1p(crp)]
    a[, log_calprotectin := log1p(calprotectin)]
    a[, anti_tnf_history := factor(anti_tnf_history)]
  }
  # Each signature therefore shares the same eligible and complete-case sample.
  stopifnot(all(complete.cases(a[, c(ids, "failure", "treatment"), with = FALSE])))
  a[, active := as.integer(treatment != "placebo")]
  a[, trial := tr]
  adjustment <- covariates[[tr]]
  if (length(adjustment)) {
    miss <- is.na(a[, ..adjustment])
    a[, keep_case := rowSums(miss) == 0]
    a[, missing_pattern := apply(miss, 1L, function(z) {
      if (any(z)) paste(adjustment[z], collapse = " + ") else "complete"
    })]
    missing_variables[[tr]] <- data.table(trial=tr, variable=adjustment,
                                          missing_n=colSums(miss), eligible_n=nrow(a))
  } else {
    a[, `:=`(keep_case = TRUE, missing_pattern = "complete")]
  }
  flow[[tr]] <- a[, .(eligible_n=.N, complete_n=sum(keep_case), excluded_n=sum(!keep_case),
                      failures_eligible=sum(failure), failures_complete=sum(failure[keep_case])),
                  by=.(trial, treatment)]
  missing_patterns[[tr]] <- a[, .(n=.N), by=.(trial, treatment, missing_pattern)]
  cohort[[tr]] <- a
}
write_table(rbindlist(flow), "participant_flow")
write_table(rbindlist(missing_patterns), "missingness_patterns")
write_table(rbindlist(missing_variables), "missingness_variables")
stopifnot(sum(vapply(cohort, nrow, integer(1))) == 673L,
          sum(vapply(cohort, function(x) sum(x$keep_case), integer(1))) == 630L)

interval <- function(b, se) {
  data.table(log_or=b, se=se, or=exp(b), ci_low=exp(b-1.96*se),
             ci_high=exp(b+1.96*se), p_value=2*pnorm(-abs(b/se)))
}
fit_one <- function(d, id, adjustment) {
  x <- copy(d)
  x[, z := as.numeric(scale(get(id)))]
  f <- reformulate(c("active * z", adjustment), response="failure")
  # Compare two IRLS stopping tolerances to guard against a drifting MLE.
  # Stop before reporting any unstable Wald estimate.
  warnings <- character()
  fits <- withCallingHandlers(lapply(c(1e-8, 1e-12), function(e)
    glm(f, x, family=binomial(), control=glm.control(epsilon=e, maxit=100))),
    warning=function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  valid <- all(vapply(fits, function(g) g$converged && !g$boundary &&
    g$rank == length(coef(g)) && all(is.finite(coef(g))), logical(1)))
  if (!valid || length(warnings) || max(abs(coef(fits[[1]]) - coef(fits[[2]]))) > 1e-4) {
    stop("Unstable interaction fit for ", id, ": ", paste(warnings, collapse="; "))
  }
  # Use the same stopping rule as the locked primary models after validation.
  g <- fits[[1]]
  term <- "active:z"
  se <- sqrt(vcov(g)[term, term])
  stopifnot(is.finite(se), se>0, nobs(g)==nrow(x))
  cbind(data.table(n=nrow(x), active_n=sum(x$active), placebo_n=sum(1-x$active)),
        interval(unname(coef(g)[term]), se))
}

models <- c("adjusted_complete", "unadjusted_complete", "unadjusted_all")
trial <- rbindlist(lapply(names(cohort), function(tr) {
  a <- cohort[[tr]]
  rbindlist(lapply(models, function(spec) {
    d <- if (spec == "unadjusted_all") a else a[keep_case == TRUE]
    vars <- if (spec == "adjusted_complete") covariates[[tr]] else character()
    rbindlist(lapply(ids, function(id) cbind(data.table(signature_id=id, trial=tr,
      cohort=names_by_trial[[tr]], model=spec), fit_one(d,id,vars))))
  }))
}))
trial[, fdr := p.adjust(p_value,"BH"), by=.(model,trial)]
# Agreement with the historical adjusted fits after tighter convergence.
audit <- merge(trial[model == "adjusted_complete"], prior, by=c("signature_id","cohort"))
stopifnot(nrow(audit)==21L, max(abs(audit$log_or.x-audit$log_or.y))<1e-5,
          max(abs(audit$se.x-audit$se.y))<1e-4, all(audit$n.x==audit$n.y))
write_table(trial, "trial_models")

sets <- list(all_trials=names(names_by_trial), anti_TNF=c("PURSUIT-SC","ACT-1"))
pools <- rbindlist(lapply(names(sets), function(pool_id) {
  trial[trial %in% sets[[pool_id]], {
    fit <- random_effects_meta(log_or, se)
    cbind(data.table(pool=pool_id,k=.N,n=sum(n)), fit)
  }, by=.(signature_id,model)]
}))
pools[, fdr := p.adjust(p_value,"BH"), by=.(model,pool)]
write_table(pools, "pooled_models")
weights <- merge(trial, pools[, .(signature_id,model,pool,tau2)],
                 by=c("signature_id","model"),allow.cartesian=TRUE)
weights <- weights[pool == "all_trials" | trial != "UNIFI"]
weights[, w := 1/(se^2+tau2)]
weights[, `:=`(weight_pct=100*w/sum(w), participant_pct=100*n/sum(n)),
        by=.(signature_id,model,pool)]
write_table(weights[, .(model,pool,signature_id,trial,n,participant_pct,tau2,weight_pct)],
            "random_effects_weights")
contrast <- merge(pools[pool == "anti_TNF"], trial[trial == "UNIFI"],
                  by=c("signature_id","model"), suffixes=c("_anti_TNF","_UNIFI"))
contrast <- contrast[, interval(log_or_anti_TNF-log_or_UNIFI, sqrt(se_anti_TNF^2+se_UNIFI^2)),
                     by=.(signature_id,model)]
contrast[, fdr := p.adjust(p_value, "BH"), by=model]
write_table(contrast, "treatment_class_contrasts")

placebo <- fread(file.path(previous,"placebo_meta.csv"))[, .(signature_id, or, ci_low, ci_high)]
placebo[, panel := "F   Placebo association (225)"]
rows <- trial[model=="adjusted_complete", .(signature_id,panel=trial,or,ci_low,ci_high)]
panels <- c("A   UNIFI (499)", "B   PURSUIT-SC (87)", "C   ACT-1 (44)",
            "D   Anti-TNF trials (131)", "E   All trials (630)", "F   Placebo association (225)")
rows[, panel := panels[match(panel, names(names_by_trial))]]
poolrows <- pools[model=="adjusted_complete",
                  .(signature_id, panel=ifelse(pool=="anti_TNF",panels[4],panels[5]),or,ci_low,ci_high)]
rows <- rbindlist(list(rows, poolrows, placebo),use.names=TRUE)
rows[, panel := factor(panel, levels=panels)]
names_sig <- c(IFX_ARJS_5="Arijs anti-TNF", OSM_OSMR="OSM–OSMR", IL22_NEUTROPHIL_6="IL-22 chemokines",
               VDZ_NR_3="Vedolizumab nonresponse stromal", VDZ_4="Vedolizumab 4-gene",
               ETRO_2="Etrolizumab GZMA–ITGAE", VDZ_R_5="Vedolizumab crypt-base")
rows[, signature := factor(signature_id, levels=rev(names(names_sig)),labels=rev(unname(names_sig)))]
fig <- ggplot(rows,aes(or,signature)) +
  geom_vline(xintercept=1,linetype="dashed",colour="#6B7880",linewidth=.35) +
  geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation="y",width=.16,linewidth=.55,colour="#17536A") +
  geom_point(size=1.9,colour="#17536A") +
  facet_wrap(~panel, ncol=3) +
  scale_x_log10(breaks=c(.125,.5,2,8,32,128),labels=c("⅛","½","2","8","32","128")) +
  labs(x="Odds ratio (log scale)",y=NULL) +
  theme_bw(base_size=10.5) +
  theme(strip.background=element_rect(fill="#EDF2F3",colour=NA),
        strip.text=element_text(hjust=0,face="bold",colour="#173D4B"),
        panel.border=element_rect(colour="#B7C3C8",fill=NA,linewidth=.35),
        panel.grid.major.y=element_line(colour="#E7EDEE",linewidth=.25),
        panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
        axis.text=element_text(colour="#263C46"))
ggsave(file.path(out,"Figure_2.tiff"),fig,device=ragg::agg_tiff,compression="lzw",
       width=10.0,height=5.6,units="in",dpi=350,bg="white")
ggsave(file.path(out,"Figure_2.png"),fig,device=ragg::agg_png,
       width=10.0,height=5.6,units="in",dpi=150,bg="white")
