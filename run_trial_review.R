# Run from the extracted R project root. With "clinical", reuse completed
# overlap-review results. The default starts at the public-source data.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L || (length(args) && !identical(args, "clinical"))) {
  stop("Usage: Rscript --vanilla run_trial_review.R [clinical]")
}
if (!file.exists("UC_MucosalRecovery_Q2_10IF_project.Rproj")) {
  stop("Run from the directory containing the .Rproj file.")
}
steps <- c("scripts/28_trial_mechanism_interactions.R", "scripts/29_unifi_residual_auc.R")
if (!length(args)) steps <- c("run_overlap_review.R", steps)
rscript <- file.path(R.home("bin"), "Rscript")
for (step in steps) {
  message("Running ", step)
  status <- system2(rscript, c("--vanilla", step))
  if (status != 0L) stop("Stopped at ", step, " (status ", status, ").")
}
target <- file.path("results", "version_audit")
dir.create(target, recursive = TRUE, showWarnings = FALSE)
capture.output(sessionInfo(), file = file.path(target, "sessionInfo_clinical_driver.txt"))
