# Full source-data analysis by default. Use "validation" to reuse a completed
# run_trial_review.R analysis in the same extracted R project.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) && !identical(args, "validation")) {
  stop("Usage: Rscript --vanilla run_validation_review.R [validation]")
}
if (!file.exists("UC_MucosalRecovery_Q2_10IF_project.Rproj")) {
  stop("Open the folder containing the extracted .Rproj file.")
}
steps <- c("scripts/30_unifi_repeated_cv.R", "scripts/31_telesco_proxy_results.R",
           "scripts/32_antitnf_pool_sensitivity.R")
if (!length(args)) steps <- c("run_trial_review.R", steps)
for (step in steps) {
  message("Running ", step)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", step))
  if (status != 0L) stop("Stopped at ", step, " (exit status ", status, ").")
}
message("Validation analyses completed. Compare results with the reference tables ",
        "using Rscript --vanilla scripts/check_review_tables.R")
