# Run from the directory containing this script and the .Rproj file.
# Each stage runs in a fresh R session so large matrices can be released.
if (!file.exists("UC_MucosalRecovery_Q2_10IF_project.Rproj")) {
  stop("Set your working directory to the extracted R project.")
}
rscript <- file.path(R.home("bin"), "Rscript")
stages <- c(
  "run_stage1_core.R", "run_stage2_longitudinal.R", "run_stage3_singlecell.R",
  "run_stage4_figures.R", "run_stage5_benchmark_upgrade.R", "run_stage6_cbm_upgrade.R",
  "scripts/audit_geo_overlap.R", "scripts/26_overlap_reanalysis.R",
  "scripts/review_pathways.R", "scripts/27_overlap_figures.R",
  "scripts/audit_recorded_versions.R"
)
for (stage in stages) {
  message("Running ", stage)
  status <- system2(rscript, c("--vanilla", stage))
  if (status != 0L) stop("Stopped at ", stage, " (exit status ", status, ").")
}
