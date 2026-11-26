# 在既有 benchmark 完成后运行：机制通路、异质性、精细单细胞、图和验收。

source("config.R")
source(file.path(PROJECT_ROOT, "scripts", "helpers.R"))
source(file.path(PROJECT_ROOT, "scripts", "benchmark_helpers.R"))

start_time <- Sys.time()
source(file.path(PROJECT_ROOT, "scripts", "16_pathway_mechanism.R"))
source(file.path(PROJECT_ROOT, "scripts", "17_heterogeneity_stability.R"))
source(file.path(PROJECT_ROOT, "scripts", "18_singlecell_fine_localization.R"))
source(file.path(PROJECT_ROOT, "scripts", "22_review_critique_sensitivity.R"))
source(file.path(PROJECT_ROOT, "scripts", "19_make_cbm_upgrade_figures.R"))
source(file.path(PROJECT_ROOT, "scripts", "23_make_review_sensitivity_figures.R"))
source(file.path(PROJECT_ROOT, "scripts", "20_validate_cbm_upgrade.R"))
source(file.path(PROJECT_ROOT, "scripts", "24_validate_review_sensitivity.R"))
elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
message("CBM 升级阶段完成，用时约 ", round(elapsed, 1), " 分钟。")
