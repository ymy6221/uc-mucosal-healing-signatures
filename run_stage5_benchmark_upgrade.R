# 高影响力升级总入口。请先完成 run_stage1_core.R；之后只需运行本文件。

source("config.R")
source(file.path(PROJECT_ROOT, "scripts", "helpers.R"))
source(file.path(PROJECT_ROOT, "scripts", "benchmark_helpers.R"))

start_time <- Sys.time()
source(file.path(PROJECT_ROOT, "scripts", "10_download_benchmark.R"))
source(file.path(PROJECT_ROOT, "scripts", "15_verify_public_downloads.R"))
source(file.path(PROJECT_ROOT, "scripts", "11_score_benchmark.R"))

# 后续统计、作图和高影响力门槛检查由下面三个脚本完成。
source(file.path(PROJECT_ROOT, "scripts", "12_analyze_benchmark.R"))
source(file.path(PROJECT_ROOT, "scripts", "13_make_benchmark_figures.R"))
source(file.path(PROJECT_ROOT, "scripts", "14_validate_benchmark.R"))

elapsed_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
message("高影响力升级阶段完成，用时约 ", round(elapsed_minutes, 1), " 分钟。")
message("请先打开 results/high_impact_claim_gate.csv，再决定是否进入投稿写作。")
