# 从公共下载开始完整复跑。已有且哈希正确的数据会跳过下载，任何原始文件都不删除。

start_time <- Sys.time()
source("run_stage1_core.R")
source("run_stage2_longitudinal.R")
source("run_stage3_singlecell.R")
source("run_stage4_figures.R")
source("run_stage5_benchmark_upgrade.R")
source("run_stage6_cbm_upgrade.R")
elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
message("全部阶段从公共数据完整复跑成功，总用时约 ", round(elapsed, 1), " 分钟。")
