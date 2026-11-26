# 第一阶段总入口：下载、评分、队列分析、随机效应 Meta。

# 载入配置和公共函数。
source("config.R")
source(file.path(PROJECT_ROOT, "scripts", "helpers.R"))

# 记录开始时间，便于评估电脑负担。
start_time <- Sys.time()

# 按顺序执行。中途报错时，修复后重新运行即可；已下载文件会自动跳过。
source(file.path(PROJECT_ROOT, "scripts", "01_download_core.R"))
source(file.path(PROJECT_ROOT, "scripts", "02_score_cohorts.R"))
source(file.path(PROJECT_ROOT, "scripts", "03_analyze_primary.R"))
source(file.path(PROJECT_ROOT, "scripts", "04_meta_analysis.R"))

# 写入运行环境，投稿补充材料可直接引用。
writeLines(capture.output(sessionInfo()), file.path(DIR_RESULTS, "sessionInfo.txt"))

# 打印耗时和下一步提示。
elapsed_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
message("第一阶段完成，用时约 ", round(elapsed_minutes, 1), " 分钟。")
message("先打开 results/quality_gate.csv，再运行 run_stage2_longitudinal.R")

