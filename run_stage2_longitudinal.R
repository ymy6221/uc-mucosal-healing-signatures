# 第二阶段入口：患者配对的纵向恢复分析和第一批投稿图。

source("config.R")
source(file.path(PROJECT_ROOT, "scripts", "helpers.R"))

# 若第一阶段未完成，直接给出明确错误。
if (!file.exists(file.path(DIR_RESULTS, "quality_gate.csv"))) {
  stop("请先成功运行 run_stage1_core.R")
}

# 运行纵向分析和作图。
source(file.path(PROJECT_ROOT, "scripts", "05_longitudinal.R"))
source(file.path(PROJECT_ROOT, "scripts", "07_make_figures.R"))

message("第二阶段完成。请查看 results/longitudinal_summary.csv 和 figures/。")

