# 第三阶段入口：低内存单细胞伪汇总和细胞定位。

source("config.R")
source(file.path(PROJECT_ROOT, "scripts", "helpers.R"))

# 核心临床证据必须先完成，避免先烧资源再发现主线不成立。
if (!file.exists(file.path(DIR_RESULTS, "quality_gate.csv"))) {
  stop("请先成功运行 run_stage1_core.R")
}

# 运行聚合和作图。
source(file.path(PROJECT_ROOT, "scripts", "06_singlecell_localization.R"))
source(file.path(PROJECT_ROOT, "scripts", "08_make_singlecell_figure.R"))

message("第三阶段完成。请查看 results/singlecell_epithelial_test.csv 和 Figure4。")

