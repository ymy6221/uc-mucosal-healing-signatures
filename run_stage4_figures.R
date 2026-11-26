# 只重新作图，不重新下载或评分。

source("config.R")
source(file.path(PROJECT_ROOT, "scripts", "helpers.R"))

if (!file.exists(file.path(DIR_RESULTS, "longitudinal_pairs.csv"))) {
  stop("请先运行 run_stage2_longitudinal.R")
}

source(file.path(PROJECT_ROOT, "scripts", "07_make_figures.R"))
message("作图完成：", DIR_FIGURES)

