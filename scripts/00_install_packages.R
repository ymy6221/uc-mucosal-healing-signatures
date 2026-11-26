# 只在第一次运行时执行本脚本；以后不用重复安装。

# CRAN 包负责读表、作图和 ROC。
cran_packages <- c("data.table", "ggplot2", "pROC", "patchwork", "Matrix",
                   "SeuratObject", "scales", "digest")

# 找出尚未安装的 CRAN 包。
missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]

# 仅安装缺失项，避免无意义地覆盖现有环境。
if (length(missing_cran) > 0L) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

# GEOquery 和 Biobase 来自 Bioconductor。
bioc_packages <- c("GEOquery", "Biobase", "AnnotationDbi", "org.Hs.eg.db")

# 若 BiocManager 不存在，先从 CRAN 安装。
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# 找出尚未安装的 Bioconductor 包。
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]

# 仅安装缺失项；ask=FALSE 防止中途等待键盘输入。
if (length(missing_bioc) > 0L) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

# 打印确认信息。全部为 TRUE 才能进入下一步。
all_packages <- c(cran_packages, bioc_packages)
print(vapply(all_packages, requireNamespace, logical(1), quietly = TRUE))
