# 只安装 CBM 升级新增且本机缺少的包；已安装的包不会被覆盖。

cran_packages <- c("msigdbr", "Seurat", "SeuratObject", "ggrepel")
missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace,
                                     logical(1), quietly = TRUE)]
if (length(missing_cran) > 0L) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
bioc_packages <- c("GSVA", "BiocParallel")
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace,
                                     logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0L) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

all_packages <- c(cran_packages, bioc_packages)
print(vapply(all_packages, requireNamespace, logical(1), quietly = TRUE))
