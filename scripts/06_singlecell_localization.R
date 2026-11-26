# 用 GSE250487 做细胞来源定位。策略是“直接伪汇总”，不重聚类、不跑 UMAP。

check_packages()

# 单细胞阶段只额外需要稀疏矩阵和 Seurat 对象读取工具。
single_cell_packages <- c("Matrix", "SeuratObject")
single_cell_ok <- vapply(single_cell_packages, requireNamespace, logical(1), quietly = TRUE)
if (!all(single_cell_ok)) {
  stop("缺少单细胞轻量依赖：", paste(single_cell_packages[!single_cell_ok], collapse = ", "),
       "。请重新运行 scripts/00_install_packages.R")
}

# 下载 GEO 提供的处理后 RDS。该文件约 389 MB。
single_cell_path <- file.path(PROJECT_ROOT, "data", "singlecell", "GSE250487_biopsy_RNA.rds.gz")
download_with_retry(SC_RDS_URL, single_cell_path, min_bytes = 300000000L)

# GEO 文件外层和 RDS 本身各有一层 gzip；用两个连接直接读取，不生成第二份 410 MB 文件。
read_double_gzip_rds <- function(path) {
  outer_connection <- gzfile(path, open = "rb")
  inner_connection <- gzcon(outer_connection)
  on.exit(close(inner_connection), add = TRUE)
  readRDS(inner_connection)
}

# 读取对象；对象约 1.9 GB，后面会尽快只保留聚合结果。
message("读取 GSE250487；普通电脑约需数秒到数十秒。")
seurat_object <- read_double_gzip_rds(single_cell_path)

# 检查作者提供的元数据字段是否仍与教程一致。
required_metadata <- c(
  "CoLabs_sample", "patient_short", "coarse_annotations_MK", "condition",
  "Endoscopic_Severity_V1", "nCount_RNA"
)
missing_metadata <- setdiff(required_metadata, colnames(seurat_object@meta.data))
if (length(missing_metadata) > 0L) {
  stop("GSE250487 缺少预期元数据列：", paste(missing_metadata, collapse = ", "))
}

# 建立“样本×粗细胞类型”分组；每个细胞只属于一个组。
cell_metadata <- data.table::as.data.table(seurat_object@meta.data, keep.rownames = "cell")
cell_metadata[, group := paste(CoLabs_sample, coarse_annotations_MK, sep = "___")]
group_factor <- factor(cell_metadata$group)
design_matrix <- Matrix::sparse.model.matrix(~ 0 + group_factor)
colnames(design_matrix) <- levels(group_factor)

# 从 RNA assay 取稀疏原始计数，并一次矩阵乘法完成伪汇总。
count_matrix <- SeuratObject::LayerData(seurat_object[["RNA"]], layer = "counts")
pseudobulk_counts <- count_matrix %*% design_matrix

# 每组总文库量来自作者给出的每细胞 nCount_RNA；用于 CPM 标准化。
pseudobulk_library <- as.numeric(Matrix::crossprod(design_matrix, as.numeric(cell_metadata$nCount_RNA)))
pseudobulk_logcpm <- log2(sweep(as.matrix(pseudobulk_counts), 2L, pseudobulk_library, "/") * 1e6 + 1)

# 用与 bulk 完全相同的样本内秩算法计算 ERI，保持尺度定义一致。
signature_table <- data.table::fread(file.path(PROJECT_ROOT, "signatures", "epithelial_restitution.csv"))
pseudobulk_scores <- calculate_eri(pseudobulk_logcpm, signature_table)
data.table::setnames(pseudobulk_scores, "geo_accession", "group")

# 每个组的样本、患者、细胞类型和严重度信息。
group_information <- unique(cell_metadata[, .(
  group, CoLabs_sample, patient_short, coarse_annotations_MK, condition, Endoscopic_Severity_V1
)])
group_information <- merge(group_information, pseudobulk_scores, by = "group", all.x = TRUE, sort = FALSE)
group_information[, n_cells := as.numeric(Matrix::colSums(design_matrix))[match(group, colnames(design_matrix))]]

# 少于 20 个细胞的伪汇总不稳定，预先规定为排除。
group_information <- group_information[n_cells >= 20L]
group_information[, severity_order := data.table::fcase(
  Endoscopic_Severity_V1 == "HC", 0,
  Endoscopic_Severity_V1 == "Severity_0_1", 1,
  Endoscopic_Severity_V1 == "Severity_2_3", 2,
  default = NA_real_
)]

# 对每个细胞类型计算重度 UC 相对健康的中位 ERI 位移。
cell_type_shift <- group_information[
  Endoscopic_Severity_V1 %in% c("HC", "Severity_2_3"),
  .(
    n_hc = sum(Endoscopic_Severity_V1 == "HC"),
    n_severe = sum(Endoscopic_Severity_V1 == "Severity_2_3"),
    median_hc = stats::median(eri[Endoscopic_Severity_V1 == "HC"]),
    median_severe = stats::median(eri[Endoscopic_Severity_V1 == "Severity_2_3"]),
    shift_severe_minus_hc = stats::median(eri[Endoscopic_Severity_V1 == "Severity_2_3"]) -
      stats::median(eri[Endoscopic_Severity_V1 == "HC"])
  ),
  by = coarse_annotations_MK
]

# 某些稀有细胞类型在健康或重度组没有达到 20 细胞，不能计算位移，明确排除。
cell_type_shift <- cell_type_shift[
  n_hc > 0L & n_severe > 0L & is.finite(shift_severe_minus_hc)
]

# 上皮伪汇总是预先指定的主定位检验。
epithelial <- group_information[coarse_annotations_MK == "10-Epithelial" & !is.na(severity_order)]
epithelial_trend <- stats::cor.test(epithelial$eri, epithelial$severity_order,
                                    method = "spearman", exact = FALSE)
epithelial_extreme <- stats::wilcox.test(
  eri ~ Endoscopic_Severity_V1,
  data = epithelial[Endoscopic_Severity_V1 %in% c("HC", "Severity_2_3")],
  exact = FALSE
)
epithelial_test <- data.table::data.table(
  n_pseudobulks = nrow(epithelial),
  spearman_rho = unname(epithelial_trend$estimate),
  spearman_p = epithelial_trend$p.value,
  severe_vs_hc_wilcoxon_p = epithelial_extreme$p.value
)

# 保存小型聚合结果，后续作图不再读取 1.9 GB 对象。
data.table::fwrite(group_information, file.path(DIR_RESULTS, "singlecell_pseudobulk.csv"))
data.table::fwrite(cell_type_shift, file.path(DIR_RESULTS, "singlecell_celltype_shift.csv"))
data.table::fwrite(epithelial_test, file.path(DIR_RESULTS, "singlecell_epithelial_test.csv"))

# 主动释放大对象，降低后续作图的内存占用。
rm(seurat_object, count_matrix, pseudobulk_counts, pseudobulk_logcpm, design_matrix)
invisible(gc())
print(epithelial_test)
