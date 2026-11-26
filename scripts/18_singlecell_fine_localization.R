# CBM 升级 3：用 GSE250487 的作者提供 51 类精细注释定位四个稳健签名。
#
# 推断使用全部 93,900 个细胞的患者×亚群伪汇总；UMAP 只用于展示，并在
# 上皮/基质细胞内分层抽样至不超过 20,000 个细胞，避免普通电脑内存过载。

single_cell_packages <- c("Matrix", "SeuratObject", "Seurat", "data.table")
single_cell_ok <- vapply(single_cell_packages, requireNamespace, logical(1), quietly = TRUE)
if (!all(single_cell_ok)) {
  stop("缺少精细单细胞分析包：",
       paste(single_cell_packages[!single_cell_ok], collapse = ", "),
       "。请运行 scripts/21_install_cbm_upgrade_packages.R")
}

DIR_CBM_RESULTS <- file.path(DIR_RESULTS, "cbm_upgrade")
DIR_CBM_FIGURES <- file.path(DIR_FIGURES, "cbm_upgrade")
dir.create(DIR_CBM_RESULTS, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_CBM_FIGURES, recursive = TRUE, showWarnings = FALSE)

single_cell_path <- file.path(PROJECT_ROOT, "data", "singlecell",
                             "GSE250487_biopsy_RNA.rds.gz")
if (!file.exists(single_cell_path) || file.info(single_cell_path)$size < 300000000L) {
  stop("缺少完整 GSE250487 RDS；请先运行 run_stage3_singlecell.R")
}

read_double_gzip_rds <- function(path) {
  outer_connection <- gzfile(path, open = "rb")
  inner_connection <- gzcon(outer_connection)
  on.exit(close(inner_connection), add = TRUE)
  readRDS(inner_connection)
}

message("读取 GSE250487（只读取一次，约需 2 GB 内存）。")
seurat_object <- read_double_gzip_rds(single_cell_path)
required_metadata <- c(
  "CoLabs_sample", "patient_short", "coarse_annotations_MK",
  "fine_annotations_MK_V4", "condition", "Endoscopic_Severity_V1"
)
missing_metadata <- setdiff(required_metadata, colnames(seurat_object@meta.data))
if (length(missing_metadata) > 0L) {
  stop("GSE250487 缺少预期字段：", paste(missing_metadata, collapse = ", "))
}

cell_metadata <- data.table::as.data.table(seurat_object@meta.data,
                                           keep.rownames = "cell")
cell_metadata[, group := paste(CoLabs_sample, fine_annotations_MK_V4, sep = "___")]
group_sizes <- cell_metadata[, .(n_cells = .N), by = group]
valid_groups <- group_sizes[n_cells >= 20L, group]
aggregation_metadata <- cell_metadata[group %in% valid_groups]
aggregation_metadata[, group_factor := factor(group, levels = valid_groups)]

design_matrix <- Matrix::sparse.model.matrix(~ 0 + group_factor,
                                              data = aggregation_metadata)
colnames(design_matrix) <- sub("^group_factor", "", colnames(design_matrix))
count_matrix <- SeuratObject::LayerData(seurat_object[["RNA"]], layer = "counts")
count_matrix <- count_matrix[, aggregation_metadata$cell, drop = FALSE]
pseudobulk_counts <- count_matrix %*% design_matrix
pseudobulk_library <- as.numeric(Matrix::colSums(pseudobulk_counts))
pseudobulk_logcpm <- log2(sweep(as.matrix(pseudobulk_counts), 2L,
                                pseudobulk_library, "/") * 1e6 + 1)

registry <- read_signature_registry()
robust_ids <- c("IFX_ARJS_5", "OSM_OSMR", "IL22_NEUTROPHIL_6", "VDZ_NR_3")
robust_registry <- registry[signature_id %in% robust_ids]
scored <- score_registered_signatures(pseudobulk_logcpm, robust_registry)
pseudobulk_scores <- data.table::copy(scored$scores)
data.table::setnames(pseudobulk_scores, "geo_accession", "group")
coverage <- data.table::copy(scored$coverage)
rm(scored)

group_information <- unique(aggregation_metadata[, .(
  group, CoLabs_sample, patient_short, coarse_annotations_MK,
  fine_annotations_MK_V4, condition, Endoscopic_Severity_V1
)])
group_information <- merge(group_information, group_sizes, by = "group",
                           all.x = TRUE, sort = FALSE)
group_information <- merge(group_information, pseudobulk_scores, by = "group",
                           all.x = TRUE, sort = FALSE)
group_information[, severity_order := data.table::fcase(
  Endoscopic_Severity_V1 == "HC", 0,
  Endoscopic_Severity_V1 == "Severity_0_1", 1,
  Endoscopic_Severity_V1 == "Severity_2_3", 2,
  default = NA_real_
)]

# 同一患者左右结肠的同一亚群取中位数；严重度采用该患者两处活检的较高等级。
patient_scores <- group_information[!is.na(severity_order), c(
  list(
    coarse_cell_type = names(sort(table(coarse_annotations_MK), decreasing = TRUE))[1L],
    n_samples = data.table::uniqueN(CoLabs_sample),
    n_cells = sum(n_cells),
    severity_order = max(severity_order),
    condition = paste(sort(unique(condition)), collapse = ";")
  ),
  lapply(.SD, stats::median, na.rm = TRUE)
), by = .(patient_short, fine_annotations_MK_V4),
.SDcols = robust_ids]
patient_scores[, severity_label := data.table::fcase(
  severity_order == 0, "Healthy control",
  severity_order == 1, "UC severity 0-1",
  severity_order == 2, "UC severity 2-3",
  default = NA_character_
)]

# 主要定位检验：在每位患者内部把各亚群分数减去该患者所有可测亚群的中位数，
# 再检验某亚群是否在患者间稳定高于 0。这样检验的是“信号来源”，而不是把细胞
# 当成独立重复；严重度趋势保留为次要、探索性分析。
localization_scores <- data.table::rbindlist(lapply(robust_ids, function(sig_current) {
  one <- patient_scores[, .(
    patient_short, fine_annotations_MK_V4, coarse_cell_type,
    signature_id = sig_current,
    score = get(sig_current)
  )]
  one[, centered_score := score - stats::median(score, na.rm = TRUE),
      by = patient_short]
  one
}))
localization_tests <- localization_scores[stats::complete.cases(centered_score), {
  if (.N < 5L || stats::sd(centered_score) == 0) {
    list(n_patients = .N, median_centered_score = stats::median(centered_score),
         localization_p = NA_real_)
  } else {
    test <- suppressWarnings(stats::wilcox.test(
      centered_score, mu = 0, alternative = "greater", exact = FALSE
    ))
    list(n_patients = .N, median_centered_score = stats::median(centered_score),
         localization_p = test$p.value)
  }
}, by = .(signature_id, fine_cell_type = fine_annotations_MK_V4,
          coarse_cell_type)]
localization_tests[, localization_fdr := stats::p.adjust(localization_p, method = "BH")]

fine_tests <- list()
test_index <- 0L
for (sig_current in robust_ids) {
  for (cell_type in unique(patient_scores$fine_annotations_MK_V4)) {
    one <- patient_scores[fine_annotations_MK_V4 == cell_type &
                            stats::complete.cases(get(sig_current), severity_order)]
    n_hc <- sum(one$severity_order == 0)
    n_severe <- sum(one$severity_order == 2)
    if (nrow(one) < 6L || data.table::uniqueN(one$severity_order) < 2L ||
        stats::sd(one[[sig_current]]) == 0) next
    trend <- suppressWarnings(stats::cor.test(one[[sig_current]], one$severity_order,
                                              method = "spearman", exact = FALSE))
    wilcoxon_p <- NA_real_
    median_shift <- NA_real_
    if (n_hc >= 2L && n_severe >= 2L) {
      extremes <- one[severity_order %in% c(0, 2)]
      wilcoxon_p <- suppressWarnings(stats::wilcox.test(
        extremes[[sig_current]] ~ extremes$severity_order, exact = FALSE
      )$p.value)
      median_shift <- stats::median(one[severity_order == 2, get(sig_current)]) -
        stats::median(one[severity_order == 0, get(sig_current)])
    }
    test_index <- test_index + 1L
    fine_tests[[test_index]] <- data.table::data.table(
      signature_id = sig_current,
      fine_cell_type = cell_type,
      coarse_cell_type = one$coarse_cell_type[1L],
      n_patients = nrow(one),
      n_hc = n_hc,
      n_severe = n_severe,
      spearman_rho = unname(trend$estimate),
      spearman_p = trend$p.value,
      severe_minus_hc = median_shift,
      severe_vs_hc_p = wilcoxon_p
    )
  }
}
fine_tests <- data.table::rbindlist(fine_tests, fill = TRUE)
fine_tests[, spearman_fdr := stats::p.adjust(spearman_p, method = "BH")]
fine_tests[, severe_vs_hc_fdr := stats::p.adjust(severe_vs_hc_p, method = "BH")]

# 保存构成四个签名的单基因 logCPM，特别用于区分 OSM 配体和 OSMR 受体。
robust_genes <- unique(robust_registry$gene)
measured_robust_genes <- intersect(robust_genes, rownames(pseudobulk_logcpm))
gene_expression <- data.table::as.data.table(
  t(pseudobulk_logcpm[measured_robust_genes, , drop = FALSE]),
  keep.rownames = "group"
)
gene_expression <- merge(group_information[, .(
  group, CoLabs_sample, patient_short, coarse_annotations_MK,
  fine_annotations_MK_V4, Endoscopic_Severity_V1, severity_order, n_cells
)], gene_expression, by = "group", all.x = TRUE, sort = FALSE)

# OSM 配体和 OSMR 受体可能来自不同细胞群，因此分开做患者内定位检验。
axis_genes <- intersect(c("OSM", "OSMR"), names(gene_expression))
axis_gene_patient <- gene_expression[, lapply(.SD, stats::median, na.rm = TRUE),
                                     by = .(patient_short, fine_annotations_MK_V4,
                                            coarse_annotations_MK),
                                     .SDcols = axis_genes]
axis_gene_long <- data.table::melt(
  axis_gene_patient,
  id.vars = c("patient_short", "fine_annotations_MK_V4", "coarse_annotations_MK"),
  variable.name = "gene", value.name = "logcpm"
)
axis_gene_long[, centered_logcpm := logcpm - stats::median(logcpm, na.rm = TRUE),
               by = .(patient_short, gene)]
axis_gene_tests <- axis_gene_long[stats::complete.cases(centered_logcpm), {
  if (.N < 5L || stats::sd(centered_logcpm) == 0) {
    list(n_patients = .N, median_centered_logcpm = stats::median(centered_logcpm),
         localization_p = NA_real_)
  } else {
    test <- suppressWarnings(stats::wilcox.test(
      centered_logcpm, mu = 0, alternative = "greater", exact = FALSE
    ))
    list(n_patients = .N, median_centered_logcpm = stats::median(centered_logcpm),
         localization_p = test$p.value)
  }
}, by = .(gene, fine_cell_type = fine_annotations_MK_V4,
          coarse_cell_type = coarse_annotations_MK)]
axis_gene_tests[, localization_fdr := stats::p.adjust(localization_p, method = "BH")]

data.table::fwrite(group_information,
                   file.path(DIR_CBM_RESULTS, "singlecell_fine_sample_pseudobulk.csv"))
data.table::fwrite(patient_scores,
                   file.path(DIR_CBM_RESULTS, "singlecell_fine_patient_scores.csv"))
data.table::fwrite(fine_tests,
                   file.path(DIR_CBM_RESULTS, "singlecell_fine_tests.csv"))
data.table::fwrite(localization_scores,
                   file.path(DIR_CBM_RESULTS, "singlecell_fine_localization_scores.csv"))
data.table::fwrite(localization_tests,
                   file.path(DIR_CBM_RESULTS, "singlecell_fine_localization_tests.csv"))
data.table::fwrite(gene_expression,
                   file.path(DIR_CBM_RESULTS, "singlecell_fine_gene_expression.csv"))
data.table::fwrite(axis_gene_long,
                   file.path(DIR_CBM_RESULTS, "singlecell_osm_axis_patient_expression.csv"))
data.table::fwrite(axis_gene_tests,
                   file.path(DIR_CBM_RESULTS, "singlecell_osm_axis_localization_tests.csv"))
data.table::fwrite(coverage,
                   file.path(DIR_CBM_RESULTS, "singlecell_fine_signature_coverage.csv"))

# ---------- 轻量 UMAP：仅上皮和基质，按作者亚群分层抽样 ----------
set.seed(SEED)
visual_metadata <- cell_metadata[
  coarse_annotations_MK %in% c("10-Epithelial", "11-Stromal-fibroblast")
]
selected_cells <- visual_metadata[, {
  take <- min(.N, 1200L)
  list(cell = sample(cell, take))
}, by = fine_annotations_MK_V4]$cell
if (length(selected_cells) > 20000L) {
  selected_cells <- sample(selected_cells, 20000L)
}
visual_metadata <- visual_metadata[match(selected_cells, cell)]
visual_counts <- SeuratObject::LayerData(seurat_object[["RNA"]], layer = "counts")
visual_counts <- visual_counts[, selected_cells, drop = FALSE]
visual_meta_frame <- as.data.frame(visual_metadata)
rownames(visual_meta_frame) <- visual_meta_frame$cell

visual_object <- SeuratObject::CreateSeuratObject(
  counts = visual_counts,
  meta.data = visual_meta_frame,
  project = "GSE250487_epithelial_stromal"
)
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")
visual_object <- Seurat::NormalizeData(visual_object, verbose = FALSE)
visual_object <- Seurat::FindVariableFeatures(
  visual_object, selection.method = "vst", nfeatures = 2000L, verbose = FALSE
)
visual_object <- Seurat::ScaleData(
  visual_object, features = SeuratObject::VariableFeatures(visual_object),
  verbose = FALSE
)
visual_object <- Seurat::RunPCA(
  visual_object, features = SeuratObject::VariableFeatures(visual_object),
  npcs = 25L, verbose = FALSE
)
visual_object <- Seurat::RunUMAP(
  visual_object, dims = 1:20, n.neighbors = 30L, min.dist = 0.30,
  seed.use = SEED, verbose = FALSE
)

module_features <- lapply(robust_ids, function(sig_current) {
  intersect(robust_registry[robust_registry$signature_id == sig_current, gene],
            rownames(visual_object))
})
names(module_features) <- robust_ids
visual_object <- Seurat::AddModuleScore(
  visual_object,
  features = unname(module_features),
  name = "robust_signature_",
  ctrl = 50L,
  seed = SEED,
  search = FALSE
)
umap <- data.table::as.data.table(Seurat::Embeddings(visual_object, "umap"),
                                  keep.rownames = "cell")
umap_metadata <- data.table::as.data.table(visual_object@meta.data,
                                           keep.rownames = "cell")
for (i in seq_along(robust_ids)) {
  data.table::setnames(umap_metadata, paste0("robust_signature_", i), robust_ids[i])
}
umap <- merge(umap, umap_metadata[, c(
  "cell", "patient_short", "CoLabs_sample", "condition",
  "Endoscopic_Severity_V1", "coarse_annotations_MK",
  "fine_annotations_MK_V4", robust_ids
), with = FALSE], by = "cell", all.x = TRUE, sort = FALSE)
data.table::fwrite(umap, file.path(DIR_CBM_RESULTS, "singlecell_fine_umap.csv"))

resource_summary <- data.table::data.table(
  source_cells = ncol(seurat_object),
  source_genes = nrow(seurat_object),
  valid_sample_by_subtype_pseudobulks = nrow(group_information),
  patient_by_subtype_units = nrow(patient_scores),
  author_fine_cell_types = data.table::uniqueN(cell_metadata$fine_annotations_MK_V4),
  umap_cells = nrow(umap),
  umap_scope = "author-annotated epithelial and stromal cells only",
  minimum_cells_per_pseudobulk = 20L
)
data.table::fwrite(resource_summary,
                   file.path(DIR_CBM_RESULTS, "singlecell_resource_summary.csv"))

rm(seurat_object, count_matrix, pseudobulk_counts, pseudobulk_logcpm,
   design_matrix, visual_counts, visual_object)
invisible(gc())
message("精细单细胞定位完成。")
print(resource_summary)
