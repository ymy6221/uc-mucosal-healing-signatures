# CBM 升级阶段的技术验收、输入数据不可变性检查和评分工具参考分布。

DIR_CBM_RESULTS <- file.path(DIR_RESULTS, "cbm_upgrade")
DIR_CBM_FIGURES <- file.path(DIR_FIGURES, "cbm_upgrade")

required_results <- c(
  "pathway_registry_msigdb.csv", "pathway_overlap_audit.csv",
  "pathway_cohort_correlations.csv", "pathway_meta_correlations.csv",
  "pathway_pc1_sample_scores.csv", "meta_prediction_and_sensitivity.csv",
  "meta_leave_one_cohort_out.csv", "meta_influence_baujat.csv",
  "singlecell_fine_sample_pseudobulk.csv", "singlecell_fine_patient_scores.csv",
  "singlecell_fine_localization_tests.csv", "singlecell_fine_tests.csv",
  "singlecell_fine_umap.csv", "singlecell_resource_summary.csv",
  "singlecell_osm_axis_localization_tests.csv"
)
required_figures <- c(
  "CBM_Figure5_pathway_mechanism.pdf",
  "CBM_Figure5_pathway_mechanism.tiff",
  "CBM_Figure6_singlecell_fine_localization.pdf",
  "CBM_Figure6_singlecell_fine_localization.tiff",
  "CBM_FigureS2_heterogeneity_stability.pdf",
  "CBM_FigureS2_heterogeneity_stability.tiff",
  "CBM_FigureS3_OSM_OSMR_localization.pdf",
  "CBM_FigureS3_OSM_OSMR_localization.tiff",
  "CBM_Graphical_Abstract.png",
  "CBM_Graphical_Abstract.tiff"
)
result_paths <- file.path(DIR_CBM_RESULTS, required_results)
figure_paths <- file.path(DIR_CBM_FIGURES, required_figures)

pathway_meta <- data.table::fread(file.path(DIR_CBM_RESULTS,
                                            "pathway_meta_correlations.csv"))
pathway_cohort <- data.table::fread(file.path(DIR_CBM_RESULTS,
                                              "pathway_cohort_correlations.csv"))
overlap <- data.table::fread(file.path(DIR_CBM_RESULTS, "pathway_overlap_audit.csv"))
stability <- data.table::fread(file.path(DIR_CBM_RESULTS,
                                         "meta_prediction_and_sensitivity.csv"))
leave_one <- data.table::fread(file.path(DIR_CBM_RESULTS,
                                         "meta_leave_one_cohort_out.csv"))
singlecell <- data.table::fread(file.path(DIR_CBM_RESULTS,
                                          "singlecell_resource_summary.csv"))
patient_scores <- data.table::fread(file.path(DIR_CBM_RESULTS,
                                               "singlecell_fine_patient_scores.csv"))
localization <- data.table::fread(file.path(DIR_CBM_RESULTS,
                                             "singlecell_fine_localization_tests.csv"))

# 根据锁定清单重新核对 12 个 bulk 输入；另锁定单细胞文件的字节数和 SHA-256。
manifest <- data.table::fread(file.path(PROJECT_ROOT, "data",
                                        "public_download_manifest.csv"))
manifest[, path := file.path(DIR_RAW, filename)]
manifest <- data.table::rbindlist(list(
  manifest[, .(filename, path, expected_bytes = verified_bytes,
               expected_sha256 = sha256)],
  data.table::data.table(
    filename = "GSE250487_biopsy_RNA.rds.gz",
    path = file.path(PROJECT_ROOT, "data", "singlecell",
                     "GSE250487_biopsy_RNA.rds.gz"),
    expected_bytes = 408316301,
    expected_sha256 = "1bca3338b5ad1482332f21782dbfe392b547cee14fe1c4934f7da1e433e39f2f"
  )
), fill = TRUE)
manifest[, actual_bytes := ifelse(file.exists(path), file.info(path)$size, NA_real_)]
manifest[, actual_sha256 := vapply(path, function(one_path) {
  if (!file.exists(one_path)) return(NA_character_)
  digest::digest(file = one_path, algo = "sha256", serialize = FALSE)
}, character(1))]
manifest[, `:=`(
  bytes_passed = actual_bytes == expected_bytes,
  sha256_passed = actual_sha256 == expected_sha256
)]
manifest[, passed := bytes_passed & sha256_passed]
data.table::fwrite(manifest,
                   file.path(DIR_CBM_RESULTS, "public_input_immutability_check.csv"))

validation <- data.table::data.table(
  check = c(
    "all CBM upgrade result tables exist",
    "all CBM figures and graphical abstract exist",
    "12 locked pathways have full and overlap-excluded estimates",
    "every pathway variant contains eight cohort correlations",
    "overlap-excluded pathways retain at least ten measured-definition genes",
    "heterogeneity table covers seven signatures and four analyses",
    "leave-one-out table contains seven by eight estimates",
    "single-cell source has exactly 93,900 cells and 51 fine annotations",
    "single-cell inference uses patient units and UMAP remains below 20,000 cells",
    "four robust signatures are present in fine-cell localization",
    "all 13 public inputs retain locked byte counts and SHA-256 hashes",
    "upstream benchmark technical validation remains PASS"
  ),
  observed = c(
    paste0(sum(file.exists(result_paths) & file.info(result_paths)$size > 100L),
           "/", length(result_paths)),
    paste0(sum(file.exists(figure_paths) & file.info(figure_paths)$size > 100L),
           "/", length(figure_paths)),
    paste(pathway_meta[, .N, by = variant][, paste0(variant, "=", N)], collapse = ";"),
    paste(pathway_cohort[, .(minimum_cohorts = data.table::uniqueN(cohort)),
                         by = .(pathway, variant)][, min(minimum_cohorts)]),
    paste0("minimum ", min(overlap$genes_after_exclusion)),
    paste0(data.table::uniqueN(stability$signature_id), " signatures; ",
           data.table::uniqueN(stability$analysis), " analyses"),
    paste0(nrow(leave_one), " rows"),
    paste0(singlecell$source_cells, " cells; ",
           singlecell$author_fine_cell_types, " types"),
    paste0(data.table::uniqueN(patient_scores$patient_short),
           " patients; UMAP n=", singlecell$umap_cells),
    paste0(data.table::uniqueN(localization$signature_id), " signatures"),
    paste0(sum(manifest$passed), "/", nrow(manifest)),
    data.table::fread(file.path(DIR_RESULTS, "submission_readiness.csv"))$technical_reproducibility
  ),
  passed = c(
    all(file.exists(result_paths) & file.info(result_paths)$size > 100L),
    all(file.exists(figure_paths) & file.info(figure_paths)$size > 100L),
    nrow(pathway_meta[variant == "full"]) == 12L &&
      nrow(pathway_meta[variant == "overlap_excluded"]) == 12L,
    all(pathway_cohort[, data.table::uniqueN(cohort), by = .(pathway, variant)]$V1 == 8L),
    min(overlap$genes_after_exclusion) >= 10L,
    data.table::uniqueN(stability$signature_id) == 7L &&
      data.table::uniqueN(stability$analysis) == 4L,
    nrow(leave_one) == 56L,
    singlecell$source_cells == 93900L && singlecell$author_fine_cell_types == 51L,
    data.table::uniqueN(patient_scores$patient_short) == 12L &&
      singlecell$umap_cells <= 20000L && max(patient_scores$severity_order) == 2L,
    data.table::uniqueN(localization$signature_id) == 4L,
    all(manifest$passed),
    data.table::fread(file.path(DIR_RESULTS, "submission_readiness.csv"))$technical_reproducibility == "PASS"
  )
)
data.table::fwrite(validation,
                   file.path(DIR_CBM_RESULTS, "cbm_upgrade_validation.csv"))
if (!all(validation$passed)) {
  print(validation)
  stop("CBM 升级技术验收失败。")
}

# 为轻量评分脚本生成不含结局的公共参考分布。
scored_files <- file.path(DIR_INTERMEDIATE, paste0(
  "benchmark_",
  c("GSE206285", "GSE92415", "GSE16879", "GSE12251",
    "GSE73661_IFX", "GSE73661_VDZ", "GSE23597", "GSE72819"),
  "_scored.rds"
))
primary_ids <- unique(read_signature_registry()[primary_set == TRUE, signature_id])
robust_ids <- c("IFX_ARJS_5", "OSM_OSMR", "IL22_NEUTROPHIL_6", "VDZ_NR_3")
reference <- data.table::rbindlist(lapply(scored_files, function(path) {
  cohort_name <- sub("_scored\\.rds$", "", sub("^benchmark_", "", basename(path)))
  one <- data.table::as.data.table(readRDS(path))
  one[, robust_consensus_rank := rowMeans(.SD, na.rm = TRUE), .SDcols = robust_ids]
  one[, cohort := cohort_name]
  one[, c("cohort", "geo_accession", primary_ids, "robust_consensus_rank"), with = FALSE]
}), fill = TRUE)
reference_long <- data.table::melt(
  reference,
  id.vars = c("cohort", "geo_accession"),
  variable.name = "score_id", value.name = "score"
)
data.table::fwrite(reference_long,
                   file.path(DIR_CBM_RESULTS, "signature_reference_distribution.csv"))

version_packages <- c("R", "data.table", "GSVA", "BiocParallel", "msigdbr",
                      "Seurat", "SeuratObject", "Matrix", "ggplot2",
                      "patchwork", "ggrepel", "digest")
versions <- data.table::data.table(
  package = version_packages,
  version = c(as.character(getRversion()), vapply(version_packages[-1L], function(package) {
    as.character(utils::packageVersion(package))
  }, character(1)))
)
data.table::fwrite(versions,
                   file.path(DIR_CBM_RESULTS, "package_versions_cbm_upgrade.csv"))
capture.output(sessionInfo(),
               file = file.path(DIR_CBM_RESULTS, "sessionInfo_cbm_upgrade.txt"))
print(validation)
message("CBM 升级技术验收全部通过。")
