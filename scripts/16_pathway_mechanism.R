# CBM 升级 1：把 outcome-free 的签名收敛轴翻译为可解释的分子通路。
#
# 设计原则：
# 1) 每个 GEO 队列单独评分，不合并原始表达矩阵；
# 2) 通路分数不使用临床结局；
# 3) 相关系数在队列层面计算，再用 Fisher z 随机效应合并；
# 4) 同时删除七个主签名的全部基因，排除机械性的基因重叠解释；
# 5) 结果写入 results/cbm_upgrade，不覆盖既有结果或公共数据。

upgrade_packages <- c("data.table", "GSVA", "BiocParallel", "msigdbr")
upgrade_ok <- vapply(upgrade_packages, requireNamespace, logical(1), quietly = TRUE)
if (!all(upgrade_ok)) {
  stop("缺少通路分析包：", paste(upgrade_packages[!upgrade_ok], collapse = ", "),
       "。请运行 scripts/21_install_cbm_upgrade_packages.R")
}

DIR_CBM_RESULTS <- file.path(DIR_RESULTS, "cbm_upgrade")
DIR_CBM_FIGURES <- file.path(DIR_FIGURES, "cbm_upgrade")
dir.create(DIR_CBM_RESULTS, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_CBM_FIGURES, recursive = TRUE, showWarnings = FALSE)

# 事先锁定 12 条与炎症、基质重塑和屏障修复直接相关的通路。
pathway_specification <- data.table::data.table(
  collection = c(rep("H", 8L), rep("C2", 4L)),
  subcollection = c(rep("", 8L), rep("CP:KEGG_LEGACY", 4L)),
  pathway = c(
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_HYPOXIA",
    "HALLMARK_APICAL_JUNCTION",
    "KEGG_ECM_RECEPTOR_INTERACTION",
    "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
    "KEGG_FOCAL_ADHESION",
    "KEGG_TIGHT_JUNCTION"
  ),
  display_name = c(
    "Inflammatory response", "TNF-alpha/NF-kappaB", "IL-6/JAK/STAT3",
    "Epithelial-mesenchymal transition", "Interferon-gamma response",
    "Complement", "Hypoxia", "Apical junction", "ECM-receptor interaction",
    "Cytokine-receptor interaction", "Focal adhesion", "Tight junction"
  )
)

hallmark <- data.table::as.data.table(msigdbr::msigdbr(
  species = "Homo sapiens", collection = "H"
))
kegg <- data.table::as.data.table(msigdbr::msigdbr(
  species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY"
))
msigdb <- data.table::rbindlist(list(hallmark, kegg), fill = TRUE)
if (!all(pathway_specification$pathway %in% msigdb$gs_name)) {
  stop("MSigDB 版本中缺少锁定通路：",
       paste(setdiff(pathway_specification$pathway, msigdb$gs_name), collapse = ", "))
}

pathway_registry <- msigdb[gs_name %in% pathway_specification$pathway,
                           .(pathway = gs_name, gene = gene_symbol,
                             database_version = db_version,
                             collection = gs_collection,
                             subcollection = gs_subcollection)]
pathway_registry <- unique(pathway_registry)
pathway_registry <- merge(pathway_registry, pathway_specification,
                          by = c("pathway", "collection", "subcollection"), all.x = TRUE)
data.table::setorder(pathway_registry, pathway, gene)
data.table::fwrite(pathway_registry,
                   file.path(DIR_CBM_RESULTS, "pathway_registry_msigdb.csv"))

registry <- read_signature_registry()
primary_ids <- unique(registry[primary_set == TRUE, signature_id])
primary_genes <- unique(registry[primary_set == TRUE, gene])
loadings <- data.table::fread(file.path(DIR_RESULTS, "signature_convergence_loadings.csv"))
loadings <- loadings[signature_id %in% primary_ids]
if (!setequal(loadings$signature_id, primary_ids)) {
  stop("PC1 载荷与七个主签名不一致；请先运行 run_stage5_benchmark_upgrade.R")
}

full_gene_sets <- split(pathway_registry$gene, pathway_registry$pathway)
excluded_gene_sets <- lapply(full_gene_sets, setdiff, y = primary_genes)
overlap_summary <- data.table::data.table(
  pathway = names(full_gene_sets),
  genes_full = lengths(full_gene_sets),
  primary_signature_overlap = lengths(Map(intersect, full_gene_sets,
                                           MoreArgs = list(y = primary_genes))),
  genes_after_exclusion = lengths(excluded_gene_sets)
)
overlap_summary <- merge(overlap_summary,
                         pathway_specification[, .(pathway, display_name)],
                         by = "pathway", all.x = TRUE)
data.table::fwrite(overlap_summary,
                   file.path(DIR_CBM_RESULTS, "pathway_overlap_audit.csv"))

source_ids <- c("GSE206285", "GSE92415", "GSE16879", "GSE12251",
                "GSE73661", "GSE23597", "GSE72819")
cohort_correlations <- list()
sample_scores <- list()
correlation_index <- 0L
sample_index <- 0L

for (source_id in source_ids) {
  message("通路评分：", source_id)
  cohort_objects <- prepare_source_cohorts(source_id)
  for (object in cohort_objects) {
    scored <- score_registered_signatures(object$expression, registry)
    signature_scores <- data.table::copy(scored$scores)
    signature_scores <- signature_scores[match(colnames(object$expression), geo_accession)]

    standardized <- sapply(primary_ids, function(signature_id) {
      values <- signature_scores[[signature_id]]
      if (sum(is.finite(values)) < 2L || stats::sd(values, na.rm = TRUE) == 0) {
        rep(0, length(values))
      } else {
        values <- as.numeric(scale(values))
        values[!is.finite(values)] <- 0
        values
      }
    })
    if (is.null(dim(standardized))) standardized <- matrix(standardized, ncol = 1L)
    colnames(standardized) <- primary_ids
    loading_vector <- loadings$pc1_loading[match(primary_ids, loadings$signature_id)]
    pc1_score <- as.numeric(standardized %*% loading_vector)

    measured_full <- lapply(full_gene_sets, intersect, y = rownames(object$expression))
    measured_excluded <- lapply(excluded_gene_sets, intersect, y = rownames(object$expression))
    combined_gene_sets <- c(
      stats::setNames(measured_full, paste0("full::", names(measured_full))),
      stats::setNames(measured_excluded, paste0("overlap_excluded::", names(measured_excluded)))
    )
    combined_gene_sets <- combined_gene_sets[lengths(combined_gene_sets) >= 10L]
    gsva_parameters <- GSVA::gsvaParam(
      exprData = object$expression,
      geneSets = combined_gene_sets,
      minSize = 10L,
      maxSize = 500L,
      kcdf = "Gaussian",
      sparse = FALSE,
      verbose = FALSE
    )
    pathway_scores <- GSVA::gsva(
      gsva_parameters,
      verbose = FALSE,
      BPPARAM = BiocParallel::SerialParam()
    )

    phenotype <- data.table::copy(object$phenotype)
    phenotype <- phenotype[match(colnames(object$expression), geo_accession)]
    sample_table <- data.table::data.table(
      cohort = object$cohort,
      geo_accession = colnames(object$expression),
      treatment = phenotype$treatment,
      pc1 = pc1_score
    )
    score_wide <- data.table::as.data.table(t(pathway_scores), keep.rownames = "geo_accession")
    score_wide <- score_wide[match(sample_table$geo_accession, geo_accession)]
    score_wide[, geo_accession := NULL]
    sample_index <- sample_index + 1L
    sample_scores[[sample_index]] <- cbind(sample_table, score_wide)

    for (score_name in rownames(pathway_scores)) {
      split_name <- strsplit(score_name, "::", fixed = TRUE)[[1L]]
      variant <- split_name[1L]
      pathway <- split_name[2L]
      values <- as.numeric(pathway_scores[score_name, ])
      keep <- stats::complete.cases(values, pc1_score)
      if (sum(keep) < 10L || stats::sd(values[keep]) == 0 || stats::sd(pc1_score[keep]) == 0) next
      test <- suppressWarnings(stats::cor.test(
        values[keep], pc1_score[keep], method = "spearman", exact = FALSE
      ))
      rho <- unname(test$estimate)
      rho <- min(0.999999, max(-0.999999, rho))
      correlation_index <- correlation_index + 1L
      cohort_correlations[[correlation_index]] <- data.table::data.table(
        cohort = object$cohort,
        pathway = pathway,
        variant = variant,
        n = sum(keep),
        rho = rho,
        fisher_z = atanh(rho),
        se = 1 / sqrt(sum(keep) - 3),
        p_value = test$p.value,
        genes_measured = length(combined_gene_sets[[score_name]])
      )
    }
    rm(scored, standardized, pathway_scores, object)
    invisible(gc())
  }
  rm(cohort_objects)
  invisible(gc())
}

cohort_correlations <- data.table::rbindlist(cohort_correlations, fill = TRUE)
sample_scores <- data.table::rbindlist(sample_scores, fill = TRUE)
cohort_correlations <- merge(
  cohort_correlations,
  pathway_specification[, .(pathway, display_name)],
  by = "pathway", all.x = TRUE
)

meta_results <- cohort_correlations[, {
  meta <- random_effects_meta(fisher_z, se)
  k <- .N
  prediction_multiplier <- if (k >= 3L) stats::qt(0.975, df = k - 2L) else NA_real_
  prediction_se <- sqrt(meta$tau2 + meta$se^2)
  .(
    cohorts = k,
    total_n = sum(n),
    fisher_z = meta$log_or,
    rho = tanh(meta$log_or),
    ci_low = tanh(meta$log_or - 1.96 * meta$se),
    ci_high = tanh(meta$log_or + 1.96 * meta$se),
    p_value = meta$p_value,
    tau2 = meta$tau2,
    i2 = meta$i2,
    prediction_low = if (k >= 3L) tanh(meta$log_or - prediction_multiplier * prediction_se) else NA_real_,
    prediction_high = if (k >= 3L) tanh(meta$log_or + prediction_multiplier * prediction_se) else NA_real_
  )
}, by = .(pathway, display_name, variant)]
meta_results[, fdr := stats::p.adjust(p_value, method = "BH"), by = variant]
data.table::setorder(meta_results, variant, -rho)

data.table::fwrite(cohort_correlations,
                   file.path(DIR_CBM_RESULTS, "pathway_cohort_correlations.csv"))
data.table::fwrite(meta_results,
                   file.path(DIR_CBM_RESULTS, "pathway_meta_correlations.csv"))
data.table::fwrite(sample_scores,
                   file.path(DIR_CBM_RESULTS, "pathway_pc1_sample_scores.csv"))

message("通路机制分析完成：", nrow(meta_results), " 个合并结果。")
print(meta_results)
