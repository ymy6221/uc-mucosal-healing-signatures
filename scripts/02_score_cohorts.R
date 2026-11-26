# 对每个队列逐一计算 ERI，保存后立即释放大矩阵，控制内存峰值。

check_packages()

# 读取预先锁定的基因表；本项目禁止根据响应结局增删基因。
signature_path <- file.path(PROJECT_ROOT, "signatures", "epithelial_restitution.csv")
signature_table <- data.table::fread(signature_path)

# 循环一次只处理一个队列。
for (i in seq_len(nrow(CORE_FILES))) {
  gse_id <- CORE_FILES$id[i]
  gpl_id <- CORE_FILES$platform[i]
  message("开始评分：", gse_id, " / ", gpl_id)
  eset <- read_geo_eset(gse_id)
  phenotype <- data.table::as.data.table(Biobase::pData(eset), keep.rownames = FALSE)
  gene_expression <- probe_to_gene_matrix(eset, gpl_id)
  score_table <- calculate_eri(gene_expression, signature_table)
  phenotype[, geo_accession := as.character(geo_accession)]
  output <- merge(phenotype, score_table, by = "geo_accession", all.x = TRUE, sort = FALSE)
  saveRDS(output, file.path(DIR_INTERMEDIATE, paste0(gse_id, "_scored.rds")), compress = "xz")
  coverage <- unique(score_table[, .(injury_genes_n, mature_genes_n)])
  coverage[, cohort := gse_id]
  data.table::fwrite(coverage, file.path(DIR_RESULTS, paste0(gse_id, "_signature_coverage.csv")))
  rm(eset, phenotype, gene_expression, score_table, output)
  invisible(gc())
  message("完成评分：", gse_id)
}

