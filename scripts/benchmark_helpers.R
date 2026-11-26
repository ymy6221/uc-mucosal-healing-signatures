# 高影响力升级阶段的专用函数。所有队列都在各自平台内转换为样本内基因百分位秩，
# 不合并原始表达值，因此避免了不同平台的批次校正和高内存消耗。

check_benchmark_packages <- function() {
  needed <- c("GEOquery", "Biobase", "data.table", "pROC", "Matrix",
              "AnnotationDbi", "org.Hs.eg.db")
  ok <- vapply(needed, requireNamespace, logical(1), quietly = TRUE)
  if (!all(ok)) {
    stop("缺少升级阶段 R 包：", paste(needed[!ok], collapse = ", "),
         "。请先运行 scripts/00_install_packages.R")
  }
  invisible(TRUE)
}

# 读入已发表签名注册表，并追加本项目事先锁定的 ERI 作为内部对照。
read_signature_registry <- function() {
  published <- data.table::fread(
    file.path(PROJECT_ROOT, "signatures", "published_signature_registry.csv")
  )
  # 主分析只允许作者原文中可逐字恢复的基因集合。缺少原始系数的代理集合、
  # 会议摘要集合和本项目重建集合仍会完整计算，但只能进入敏感性分析。
  published[definition_status != "exact_text", primary_set := FALSE]
  eri <- data.table::fread(file.path(PROJECT_ROOT, "signatures", "epithelial_restitution.csv"))
  eri_rows <- eri[, .(
    signature_id = "ERI_LOCKED_INTERNAL",
    display_name = "Locked epithelial restitution index",
    gene,
    component_weight = ifelse(module == "injury", 1, -1),
    risk_orientation = 1,
    primary_set = FALSE,
    source_doi = NA_character_,
    source_pmid = NA_character_,
    source_cohorts = "none_outcome_selected",
    definition_status = "locked_internal",
    note = "Genes locked before outcome analysis; internal comparator rather than published classifier"
  )]
  data.table::rbindlist(list(published, eri_rows), fill = TRUE)
}

# 读取一个微阵列队列，并用 GEO 官方注释转为基因符号矩阵。
load_microarray_source <- function(gse_id, gpl_id) {
  eset <- read_geo_eset(gse_id)
  list(
    expression = probe_to_gene_matrix(eset, gpl_id),
    phenotype = data.table::as.data.table(Biobase::pData(eset), keep.rownames = FALSE)
  )
}

# 解析 GSE72819 官方 TAR 中 73 个逐样本 RPKM 文件。临时解压目录会在函数退出时删除。
load_gse72819_source <- function() {
  metadata_eset <- read_geo_eset("GSE72819")
  phenotype <- data.table::as.data.table(Biobase::pData(metadata_eset), keep.rownames = FALSE)
  tar_path <- file.path(DIR_RAW, "GSE72819_RAW.tar")
  if (!file.exists(tar_path)) {
    stop("缺少 GSE72819_RAW.tar；请先运行 scripts/10_download_benchmark.R")
  }

  extraction_dir <- tempfile(pattern = "GSE72819_extract_", tmpdir = tempdir())
  dir.create(extraction_dir, recursive = TRUE)
  on.exit(unlink(extraction_dir, recursive = TRUE, force = TRUE), add = TRUE)
  utils::untar(tar_path, exdir = extraction_dir)
  files <- list.files(extraction_dir, pattern = "\\.tab\\.gz$", full.names = TRUE)
  if (length(files) != 73L) {
    stop("GSE72819 预期 73 个处理后样本文件，实际为 ", length(files))
  }
  accessions <- sub("_.*$", "", basename(files))

  first <- data.table::fread(files[1L], select = c("name", "rpkm"), showProgress = FALSE)
  gene_id <- sub("^GeneID:", "", first$name)
  expression <- matrix(NA_real_, nrow = nrow(first), ncol = length(files),
                       dimnames = list(gene_id, accessions))
  expression[, 1L] <- first$rpkm
  if (length(files) > 1L) {
    for (i in 2:length(files)) {
      one <- data.table::fread(files[i], select = c("name", "rpkm"), showProgress = FALSE)
      one_id <- sub("^GeneID:", "", one$name)
      if (!identical(one_id, gene_id)) {
        expression[, i] <- one$rpkm[match(gene_id, one_id)]
      } else {
        expression[, i] <- one$rpkm
      }
    }
  }

  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = rownames(expression),
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  )
  keep <- !is.na(symbols) & symbols != ""
  expression <- log2(expression[keep, , drop = FALSE] + 1)
  symbols <- as.character(symbols[keep])
  row_iqr <- apply(expression, 1L, stats::IQR, na.rm = TRUE)
  order_table <- data.table::data.table(row_number = seq_along(symbols), symbol = symbols, iqr = row_iqr)
  data.table::setorder(order_table, symbol, -iqr)
  order_table <- order_table[!duplicated(symbol)]
  expression <- expression[order_table$row_number, , drop = FALSE]
  rownames(expression) <- order_table$symbol
  list(expression = expression, phenotype = phenotype)
}

# 保证表达矩阵和表型表样本一一对应，并且只保留有明确结局的基线样本。
make_cohort_object <- function(expression, phenotype, cohort, active_treatment) {
  phenotype <- phenotype[!duplicated(geo_accession)]
  phenotype <- phenotype[stats::complete.cases(geo_accession, treatment, failure)]
  phenotype <- phenotype[geo_accession %in% colnames(expression)]
  phenotype <- phenotype[match(intersect(colnames(expression), phenotype$geo_accession), geo_accession)]
  expression <- expression[, phenotype$geo_accession, drop = FALSE]
  if (!identical(colnames(expression), phenotype$geo_accession)) {
    stop(cohort, " 的表达矩阵和表型表顺序不一致。")
  }
  list(cohort = cohort, active_treatment = active_treatment,
       expression = expression, phenotype = phenotype)
}

# 将各 GEO 的原始字段统一为 treatment/failure。failure=1 始终表示未应答或未缓解。
prepare_source_cohorts <- function(source_id) {
  if (source_id == "GSE206285") {
    source <- load_microarray_source("GSE206285", "GPL13158")
    p <- source$phenotype
    p[, geo_accession := as.character(geo_accession)]
    p[, treatment := data.table::fcase(
      grepl("Ustekinumab", get("treatment:ch1")), "ustekinumab",
      grepl("Placebo", get("treatment:ch1")), "placebo",
      default = NA_character_
    )]
    p[, failure := data.table::fcase(
      get("mucosal healing at week 8:ch1") == "N", 1L,
      get("mucosal healing at week 8:ch1") == "Y", 0L,
      default = NA_integer_
    )]
    clinical <- data.table::fread(file.path(DIR_RAW, "GSE206285_additional_sample_metadata.tsv.gz"))
    data.table::setnames(
      clinical,
      c("Accession", "Field[CRP]", "Field[Calprotectin]", "Field[Mayo_score]",
        "Field[Disease_duration]", "Field[Anti_TNF_history]", "Field[Neutrophils_LP]",
        "Field[Neutrophils_Epithelium]"),
      c("geo_accession", "crp", "calprotectin", "mayo", "disease_duration",
        "anti_tnf_history", "neutrophils_lp", "neutrophils_epithelium")
    )
    p <- merge(p, clinical, by = "geo_accession", all.x = TRUE, sort = FALSE)
    for (column in c("crp", "calprotectin", "mayo", "disease_duration",
                     "neutrophils_lp", "neutrophils_epithelium")) {
      data.table::set(p, j = column, value = as_numeric_safely(p[[column]]))
    }
    return(list(make_cohort_object(source$expression, p, "GSE206285", "ustekinumab")))
  }

  if (source_id == "GSE92415") {
    source <- load_microarray_source("GSE92415", "GPL13158")
    p <- source$phenotype[get("visit:ch1") == "Week 0" &
                            get("treatment:ch1") %in% c("golimumab", "Placebo")]
    p[, geo_accession := as.character(geo_accession)]
    p[, treatment := tolower(get("treatment:ch1"))]
    p[, failure := data.table::fcase(
      get("wk6response:ch1") == "No", 1L,
      get("wk6response:ch1") == "Yes", 0L,
      default = NA_integer_
    )]
    p[, mayo := as_numeric_safely(get("mayo score:ch1"))]
    p[, age := as_numeric_safely(get("age:ch1"))]
    return(list(make_cohort_object(source$expression, p, "GSE92415", "golimumab")))
  }

  if (source_id == "GSE16879") {
    source <- load_microarray_source("GSE16879", "GPL570")
    p <- source$phenotype[
      get("disease:ch1") == "UC" &
        grepl("Before", get("before or after first infliximab treatment:ch1"))
    ]
    p[, geo_accession := as.character(geo_accession)]
    p[, treatment := "infliximab"]
    p[, failure := data.table::fcase(
      get("response to infliximab:ch1") == "No", 1L,
      get("response to infliximab:ch1") == "Yes", 0L,
      default = NA_integer_
    )]
    return(list(make_cohort_object(source$expression, p, "GSE16879", "infliximab")))
  }

  if (source_id == "GSE12251") {
    source <- load_microarray_source("GSE12251", "GPL570")
    p <- source$phenotype
    p[, geo_accession := as.character(geo_accession)]
    p[, patient_id := sub("/.*$", "", title)]
    p <- p[!duplicated(patient_id)]
    p[, treatment := "infliximab"]
    p[, failure := data.table::fcase(
      get("WK8RSPHM:ch1") == "No", 1L,
      get("WK8RSPHM:ch1") == "Yes", 0L,
      default = NA_integer_
    )]
    return(list(make_cohort_object(source$expression, p, "GSE12251", "infliximab")))
  }

  if (source_id == "GSE73661") {
    source <- load_microarray_source("GSE73661", "GPL6244")
    p <- source$phenotype
    p[, geo_accession := as.character(geo_accession)]
    p[, patient_id := trimws(as.character(get("study individual number:ch1")))]
    p[, regimen := as.character(get("induction therapy_maintenance therapy:ch1"))]
    p[, week := as.character(get("week (w):ch1"))]
    p[, response_marker := data.table::fcase(
      grepl("_UC NR ", title), "NR",
      grepl("_UC R ", title), "R",
      default = NA_character_
    )]
    p[, week_number := as_numeric_safely(gsub("[^0-9]", "", week))]
    p[, key := paste(regimen, patient_id, sep = "__")]
    response_map <- p[!is.na(response_marker)][order(week_number), .SD[1L], by = key]
    response_map <- response_map[, .(key, assigned_response = response_marker)]
    p <- merge(p, response_map, by = "key", all.x = TRUE, sort = FALSE)
    p[, treatment := data.table::fcase(
      regimen == "IFX", "infliximab",
      grepl("vdz", regimen) & regimen != "plac_plac", "vedolizumab",
      default = NA_character_
    )]
    p[, failure := data.table::fcase(
      assigned_response == "NR", 1L,
      assigned_response == "R", 0L,
      default = NA_integer_
    )]
    p <- p[week == "W0" & !is.na(treatment) & !is.na(failure)]
    return(list(
      make_cohort_object(source$expression, p[treatment == "infliximab"],
                         "GSE73661_IFX", "infliximab"),
      make_cohort_object(source$expression, p[treatment == "vedolizumab"],
                         "GSE73661_VDZ", "vedolizumab")
    ))
  }

  if (source_id == "GSE23597") {
    source <- load_microarray_source("GSE23597", "GPL570")
    p <- source$phenotype[get("time:ch1") == "W0"]
    p[, geo_accession := as.character(geo_accession)]
    p[, treatment := ifelse(tolower(get("dose:ch1")) == "placebo", "placebo", "infliximab")]
    p[, failure := data.table::fcase(
      get("wk8 response:ch1") == "No", 1L,
      get("wk8 response:ch1") == "Yes", 0L,
      default = NA_integer_
    )]
    p[, patient_id := as.character(get("subject:ch1"))]
    return(list(make_cohort_object(source$expression, p, "GSE23597", "infliximab")))
  }

  if (source_id == "GSE72819") {
    source <- load_gse72819_source()
    p <- source$phenotype
    p[, geo_accession := as.character(geo_accession)]
    p[, treatment := "etrolizumab"]
    p[, failure := data.table::fcase(
      get("remission at week 10:ch1") == "Non-remitter", 1L,
      get("remission at week 10:ch1") == "Remitter", 0L,
      default = NA_integer_
    )]
    p[, prior_anti_tnf := as.character(get("prior anti.tnf:ch1"))]
    return(list(make_cohort_object(source$expression, p, "GSE72819", "etrolizumab")))
  }

  stop("未定义的 source_id：", source_id)
}

# 根据注册表计算全部签名。risk_score 越高始终表示失败风险越高。
score_registered_signatures <- function(expression, registry) {
  rank_matrix <- rank_transform_columns(expression)
  signature_ids <- unique(registry$signature_id)
  scores <- data.table::data.table(geo_accession = colnames(expression))
  coverage <- vector("list", length(signature_ids))

  for (i in seq_along(signature_ids)) {
    sig_id <- signature_ids[i]
    definition <- registry[signature_id == sig_id]
    present <- definition[gene %in% rownames(rank_matrix)]
    required <- if (nrow(definition) <= 3L) nrow(definition) else ceiling(0.60 * nrow(definition))
    if (nrow(present) < required) {
      risk_score <- rep(NA_real_, ncol(rank_matrix))
    } else if (any(present$component_weight < 0)) {
      positive <- present[component_weight > 0, gene]
      negative <- present[component_weight < 0, gene]
      risk_score <- colMeans(rank_matrix[positive, , drop = FALSE], na.rm = TRUE) -
        colMeans(rank_matrix[negative, , drop = FALSE], na.rm = TRUE)
      risk_score <- risk_score * unique(present$risk_orientation)
    } else {
      risk_score <- colMeans(rank_matrix[present$gene, , drop = FALSE], na.rm = TRUE) *
        unique(present$risk_orientation)
    }
    scores[, (sig_id) := as.numeric(risk_score)]
    coverage[[i]] <- data.table::data.table(
      signature_id = sig_id,
      genes_total = nrow(definition),
      genes_present = nrow(present),
      coverage = nrow(present) / nrow(definition),
      genes_missing = paste(setdiff(definition$gene, present$gene), collapse = ";")
    )
  }
  list(scores = scores, coverage = data.table::rbindlist(coverage), ranks = rank_matrix)
}

# 为匹配随机基因集生成 5x5 的平均表达秩×变异度分层。
make_gene_matching_pools <- function(rank_matrix) {
  means <- rowMeans(rank_matrix, na.rm = TRUE)
  iqrs <- apply(rank_matrix, 1L, stats::IQR, na.rm = TRUE)
  mean_bin <- pmin(5L, pmax(1L, ceiling(rank(means, ties.method = "average") / length(means) * 5)))
  iqr_bin <- pmin(5L, pmax(1L, ceiling(rank(iqrs, ties.method = "average") / length(iqrs) * 5)))
  bins <- paste(mean_bin, iqr_bin, sep = "_")
  pools <- split(seq_along(bins), bins)
  list(bins = bins, pools = pools)
}

# 将连续得分与二分结局的相关转为 Fisher Z，作为匹配随机零分布的效应统计量。
score_matrix_to_fisher_z <- function(score_matrix, outcome) {
  y <- as.numeric(outcome)
  y_centered <- y - mean(y)
  score_matrix <- as.matrix(score_matrix)
  score_centered <- score_matrix - rowMeans(score_matrix)
  numerator <- as.numeric(score_centered %*% y_centered)
  denominator <- sqrt(rowSums(score_centered^2) * sum(y_centered^2))
  correlation <- numerator / denominator
  correlation <- pmin(0.999999, pmax(-0.999999, correlation))
  atanh(correlation) * sqrt(length(y) - 3)
}

# 每个真实签名在同一队列中生成 n_random 个大小、平均表达和 IQR 匹配的随机签名。
# 使用稀疏矩阵乘法，一次只保留 500 个随机签名，控制内存峰值。
matched_random_null <- function(rank_matrix, outcome, registry, n_random, seed_offset = 0L) {
  set.seed(SEED + as.integer(seed_offset))
  pools_object <- make_gene_matching_pools(rank_matrix)
  signature_ids <- unique(registry$signature_id)
  output <- matrix(NA_real_, nrow = n_random, ncol = length(signature_ids),
                   dimnames = list(NULL, signature_ids))
  block_size <- 500L

  for (sig_id in signature_ids) {
    definition <- registry[signature_id == sig_id]
    present <- definition[gene %in% rownames(rank_matrix)]
    required <- if (nrow(definition) <= 3L) nrow(definition) else ceiling(0.60 * nrow(definition))
    if (nrow(present) < required) next

    gene_index <- match(present$gene, rownames(rank_matrix))
    forbidden <- gene_index
    matched_pools <- lapply(gene_index, function(index) {
      pool <- pools_object$pools[[pools_object$bins[index]]]
      candidate <- setdiff(pool, forbidden)
      if (length(candidate) < 20L) candidate <- setdiff(seq_len(nrow(rank_matrix)), forbidden)
      candidate
    })
    if (any(present$component_weight < 0)) {
      positive_n <- sum(present$component_weight > 0)
      negative_n <- sum(present$component_weight < 0)
      weights <- ifelse(present$component_weight > 0, 1 / positive_n, -1 / negative_n)
    } else {
      weights <- rep(1 / nrow(present), nrow(present))
    }
    weights <- weights * unique(present$risk_orientation)

    starts <- seq.int(1L, n_random, by = block_size)
    for (start in starts) {
      end <- min(n_random, start + block_size - 1L)
      size <- end - start + 1L
      sampled <- do.call(cbind, lapply(matched_pools, sample, size = size, replace = TRUE))
      incidence <- Matrix::sparseMatrix(
        i = rep(seq_len(size), each = nrow(present)),
        j = as.vector(t(sampled)),
        x = rep(weights, times = size),
        dims = c(size, nrow(rank_matrix))
      )
      random_scores <- incidence %*% rank_matrix
      output[start:end, sig_id] <- score_matrix_to_fisher_z(random_scores, outcome)
    }
  }
  output
}

# 对一个治疗臂拟合所有签名的单变量 logistic 模型并计算 AUC。
fit_all_signature_effects <- function(scored_phenotype, cohort, treatment_name, registry) {
  signature_ids <- unique(registry$signature_id)
  results <- lapply(signature_ids, function(sig_id) {
    d <- scored_phenotype[treatment == treatment_name &
                            stats::complete.cases(failure, get(sig_id))]
    if (nrow(d) < 10L || length(unique(d$failure)) < 2L || stats::sd(d[[sig_id]]) == 0) {
      return(NULL)
    }
    d[, score_z := standardize_score(get(sig_id))]
    model <- stats::glm(failure ~ score_z, family = stats::binomial(), data = d)
    effect <- extract_glm_term(model, "score_z", cohort, treatment_name, "signature_unadjusted")
    auc <- calculate_auc(d$failure, d$score_z)
    cbind(data.table::data.table(signature_id = sig_id), effect, auc)
  })
  data.table::rbindlist(results, fill = TRUE)
}

# 安全抽取交互项；完全分离时保留结果但打上警告字段。
fit_interaction_term <- function(model, term, signature_id, cohort) {
  table <- summary(model)$coefficients
  if (!term %in% rownames(table)) return(NULL)
  estimate <- table[term, "Estimate"]
  se <- table[term, "Std. Error"]
  data.table::data.table(
    signature_id = signature_id,
    cohort = cohort,
    term = term,
    log_or = estimate,
    se = se,
    or = exp(estimate),
    ci_low = exp(estimate - 1.96 * se),
    ci_high = exp(estimate + 1.96 * se),
    p_value = table[term, "Pr(>|z|)"],
    n = stats::nobs(model),
    convergence_warning = !isTRUE(model$converged) || abs(estimate) > 10
  )
}
