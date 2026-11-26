# 本文件集中放可复用函数。每个函数只做一件事，便于排错。

# 检查核心 R 包；缺包时给出明确提示而不是运行到一半才报错。
check_packages <- function() {
  needed <- c("GEOquery", "Biobase", "data.table", "ggplot2", "pROC")
  ok <- vapply(needed, requireNamespace, logical(1), quietly = TRUE)
  if (!all(ok)) {
    stop("缺少 R 包：", paste(needed[!ok], collapse = ", "),
         "。请先运行 scripts/00_install_packages.R")
  }
  invisible(TRUE)
}

# 带重试的下载函数。已有且大小正常的文件会跳过，不重复下载。
download_with_retry <- function(url, destination, min_bytes = 1024L, attempts = 3L) {
  if (file.exists(destination) && file.info(destination)$size >= min_bytes) {
    message("已存在，跳过下载：", basename(destination))
    return(invisible(destination))
  }
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(attempts)) {
    message("下载第 ", i, "/", attempts, " 次：", basename(destination))
    status <- try(
      utils::download.file(url, destination, mode = "wb", method = "libcurl", quiet = FALSE),
      silent = TRUE
    )
    valid <- !inherits(status, "try-error") && file.exists(destination) &&
      file.info(destination)$size >= min_bytes
    if (valid) {
      return(invisible(destination))
    }
    if (file.exists(destination)) {
      unlink(destination)
    }
  }
  stop("下载失败：", url, "\n请检查网络；也可用浏览器下载后放入：", destination)
}

# 读取 GEO Series Matrix，但禁止 GEOquery 再联网下载平台文件。
read_geo_eset <- function(gse_id) {
  path <- file.path(DIR_RAW, paste0(gse_id, "_series_matrix.txt.gz"))
  if (!file.exists(path)) {
    stop("找不到 ", path, "；请先运行 scripts/01_download_core.R")
  }
  GEOquery::getGEO(filename = path, getGPL = FALSE)
}

# 读取 GEO 官方平台注释，只保留探针 ID 和基因符号两列。
read_geo_annotation <- function(gpl_id) {
  path <- file.path(DIR_RAW, paste0(gpl_id, ".annot.gz"))
  if (!file.exists(path)) {
    stop("找不到 ", path, "；请先运行 scripts/01_download_core.R")
  }
  ann <- data.table::fread(
    path,
    skip = 27L,
    select = c("ID", "Gene symbol"),
    showProgress = FALSE
  )
  data.table::setnames(ann, c("probe", "symbol"))
  ann[, probe := as.character(probe)]
  ann[, symbol := sub(" ///.*$", "", symbol)]
  ann <- ann[!is.na(symbol) & symbol != "" & symbol != "---"]
  ann
}

# 把探针矩阵转成基因矩阵。多探针对应同一基因时，选 IQR 最大的探针。
# 该选择完全不使用临床结局，因此不会发生标签泄漏。
probe_to_gene_matrix <- function(eset, gpl_id) {
  probe_expression <- Biobase::exprs(eset)
  annotation <- read_geo_annotation(gpl_id)
  annotation <- annotation[probe %in% rownames(probe_expression)]
  if (nrow(annotation) == 0L) {
    stop(gpl_id, " 的注释与表达矩阵没有匹配探针。")
  }
  matched_expression <- probe_expression[annotation$probe, , drop = FALSE]
  probe_iqr <- apply(matched_expression, 1L, stats::IQR, na.rm = TRUE)
  annotation[, probe_iqr := probe_iqr]
  data.table::setorder(annotation, symbol, -probe_iqr)
  annotation <- annotation[!duplicated(symbol)]
  gene_expression <- probe_expression[match(annotation$probe, rownames(probe_expression)), , drop = FALSE]
  rownames(gene_expression) <- annotation$symbol
  gene_expression
}

# 计算每个样本内部的百分位秩；这样不同芯片平台无需批次合并。
rank_transform_columns <- function(gene_expression) {
  ranks <- apply(
    gene_expression,
    2L,
    rank,
    ties.method = "average",
    na.last = "keep"
  )
  if (is.null(dim(ranks))) {
    ranks <- matrix(ranks, ncol = 1L)
  }
  ranks / (nrow(gene_expression) + 1)
}

# 从固定基因表计算 ERI：损伤模块秩均值减去成熟模块秩均值。
calculate_eri <- function(gene_expression, signature_table) {
  rank_matrix <- rank_transform_columns(gene_expression)
  injury_genes <- intersect(signature_table[module == "injury", gene], rownames(rank_matrix))
  mature_genes <- intersect(signature_table[module == "mature", gene], rownames(rank_matrix))
  if (length(injury_genes) < 10L || length(mature_genes) < 10L) {
    stop("签名基因覆盖不足：injury=", length(injury_genes), ", mature=", length(mature_genes))
  }
  injury_score <- colMeans(rank_matrix[injury_genes, , drop = FALSE], na.rm = TRUE)
  mature_score <- colMeans(rank_matrix[mature_genes, , drop = FALSE], na.rm = TRUE)
  data.table::data.table(
    geo_accession = colnames(gene_expression),
    injury_score = as.numeric(injury_score),
    mature_score = as.numeric(mature_score),
    eri = as.numeric(injury_score - mature_score),
    injury_genes_n = length(injury_genes),
    mature_genes_n = length(mature_genes)
  )
}

# 安全地把字符转为数字；无法转换的值变成 NA。
as_numeric_safely <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

# 在一个队列内部把 ERI 标准化为每 1 个标准差，便于合并 OR。
standardize_score <- function(x) {
  as.numeric(scale(as.numeric(x)))
}

# 从 glm 中抽取指定项的 OR、标准误、置信区间和 P 值。
extract_glm_term <- function(model, term, cohort, treatment, model_name) {
  coefficient_table <- summary(model)$coefficients
  if (!term %in% rownames(coefficient_table)) {
    stop("模型中不存在项：", term)
  }
  estimate <- coefficient_table[term, "Estimate"]
  standard_error <- coefficient_table[term, "Std. Error"]
  data.table::data.table(
    cohort = cohort,
    treatment = treatment,
    model = model_name,
    term = term,
    log_or = estimate,
    se = standard_error,
    or = exp(estimate),
    ci_low = exp(estimate - 1.96 * standard_error),
    ci_high = exp(estimate + 1.96 * standard_error),
    p_value = coefficient_table[term, "Pr(>|z|)"],
    n = stats::nobs(model),
    events = sum(stats::model.response(stats::model.frame(model)) == 1L)
  )
}

# 计算 AUC 及 DeLong 置信区间。结局 1 表示失败，ERI 越高风险越高。
calculate_auc <- function(outcome, score) {
  keep <- stats::complete.cases(outcome, score)
  roc_object <- pROC::roc(outcome[keep], score[keep], direction = "<", quiet = TRUE)
  interval <- as.numeric(pROC::ci.auc(roc_object, method = "delong"))
  data.table::data.table(auc = as.numeric(pROC::auc(roc_object)),
                         auc_low = interval[1L], auc_high = interval[3L])
}

# 单队列基础模型：failure ~ ERI_z，并同时返回效应量和 AUC。
fit_simple_cohort <- function(data, cohort, treatment) {
  analysis_data <- data[stats::complete.cases(failure, eri)]
  analysis_data[, eri_z := standardize_score(eri)]
  model <- stats::glm(failure ~ eri_z, family = stats::binomial(), data = analysis_data)
  result <- extract_glm_term(model, "eri_z", cohort, treatment, "unadjusted")
  cbind(result, calculate_auc(analysis_data$failure, analysis_data$eri_z))
}

# DerSimonian-Laird 随机效应 Meta；不额外安装 metafor，降低环境负担。
random_effects_meta <- function(log_or, se) {
  variance <- se^2
  fixed_weight <- 1 / variance
  fixed_mean <- sum(fixed_weight * log_or) / sum(fixed_weight)
  q <- sum(fixed_weight * (log_or - fixed_mean)^2)
  df <- length(log_or) - 1L
  c_value <- sum(fixed_weight) - sum(fixed_weight^2) / sum(fixed_weight)
  tau2 <- max(0, (q - df) / c_value)
  random_weight <- 1 / (variance + tau2)
  pooled <- sum(random_weight * log_or) / sum(random_weight)
  pooled_se <- sqrt(1 / sum(random_weight))
  i2 <- if (q > 0) max(0, (q - df) / q) * 100 else 0
  data.table::data.table(
    log_or = pooled,
    se = pooled_se,
    or = exp(pooled),
    ci_low = exp(pooled - 1.96 * pooled_se),
    ci_high = exp(pooled + 1.96 * pooled_se),
    p_value = 2 * stats::pnorm(-abs(pooled / pooled_se)),
    tau2 = tau2,
    q = q,
    df = df,
    i2 = i2
  )
}

