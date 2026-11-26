# Run from the project root after Stage 5. Only the audit directory is written.
source("config.R")
source("scripts/helpers.R")
source("scripts/benchmark_helpers.R")
suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)
out <- file.path(DIR_RESULTS, "overlap_audit")
dir.create(out, showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(out, name), na = "NA")
read_arm <- function(k) as.data.table(readRDS(file.path(
  DIR_INTERMEDIATE, paste0("benchmark_", k, "_scored.rds"))))

# Stage 5 defines the eligible GSMs. GEO files supply expression and the patient keys.
arms <- c("GSE12251", "GSE23597", "GSE16879", "GSE73661_IFX", "GSE73661_VDZ",
          "GSE206285", "GSE92415", "GSE72819")
ph <- setNames(lapply(arms, read_arm), arms)
stopifnot(all(vapply(ph, function(p) !anyDuplicated(p$geo_accession), logical(1))))
compare_gsm <- rbindlist(lapply(arms, function(k) {
  p <- ph[[k]]
  data.table(arm = k, gsm = p$geo_accession, treatment = p$treatment)
}))
save_table(compare_gsm, "baseline_gsm_registry.csv")
stopifnot(!anyDuplicated(compare_gsm$gsm))
within_ids <- rbindlist(lapply(arms, function(k) {
  p <- ph[[k]]
  field <- intersect(c("donor id:ch1", "subject:ch1", "patient_id", "title"), names(p))[1]
  p[, .(arm = k, key = field, subject = get(field), treatment, geo_accession)]
}))
stopifnot(!anyNA(within_ids$subject), all(nzchar(within_ids$subject)))
save_table(within_ids[, .(baseline_arrays = .N, unique_key = uniqueN(subject)),
                      by = .(arm, treatment, key)], "within_arm_subject_counts.csv")
act1 <- ph$GSE12251
act2 <- ph$GSE23597[treatment == "infliximab"]
ids_1 <- sub("/.*", "", act1$title)
ids_2 <- sub("/.*", "", act2$title)
key <- function(title) sub(" .*", "", title)
act_ids <- data.table(gsm_a = act1$geo_accession, subject = ids_1,
                      dose_visit_a = key(act1$title))
act_ids <- merge(act_ids, data.table(gsm_b = act2$geo_accession, subject = ids_2,
                                    dose_visit_b = key(act2$title)), by = "subject", all.x = TRUE)
stopifnot(!anyDuplicated(ids_1), identical(unique(ids_2[duplicated(ids_2)]), "P13"))
discovery_samples <- paste0("GSM", 364633:364656)
# Source: GSE14580 Samples list, 24 UC baseline arrays excluding its 6 controls.
discovery <- data.table(gsm = discovery_samples,
                        present_in_GSE16879_baseline = discovery_samples %in% ph$GSE16879$geo_accession)
save_table(discovery, "discovery_gse14580_gsm_check.csv")
save_table(data.table(
  comparison = "GSE73661_IFX vs GSE73661_VDZ at baseline",
  n_ifx = nrow(ph$GSE73661_IFX), n_vdz = nrow(ph$GSE73661_VDZ),
  shared_study_individual_numbers = length(intersect(ph$GSE73661_IFX$patient_id,
                                                     ph$GSE73661_VDZ$patient_id))),
  "within_series_id_check.csv")

# Same-platform audit. Include P13's second pretreatment array as a documented sensitivity.
probe_a <- read_geo_eset("GSE12251")
probe_b <- read_geo_eset("GSE23597")
meta_a <- Biobase::pData(probe_a)
a <- Biobase::exprs(probe_a)
b <- Biobase::exprs(probe_b)[, act2$geo_accession, drop = FALSE]
features <- intersect(rownames(a), rownames(b))
a <- a[features, , drop = FALSE]
b <- b[features, , drop = FALSE]
stopifnot(all(is.finite(a)), all(is.finite(b)), min(a) >= 0, min(b) >= 0)
pc <- cor(log2(1 + a), log2(1 + b))
sc <- cor(a, b, method = "spearman")
probe_pairs <- as.data.table(as.table(pc))
setnames(probe_pairs, c("gsm_a", "gsm_b", "pearson_log2p1"))
probe_pairs[, spearman := as.vector(sc)]
probe_pairs[, subject_a := sub("/.*", "", meta_a[as.character(gsm_a), "title"])]
probe_pairs[, subject_b := ids_2[match(gsm_b, act2$geo_accession)]]
probe_pairs[, same_subject := subject_a == subject_b]
probe_pairs[, in_stage5_a := gsm_a %in% act1$geo_accession]
probe_pairs[, exactly_equal_probes := vapply(seq_len(.N), function(i)
  sum(a[, as.character(gsm_a[i])] == b[, as.character(gsm_b[i])]), integer(1))]
probe_pairs[, shared_probes := length(features)]
save_table(probe_pairs, "act1_all_probe_pairs.csv")
act_ids <- merge(act_ids, probe_pairs[same_subject & in_stage5_a,
    .(gsm_a, gsm_b, pearson_log2p1, spearman, exactly_equal_probes, shared_probes)],
    by = c("gsm_a", "gsm_b"), all.x = TRUE)
save_table(act_ids, "act1_baseline_crosswalk.csv")
rm(a, b, probe_a, probe_b)
invisible(gc())

# Across platforms, use shared gene symbols and log-scale expression. Within each study
# choose a unique, unambiguously mapped probe by baseline IQR without using outcomes.
# Center and scale each gene across its own arm before comparing sample profiles.
platforms <- c(GSE12251 = "GPL570", GSE23597 = "GPL570", GSE16879 = "GPL570",
               GSE73661 = "GPL6244", GSE206285 = "GPL13158", GSE92415 = "GPL13158")
read_genes <- function(id, gsms) {
  if (id == "GSE72819") return(load_gse72819_source()$expression[, gsms, drop = FALSE])
  eset <- read_geo_eset(id)
  e <- Biobase::exprs(eset)[, gsms, drop = FALSE]
  if (platforms[[id]] == "GPL570") e <- log2(1 + e)
  ann <- fread(file.path(DIR_RAW, paste0(platforms[[id]], ".annot.gz")),
               skip = 27L, select = c("ID", "Gene symbol"), showProgress = FALSE)
  setnames(ann, c("probe", "gene"))
  ann[, probe := as.character(probe)]
  ann <- ann[probe %in% rownames(e) & !is.na(gene) & !gene %in% c("", "---") &
             !grepl("///", gene, fixed = TRUE)]
  e <- e[ann$probe, , drop = FALSE]
  ann[, iqr := matrixStats::rowIQRs(e)]
  setorder(ann, gene, -iqr, probe)
  winners <- ann[!duplicated(gene)]
  e <- e[winners$probe, , drop = FALSE]
  rownames(e) <- winners$gene
  e
}
expression <- list()
for (id in unique(sub("_.*", "", arms))) {
  keys <- arms[sub("_.*", "", arms) == id]
  gsms <- unlist(lapply(ph[keys], "[[", "geo_accession"), use.names = FALSE)
  e <- read_genes(id, gsms)
  for (k in keys) expression[[k]] <- e[, ph[[k]]$geo_accession, drop = FALSE]
  rm(e)
  invisible(gc())
}

mutual_max <- function(x) {
  by_row <- sweep(x, 1L, apply(x, 1L, max), "==")
  by_col <- sweep(x, 2L, apply(x, 2L, max), "==")
  by_row & by_col
}
scan <- function(arm_a, arm_b, detail = FALSE) {
  x <- expression[[arm_a]]
  y <- expression[[arm_b]]
  if (arm_a == "GSE12251" && arm_b == "GSE23597")
    y <- y[, act2$geo_accession, drop = FALSE]
  genes <- sort(intersect(rownames(x), rownames(y)), method = "radix")
  x <- x[genes, , drop = FALSE]
  y <- y[genes, , drop = FALSE]
  sx <- matrixStats::rowSds(x)
  sy <- matrixStats::rowSds(y)
  keep <- is.finite(sx) & is.finite(sy) & sx > 0 & sy > 0
  x <- x[keep, , drop = FALSE]
  y <- y[keep, , drop = FALSE]
  zx <- (x - rowMeans(x)) / sx[keep]
  zy <- (y - rowMeans(y)) / sy[keep]
  corr <- cor(zx, zy)
  # Requiring both fixed, non-overlapping gene halves to agree is descriptive, not identity proof.
  odd <- seq_len(nrow(zx)) %% 2L == 1L
  c1 <- cor(zx[odd, , drop = FALSE], zy[odd, , drop = FALSE])
  c2 <- cor(zx[!odd, , drop = FALSE], zy[!odd, , drop = FALSE])
  stable <- mutual_max(corr) & mutual_max(c1) & mutual_max(c2)
  dat <- as.data.table(as.table(corr))
  setnames(dat, c("gsm_a", "gsm_b", "standardised_r"))
  dat[, `:=`(arm_a = arm_a, arm_b = arm_b, shared_genes = nrow(zx),
             stable_reciprocal_best = as.vector(stable))]
  if (detail) {
    dat[, `:=`(pearson_log_expression = as.vector(cor(x, y)),
               spearman_expression = as.vector(cor(x, y, method = "spearman")),
               standardised_r_half1 = as.vector(c1), standardised_r_half2 = as.vector(c2))]
    save_table(dat, paste0("pair_matrix_", arm_a, "__", arm_b, ".csv"))
  }
  list(summary = dat[, .(arm_a = arm_a[1], arm_b = arm_b[1], shared_genes = shared_genes[1],
    pairs = .N, median_r = median(standardised_r), p99_r = quantile(standardised_r, .99),
    max_r = max(standardised_r), stable_reciprocal_best = sum(stable_reciprocal_best),
    stable_r_over_0.5 = sum(stable_reciprocal_best & standardised_r > .5))],
    hits = dat[stable_reciprocal_best == TRUE, .(arm_a, arm_b, gsm_a, gsm_b, shared_genes, standardised_r)])
}
pairs <- combn(arms, 2L)
results <- lapply(seq_len(ncol(pairs)), function(i) {
  pair <- pairs[, i]
  detail <- all(pair %in% c("GSE12251", "GSE23597")) ||
            all(pair %in% c("GSE16879", "GSE73661_IFX", "GSE73661_VDZ"))
  scan(pair[1], pair[2], detail)
})
save_table(rbindlist(lapply(results, "[[", "summary")), "all_arm_pair_summary.csv")
save_table(rbindlist(lapply(results, "[[", "hits")), "all_stable_reciprocal_pairs.csv")

# Independently repeat the Leuven screen at weeks 4–6. Relate visits only with
# keys supplied within each GEO series. Disease activity may influence both visits.
meta16 <- as.data.table(Biobase::pData(read_geo_eset("GSE16879")))
meta73 <- as.data.table(Biobase::pData(read_geo_eset("GSE73661")))
follow16 <- meta16[grepl("^UC(N?R)[0-9]+_afterT$", title),
                   .(gsm = geo_accession, subject = sub("_afterT$", "", title))]
follow73 <- meta73[get("induction therapy_maintenance therapy:ch1") == "IFX" &
                   get("week (w):ch1") == "W4_W6",
                   .(gsm = geo_accession, subject = as.character(get("study individual number:ch1")))]
stopifnot(!anyDuplicated(follow16$subject), !anyDuplicated(follow73$subject),
          setequal(follow16$subject, sub("_beforeT$", "", ph$GSE16879$title)),
          setequal(follow73$subject, ph$GSE73661_IFX$patient_id))
expression <- list(GSE16879_after = read_genes("GSE16879", follow16$gsm),
                   GSE73661_IFX_after = read_genes("GSE73661", follow73$gsm))
follow_pairs <- scan("GSE16879_after", "GSE73661_IFX_after")$hits
follow_pairs[, subject_a := follow16$subject[match(gsm_a, follow16$gsm)]]
follow_pairs[, subject_b := follow73$subject[match(gsm_b, follow73$gsm)]]
base_pairs <- fread(file.path(out, "pair_matrix_GSE16879__GSE73661_IFX.csv"))
base_pairs[, subject_a := sub("_beforeT$", "", ph$GSE16879$title[match(gsm_a, ph$GSE16879$geo_accession)])]
base_pairs[, subject_b := ph$GSE73661_IFX$patient_id[match(gsm_b, ph$GSE73661_IFX$geo_accession)]]
base_pairs <- base_pairs[stable_reciprocal_best == TRUE]
cross_visit <- merge(base_pairs[, .(subject_a, subject_b, baseline_gsm_a = gsm_a,
                             baseline_gsm_b = gsm_b, baseline_r = standardised_r)],
                     follow_pairs[, .(subject_a, subject_b, followup_gsm_a = gsm_a,
                             followup_gsm_b = gsm_b, followup_r = standardised_r)],
                     by = c("subject_a", "subject_b"), all = TRUE)
save_table(cross_visit, "leuven_cross_visit_candidates.csv")
single_act <- act_ids[exactly_equal_probes == shared_probes]
matches_a <- match(single_act$gsm_a, act1$geo_accession)
matches_b <- match(single_act$gsm_b, act2$geo_accession)
summary <- data.table(
  metric = c("ACT1 distinct subjects matched across the two series",
    "ACT1 included subjects with identical profiles at all shared probes",
    "ACT1 probe sets compared", "ACT1 matched subjects with discordant endpoint labels",
    "GSE23597 extra baseline records for subjects represented twice", "GSE14580 discovery GSMs reused by GSE16879",
    "Leuven cross-platform IFX baseline comparisons", "Leuven shared baseline genes",
    "Leuven stable reciprocal pairs at baseline", "Leuven stable reciprocal pairs at weeks 4-6",
    "Leuven stable reciprocal pairs at both visits"),
  value = c(uniqueN(act_ids$subject), nrow(single_act), unique(single_act$shared_probes),
    sum(act1$failure[matches_a] != act2$failure[matches_b]), sum(duplicated(ids_2)),
    sum(discovery$present_in_GSE16879_baseline), nrow(ph$GSE16879) * nrow(ph$GSE73661_IFX),
    unique(base_pairs$shared_genes), nrow(base_pairs), nrow(follow_pairs),
    sum(complete.cases(cross_visit))))
save_table(summary, "audit_summary.csv")
print(summary)
