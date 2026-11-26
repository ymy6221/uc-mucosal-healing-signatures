# Run after stage 6. Refit GSVA wherever the eligible set differs.
library(data.table)
source("config.R")
source("scripts/helpers.R")
source("scripts/benchmark_helpers.R")
dest <- "results/overlap_pathways"
dir.create(dest, showWarnings = FALSE, recursive = TRUE)
sets <- fread("results/signature_convergence_loadings.csv")$signature_id
path_scores <- fread("results/cbm_upgrade/pathway_pc1_sample_scores.csv")
labels <- unique(fread("results/cbm_upgrade/pathway_overlap_audit.csv")[, .(pathway, display_name)])
common <- c("GSE206285", "GSE92415", "GSE73661_VDZ", "GSE72819")
# Refit on both active and placebo participants, as in the original PCA code.
cohorts <- c(common, "GSE23597", "GSE16879")
scores <- rbindlist(lapply(cohorts, function(id) {
  x <- as.data.table(readRDS(sprintf("data/intermediate/benchmark_%s_scored.rds", id)))
  if (id == "GSE23597") x <- x[order(geo_accession)][!duplicated(patient_id)]
  x[, cohort := id]
  x[, c("cohort","geo_accession", sets), with = FALSE]
}))
scores[, (sets) := lapply(.SD, function(x) as.numeric(scale(x))), by = cohort, .SDcols = sets]
X <- as.matrix(scores[, ..sets])
stopifnot(all(is.finite(X)), nrow(X)==808L, !anyDuplicated(scores$geo_accession))
fit <- prcomp(X, center = TRUE, scale. = TRUE)
orient <- if (sum(fit$rotation[, 1]) < 0) -1 else 1
scores[, pc1 := orient * fit$x[, 1]]
vars <- fit$sdev^2 / sum(fit$sdev^2)
fwrite(data.table(component = seq_along(vars), proportion = vars), file.path(dest, "variance.csv"))
fwrite(data.table(signature_id = sets, loading = orient * fit$rotation[, 1]), file.path(dest, "loadings.csv"))
path_scores[, pc1 := NULL]
act <- prepare_source_cohorts("GSE23597")[[1]]
gsm <- scores[cohort == "GSE23597", geo_accession]
expr <- act$expression[, gsm, drop = FALSE]
reg <- fread("results/cbm_upgrade/pathway_registry_msigdb.csv")
genes <- split(reg$gene, reg$pathway)
exclude <- read_signature_registry()[primary_set == TRUE, unique(gene)]
variants <- c(setNames(genes, paste0("full::", names(genes))),
              setNames(lapply(genes, setdiff, y = exclude), paste0("overlap_excluded::", names(genes))))
params <- GSVA::gsvaParam(expr, lapply(variants, intersect, y = rownames(expr)),
                          kcdf = "Gaussian", minSize = 10L, maxSize = 500L, sparse = FALSE, verbose = FALSE)
new_gsva <- GSVA::gsva(params, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
act_rows <- path_scores[cohort == "GSE23597" & geo_accession %in% gsm, .(cohort, geo_accession, treatment)]
new_rows <- merge(act_rows, as.data.table(t(new_gsva), keep.rownames = "geo_accession"), by = "geo_accession")
stopifnot(nrow(new_rows) == 44L, nrow(new_gsva) == 24L)
path_scores <- rbindlist(list(path_scores[cohort != "GSE23597"], new_rows), use.names = TRUE)
j <- merge(scores[, .(cohort, geo_accession, pc1)], path_scores, by = c("cohort","geo_accession"))
stopifnot(nrow(j) == nrow(X), j[treatment == "placebo", .N] == 225)
long <- melt(j, id.vars = c("cohort","geo_accession","treatment","pc1"),
             variable.name = "key", value.name = "score")
long[, c("variant", "pathway") := tstrsplit(key, "::", fixed = TRUE)]
est <- long[, .(rho = cor(pc1, score, method = "spearman"), n = .N), by = .(cohort, variant, pathway)]
tab <- est[, {
  m <- random_effects_meta(atanh(rho), 1/sqrt(n-3))
  .(k = .N, n = sum(n), rho = tanh(m$log_or), low = tanh(m$log_or - 1.96*m$se),
    high = tanh(m$log_or + 1.96*m$se), p = m$p_value)
}, by = .(variant, pathway)]
tab[, fdr := p.adjust(p, "BH"), by = variant]
fwrite(est, file.path(dest, "cohort_correlations.csv"))
fwrite(merge(tab, labels, by = "pathway"), file.path(dest, "pathways.csv"))
fwrite(j[, .(cohort, geo_accession, pc1)], file.path(dest, "projected_scores.csv"))

figure <- merge(tab, labels, by = "pathway")
order_names <- figure[variant == "full"][order(rho), display_name]
draw <- function() {
  par(mar = c(4, 12.5, 2, 1))
  plot(NA, xlim = c(-0.55, 1), ylim = c(0.3, 13), axes = FALSE,
       xlab = "Pooled Spearman correlation with PC1", ylab = "")
  axis(1, at = seq(-0.5, 1, 0.25))
  axis(2, at = seq_along(order_names), labels = order_names, las = 1, tick = FALSE, cex.axis = 0.80)
  abline(v = 0, col = "grey65", lty = 2)
  for (i in 1:2) {
    subset <- figure[variant == c("full", "overlap_excluded")[i]]
    y <- match(subset$display_name, order_names) + c(0.13,-0.13)[i]
    colour <- c("#18557A", "#B25827")[i]
    segments(subset$low, y, subset$high, y, col = colour, lwd = 1.9)
    points(subset$rho, y, pch = c(16,17)[i], col = colour)
  }
  legend("topleft", c("Full pathway", "Signature genes excluded"), pch = c(16,17),
         col = c("#18557A","#B25827"), bty = "n", cex = 0.8)
}
ragg::agg_tiff(file.path(dest, "pathway_correlations.tiff"), width = 8.5, height = 6.8, units = "in",
               res = 320, compression = "lzw")
draw()
invisible(dev.off())
