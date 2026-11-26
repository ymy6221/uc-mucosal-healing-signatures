# CBM 升级图：通路机制、精细单细胞定位、异质性稳健性和图文摘要。

figure_packages <- c("data.table", "ggplot2", "patchwork", "scales", "ggrepel")
figure_ok <- vapply(figure_packages, requireNamespace, logical(1), quietly = TRUE)
if (!all(figure_ok)) {
  stop("缺少作图包：", paste(figure_packages[!figure_ok], collapse = ", "))
}

DIR_CBM_RESULTS <- file.path(DIR_RESULTS, "cbm_upgrade")
DIR_CBM_FIGURES <- file.path(DIR_FIGURES, "cbm_upgrade")
dir.create(DIR_CBM_FIGURES, recursive = TRUE, showWarnings = FALSE)

robust_ids <- c("IFX_ARJS_5", "OSM_OSMR", "IL22_NEUTROPHIL_6", "VDZ_NR_3")
signature_labels <- c(
  IFX_ARJS_5 = "Arijs 5-gene",
  OSM_OSMR = "OSM-OSMR",
  IL22_NEUTROPHIL_6 = "IL-22 chemokines",
  VDZ_NR_3 = "VDZ stromal NR"
)

theme_cbm <- ggplot2::theme_minimal(base_size = 9, base_family = "sans") +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", size = 10),
    plot.subtitle = ggplot2::element_text(size = 8, colour = "#4D4D4D"),
    axis.title = ggplot2::element_text(face = "bold"),
    legend.title = ggplot2::element_text(face = "bold"),
    plot.margin = ggplot2::margin(6, 7, 6, 7)
  )

save_publication_figure <- function(plot, stem, width, height) {
  ggplot2::ggsave(file.path(DIR_CBM_FIGURES, paste0(stem, ".pdf")), plot,
                  width = width, height = height, units = "in", device = "pdf")
  ggplot2::ggsave(file.path(DIR_CBM_FIGURES, paste0(stem, ".tiff")), plot,
                  width = width, height = height, units = "in", dpi = 600,
                  compression = "lzw")
}

# ---------- Figure 5: pathway mechanism ----------
pathway_meta <- data.table::fread(file.path(DIR_CBM_RESULTS,
                                            "pathway_meta_correlations.csv"))
pathway_cohort <- data.table::fread(file.path(DIR_CBM_RESULTS,
                                              "pathway_cohort_correlations.csv"))
pathway_samples <- data.table::fread(file.path(DIR_CBM_RESULTS,
                                               "pathway_pc1_sample_scores.csv"))

pathway_order <- pathway_meta[variant == "full"][order(rho), display_name]
pathway_meta[, display_factor := factor(display_name, levels = pathway_order)]
pathway_meta[, variant_label := data.table::fcase(
  variant == "full", "Full pathway",
  variant == "overlap_excluded", "Signature genes removed"
)]
panel_5a <- ggplot2::ggplot(pathway_meta,
                            ggplot2::aes(rho, display_factor,
                                         colour = variant_label)) +
  ggplot2::geom_vline(xintercept = 0, colour = "#888888", linewidth = 0.35) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_low, xmax = ci_high),
                          height = 0.18, linewidth = 0.55,
                          position = ggplot2::position_dodge(width = 0.45)) +
  ggplot2::geom_point(size = 2.1, position = ggplot2::position_dodge(width = 0.45)) +
  ggplot2::scale_colour_manual(values = c("Full pathway" = "#2166AC",
                                          "Signature genes removed" = "#B2182B")) +
  ggplot2::labs(title = "A  Cross-cohort pathway correlation with outcome-free PC1",
                subtitle = "Random-effects Fisher-z meta-analysis; 8 cohorts, n=854",
                x = "Meta-analytic Spearman rho (95% CI)", y = NULL,
                colour = NULL) +
  ggplot2::coord_cartesian(xlim = c(-0.45, 0.95)) + theme_cbm +
  ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_text(size = 7.2))

full_meta <- pathway_meta[variant == "full",
                          .(pathway, display_name, rho_full = rho)]
excluded_meta <- pathway_meta[variant == "overlap_excluded",
                              .(pathway, rho_excluded = rho)]
comparison <- merge(full_meta, excluded_meta, by = "pathway")
comparison[, short_label := data.table::fcase(
  grepl("Inflammatory", display_name), "Inflammation",
  grepl("Cytokine", display_name), "Cytokine-R",
  grepl("mesenchymal", display_name), "EMT",
  grepl("TNF", display_name), "TNF/NF-kB",
  grepl("IL-6", display_name), "IL6/STAT3",
  grepl("ECM", display_name), "ECM-R",
  grepl("Focal", display_name), "Focal adhesion",
  grepl("Tight", display_name), "Tight junction",
  default = ""
)]
panel_5b <- ggplot2::ggplot(comparison,
                            ggplot2::aes(rho_full, rho_excluded)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#777777") +
  ggplot2::geom_point(size = 2.4, colour = "#542788") +
  ggplot2::geom_text(ggplot2::aes(label = short_label), nudge_y = 0.035,
                     size = 2.5, check_overlap = TRUE) +
  ggplot2::labs(title = "B  Gene-overlap exclusion audit",
                subtitle = "Persistence near the diagonal argues against circular overlap",
                x = "Full-pathway rho", y = "Rho after removing all primary-signature genes") +
  ggplot2::coord_equal(xlim = c(-0.4, 0.95), ylim = c(-0.4, 0.95)) + theme_cbm

heatmap_data <- pathway_cohort[variant == "overlap_excluded"]
heatmap_data[, display_factor := factor(display_name, levels = rev(pathway_order))]
heatmap_data[, cohort := factor(cohort, levels = c(
  "GSE206285", "GSE92415", "GSE23597", "GSE16879", "GSE12251",
  "GSE73661_IFX", "GSE73661_VDZ", "GSE72819"
))]
panel_5c <- ggplot2::ggplot(heatmap_data,
                            ggplot2::aes(cohort, display_factor, fill = rho)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                midpoint = 0, limits = c(-1, 1)) +
  ggplot2::labs(title = "C  Replication after overlap exclusion",
                x = NULL, y = NULL, fill = "Spearman\nrho") + theme_cbm +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
                 axis.text.y = ggplot2::element_text(size = 7))

top_column <- "overlap_excluded::HALLMARK_INFLAMMATORY_RESPONSE"
scatter_data <- pathway_samples[cohort == "GSE206285",
                                .(pc1, pathway_score = get(top_column))]
scatter_rho <- suppressWarnings(stats::cor(scatter_data$pc1,
                                           scatter_data$pathway_score,
                                           method = "spearman"))
panel_5d <- ggplot2::ggplot(scatter_data, ggplot2::aes(pc1, pathway_score)) +
  ggplot2::geom_point(alpha = 0.48, size = 1.25, colour = "#2166AC") +
  ggplot2::geom_smooth(method = "lm", se = TRUE, linewidth = 0.7,
                       colour = "#B2182B", fill = "#F4A582") +
  ggplot2::annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
                    label = sprintf("Spearman rho = %.3f", scatter_rho), size = 2.7) +
  ggplot2::labs(title = "D  Largest-cohort illustration",
                subtitle = "GSE206285; signature genes removed from pathway",
                x = "Outcome-free signature PC1", y = "Inflammatory-response GSVA score") +
  theme_cbm

figure_5 <- (panel_5a | panel_5b) / (panel_5c | panel_5d) +
  patchwork::plot_annotation(tag_levels = NULL)
save_publication_figure(figure_5, "CBM_Figure5_pathway_mechanism", 13.2, 9.0)

# ---------- Figure 6: cell-state localization ----------
umap <- data.table::fread(file.path(DIR_CBM_RESULTS, "singlecell_fine_umap.csv"))
umap_columns <- grep("^umap", names(umap), value = TRUE, ignore.case = TRUE)
if (length(umap_columns) < 2L) stop("UMAP 坐标列缺失。")
data.table::setnames(umap, umap_columns[1:2], c("UMAP_1", "UMAP_2"))
umap[, short_cell_type := sub("^[0-9]+[A-Z]?-", "", fine_annotations_MK_V4)]
umap[, short_cell_type := data.table::fcase(
  coarse_annotations_MK == "10-Epithelial",
  paste0("E: ", sub("^Epithelial[ -]?", "", short_cell_type)),
  coarse_annotations_MK == "11-Stromal-fibroblast",
  paste0("S: ", sub("^Stromal[ -]?", "", short_cell_type)),
  default = short_cell_type
)]
centres <- umap[, .(UMAP_1 = stats::median(UMAP_1),
                    UMAP_2 = stats::median(UMAP_2)), by = short_cell_type]

panel_6a <- ggplot2::ggplot(umap, ggplot2::aes(UMAP_1, UMAP_2,
                                               colour = short_cell_type)) +
  ggplot2::geom_point(size = 0.25, alpha = 0.55) +
  ggrepel::geom_label_repel(data = centres,
                            ggplot2::aes(label = short_cell_type),
                            size = 1.65, colour = "black",
                            fill = scales::alpha("white", 0.82),
                            linewidth = 0.12, box.padding = 0.25,
                            point.padding = 0.08, min.segment.length = 0,
                            max.overlaps = Inf, seed = SEED,
                            show.legend = FALSE) +
  ggplot2::labs(title = "A  Author-annotated epithelial and stromal cell states",
                subtitle = paste0(format(nrow(umap), big.mark = ","),
                                  " stratified cells; UMAP is descriptive"),
                x = "UMAP 1", y = "UMAP 2") + theme_cbm +
  ggplot2::theme(legend.position = "none")

make_feature_panel <- function(signature_id, panel_letter) {
  one <- data.table::copy(umap)
  limits <- stats::quantile(one[[signature_id]], c(0.01, 0.99), na.rm = TRUE)
  data.table::setorderv(one, signature_id)
  ggplot2::ggplot(one, ggplot2::aes(UMAP_1, UMAP_2, colour = get(signature_id))) +
    ggplot2::geom_point(size = 0.30, alpha = 0.72) +
    ggplot2::scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7",
                                    high = "#B2182B", midpoint = 0,
                                    limits = limits, oob = scales::squish) +
    ggplot2::labs(title = paste0(panel_letter, "  ", signature_labels[[signature_id]]),
                  x = "UMAP 1", y = "UMAP 2", colour = "Module\nscore") + theme_cbm +
    ggplot2::theme(legend.position = "right")
}

panel_6b <- make_feature_panel("IFX_ARJS_5", "B")
panel_6c <- make_feature_panel("OSM_OSMR", "C")
panel_6d <- make_feature_panel("IL22_NEUTROPHIL_6", "D")
panel_6e <- make_feature_panel("VDZ_NR_3", "E")

localization <- data.table::fread(file.path(
  DIR_CBM_RESULTS, "singlecell_fine_localization_tests.csv"
))
localization <- localization[grepl("^(10|11)", fine_cell_type)]
selected_types <- localization[, .(max_signal = max(median_centered_score, na.rm = TRUE),
                                   min_fdr = min(localization_fdr, na.rm = TRUE)),
                               by = fine_cell_type][order(min_fdr, -max_signal)][1:min(.N, 14L),
                                                                                 fine_cell_type]
heat_local <- localization[fine_cell_type %in% selected_types]
heat_local[, signature_label := factor(signature_labels[signature_id],
                                        levels = unname(signature_labels[robust_ids]))]
type_order <- heat_local[, .(max_signal = max(median_centered_score, na.rm = TRUE)),
                         by = fine_cell_type][order(max_signal), fine_cell_type]
heat_local[, fine_factor := factor(sub("^[0-9]+[A-Z]?-", "", fine_cell_type),
                                   levels = sub("^[0-9]+[A-Z]?-", "", type_order))]
heat_local[, significance := ifelse(localization_fdr < 0.05, "*", "")]
panel_6f <- ggplot2::ggplot(heat_local,
                            ggplot2::aes(signature_label, fine_factor,
                                         fill = median_centered_score)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
  ggplot2::geom_text(ggplot2::aes(label = significance), size = 3) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                midpoint = 0) +
  ggplot2::labs(title = "F  Patient-level pseudobulk localization",
                subtitle = "Within-patient centering; * FDR < 0.05",
                x = NULL, y = NULL, fill = "Median\ncentered score") + theme_cbm +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, size = 7),
                 axis.text.y = ggplot2::element_text(size = 7))

figure_6 <- (panel_6a | panel_6b | panel_6c) /
  (panel_6d | panel_6e | panel_6f)
save_publication_figure(figure_6, "CBM_Figure6_singlecell_fine_localization", 13.2, 8.2)

# ---------- Supplementary heterogeneity figure ----------
stability <- data.table::fread(file.path(DIR_CBM_RESULTS,
                                         "meta_prediction_and_sensitivity.csv"))
stability <- stability[signature_id %in% robust_ids]
stability[, signature_label := factor(signature_labels[signature_id],
                                       levels = rev(unname(signature_labels[robust_ids])))]
stability[, analysis_label := factor(analysis,
  levels = c("all_active", "microarray_only", "randomized_active_arms",
             "excluding_largest_cohort"),
  labels = c("All active cohorts", "Microarray only", "Randomized active arms",
             "Excluding largest cohort")
)]
panel_s2a <- ggplot2::ggplot(stability,
                             ggplot2::aes(or, signature_label,
                                          colour = analysis_label)) +
  ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "#777777") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_low, xmax = ci_high),
                          height = 0.16, position = ggplot2::position_dodge(width = 0.55)) +
  ggplot2::geom_point(size = 2, position = ggplot2::position_dodge(width = 0.55)) +
  ggplot2::scale_x_log10() +
  ggplot2::labs(title = "A  Prespecified sensitivity analyses",
                x = "Odds ratio for failure per SD (95% CI)", y = NULL,
                colour = NULL) + theme_cbm +
  ggplot2::theme(legend.position = "bottom")

all_active <- stability[analysis == "all_active"]
panel_s2b <- ggplot2::ggplot(all_active,
                             ggplot2::aes(or, signature_label)) +
  ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "#777777") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = prediction_low,
                                       xmax = prediction_high),
                          height = 0, linewidth = 1.1, colour = "#BDBDBD") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_low, xmax = ci_high),
                          height = 0.16, linewidth = 0.75, colour = "#2166AC") +
  ggplot2::geom_point(size = 2.4, colour = "#B2182B") +
  ggplot2::scale_x_log10() +
  ggplot2::labs(title = "B  Confidence versus prediction intervals",
                subtitle = "Grey: 95% prediction interval; blue: 95% CI",
                x = "Odds ratio for failure", y = NULL) + theme_cbm

figure_s2 <- panel_s2a | panel_s2b
save_publication_figure(figure_s2, "CBM_FigureS2_heterogeneity_stability", 12.4, 5.2)

# ---------- Supplementary OSM ligand-receptor localization ----------
osm_tests_path <- file.path(DIR_CBM_RESULTS, "singlecell_osm_axis_localization_tests.csv")
if (file.exists(osm_tests_path)) {
  osm_tests <- data.table::fread(osm_tests_path)
  top_osm_types <- osm_tests[, .(min_fdr = min(localization_fdr, na.rm = TRUE),
                                max_signal = max(median_centered_logcpm, na.rm = TRUE)),
                            by = fine_cell_type][order(min_fdr, -max_signal)][1:min(.N, 18L),
                                                                              fine_cell_type]
  osm_plot_data <- osm_tests[fine_cell_type %in% top_osm_types]
  osm_order <- osm_plot_data[, .(signal = max(median_centered_logcpm, na.rm = TRUE)),
                             by = fine_cell_type][order(signal), fine_cell_type]
  osm_plot_data[, fine_factor := factor(sub("^[0-9]+[A-Z]?-", "", fine_cell_type),
                                        levels = sub("^[0-9]+[A-Z]?-", "", osm_order))]
  osm_plot_data[, significance := ifelse(localization_fdr < 0.05, "*", "")]
  figure_s3 <- ggplot2::ggplot(osm_plot_data,
                               ggplot2::aes(gene, fine_factor,
                                            colour = median_centered_logcpm,
                                            size = n_patients)) +
    ggplot2::geom_point() +
    ggplot2::geom_text(ggplot2::aes(label = significance), colour = "black",
                       size = 2.5, vjust = 0.35) +
    ggplot2::scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                    midpoint = 0) +
    ggplot2::labs(title = "OSM ligand and OSMR receptor localize to distinct fine cell states",
                  subtitle = "Patient-level pseudobulk; within-patient centered logCPM; * FDR < 0.05",
                  x = NULL, y = NULL, colour = "Median centered\nlogCPM",
                  size = "Patients") + theme_cbm
  save_publication_figure(figure_s3, "CBM_FigureS3_OSM_OSMR_localization", 7.5, 6.0)
}

# ---------- Graphical abstract (programmatic scientific artwork) ----------
box_data <- data.table::data.table(
  xmin = c(0.3, 3.5, 6.7, 9.9), xmax = c(2.8, 6.0, 9.2, 12.7),
  ymin = 1.4, ymax = 4.5,
  fill = c("#E8F1F8", "#F3EAF7", "#FDE9E7", "#EAF4EA"),
  title = c("PUBLIC DATA", "LOCKED COMPUTATION", "RESULT", "BIOLOGICAL MAP"),
  body = c("8 active cohorts (n=629)\n3 placebo arms (n=225)\n+ 93,900 single cells",
           "7 exact published modules\nwithin-sample ranks\n10,000 matched controls",
           "Common PC1: 47.2%\n0/7 treatment interactions\npassed FDR",
           "Inflammation / TNF / EMT\nECM and cell-state\nlocalization")
)
ga <- ggplot2::ggplot() +
  ggplot2::geom_rect(data = box_data,
                     ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin,
                                  ymax = ymax, fill = fill),
                     colour = "#334E68", linewidth = 0.55, show.legend = FALSE) +
  ggplot2::scale_fill_identity() +
  ggplot2::geom_text(data = box_data,
                     ggplot2::aes(x = (xmin + xmax) / 2, y = 3.85, label = title),
                     fontface = "bold", size = 3.3, colour = "#102A43") +
  ggplot2::geom_text(data = box_data,
                     ggplot2::aes(x = (xmin + xmax) / 2, y = 2.75, label = body),
                     size = 3.0, lineheight = 1.15, colour = "#243B53") +
  ggplot2::annotate("segment", x = 2.85, xend = 3.4, y = 2.95, yend = 2.95,
                    arrow = grid::arrow(length = grid::unit(0.13, "in")),
                    colour = "#486581", linewidth = 0.8) +
  ggplot2::annotate("segment", x = 6.05, xend = 6.6, y = 2.95, yend = 2.95,
                    arrow = grid::arrow(length = grid::unit(0.13, "in")),
                    colour = "#486581", linewidth = 0.8) +
  ggplot2::annotate("segment", x = 9.25, xend = 9.8, y = 2.95, yend = 2.95,
                    arrow = grid::arrow(length = grid::unit(0.13, "in")),
                    colour = "#486581", linewidth = 0.8) +
  ggplot2::annotate("label", x = 6.5, y = 0.55,
                    label = "Portable prognostic biology does not equal drug-specific treatment selection",
                    size = 3.6, fontface = "bold", colour = "#7F0000",
                    fill = "white", label.size = 0.35) +
  ggplot2::coord_cartesian(xlim = c(0, 13), ylim = c(0.1, 4.9), expand = FALSE) +
  ggplot2::theme_void(base_family = "sans")
ggplot2::ggsave(file.path(DIR_CBM_FIGURES, "CBM_Graphical_Abstract.png"), ga,
                width = 13.28, height = 5.31, units = "in", dpi = 200,
                bg = "white")
ggplot2::ggsave(file.path(DIR_CBM_FIGURES, "CBM_Graphical_Abstract.tiff"), ga,
                width = 13.28, height = 5.31, units = "in", dpi = 300,
                compression = "lzw", bg = "white")

message("CBM 升级图全部生成。")
