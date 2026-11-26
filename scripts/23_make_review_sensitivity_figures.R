# 审稿意见定向升级图：不重算原始数据，只可视化 scripts/22 的派生结果。

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("缺少 ggplot2；请先运行 scripts/21_install_cbm_upgrade_packages.R")
}

DIR_CBM_RESULTS <- file.path(DIR_RESULTS, "cbm_upgrade")
DIR_CBM_FIGURES <- file.path(DIR_FIGURES, "cbm_upgrade")
dir.create(DIR_CBM_FIGURES, recursive = TRUE, showWarnings = FALSE)

signature_labels <- c(
  IFX_ARJS_5 = "Arijs anti-TNF",
  VDZ_4 = "Vedolizumab 4-gene",
  OSM_OSMR = "OSM-OSMR",
  ETRO_2 = "GZMA-ITGAE",
  IL22_NEUTROPHIL_6 = "IL-22 chemokines",
  VDZ_NR_3 = "VDZ nonresponse stromal",
  VDZ_R_5 = "VDZ response crypt-base"
)
robust_ids <- c("IFX_ARJS_5", "OSM_OSMR", "IL22_NEUTROPHIL_6", "VDZ_NR_3")

theme_submission <- ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "#F0F3F5", colour = "#4A4A4A"),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold", size = 11),
    plot.subtitle = ggplot2::element_text(size = 9),
    axis.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )

save_two_formats <- function(plot, stem, width, height) {
  ggplot2::ggsave(
    file.path(DIR_CBM_FIGURES, paste0(stem, ".pdf")), plot,
    width = width, height = height, units = "in", device = grDevices::pdf,
    useDingbats = FALSE
  )
  ggplot2::ggsave(
    file.path(DIR_CBM_FIGURES, paste0(stem, ".tiff")), plot,
    width = width, height = height, units = "in", dpi = 300,
    compression = "lzw", bg = "white"
  )
}

# Supplementary Figure S3: 交互效应精度。
interaction <- data.table::fread(file.path(
  DIR_CBM_RESULTS, "interaction_precision_audit.csv"
))
interaction[, label := unname(signature_labels[signature_id])]
interaction[, label := factor(label, levels = rev(unname(signature_labels)))]

figure_s3 <- ggplot2::ggplot(interaction, ggplot2::aes(y = label, x = or)) +
  ggplot2::annotate(
    "rect", xmin = 2 / 3, xmax = 1.5, ymin = -Inf, ymax = Inf,
    fill = "#DCEFE3", alpha = 0.70
  ) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "#333333") +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.18,
    linewidth = 0.65, colour = "#24557A"
  ) +
  ggplot2::geom_point(size = 2.6, shape = 21, fill = "#E07A5F", colour = "#7A2E1D") +
  ggplot2::scale_x_log10(
    breaks = c(0.5, 2 / 3, 1, 1.5, 2, 3, 4),
    labels = c("0.50", "0.67", "1.00", "1.50", "2.00", "3.00", "4.00")
  ) +
  ggplot2::coord_cartesian(xlim = c(0.48, 4.1)) +
  ggplot2::labs(
    title = "Randomised treatment-by-signature interaction precision",
    subtitle = paste0("Three placebo-controlled trials (n=",
                      unique(interaction$total_randomized_n),
                      "); shaded band is exploratory"),
    x = "Interaction odds ratio (95% CI; logarithmic scale)", y = NULL,
    caption = paste0(
      "The 0.67-1.50 margin is not clinically validated. No 95% CI was contained within it; ",
      "non-significance is not equivalence."
    )
  ) +
  theme_submission +
  ggplot2::theme(legend.position = "none")
save_two_formats(figure_s3, "APT_FigureS3_interaction_precision", 7.5, 5.2)

# Supplementary Figure S4: 客观/混合黏膜终点 vs 纯临床终点。
endpoint <- data.table::fread(file.path(
  DIR_CBM_RESULTS, "endpoint_stratified_meta.csv"
))
endpoint <- endpoint[signature_id %in% robust_ids & stratum %in% c(
  "all_active", "objective_or_composite_mucosal", "clinical_only"
)]
endpoint[, signature_label := factor(
  unname(signature_labels[signature_id]),
  levels = unname(signature_labels[robust_ids])
)]
endpoint[, stratum_label := factor(
  stratum,
  levels = c("clinical_only", "objective_or_composite_mucosal", "all_active"),
  labels = c("Clinical only (3 cohorts)",
             "Objective/composite mucosal (5 cohorts)",
             "All active-treatment cohorts (8 cohorts)")
)]

figure_s4 <- ggplot2::ggplot(endpoint, ggplot2::aes(y = stratum_label, x = or)) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "#555555") +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.16,
    linewidth = 0.55, colour = "#24557A"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(size = total_n, fill = stratum_label),
    shape = 21, colour = "#333333", stroke = 0.4
  ) +
  ggplot2::facet_wrap(~ signature_label, ncol = 2) +
  ggplot2::scale_x_log10(breaks = c(1, 1.5, 2, 3, 4, 6)) +
  ggplot2::scale_fill_manual(values = c(
    "Clinical only (3 cohorts)" = "#F2CC8F",
    "Objective/composite mucosal (5 cohorts)" = "#81B29A",
    "All active-treatment cohorts (8 cohorts)" = "#3D5A80"
  )) +
  ggplot2::scale_size_continuous(range = c(2.2, 4.6)) +
  ggplot2::labs(
    title = "Endpoint-stratified random-effects sensitivity analysis",
    subtitle = "Failure-risk association per 1-SD higher baseline signature score",
    x = "Odds ratio (95% CI; logarithmic scale)", y = NULL,
    size = "Participants", fill = "Endpoint stratum",
    caption = "Strata retain each study's locked endpoint definition; between-stratum tests are exploratory."
  ) +
  theme_submission +
  ggplot2::guides(fill = "none", size = "none") +
  ggplot2::theme(legend.position = "none")
save_two_formats(figure_s4, "APT_FigureS4_endpoint_stratified", 8.2, 6.5)

# Supplementary Figure S5: 单细胞定位逐患者剔除稳定性。
lopo <- data.table::fread(file.path(
  DIR_CBM_RESULTS, "singlecell_leave_one_patient_out_summary.csv"
))
lopo <- lopo[order(signature_id, original_localization_fdr,
                   -original_median_centered_score)]
lopo[, cell_label := factor(fine_cell_type, levels = rev(unique(fine_cell_type)))]
lopo[, signature_label := factor(
  unname(signature_labels[signature_id]),
  levels = unname(signature_labels[robust_ids])
)]
lopo[, stability_percent := 100 * fdr_significance_fraction]

figure_s5 <- ggplot2::ggplot(
  lopo, ggplot2::aes(y = cell_label, x = stability_percent)
) +
  ggplot2::geom_vline(xintercept = 100, colour = "#999999", linewidth = 0.4) +
  ggplot2::geom_col(width = 0.68, fill = "#4C78A8") +
  ggplot2::geom_point(
    ggplot2::aes(shape = sign_stable_all_omissions),
    x = 103, size = 2.2, colour = "#B23A48"
  ) +
  ggplot2::facet_wrap(~ signature_label, scales = "free_y", ncol = 2) +
  ggplot2::scale_x_continuous(
    limits = c(0, 108), breaks = c(0, 25, 50, 75, 100),
    labels = function(x) paste0(x, "%")
  ) +
  ggplot2::scale_shape_manual(
    values = c(`TRUE` = 8, `FALSE` = 4),
    labels = c(`TRUE` = "positive after every omission",
               `FALSE` = "sign changed after at least one omission")
  ) +
  ggplot2::labs(
    title = "Leave-one-patient-out stability of fine-cell localisation",
    subtitle = paste0(
      "Originally FDR-significant tests; all tests re-estimated and BH-adjusted after each of ",
      unique(lopo$omissions), " patient omissions"
    ),
    x = "Patient omissions retaining localisation FDR <0.05", y = NULL,
    shape = "Direction stability",
    caption = "Inference uses patients, not cells; 93,900 cells originated from 12 patients."
  ) +
  theme_submission +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 7.5),
    legend.position = "bottom"
  )
save_two_formats(figure_s5, "APT_FigureS5_singlecell_LOPO", 9.0, 8.0)

message("审稿意见定向补充图 S3-S5 已生成。")
