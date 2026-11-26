# 为系统基准部分生成三张主图。图形只读取已经锁定的统计结果，不重新拟合模型。

needed_plot_packages <- c("ggplot2", "patchwork", "scales")
plot_ok <- vapply(needed_plot_packages, requireNamespace, logical(1), quietly = TRUE)
if (!all(plot_ok)) {
  stop("缺少作图包：", paste(needed_plot_packages[!plot_ok], collapse = ", "),
       "。请先运行 scripts/00_install_packages.R")
}

paper_theme_benchmark <- ggplot2::theme_classic(base_size = 10.5) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 11.5),
    plot.subtitle = ggplot2::element_text(size = 9.5, color = "grey30"),
    axis.title = ggplot2::element_text(face = "bold"),
    legend.position = "top",
    strip.background = ggplot2::element_rect(fill = "grey94", color = "grey70"),
    strip.text = ggplot2::element_text(face = "bold")
  )

save_benchmark_figure <- function(plot_object, filename, width, height) {
  ggplot2::ggsave(
    file.path(DIR_FIGURES, paste0(filename, ".pdf")), plot_object,
    width = width, height = height, units = "in", bg = "white"
  )
  ggplot2::ggsave(
    file.path(DIR_FIGURES, paste0(filename, ".tiff")), plot_object,
    width = width, height = height, units = "in", dpi = 600,
    compression = "lzw", bg = "white"
  )
}

registry_plot <- read_signature_registry()
signature_labels <- unique(registry_plot[, .(signature_id, display_name, primary_set)])
primary_ids_plot <- signature_labels[primary_set == TRUE, signature_id]

# ---------- Main Figure 1：升级后的系统基准研究设计 ----------
design_boxes <- data.table::data.table(
  xmin = c(0.2, 0.2, 3.35, 3.35, 3.35, 6.5, 6.5, 6.5),
  xmax = c(2.65, 2.65, 5.8, 5.8, 5.8, 8.95, 8.95, 8.95),
  ymin = c(5.25, 2.15, 6.25, 3.7, 1.15, 6.25, 3.7, 1.15),
  ymax = c(7.35, 4.25, 8.0, 5.45, 2.9, 8.0, 5.45, 2.9),
  label = c(
    "Public trial transcriptomes\n8 active cohorts (n=629)\n5 agents / 4 target classes",
    "Randomized controls\n3 placebo arms (n=225)\nUNIFI / PURSUIT-SC / ACT-1",
    "Outcome-blind registry\n7 exact published modules\n4 sensitivity comparators",
    "Portable scoring\nwithin-sample gene ranks\nno cross-platform merging",
    "Source-cohort exclusion\n10,000 matched random sets\nrandom-effects Meta-analysis",
    "Prognostic evidence\nactive + placebo effects\nseverity residualization",
    "Predictive evidence\n3 treatment-by-signature tests\nand interaction Meta-analysis",
    "Convergence and calibration\noutcome-free PCA\ncontext-restricted recovery"
  ),
  group = c("Data", "Data", "Locked methods", "Locked methods", "Locked methods",
            "Inference", "Inference", "Inference")
)
design_arrows <- data.table::data.table(
  x = c(2.72, 2.72, 5.87, 5.87, 5.87),
  xend = c(3.27, 3.27, 6.42, 6.42, 6.42),
  y = c(6.3, 3.2, 7.1, 4.58, 2.03),
  yend = c(7.1, 2.03, 7.1, 4.58, 2.03)
)
benchmark_design <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = design_boxes,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = group),
    color = "#333333", linewidth = 0.45
  ) +
  ggplot2::geom_text(
    data = design_boxes,
    ggplot2::aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    size = 3.15, lineheight = 0.95
  ) +
  ggplot2::geom_segment(
    data = design_arrows,
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
    arrow = grid::arrow(length = grid::unit(0.11, "in")), linewidth = 0.65
  ) +
  ggplot2::scale_fill_manual(
    values = c("Data" = "#D9EAF7", "Locked methods" = "#F4D6A0", "Inference" = "#DDEED8"),
    guide = "none"
  ) +
  ggplot2::coord_cartesian(xlim = c(0, 9.15), ylim = c(0.9, 8.2), clip = "off") +
  ggplot2::labs(
    title = "Placebo-aware benchmark of eligible published mucosal response signatures",
    subtitle = "Primary inference was locked to exact-text signatures before outcome analysis"
  ) +
  ggplot2::theme_void(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
    plot.margin = ggplot2::margin(8, 8, 8, 8)
  )
save_benchmark_figure(benchmark_design, "MainFigure1_upgraded_study_design", 11.5, 7.2)

# ---------- Figure 5：跨 8 个治疗队列的效应及 Meta ----------
cohort_effects <- data.table::fread(file.path(DIR_RESULTS, "benchmark_cohort_effects.csv"))
cohort_effects <- cohort_effects[
  signature_id %in% primary_ids_plot &
    treatment %in% c("ustekinumab", "golimumab", "infliximab", "vedolizumab", "etrolizumab")
]
cohort_effects <- merge(cohort_effects, signature_labels, by = "signature_id", all.x = TRUE)
cohort_order <- c("GSE206285", "GSE92415", "GSE16879", "GSE12251",
                  "GSE73661_IFX", "GSE73661_VDZ", "GSE23597", "GSE72819")
cohort_effects[, cohort := factor(cohort, levels = cohort_order)]
cohort_effects[, display_name := factor(display_name, levels = rev(unique(signature_labels[
  signature_id %in% primary_ids_plot, display_name
])))]
cohort_effects[, significance := ifelse(p_value < 0.05, "*", "")]

effect_heatmap <- ggplot2::ggplot(
  cohort_effects, ggplot2::aes(x = cohort, y = display_name, fill = log_or)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.35) +
  ggplot2::geom_text(ggplot2::aes(label = significance), size = 4.1) +
  ggplot2::scale_fill_gradient2(
    low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
    limits = c(-1.5, 1.5), oob = scales::squish,
    name = "log OR\nfor failure"
  ) +
  ggplot2::labs(
    x = NULL, y = NULL,
    title = "A. Published signatures across eight active-treatment cohorts",
    subtitle = "Positive values indicate higher odds of treatment failure; * nominal P < 0.05"
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold")
  )

active_meta_plot <- data.table::fread(file.path(DIR_RESULTS, "benchmark_active_meta.csv"))
external_meta_plot <- data.table::fread(file.path(DIR_RESULTS, "benchmark_leave_source_out_meta.csv"))
forest_data <- data.table::rbindlist(list(
  active_meta_plot[signature_id %in% primary_ids_plot, .(
    signature_id, estimate = or, low = ci_low, high = ci_high, analysis = "All treated cohorts"
  )],
  external_meta_plot[signature_id %in% primary_ids_plot, .(
    signature_id, estimate = or, low = ci_low, high = ci_high, analysis = "Source cohorts excluded"
  )]
), fill = TRUE)
forest_data <- merge(forest_data, signature_labels, by = "signature_id", all.x = TRUE)
forest_data[, display_name := factor(display_name, levels = levels(cohort_effects$display_name))]

meta_forest <- ggplot2::ggplot(
  forest_data,
  ggplot2::aes(x = estimate, y = display_name, xmin = low, xmax = high, color = analysis)
) +
  ggplot2::geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
  ggplot2::geom_errorbar(orientation = "y", width = 0.14,
                         position = ggplot2::position_dodge(width = 0.48)) +
  ggplot2::geom_point(size = 2.2, position = ggplot2::position_dodge(width = 0.48)) +
  ggplot2::scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8, 16)) +
  ggplot2::scale_color_manual(values = c(
    "All treated cohorts" = "#A51C30", "Source cohorts excluded" = "#2166AC"
  )) +
  ggplot2::labs(
    x = "Random-effects odds ratio per 1-SD higher score", y = NULL, color = NULL,
    title = "B. Full and source-cohort-excluded Meta-analysis"
  ) + paper_theme_benchmark

figure5 <- effect_heatmap / meta_forest + patchwork::plot_layout(heights = c(1.18, 1))
save_benchmark_figure(figure5, "Figure5_published_signature_benchmark", 11.5, 10.0)
save_benchmark_figure(figure5, "MainFigure2_published_signature_benchmark", 11.5, 10.0)

# ---------- Figure 6：安慰剂效应与药物特异性交互 ----------
evidence_plot <- data.table::fread(file.path(DIR_RESULTS, "published_signature_evidence_matrix.csv"))
evidence_plot[, short_label := data.table::fcase(
  signature_id == "IFX_ARJS_5", "Arijs",
  signature_id == "VDZ_4", "VDZ4",
  signature_id == "OSM_OSMR", "OSM",
  signature_id == "ETRO_2", "ETRO",
  signature_id == "IL22_NEUTROPHIL_6", "IL-22",
  signature_id == "VDZ_NR_3", "VDZ-NR",
  signature_id == "VDZ_R_5", "VDZ-R",
  default = signature_id
)]
specificity_scatter <- ggplot2::ggplot(
  evidence_plot,
  ggplot2::aes(x = log(or), y = log(placebo_or), color = evidence_class, label = short_label)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey65") +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey65") +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 3, color = "grey55") +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(nudge_x = 0.018, nudge_y = 0.045, check_overlap = FALSE, size = 3) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.16))) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.18))) +
  ggplot2::scale_color_manual(values = c(
    general_healing_propensity_supported = "#D55E00",
    treated_association_without_specificity = "#0072B2",
    drug_specific_interaction_supported = "#009E73",
    insufficient_or_inconsistent = "grey55"
  )) +
  ggplot2::labs(
    x = "log OR in active-treatment cohorts",
    y = "log OR in placebo arms", color = "Evidence class",
    title = "A. Active-treatment association versus placebo association",
    subtitle = "A positive placebo effect is incompatible with a purely drug-specific predictor"
  ) + paper_theme_benchmark

interaction_plot <- data.table::fread(file.path(DIR_RESULTS, "benchmark_interaction_meta.csv"))
interaction_plot <- interaction_plot[signature_id %in% primary_ids_plot]
interaction_plot <- merge(interaction_plot, signature_labels, by = "signature_id", all.x = TRUE)
interaction_plot[, display_name := factor(display_name, levels = levels(cohort_effects$display_name))]
interaction_forest <- ggplot2::ggplot(
  interaction_plot,
  ggplot2::aes(x = or, y = display_name, xmin = ci_low, xmax = ci_high)
) +
  ggplot2::geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
  ggplot2::geom_errorbar(orientation = "y", width = 0.15, color = "#4D4D4D") +
  ggplot2::geom_point(ggplot2::aes(color = fdr < 0.05), size = 2.5) +
  ggplot2::scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
  ggplot2::scale_color_manual(values = c(`FALSE` = "#4D4D4D", `TRUE` = "#009E73"),
                              guide = "none") +
  ggplot2::labs(
    x = "Pooled treatment × signature interaction OR", y = NULL,
    title = "B. Drug-specific interaction across three placebo-controlled trials",
    subtitle = "Interaction Meta-analysis: UNIFI, PURSUIT-SC and ACT-1"
  ) + paper_theme_benchmark

figure6 <- specificity_scatter + interaction_forest + patchwork::plot_layout(widths = c(1.08, 1))
save_benchmark_figure(figure6, "Figure6_placebo_and_specificity", 12.0, 6.2)
save_benchmark_figure(figure6, "MainFigure3_placebo_and_specificity", 12.0, 6.2)

# ---------- Figure 7：匹配随机基因集与无监督收敛 ----------
empirical_plot <- data.table::fread(file.path(DIR_RESULTS, "matched_random_empirical_tests.csv"))
empirical_plot <- empirical_plot[signature_id %in% primary_ids_plot]
empirical_plot <- merge(empirical_plot, signature_labels, by = "signature_id", all.x = TRUE)
empirical_plot[, display_name := factor(display_name, levels = levels(cohort_effects$display_name))]
empirical_panel <- ggplot2::ggplot(
  empirical_plot,
  ggplot2::aes(x = -log10(empirical_p_one_sided), y = display_name,
               color = empirical_fdr < 0.05)
) +
  ggplot2::geom_vline(xintercept = -log10(0.05), linetype = 2, color = "grey60") +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = -log10(empirical_p_one_sided), yend = display_name),
                        color = "grey75") +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::scale_color_manual(values = c(`FALSE` = "grey55", `TRUE` = "#A51C30"),
                              guide = "none") +
  ggplot2::labs(
    x = expression(-log[10]("matched-random empirical P")), y = NULL,
    title = "A. Matched-random gene-set benchmark"
  ) + paper_theme_benchmark

correlation_plot <- data.table::fread(file.path(DIR_RESULTS, "signature_score_correlations.csv"))
correlation_plot <- correlation_plot[
  signature_1 %in% primary_ids_plot & signature_2 %in% primary_ids_plot
]
correlation_panel <- ggplot2::ggplot(
  correlation_plot,
  ggplot2::aes(x = signature_1, y = signature_2, fill = spearman_rho)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.25) +
  ggplot2::scale_fill_gradient2(
    low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
    limits = c(-1, 1), name = "Spearman\ncorrelation"
  ) +
    ggplot2::labs(x = NULL, y = NULL, title = "B. Outcome-free convergence") +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold")
  )

loadings_plot <- data.table::fread(file.path(DIR_RESULTS, "signature_convergence_loadings.csv"))
variance_plot <- data.table::fread(file.path(DIR_RESULTS, "signature_convergence_variance.csv"))
loadings_plot <- merge(loadings_plot, signature_labels, by = "signature_id", all.x = TRUE)
loadings_plot[, display_name := factor(display_name, levels = levels(cohort_effects$display_name))]
pc1_percent <- 100 * variance_plot[component == "PC1", variance_explained][1L]
loadings_panel <- ggplot2::ggplot(
  loadings_plot, ggplot2::aes(x = pc1_loading, y = display_name, fill = pc1_loading > 0)
) +
  ggplot2::geom_vline(xintercept = 0, color = "grey65") +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::scale_fill_manual(values = c(`FALSE` = "#2C7BB6", `TRUE` = "#D7191C"), guide = "none") +
  ggplot2::labs(
    x = "PC1 loading", y = NULL,
    title = paste0("C. PC1 loadings (", sprintf("%.1f", pc1_percent), "% variance)")
  ) + paper_theme_benchmark

figure7 <- empirical_panel + correlation_panel + loadings_panel +
  patchwork::plot_layout(widths = c(1.0, 1.08, 0.9))
save_benchmark_figure(figure7, "Figure7_random_controls_and_convergence", 15.0, 6.2)
save_benchmark_figure(figure7, "MainFigure4_random_controls_and_convergence", 15.0, 6.2)

message("最终投稿 Main Figure 1-4 和系统基准 Figure 5-7 已保存（PDF + 600 dpi TIFF）。")
