# 生成核心投稿图。所有图同时保存 PDF 和 600 dpi TIFF。

check_packages()

# 统一期刊风格。
paper_theme <- ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    axis.title = ggplot2::element_text(face = "bold"),
    legend.position = "top"
  )

# 一个保存函数，避免每幅图重复写参数。
save_paper_figure <- function(plot_object, filename, width, height) {
  ggplot2::ggsave(file.path(DIR_FIGURES, paste0(filename, ".pdf")),
                  plot_object, width = width, height = height, units = "in")
  ggplot2::ggsave(file.path(DIR_FIGURES, paste0(filename, ".tiff")),
                  plot_object, width = width, height = height, units = "in",
                  dpi = 600, compression = "lzw")
}

# ---------- Figure 1：研究设计和逻辑链 ----------
box_data <- data.table::data.table(
  xmin = c(0.2, 0.2, 0.2, 0.2, 3.4, 6.6, 6.6, 6.6, 6.6),
  xmax = c(2.5, 2.5, 2.5, 2.5, 5.7, 8.9, 8.9, 8.9, 8.9),
  ymin = c(7.2, 5.3, 3.4, 1.5, 3.4, 7.2, 5.3, 3.4, 1.5),
  ymax = c(8.5, 6.6, 4.7, 2.8, 6.6, 8.5, 6.6, 4.7, 2.8),
  label = c(
    "UNIFI\nUST + placebo",
    "PURSUIT-SC\nGLM + placebo",
    "Three IFX cohorts\nindependent validation",
    "VDZ cohort + scRNA\ncross-compartment evidence",
    "Locked 32-gene ERI\nwithin-sample ranks\ninjury - maturation\n(no outcome-based weights)",
    "Cross-biologic meta-analysis",
    "Treatment × ERI interaction",
    "Patient-paired longitudinal change",
    "Single-cell pseudobulk localization"
  ),
  type = c(rep("data", 4), "score", rep("analysis", 4))
)
workflow <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = box_data,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = type),
    color = "#333333", linewidth = 0.45, alpha = 0.9
  ) +
  ggplot2::geom_text(
    data = box_data,
    ggplot2::aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    size = 3.2, lineheight = 0.95
  ) +
  ggplot2::annotate("segment", x = 2.55, xend = 3.3, y = 4.95, yend = 4.95,
                    arrow = grid::arrow(length = grid::unit(0.12, "in")), linewidth = 0.7) +
  ggplot2::annotate("segment", x = 5.8, xend = 6.5, y = 4.95, yend = 4.95,
                    arrow = grid::arrow(length = grid::unit(0.12, "in")), linewidth = 0.7) +
  ggplot2::scale_fill_manual(values = c(data = "#D9EAF7", score = "#F4D6A0", analysis = "#DDEED8"),
                             guide = "none") +
  ggplot2::coord_cartesian(xlim = c(0, 9.1), ylim = c(1.2, 8.8), clip = "off") +
  ggplot2::labs(title = "Placebo-aware, cross-trial study design") +
  ggplot2::theme_void(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13),
                 plot.background = ggplot2::element_rect(fill = "white", color = NA),
                 panel.background = ggplot2::element_rect(fill = "white", color = NA),
                 plot.margin = ggplot2::margin(8, 8, 8, 8))
save_paper_figure(workflow, "Figure1_study_design", 11, 6.6)

# ---------- Figure 2A：六个治疗队列随机效应森林图 ----------
meta <- data.table::fread(file.path(DIR_RESULTS, "meta_analysis.csv"))
meta[, label := ifelse(cohort == "Random-effects pooled", "Pooled", paste0(cohort, " / ", treatment))]
meta[, label := factor(label, levels = rev(label))]
forest <- ggplot2::ggplot(meta, ggplot2::aes(x = or, y = label)) +
  ggplot2::geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low, xmax = ci_high),
                         orientation = "y", width = 0.15, color = "#234E70") +
  ggplot2::geom_point(ggplot2::aes(shape = cohort == "Random-effects pooled"),
                      size = 2.8, color = "#B33A3A") +
  ggplot2::scale_x_log10(breaks = c(0.5, 1, 2, 4, 8, 16)) +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
  ggplot2::labs(x = "Odds ratio for failure per 1-SD higher ERI", y = NULL,
                title = "A. Cross-biologic association") +
  paper_theme

# ---------- Figure 2B：治疗组和安慰剂组效应并列 ----------
all_effects <- data.table::fread(file.path(DIR_RESULTS, "cohort_effects.csv"))
randomized <- all_effects[cohort %in% c("GSE206285", "GSE92415")]
randomized[, trial := ifelse(cohort == "GSE206285", "UNIFI", "PURSUIT-SC")]
randomized[, treatment := factor(treatment,
                                 levels = c("placebo", "ustekinumab", "golimumab"))]
placebo_aware <- ggplot2::ggplot(randomized,
  ggplot2::aes(x = treatment, y = or, ymin = ci_low, ymax = ci_high, color = treatment)) +
  ggplot2::geom_hline(yintercept = 1, linetype = 2, color = "grey55") +
  ggplot2::geom_pointrange(size = 0.5) +
  ggplot2::facet_wrap(~trial, scales = "free_x") +
  ggplot2::scale_y_log10(breaks = c(0.5, 1, 2, 4)) +
  ggplot2::scale_color_manual(values = c(placebo = "#777777", ustekinumab = "#0072B2", golimumab = "#D55E00")) +
  ggplot2::labs(x = NULL, y = "Odds ratio for failure", title = "B. Placebo-aware trial analysis") +
  paper_theme + ggplot2::theme(legend.position = "none")

# 如果 patchwork 存在，拼成一幅图；否则分别保存也不影响分析。
if (requireNamespace("patchwork", quietly = TRUE)) {
  figure2 <- forest + placebo_aware + patchwork::plot_layout(widths = c(1.35, 1))
  save_paper_figure(figure2, "Figure2_forest_and_placebo", 10.5, 5.5)
} else {
  save_paper_figure(forest, "Figure2A_forest", 6.2, 5.5)
  save_paper_figure(placebo_aware, "Figure2B_placebo", 5.2, 4.8)
}

# ---------- Figure 3：不同响应组的患者配对 ERI 变化 ----------
pairs <- data.table::fread(file.path(DIR_RESULTS, "longitudinal_pairs.csv"))
pairs[, pair_id := paste(cohort, treatment, patient_id, sep = "__")]
long_plot_data <- data.table::melt(
  pairs,
  id.vars = c("pair_id", "cohort", "treatment", "response"),
  measure.vars = c("eri_baseline", "eri_followup"),
  variable.name = "timepoint",
  value.name = "eri"
)
long_plot_data[, timepoint := factor(timepoint,
                                     levels = c("eri_baseline", "eri_followup"),
                                     labels = c("Baseline", "Follow-up"))]
longitudinal_plot <- ggplot2::ggplot(
  long_plot_data,
  ggplot2::aes(x = timepoint, y = eri, group = pair_id, color = response)
) +
  ggplot2::geom_line(alpha = 0.28, linewidth = 0.35) +
  ggplot2::geom_point(alpha = 0.65, size = 1.2) +
  ggplot2::stat_summary(ggplot2::aes(group = response), fun = stats::median,
                        geom = "line", linewidth = 1.2) +
  ggplot2::facet_wrap(~cohort + treatment, scales = "free_y") +
  ggplot2::scale_color_manual(values = c(R = "#0072B2", NR = "#D55E00")) +
  ggplot2::labs(x = NULL, y = "Epithelial restitution imbalance (ERI)",
                color = "Response", title = "Longitudinal molecular recovery") +
  paper_theme
save_paper_figure(longitudinal_plot, "Figure3_longitudinal_ERI", 10.5, 6.4)

# ---------- Supplementary Figure：各基线队列响应/不响应分布 ----------
samples <- data.table::fread(file.path(DIR_RESULTS, "baseline_analysis_samples.csv"))
samples[, outcome := factor(failure, levels = c(0, 1), labels = c("Response/healing", "Failure"))]
baseline_plot <- ggplot2::ggplot(samples, ggplot2::aes(x = outcome, y = eri_z, fill = outcome)) +
  ggplot2::geom_violin(trim = FALSE, alpha = 0.55, color = NA) +
  ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.8) +
  ggplot2::facet_wrap(~cohort + treatment, scales = "free_x") +
  ggplot2::scale_fill_manual(values = c("Response/healing" = "#56B4E9", "Failure" = "#E69F00")) +
  ggplot2::labs(x = NULL, y = "ERI (within-cohort SD)", title = "Baseline ERI by outcome") +
  paper_theme + ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
save_paper_figure(baseline_plot, "FigureS1_baseline_distributions", 11.5, 7.0)
