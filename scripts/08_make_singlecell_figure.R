# 用已经聚合的小表作单细胞定位图，不再打开 Seurat 对象。

paper_theme <- ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    axis.title = ggplot2::element_text(face = "bold"),
    legend.position = "top"
  )

# 与核心作图脚本相同的双格式保存函数。
save_singlecell_figure <- function(plot_object, filename, width, height) {
  ggplot2::ggsave(file.path(DIR_FIGURES, paste0(filename, ".pdf")),
                  plot_object, width = width, height = height, units = "in")
  ggplot2::ggsave(file.path(DIR_FIGURES, paste0(filename, ".tiff")),
                  plot_object, width = width, height = height, units = "in",
                  dpi = 600, compression = "lzw")
}

# ---------- A：不同细胞类型的重度 UC-健康位移 ----------
shift <- data.table::fread(file.path(DIR_RESULTS, "singlecell_celltype_shift.csv"))
shift[, clean_cell_type := sub("^[0-9]+-", "", coarse_annotations_MK)]
shift[, clean_cell_type := factor(clean_cell_type,
                                  levels = clean_cell_type[order(shift_severe_minus_hc)])]
shift_plot <- ggplot2::ggplot(
  shift,
  ggplot2::aes(x = shift_severe_minus_hc, y = clean_cell_type,
               color = coarse_annotations_MK == "10-Epithelial")
) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = shift_severe_minus_hc,
                                     yend = clean_cell_type), linewidth = 0.7) +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::scale_color_manual(values = c(`FALSE` = "#777777", `TRUE` = "#B33A3A"), guide = "none") +
  ggplot2::labs(x = "Median ERI shift: severe UC minus healthy", y = NULL,
                title = "A. Cell-compartment localization") +
  paper_theme

# ---------- B：上皮伪汇总随内镜严重度递增 ----------
pseudobulk <- data.table::fread(file.path(DIR_RESULTS, "singlecell_pseudobulk.csv"))
epithelial <- pseudobulk[coarse_annotations_MK == "10-Epithelial"]
epithelial[, severity := factor(
  Endoscopic_Severity_V1,
  levels = c("HC", "Severity_0_1", "Severity_2_3"),
  labels = c("Healthy", "Mayo 0-1", "Mayo 2-3")
)]
epithelial_plot <- ggplot2::ggplot(epithelial, ggplot2::aes(x = severity, y = eri, fill = severity)) +
  ggplot2::geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.65) +
  ggplot2::geom_jitter(width = 0.09, size = 2, alpha = 0.8) +
  ggplot2::scale_fill_manual(values = c("#56B4E9", "#E69F00", "#D55E00"), guide = "none") +
  ggplot2::labs(x = NULL, y = "Epithelial pseudobulk ERI",
                title = "B. Epithelial ERI tracks severity") +
  paper_theme

# 拼图并保存。
if (requireNamespace("patchwork", quietly = TRUE)) {
  figure4 <- shift_plot + epithelial_plot + patchwork::plot_layout(widths = c(1.25, 1))
  save_singlecell_figure(figure4, "Figure4_singlecell_localization", 10, 5.4)
} else {
  save_singlecell_figure(shift_plot, "Figure4A_celltype_shift", 6, 5.4)
  save_singlecell_figure(epithelial_plot, "Figure4B_epithelial_severity", 4.8, 4.8)
}

