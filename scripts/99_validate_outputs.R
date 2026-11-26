# 自动核对关键文件和结论，避免手工复制数字时漏项。

required_files <- c(
  file.path(DIR_RESULTS, "quality_gate.csv"),
  file.path(DIR_RESULTS, "treatment_interactions.csv"),
  file.path(DIR_RESULTS, "clinical_adjusted_effects.csv"),
  file.path(DIR_RESULTS, "longitudinal_summary.csv"),
  file.path(DIR_RESULTS, "singlecell_epithelial_test.csv"),
  file.path(DIR_FIGURES, "Figure1_study_design.pdf"),
  file.path(DIR_FIGURES, "Figure2_forest_and_placebo.pdf"),
  file.path(DIR_FIGURES, "Figure3_longitudinal_ERI.pdf"),
  file.path(DIR_FIGURES, "Figure4_singlecell_localization.pdf")
)

# 读取结果表。
quality <- data.table::fread(file.path(DIR_RESULTS, "quality_gate.csv"))
interactions <- data.table::fread(file.path(DIR_RESULTS, "treatment_interactions.csv"))
adjusted <- data.table::fread(file.path(DIR_RESULTS, "clinical_adjusted_effects.csv"))
longitudinal <- data.table::fread(file.path(DIR_RESULTS, "longitudinal_summary.csv"))
singlecell <- data.table::fread(file.path(DIR_RESULTS, "singlecell_epithelial_test.csv"))
leave_one_out <- data.table::fread(file.path(DIR_RESULTS, "meta_leave_one_out.csv"))

# 每一行是一个可复现检查。
validation <- data.table::data.table(
  check = c(
    "all required files exist",
    "six treated cohorts",
    "all treated directions positive",
    "pooled confidence interval excludes one",
    "leave-one-out confidence intervals exclude one",
    "two treatment interactions are non-significant",
    "clinically adjusted UNIFI effect is positive",
    "four active longitudinal delta tests are below 0.05",
    "single-cell epithelial severity trend is below 0.05"
  ),
  passed = c(
    all(file.exists(required_files) & file.info(required_files)$size > 100L),
    quality$treated_cohorts == 6L,
    quality$positive_direction == 6L && quality$all_same_direction,
    quality$pooled_ci_low > 1,
    all(leave_one_out$ci_low > 1),
    nrow(interactions) == 2L && all(interactions$p_value > 0.05),
    adjusted[model == "clinical_adjusted", or] > 1,
    length(unique(longitudinal[treatment != "placebo" & delta_between_response_p < 0.05, cohort])) == 4L,
    singlecell$spearman_p < 0.05
  )
)

# 保存并打印检查表；任一失败就停止。
data.table::fwrite(validation, file.path(DIR_RESULTS, "validation_report.csv"))
print(validation)
if (!all(validation$passed)) {
  stop("输出验证未全部通过，请检查 results/validation_report.csv")
}
message("全部自动检查通过。")

