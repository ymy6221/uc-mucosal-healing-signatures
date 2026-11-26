# 审稿意见定向升级的技术验收，包括公共输入再次 SHA-256 核对。

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("缺少 digest 包，无法验证公共输入 SHA-256。")
}

DIR_CBM_RESULTS <- file.path(DIR_RESULTS, "cbm_upgrade")
DIR_CBM_FIGURES <- file.path(DIR_FIGURES, "cbm_upgrade")

required_results <- c(
  "interaction_precision_audit.csv",
  "endpoint_harmonization_registry.csv",
  "endpoint_stratified_meta.csv",
  "endpoint_subgroup_difference.csv",
  "singlecell_leave_one_patient_out.csv",
  "singlecell_leave_one_patient_out_summary.csv"
)
required_figures <- unlist(lapply(c(
  "APT_FigureS3_interaction_precision",
  "APT_FigureS4_endpoint_stratified",
  "APT_FigureS5_singlecell_LOPO"
), function(stem) paste0(stem, c(".pdf", ".tiff"))))

result_paths <- file.path(DIR_CBM_RESULTS, required_results)
figure_paths <- file.path(DIR_CBM_FIGURES, required_figures)
interaction <- data.table::fread(file.path(
  DIR_CBM_RESULTS, "interaction_precision_audit.csv"
))
endpoint_registry <- data.table::fread(file.path(
  DIR_CBM_RESULTS, "endpoint_harmonization_registry.csv"
))
endpoint_meta <- data.table::fread(file.path(
  DIR_CBM_RESULTS, "endpoint_stratified_meta.csv"
))
lopo <- data.table::fread(file.path(
  DIR_CBM_RESULTS, "singlecell_leave_one_patient_out.csv"
))
lopo_summary <- data.table::fread(file.path(
  DIR_CBM_RESULTS, "singlecell_leave_one_patient_out_summary.csv"
))

# 按原锁定清单现场重算所有 13 个公共输入的 SHA-256。
locked <- data.table::fread(file.path(
  DIR_CBM_RESULTS, "public_input_immutability_check.csv"
))
locked[, current_bytes := ifelse(file.exists(path), file.info(path)$size, NA_real_)]
locked[, current_sha256 := vapply(path, function(one_path) {
  if (!file.exists(one_path)) return(NA_character_)
  digest::digest(file = one_path, algo = "sha256", serialize = FALSE)
}, character(1))]
locked[, `:=`(
  current_bytes_passed = current_bytes == expected_bytes,
  current_sha256_passed = current_sha256 == expected_sha256
)]
locked[, current_passed := current_bytes_passed & current_sha256_passed]
data.table::fwrite(locked, file.path(
  DIR_CBM_RESULTS, "review_upgrade_public_input_immutability_check.csv"
))

primary_ids <- unique(read_signature_registry()[primary_set == TRUE, signature_id])
validation <- data.table::data.table(
  check = c(
    "all review-sensitivity result tables exist",
    "all review-sensitivity figures exist in PDF and 300-dpi TIFF",
    "interaction precision audit contains seven primary signatures",
    "interaction audit uses three trials and 674 randomized participants",
    "no interaction is incorrectly declared equivalent within the exploratory margin",
    "endpoint registry contains all eight active-treatment cohorts",
    "endpoint strata contain five objective/composite and three clinical cohorts",
    "endpoint meta-analysis covers seven signatures and five strata",
    "LOPO analysis includes every one of 12 patients",
    "LOPO summary covers all originally FDR-significant localization tests",
    "all 13 public inputs retain locked byte counts and SHA-256 hashes",
    "upstream upgrade validation remains PASS"
  ),
  observed = c(
    paste0(sum(file.exists(result_paths) & file.info(result_paths)$size > 100L),
           "/", length(result_paths)),
    paste0(sum(file.exists(figure_paths) & file.info(figure_paths)$size > 100L),
           "/", length(figure_paths)),
    paste0(data.table::uniqueN(interaction$signature_id), "/7"),
    paste0(unique(interaction$trials), " trials; n=",
           unique(interaction$total_randomized_n)),
    paste0(sum(interaction$equivalence_supported_05), "/7 equivalence calls"),
    paste0(data.table::uniqueN(endpoint_registry$cohort), "/8"),
    paste(endpoint_registry[, .N, by = endpoint_family][,
      paste0(endpoint_family, "=", N)], collapse = "; "),
    paste0(data.table::uniqueN(endpoint_meta$signature_id), " signatures; ",
           data.table::uniqueN(endpoint_meta$stratum), " strata"),
    paste0(data.table::uniqueN(lopo$omitted_patient), "/12"),
    paste0(nrow(lopo_summary), " tests"),
    paste0(sum(locked$current_passed), "/", nrow(locked)),
    paste0(sum(data.table::fread(file.path(
      DIR_CBM_RESULTS, "cbm_upgrade_validation.csv"
    ))$passed), "/12")
  ),
  passed = c(
    all(file.exists(result_paths) & file.info(result_paths)$size > 100L),
    all(file.exists(figure_paths) & file.info(figure_paths)$size > 100L),
    setequal(interaction$signature_id, primary_ids) && nrow(interaction) == 7L,
    all(interaction$trials == 3L) &&
      all(interaction$total_randomized_n == 674L) &&
      all(interaction$active_randomized_n == 449L) &&
      all(interaction$placebo_randomized_n == 225L),
    !any(interaction$equivalence_supported_05),
    data.table::uniqueN(endpoint_registry$cohort) == 8L,
    nrow(endpoint_registry[endpoint_family ==
                             "objective_or_composite_mucosal"]) == 5L &&
      nrow(endpoint_registry[endpoint_family == "clinical_only"]) == 3L,
    setequal(endpoint_meta$signature_id, primary_ids) &&
      data.table::uniqueN(endpoint_meta$stratum) == 5L,
    data.table::uniqueN(lopo$omitted_patient) == 12L,
    nrow(lopo_summary) == nrow(data.table::fread(file.path(
      DIR_CBM_RESULTS, "singlecell_fine_localization_tests.csv"
    ))[is.finite(localization_fdr) & localization_fdr < 0.05]),
    all(locked$current_passed),
    all(data.table::fread(file.path(
      DIR_CBM_RESULTS, "cbm_upgrade_validation.csv"
    ))$passed)
  )
)

data.table::fwrite(validation, file.path(
  DIR_CBM_RESULTS, "review_sensitivity_validation.csv"
))
if (!all(validation$passed)) {
  print(validation)
  stop("审稿意见定向升级验收失败。")
}
print(validation)
message("审稿意见定向升级验收全部通过。")
