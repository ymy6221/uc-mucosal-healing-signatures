# 核对 12 个主分析输入文件。字节数和 SHA-256 均一致才算下载完整。

if (!requireNamespace("data.table", quietly = TRUE) ||
    !requireNamespace("digest", quietly = TRUE)) {
  stop("缺少 data.table 或 digest；请先运行 scripts/00_install_packages.R")
}

manifest <- data.table::fread(file.path(PROJECT_ROOT, "data", "public_download_manifest.csv"))
manifest[, path := file.path(DIR_RAW, filename)]
manifest[, exists := file.exists(path)]
manifest[, actual_bytes := ifelse(exists, file.info(path)$size, NA_real_)]
manifest[, actual_sha256 := vapply(path, function(one_path) {
  if (!file.exists(one_path)) return(NA_character_)
  digest::digest(file = one_path, algo = "sha256", serialize = FALSE)
}, character(1))]
manifest[, size_passed := exists & actual_bytes == verified_bytes]
manifest[, sha256_passed := exists & actual_sha256 == sha256]
manifest[, passed := size_passed & sha256_passed]

report <- manifest[, .(
  filename, verified_bytes, actual_bytes, size_passed,
  expected_sha256 = sha256, actual_sha256, sha256_passed, passed
)]
data.table::fwrite(report, file.path(DIR_RESULTS, "public_download_verification.csv"))
print(report[, .(filename, actual_bytes, size_passed, sha256_passed, passed)])

if (!all(report$passed)) {
  stop("公共数据完整性核对失败；删除失败文件后重新运行下载脚本。")
}
message("12/12 个公共数据输入的字节数和 SHA-256 全部一致。")
