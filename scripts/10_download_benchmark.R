# 下载升级阶段的两个新公共队列。已存在且大小达标的文件会自动跳过。

check_benchmark_packages()

for (i in seq_len(nrow(BENCHMARK_FILES))) {
  destination <- file.path(DIR_RAW, BENCHMARK_FILES$filename[i])
  download_with_retry(
    BENCHMARK_FILES$url[i],
    destination,
    min_bytes = BENCHMARK_FILES$min_bytes[i]
  )
}

# 通过文件大小和 gzip/tar 结构检查阻止 HTML 错误页被当成数据。
paths <- file.path(DIR_RAW, BENCHMARK_FILES$filename)
validation <- data.table::data.table(
  file = BENCHMARK_FILES$filename,
  bytes = file.info(paths)$size,
  minimum_bytes = BENCHMARK_FILES$min_bytes,
  size_passed = file.info(paths)$size >= BENCHMARK_FILES$min_bytes
)
tar_members <- utils::untar(file.path(DIR_RAW, "GSE72819_RAW.tar"), list = TRUE)
validation[file == "GSE72819_RAW.tar", structure_passed :=
             length(tar_members) == 73L && all(grepl("\\.tab\\.gz$", tar_members))]
validation[file != "GSE72819_RAW.tar", structure_passed := TRUE]

data.table::fwrite(validation, file.path(DIR_RESULTS, "benchmark_download_validation.csv"))
print(validation)
if (!all(validation$size_passed & validation$structure_passed)) {
  stop("新队列下载完整性检查失败。")
}

