# 下载五个处理后表达矩阵、三个官方注释和一个临床变量表。

check_packages()

# 逐个下载 Series Matrix；50 KB 是防止错误网页被误认为数据文件的下限。
for (i in seq_len(nrow(CORE_FILES))) {
  destination <- file.path(DIR_RAW, paste0(CORE_FILES$id[i], "_series_matrix.txt.gz"))
  download_with_retry(CORE_FILES$url[i], destination, min_bytes = 50000L)
}

# 逐个下载平台注释；每个官方注释都应大于 100 KB。
for (i in seq_len(nrow(ANNOTATION_FILES))) {
  destination <- file.path(DIR_RAW, paste0(ANNOTATION_FILES$id[i], ".annot.gz"))
  download_with_retry(ANNOTATION_FILES$url[i], destination, min_bytes = 100000L)
}

# 下载 UNIFI 附加临床变量表。
download_with_retry(
  UNIFI_METADATA_URL,
  file.path(DIR_RAW, "GSE206285_additional_sample_metadata.tsv.gz"),
  min_bytes = 1000L
)

# 写出文件清单，便于检查数据是否完整。
downloaded <- data.table::data.table(
  file = list.files(DIR_RAW, full.names = TRUE),
  bytes = file.info(list.files(DIR_RAW, full.names = TRUE))$size
)
data.table::fwrite(downloaded, file.path(DIR_RESULTS, "downloaded_files.csv"))
print(downloaded)

