# 这个文件只放路径、随机种子和数据网址；分析逻辑不放在这里。

# 要求从项目根目录运行。双击 .Rproj 后，RStudio 会自动满足这一点。
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# 如果工作目录不对，立即停止，避免把数据下载到莫名其妙的位置。
if (!file.exists(file.path(PROJECT_ROOT, "UC_MucosalRecovery_Q2_10IF_project.Rproj"))) {
  stop("工作目录不对：请先双击 UC_MucosalRecovery_Q2_10IF_project.Rproj，再运行脚本。")
}

# 统一定义所有目录，后续脚本只引用这些变量。
DIR_RAW <- file.path(PROJECT_ROOT, "data", "raw")
DIR_INTERMEDIATE <- file.path(PROJECT_ROOT, "data", "intermediate")
DIR_RESULTS <- file.path(PROJECT_ROOT, "results")
DIR_FIGURES <- file.path(PROJECT_ROOT, "figures")
DIR_LOGS <- file.path(PROJECT_ROOT, "logs")

# 自动创建不存在的目录；recursive=TRUE 表示父目录不存在时一并创建。
for (d in c(DIR_RAW, DIR_INTERMEDIATE, DIR_RESULTS, DIR_FIGURES, DIR_LOGS)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# 固定随机种子，使 bootstrap 和抽样结果可复现。
SEED <- 20260826L
set.seed(SEED)

# 核心 GEO 处理后矩阵。总下载量约 153 MB，不下载 CEL 原始文件。
CORE_FILES <- data.frame(
  id = c("GSE206285", "GSE92415", "GSE16879", "GSE12251", "GSE73661"),
  platform = c("GPL13158", "GPL13158", "GPL570", "GPL570", "GPL6244"),
  url = c(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206285/matrix/GSE206285_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92415/matrix/GSE92415_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE16nnn/GSE16879/matrix/GSE16879_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE12nnn/GSE12251/matrix/GSE12251_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE73nnn/GSE73661/matrix/GSE73661_series_matrix.txt.gz"
  ),
  stringsAsFactors = FALSE
)

# 三个平台的 GEO 官方探针注释，避免依赖芯片注释包版本。
ANNOTATION_FILES <- data.frame(
  id = c("GPL13158", "GPL570", "GPL6244"),
  url = c(
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13158/annot/GPL13158.annot.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz"
  ),
  stringsAsFactors = FALSE
)

# UNIFI 试验额外临床变量：CRP、粪钙卫蛋白、Mayo 评分和既往抗 TNF 暴露。
UNIFI_METADATA_URL <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206285/suppl/",
  "GSE206285_additional_sample_metadata.tsv.gz"
)

# 第二阶段单细胞定位文件。核心结果通过后才下载，避免第一次运行就占用资源。
SC_RDS_URL <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250487/suppl/",
  "GSE250487_biopsy_RNA.rds.gz"
)

# 高影响力升级阶段的额外公共数据。GSE23597 是 ACT-1 英夫利昔单抗/安慰剂队列；
# GSE72819 是 EUCALYPTUS 依特珠单抗 RNA-seq 队列。只下载处理后文件，不下载 CEL/FASTQ。
BENCHMARK_FILES <- data.frame(
  id = c("GSE23597_matrix", "GSE72819_metadata", "GSE72819_counts"),
  filename = c(
    "GSE23597_series_matrix.txt.gz",
    "GSE72819_series_matrix.txt.gz",
    "GSE72819_RAW.tar"
  ),
  url = c(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE23nnn/GSE23597/matrix/GSE23597_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72819/matrix/GSE72819_series_matrix.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72819/suppl/GSE72819_RAW.tar"
  ),
  min_bytes = c(20000000L, 5000L, 20000000L),
  stringsAsFactors = FALSE
)

# 匹配随机基因集次数。10000 是主分析值；若只想调试代码，可临时改为 100。
N_MATCHED_RANDOM <- 10000L
