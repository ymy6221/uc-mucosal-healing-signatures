# From the project root: Rscript --vanilla scripts/audit_recorded_versions.R
# Compare the historical environment and pathway registry with today's local installation.
historical_dir <- file.path("results", "cbm_upgrade")
software_dir <- if (dir.exists("reference/software")) "reference/software" else historical_dir
audit_dir <- file.path("results", "version_audit")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
historical <- read.csv(file.path(software_dir, "package_versions_cbm_upgrade.csv"),
                       colClasses = "character")
historical$current <- vapply(historical$package, function(pkg) {
  if (pkg == "R") return(as.character(getRversion()))
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(packageVersion(pkg))
}, character(1))
historical$identical_version <- historical$version == historical$current
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "sessionInfo_current.txt"))

registry_path <- file.path(historical_dir, "pathway_registry_msigdb.csv")
registry <- read.csv(registry_path, colClasses = "character")
original <- unique(registry[c("pathway", "gene", "database_version")])
name_list <- unique(original$pathway)
get_pathway_members <- function(collection, subcollection = NULL) {
  db <- msigdbr::msigdbr(species = "Homo sapiens", collection = collection,
                         subcollection = subcollection)
  db <- db[db$gs_name %in% name_list, c("gs_name", "gene_symbol", "db_version")]
  names(db) <- names(original)
  unique(as.data.frame(db))
}
observed <- rbind(get_pathway_members("H"), get_pathway_members("C2", "CP:KEGG_LEGACY"))
ordered_rows <- function(x) {
  x <- x[order(x$pathway, x$gene, x$database_version), , drop = FALSE]
  rownames(x) <- NULL
  x
}
matches <- identical(ordered_rows(observed), ordered_rows(original))
registry_sha256 <- digest::digest(file = registry_path, algo = "sha256")
registry_matches_record <- registry_sha256 ==
  "b35f4806af3fa50939bd43281ceef73da598c10365b2c90fa5471360938bccae"
report <- c(
  paste("Audit date", Sys.Date()),
  paste("Historical R session", readLines(file.path(software_dir, "sessionInfo_cbm_upgrade.txt"), n = 1L)),
  paste("Current R session", R.version.string),
  "",
  capture.output(print(historical, row.names = FALSE)),
  "",
  paste("Historical database_version", paste(unique(original$database_version), collapse = ", ")),
  paste("Current db_version", paste(unique(observed$database_version), collapse = ", ")),
  paste("Pathways", length(name_list)),
  paste("Historical unique pathway/gene/version rows", nrow(original)),
  paste("Current unique pathway/gene/version rows", nrow(observed)),
  paste("Exact equality of all pathway/gene/version rows after sorting", matches),
  paste("Registry SHA-256", registry_sha256),
  paste("Registry matches original study file SHA-256", registry_matches_record)
)
writeLines(report, file.path(audit_dir, "version_comparison.txt"))
cat(report, sep = "\n")
if (!all(historical$identical_version %in% TRUE) || !matches || !registry_matches_record) {
  stop("The current installation differs from the recorded analysis environment. See the audit report.")
}
