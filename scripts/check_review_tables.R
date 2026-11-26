# Run after the full pipeline in an extracted Supporting Code S1 folder.
# Preserve both sets of files and report which reference table differs, if any.
expected <- list.files("reference", pattern = "\\.csv$", recursive = TRUE)
expected <- expected[!startsWith(expected, "software/")]
if (!length(expected)) stop("Run from the extracted S1 directory containing reference/.")
checks <- lapply(expected, function(name) {
  path <- file.path("results", name)
  if (!file.exists(path)) return(paste("missing", name))
  read_table <- function(file) read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  a <- read_table(file.path("reference", name))
  b <- read_table(path)
  difference <- all.equal(a, b, tolerance = 1e-7)
  if (isTRUE(difference)) return(NULL)
  paste(name, paste(difference, collapse = ", "))
})
failures <- unlist(checks)
if (length(failures)) stop(paste(failures, collapse = "\n"))
message("Matched ", length(expected), " reference tables at tolerance 1e-7.")
