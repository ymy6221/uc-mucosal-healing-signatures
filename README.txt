Supporting Code S1
Public-data benchmarking of mucosal transcriptomic signatures in ulcerative colitis

Inspect the tables

The reference directory contains aggregate CSV files that can be opened in a
spreadsheet program. Individual participant expression data are available from
the source GEO accessions. The R scripts regenerate results under results/ and
leave the reference files intact.

  review_repeated_cv/      Repeated, outcome-stratified internal validation
                          and full-cohort likelihood-ratio tests in UNIFI
  review_telesco/          Telesco proxy scores by arm, pooled scores and the
                          treatment interaction in PURSUIT-SC
  review_pool_sensitivity/ Anti-TNF interactions under inverse-variance common-
                          effect and DerSimonian–Laird random-effects models
  review_clinical_trials/  Trial interactions, complete-case selection, trial
                          weights and drug-class contrasts
  review_clinical_auc/     Clinical residual models, apparent AUC and the
                          original single-split analysis in UNIFI
  overlap_reanalysis/      Arm selection and treatment-outcome analyses
  overlap_audit/           Summaries of overlap across public accessions
  overlap_pathways/        Principal component and pathway analyses
  cbm_upgrade/             Pathway overlap and single-cell results
  software/                Historical and checked software versions

Telesco et al. published the panel genes and classifier thresholds
(doi:10.1053/j.gastro.2018.06.077). The historical screening-log entry about
a missing decision threshold was incorrect. The signature registry's historical
description of the original classifier as unavailable records a limitation of
the earlier implementation. It gives no basis for concluding that the authors
withheld fitted parameters or thresholds. This sensitivity analysis uses the
public 13-gene panel to calculate an unweighted rank score on the deposited
expression data. Its estimates describe that proxy, which uses a different
scoring rule from the published trained classifier.

Full rerun in an empty folder

1. Extract the entire ZIP to a writable directory. In RStudio, open the
   .Rproj file, then the Terminal tab. Run the commands below from that folder.
   Rscript must be available. Check reference/software/ for the original R
   and package versions. The installers obtain packages from CRAN and
   Bioconductor. If current packages differ, numerical results can vary.

   Rscript --vanilla scripts/00_install_packages.R
   Rscript --vanilla scripts/21_install_cbm_upgrade_packages.R
   Rscript --vanilla -e 'for (p in c("R.utils","ragg","pROC","matrixStats")) if (!requireNamespace(p,quietly=TRUE)) install.packages(p,repos="https://cloud.r-project.org")'

2. Start the entire workflow and compare its output with the reference tables.

   Rscript --vanilla run_validation_review.R
   Rscript --vanilla scripts/check_review_tables.R

   The first command runs the historical stages, the overlap review, the
   clinical revision, and scripts 30–32. The second compares all archived
   reference CSVs except the software inventory with newly generated results,
   at a numerical tolerance of 1e-7. Python is unnecessary for these analyses.

3. GEO supplies 12 bulk data, metadata and annotation files (241,537,819
   bytes) in data/raw/ and one single-cell file (408,316,301 bytes),
   GSE250487_biopsy_RNA.rds.gz, in data/singlecell/. Existing copies can be
   placed at these paths before running. See data/public_download_manifest.csv
   for the bulk inputs and config.R for the single-cell URL. Check the byte
   counts and SHA-256 hashes of all 13 files in data/expected_input_checksums.csv.
   The pathway stage obtains MSigDB from the msigdbr package cache or, on
   first use, the network. See reference/software/ for the database release.
   Set aside disk space for the analysis products and installed packages.
   The single-cell object alone requires about 2 GB of RAM when loaded, with
   additional memory used for aggregation and clustering.

Rerun just this validation revision

If the entire preceding trial-review workflow has already run in the same
working directory, including its intermediate data and model outputs, use:

   Rscript --vanilla run_validation_review.R validation
   Rscript --vanilla scripts/check_review_tables.R

Script 30 evaluates the fixed clinical predictor set with and without each
additional signature in the 330 ustekinumab-treated complete cases. Script 31
calculates the Telesco proxy on the overlap-reviewed arms. Script 32 estimates
the common-effect anti-TNF sensitivity on the trial-level estimates in
results/review_clinical_trials/trial_models.csv. These stages do not require
downloading the GEO expression files again.

Archive integrity

On macOS or Linux, run this command in the extracted directory to check each
distributed file against the package manifest:

   shasum -a 256 -c SHA256SUMS.txt
