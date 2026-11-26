# Mucosal transcriptomic signatures in ulcerative colitis

This repository accompanies Meiying Yan’s analysis of published mucosal transcriptomic signatures and induction treatment outcomes in ulcerative colitis. It contains the R scripts, public input manifest and aggregate reference results used in the analysis.

## Reproduce the analysis

See [README.txt](README.txt) for software requirements, public data sources, memory requirements and step-by-step instructions. After installing the dependencies listed there, run these commands from the project directory:

```sh
Rscript --vanilla run_validation_review.R
Rscript --vanilla scripts/check_review_tables.R
```

## Files

- [scripts/](scripts/) contains the R analysis and verification scripts.
- [data/](data/) contains the download manifest and expected checksums for the public input files.
- [signatures/](signatures/) contains the signature definitions and study registry files.
- [reference/](reference/) contains the aggregate reference tables and software inventory.

The downloadable `Supporting_Code_S1.zip` under **Releases → v1.0.0 → Assets** is the same archive supplied as Supporting Code S1 with the manuscript.
