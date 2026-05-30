# Ageing

R code and example inputs for the GBD 2021 equivalent-age analysis.

## Workflow

1. `scripts/01_prepare_gbd_inputs.R` checks processed GBD input files.
2. `scripts/02_define_age_related_causes.R` screens age-related causes and aggregates age-related DALY burden.
3. `scripts/03_equivalent_age.R` computes overall and disease-system-specific equivalent ages against the 2021 Global Both-sex age-65 benchmark.
4. `scripts/04_projection_nordpred.R` projects age-specific DALY rates, derives population-weighted aggregate regions, and maps projected burden to equivalent ages.
5. `scripts/05_aapc.R` contains the AAPC analysis.
6. `scripts/06_decomposition_shapley.R` performs Shapley decomposition of the global equivalent-age change.

## Requirements

The analysis was prepared for R 4.3.1. The scripts use:

`dplyr`, `readr`, `tidyr`, `stringr`, `tibble`, `data.table`, and `reshape2`.

`scripts/05_aapc.R` also requires the NCI Joinpoint command-line executable (`jpCommand.exe`). Set `JOINPOINT_CMD` to the full path of `jpCommand.exe` before running the full analysis.

## Installation

Install R and the packages above, then set the working directory to the repository root:

```r
install.packages(c("dplyr", "readr", "tidyr", "stringr", "tibble", "data.table", "reshape2"))
```

## Demo

The demo runs on the included `data/demo/` inputs and writes CSV outputs to `outputs/demo/`.

```bash
Rscript scripts/run_demo.R
```

Reference outputs are provided in `outputs/expected_demo/`.

## Full Analysis

Add the processed GBD disease-burden inputs described in `data/README.md` to `data/processed/`, then run:

```bash
Rscript scripts/01_prepare_gbd_inputs.R
Rscript scripts/02_define_age_related_causes.R
Rscript scripts/03_equivalent_age.R
Rscript scripts/04_projection_nordpred.R
Rscript scripts/05_aapc.R
Rscript scripts/06_decomposition_shapley.R
```

The processed historical and projected population inputs required by `scripts/04_projection_nordpred.R` are included in `data/processed/`. See `data/README.md` for file names, columns and age ranges.

## License

This repository is released under the MIT License. See `LICENSE`.
