# Manifest

## Repository Documentation

- `README.md`: overview, requirements and run instructions.
- `data/README.md`: GBD data sources and processed-input requirements.
- `LICENSE`: MIT License for the analysis code.

## Code

- `scripts/01_prepare_gbd_inputs.R`: checks processed GBD input files.
- `scripts/02_define_age_related_causes.R`: screens age-related causes and aggregates DALY burden.
- `scripts/03_equivalent_age.R`: computes overall and disease-system-specific equivalent age using the 2021 Global Both-sex age-65 benchmark.
- `scripts/04_projection_nordpred.R`: runs the Nordpred power-5 age-period-cohort projection analysis, derives population-weighted aggregate regions, and maps projected burden to equivalent age.
- `scripts/05_aapc.R`: AAPC analysis using the Joinpoint command-line executable.
- `scripts/06_decomposition_shapley.R`: Shapley decomposition of the global equivalent-age change from the aggregated historical burden output.
- `scripts/create_demo_data.R` and `scripts/run_demo.R`: demo workflow.

## Data Included

- `data/demo/`: compact inputs for the demo workflow.
- `data/processed/population_history_2021.csv.gz`: processed historical population input for projection.
- `data/processed/population_projection_who.csv.gz`: processed projected population input for projection.
- `outputs/expected_demo/`: expected demo outputs.
- `data/metadata/std_GBD2021.csv`: standard population table used by projection helper.

## Data Not Included

- Full GBD 2021 raw downloads.
