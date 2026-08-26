# Sampler Platter transfer entropy v2

This is the second, generalized Sampler Platter TE repository. It preserves the
original `Sample_Platter_PQ` repository and follows the current Disturbance
Synchrony organization and plotting conventions.

## What is included

- 25 completed site/pair analyses across Angelo, Atlanta, Baltimore, Hubbard
  Brook, Konza, Starkey, and H.J. Andrews.
- Hubbard Brook W3 and W8 analyses for both the complete record and the April-
  October growing season.
- Dataset-specific timestamp parsing and modal resolution detection. Every
  analysis tests 72 native time steps, while tables and x-axes report hours.
- Atlanta Q -> stream temperature, P -> conductivity, and conductivity ->
  stream temperature for all three sites, in addition to P -> Q. No reverse
  stream-temperature -> Q analysis is included.
- Gap-aware TE windows, so lags do not bridge missing periods or excluded
  seasons.
- Q -> stream-temperature analyses exclude missing-Q and Q=0 observations
  without collapsing their timestamps. Figures shade Q=0 periods across both
  the Q and stream-temperature time-series panels.
- The project's original histogram TE estimator, generalized into shared
  loading, analysis, plotting, and verification functions.
- The latest Disturbance-style diagnostic figure: source and sink series,
  histograms, TE, normalized TE, MI, and correlation with shuffled thresholds
  and detected-lag annotations.

Casper is intentionally pending because the supplied files do not define how
the two precipitation gauges map to the eight discharge watersheds. The exact
choice needed is documented in `00_Data/Casper_data_decision.md`.

## Repository layout

```text
00_Data/       analysis metadata and data-selection decisions
01_Functions/  generalized loaders, TE estimator, and plotting functions
02_Analysis/   audit, run, summary, figure, and verification scripts
03_Reports/    concise cross-site analysis report
04_Results/    TE tables, figures, data audit, and verification manifest
```

Raw data remain in the project-level `Data` folder and are never copied or
modified by this repository.

## Reproduce

Run from the repository root in this order:

```powershell
Rscript 02_Analysis/01_prepare_inputs.R
Rscript 02_Analysis/02_run_all_TE.R
Rscript 02_Analysis/06_compile_summary.R
Rscript 02_Analysis/04_verify_outputs.R
```

To rerun one analysis:

```powershell
Rscript 02_Analysis/03_run_one_TE.R Baltimore_P_Q
Rscript 02_Analysis/06_compile_summary.R
Rscript 02_Analysis/04_verify_outputs.R
```

The local `R_library` directory is ignored by Git. Required packages are
`data.table`, `ggplot2`, `cowplot`, `future`, and `future.apply`.

## Main outputs

- `04_Results/Tables/All_sites_TE_summary.csv`: native resolution, record
  length, dates, pairs, significant peak TE/lag, mean significant normalized TE,
  TE memory, and data notes.
- `04_Results/Tables/All_sites_peak_TE.csv`: compatibility copy of that expanded
  summary.
- `04_Results/Tables/TE_df_<Analysis_ID>.csv`: full lag tables.
- `04_Results/Figures/TE_lag_<Analysis_ID>.png` and `.pdf`: diagnostic figures.
- `03_Reports/TE_analysis_report.html`: complete styled HTML results report
  with links to every lag table and figure.
- `04_Results/verification_manifest.csv`: file and lag-unit verification.

## GitHub handoff

The repository contains no copied raw data, local R packages, RStudio state, or
analysis logs. Results and documentation are intentionally tracked. Review with
`git status`, then create the first commit and push it from your own GitHub
account when ready; this workflow does not commit or push automatically.
