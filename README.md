# The Binomial Approximation for Polytomous IRT Models

This repository contains the code for the paper *"The role of threshold
structure and the link function in the performance of the binomial
approximation for modeling polytomous responses"* by Benjamin W. Domingue,
Yiqing Liu, Sanford Student, and Klint Kanopka. The paper examines when and
why approximating the GPCM and GRM with a binomial response model succeeds or
fails, decomposing the approximation error into misfit in the link function
and misfit in threshold structure.

## Interactive supplement

An interactive companion site is available at
**https://ben-domingue.github.io/binomial-model/**. It lets you explore the
models directly in the browser — adjusting item parameters and watching how
category response functions (CRFs) and expected response functions (ERFs)
change — without installing anything. Source for the site lives in
[`docs/`](docs/); see that folder for details on adding new widgets.

## Repository layout

```
src/
  make_figures.R                 # generates most simulation + empirical figures
  isoerf_distance2.R             # isoERF distance landscape figure
  isoerf_distance2_altcontour.R  # alternate-contour variant
  isoerf_neighborhood_vs_a.R     # SI: near-isoERF neighborhood size vs. discrimination
  npspline.R                     # nonparametric spline figure
  fig_single_item_example.R      # single-item illustrative figure
  fig_np_boundary_mixture.R      # nonparametric boundary mixture figure
  fig_np_extreme_categories.R    # extreme-category behavior figure (SI)
  make_si_dataset_table.R        # supplementary dataset table
  fig_threshold_stripchart.R     # §4.3 two illustrative datasets (refits from IRW)
  fig_empirical_illustrative.R   # earlier §4.3 CRF/ERF version (not in manuscript)
  irw_fetch_cached.R             # helper: irw::irw_fetch() with a /tmp cache
  sim/                           # simulation code and precomputed results
    gpcm_binomial_report.qmd     # full simulation report (renders to HTML)
    gpcm_binomial_sweep.R        # GPCM data-generating-process simulation
    gpcm_binomial_sweep_np.R     # nonparametric data-generating-process simulation
    *.rds                        # simulation output, read by make_figures.R
  irw/                           # empirical application code
    poly_compare.qmd
    poly_compare_multi_compute.R # fits models to real datasets
    poly_compare_data/           # per-dataset results, multi_<dataset>.rds (one per IRW dataset)
  docs/                          # GitHub Pages source for the interactive supplement
  scratch/                       # diagnostics, not needed for reproduction
```

**Paths.** Output locations are currently hardcoded as absolute paths to the
author's machine — an `OUT` constant near the top of `make_figures.R`,
`make_si_dataset_table.R`, `fig_threshold_stripchart.R` and
`fig_empirical_illustrative.R`, and an absolute path inside the final
`ggsave()` call of the other figure scripts. `make_figures.R` likewise reads
inputs via `SIM`/`VIG` constants. Edit these for your own environment before
running.

## Requirements

R (tested with 4.6.1) and the CRAN packages `ggplot2`, `dplyr`, `tidyr`,
`scales`, `stringr`, `RColorBrewer`, `ggrepel`, `ggpubr`, `patchwork`,
`mgcv`, `scam`, `mirt`, `plotly`, plus `imv` and `irw`. Steps 3 and the two
§4.3 illustrative scripts fetch data from the Item Response Warehouse via
`irw::irw_fetch()`, which requires a (free) Redivis account; everything else
runs offline from the files in this repository.

## Reproducing the figures

Steps 1–3 are computationally expensive and only need to be rerun if
simulation parameters or the empirical dataset collection change. Their
outputs are precomputed and included in the repository (`sim/*.rds`,
`irw/poly_compare_data/multi_*.rds`). Steps 4–14 can run in any order and
should be run from this `src/` directory; all are fast except 13–14, which
refit two empirical datasets.

| Step | Command | Output |
|------|---------|--------|
| 1 | `Rscript sim/gpcm_binomial_sweep.R` | `sim/gpcm_binomial_results.rds` |
| 2 | `Rscript sim/gpcm_binomial_sweep_np.R` | `sim/gpcm_binomial_results_np.rds`, `sim/gpcm_binomial_results_np_oracle.rds` |
| 3 | `Rscript irw/poly_compare_multi_compute.R` | `irw/poly_compare_data/multi_*.rds` |
| 4 | `Rscript make_figures.R` | simulation + empirical figures (see inventory below) |
| 5 | `Rscript npspline.R` | nonparametric spline figure |
| 6 | `Rscript fig_single_item_example.R` | single-item illustrative figure |
| 7 | `Rscript isoerf_distance2.R` | isoERF distance landscape figure |
| 8 | `Rscript isoerf_distance2_altcontour.R` | alt-contour variant of the above |
| 9 | `Rscript isoerf_neighborhood_vs_a.R` | near-isoERF neighborhood size vs. discrimination (SI) |
| 10 | `Rscript fig_np_boundary_mixture.R` | nonparametric boundary mixture figure |
| 11 | `Rscript fig_np_extreme_categories.R` | extreme-category behavior figure (SI) |
| 12 | `Rscript make_si_dataset_table.R` | supplementary dataset table |
| 13 | `Rscript fig_threshold_stripchart.R` | §4.3 illustrative-datasets figure |
| 14 | `Rscript fig_empirical_illustrative.R` | earlier §4.3 figure (optional; not in manuscript) |

## Figure inventory

### Self-contained scripts (no data dependencies)

| Script | Description |
|--------|-------------|
| `isoerf_distance2.R` | isoERF distance landscape (Figure 1) |
| `isoerf_distance2_altcontour.R` | Figure 1, alternate-contour variant |
| `isoerf_neighborhood_vs_a.R` | Near-isoERF neighborhood size vs. discrimination, continuous sweep extending Figure 1 (SI Figure S2) |
| `npspline.R` | KL-divergence-by-category figure |
| `fig_single_item_example.R` | Single-item illustrative example |
| `fig_np_boundary_mixture.R` | Nonparametric boundary mixture illustration |

### `make_si_dataset_table.R` (reads `irw/poly_compare_data/`)

Generates the supplementary longtable summarizing each empirical dataset used.

### `fig_threshold_stripchart.R` (fetches two IRW datasets)

Section 4.3 figure: fitted GPCM thresholds for every item in a good-fit
(`assessment_time_fournier_2025_bscs`) and a poor-fit (`LSHS-E_Quidaja_2022`)
dataset. Refits both from scratch with the same pipeline as step 3 (seeded).
`fig_empirical_illustrative.R` is the earlier CRF/ERF version of this figure,
kept for reference; its header documents why other candidate datasets were
rejected.

### `fig_np_extreme_categories.R` (self-contained; simulates its own draws)

Supplementary figure showing extreme- vs. interior-category response share as
the nonparametric data-generating process becomes more flexible.

### `make_figures.R` — simulation figures

Reads `sim/gpcm_binomial_results.rds` (GPCM data-generating process):
ERF overview and grid, IMV under the GPCM DGP, IMV floor behavior, and
simulation RMSE. Also draws the SI threshold-asymmetry illustration
(`fig_threshold_asym`), which uses only the simulation's threshold function.

Reads `sim/gpcm_binomial_results_np.rds` (nonparametric DGP):
IMV vs. GPCM, implied parameter bounds, and cross-section of the fit surface.

Reads `sim/gpcm_binomial_results_np_oracle.rds` (nonparametric DGP, oracle
reference): IMV vs. the oracle model (true DGP, true θ) — a supplementary
figure.

Also generates a supplementary isoERF injectivity proof figure, which is
self-contained and has no data dependency.

### `make_figures.R` — empirical figures

Reads `irw/poly_compare_data/multi_*.rds`: CV RMSE distribution across
models, mean rank across datasets, mean IMV vs. GPCM, pairwise IMV, scatter
of model pairs, IMV gap by number of categories, IMV for 1PL-constrained
models, GRM vs. GPCM comparison, the distribution of threshold asymmetry
(`fig_asym`), and pairwise ω_c vs. threshold asymmetry (`fig_asym_imv`).

## `scratch/`

Diagnostic scripts used to check specific patterns in the results during
analysis. None produce a manuscript figure or table — their conclusions are
already reflected in the paper's text. Kept for transparency, not required
for reproduction.
