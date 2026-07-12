# Binomial Approximations to Polytomous IRT Models

This repository contains the code for the paper *"Binomial Approximations to
Polytomous Item Response Models"* by Ben Domingue et al. The paper examines
when and why approximating the GPCM and GRM with a binomial response model
succeeds or fails, focusing on the role of threshold structure, link function
choice, and the shape of the expected response function (ERF).

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
  npspline.R                     # nonparametric spline figure
  fig_single_item_example.R      # single-item illustrative figure
  fig_np_boundary_mixture.R      # nonparametric boundary mixture figure
  fig_np_extreme_categories.R    # extreme-category behavior figure (SI)
  make_si_dataset_table.R        # supplementary dataset table
  sim/                           # simulation code and precomputed results
    gpcm_binomial_report.qmd     # full simulation report (renders to HTML)
    gpcm_binomial_sweep.R        # GPCM data-generating-process simulation
    gpcm_binomial_sweep_np.R     # nonparametric data-generating-process simulation
    *.rds                        # simulation output, read by make_figures.R
  irw/                           # empirical application code
    poly_compare.qmd
    poly_compare_multi_compute.R # fits models to real datasets
    poly_compare_data/           # per-dataset results (not tracked in git; ask author)
  docs/                          # GitHub Pages source for the interactive supplement
```

Each script writes figures to a configurable output directory (an `OUT`
constant near the top of the script) and reads simulation/empirical data from
`sim/` and `irw/poly_compare_data/` respectively via similar path constants.
Update these constants for your own environment before running.

## Reproducing the figures

Steps 1–3 are computationally expensive and only need to be rerun if
simulation parameters or the empirical dataset collection change. Their
outputs are precomputed and included (`sim/*.rds`) or available from the
author on request (`irw/poly_compare_data/`, gitignored due to size). Steps
4–11 are fast, can run in any order, and all should be run from this `src/`
directory.

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
| 9 | `Rscript fig_np_boundary_mixture.R` | nonparametric boundary mixture figure |
| 10 | `Rscript fig_np_extreme_categories.R` | extreme-category behavior figure (SI) |
| 11 | `Rscript make_si_dataset_table.R` | supplementary dataset table |

## Figure inventory

### Self-contained scripts (no data dependencies)

| Script | Description |
|--------|-------------|
| `isoerf_distance2.R` | isoERF distance landscape (Figure 1) |
| `isoerf_distance2_altcontour.R` | Figure 1, alternate-contour variant |
| `npspline.R` | KL-divergence-by-category figure |
| `fig_single_item_example.R` | Single-item illustrative example |
| `fig_np_boundary_mixture.R` | Nonparametric boundary mixture illustration |

### `make_si_dataset_table.R` (reads `irw/poly_compare_data/`)

Generates the supplementary longtable summarizing each empirical dataset used.

### `fig_np_extreme_categories.R` (self-contained; simulates its own draws)

Supplementary figure showing extreme- vs. interior-category response share as
the nonparametric data-generating process becomes more flexible.

### `make_figures.R` — simulation figures

Reads `sim/gpcm_binomial_results.rds` (GPCM data-generating process):
ERF overview and grid, IMV under the GPCM DGP, IMV floor behavior, and
simulation RMSE.

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
models, GRM vs. GPCM comparison, and an asymmetry analysis.

## `scratch/`

Diagnostic scripts used to check specific patterns in the results during
analysis. None produce a manuscript figure or table — their conclusions are
already reflected in the paper's text. Kept for transparency, not required
for reproduction.
