# Reproducing figures — polytomous IRT / binomial approximation paper

This directory contains all code needed to reproduce the manuscript figures.
Figures are written to the Overleaf working directory at
`/home/ben/Dropbox/Apps/Overleaf/binomial_model/`. Three path constants near
the top of `make_figures.R` control all input/output locations:

```r
OUT <- "/home/ben/Dropbox/Apps/Overleaf/binomial_model"   # figure output
SIM <- "/home/ben/Dropbox/projects/binomial_model/src/sim" # simulation data
VIG <- "/home/ben/Dropbox/projects/binomial_model/src/irw" # empirical data
```

## Directory layout

```
src/
  make_figures.R               # generates all simulation + empirical figures
  isoerf_distance2.R           # generates Figure 1 (isoERF distance landscape)
  npspline.R                   # generates fig_kl_category.pdf
  fig_single_item_example.R    # generates fig_single_item_example.pdf
  sim/                         # simulation source code and output data
    gpcm_binomial_report.qmd   # main simulation report (renders to HTML)
    gpcm_binomial_sweep.R      # GPCM DGP simulation runner
    gpcm_binomial_sweep_np.R   # NP DGP simulation runner
    gpcm_binomial_results.rds        # GPCM DGP results (read by make_figures.R)
    gpcm_binomial_results_np.rds     # NP DGP results (read by make_figures.R)
    gpcm_binomial_results_np_oracle.rds
    gpcm_binomial_add_tutz.R
    gpcm_binomial_extend.R
    gpcm_binomial_known_theta.R
    gpcm_binomial_viz.R
  irw/                         # empirical vignette code
    poly_compare.qmd
    poly_compare_multi_compute.R  # runs empirical fits, writes multi_*.rds
    poly_compare_sim.R
    poly_compare_data/            # per-dataset RDS files (read by make_figures.R)
```

## Run order

Steps 1–3 are computationally expensive and only need to be rerun if simulation
parameters or the empirical dataset collection changes. The outputs of steps 1–2
are tracked in git (`sim/*.rds`); the outputs of step 3 (`irw/poly_compare_data/`)
are gitignored due to size but are available from the author. Steps 4–7 are fast
and can run in any order or in parallel.

| Step | Command | Output |
|------|---------|--------|
| 1 | `Rscript sim/gpcm_binomial_sweep.R` | `sim/gpcm_binomial_results.rds` |
| 2 | `Rscript sim/gpcm_binomial_sweep_np.R` | `sim/gpcm_binomial_results_np.rds` |
| 3 | `Rscript irw/poly_compare_multi_compute.R` | `irw/poly_compare_data/multi_*.rds` |
| 4 | `Rscript make_figures.R` | simulation + empirical figures → OUT |
| 5 | `Rscript npspline.R` | `fig_kl_category.pdf` → OUT |
| 6 | `Rscript fig_single_item_example.R` | `fig_single_item_example.pdf` → OUT |
| 7 | `Rscript isoerf_distance2.R` | `fig_isoerf_distance2.pdf` → OUT |

All scripts in steps 4–7 should be run from this `src/` directory.

## Figure inventory

### Self-contained scripts (no data dependencies)

| Script | Output PDF | Manuscript figure |
|--------|-----------|-------------------|
| `isoerf_distance2.R` | `fig_isoerf_distance2.pdf` | Figure 1 — isoERF distance landscape |
| `npspline.R` | `fig_kl_category.pdf` | |
| `fig_single_item_example.R` | `fig_single_item_example.pdf` | |

### `make_figures.R` — simulation figures

Reads `sim/gpcm_binomial_results.rds` (GPCM DGP):

| Output PDF | Description |
|-----------|-------------|
| `fig_erf_overview.pdf` | ERF overview |
| `fig_erf_grid.pdf` | ERF grid |
| `fig_imv_gpcm.pdf` | IMV under GPCM DGP |
| `fig_imv_floor.pdf` | IMV floor behavior |
| `fig_rmse_sim.pdf` | RMSE across simulation conditions |

Reads `sim/gpcm_binomial_results_np.rds` (NP DGP):

| Output PDF | Description |
|-----------|-------------|
| `fig_np_imv.pdf` | IMV under NP DGP |
| `fig_implied_bounds.pdf` | Implied parameter bounds |
| `fig_cross_section.pdf` | Cross-section of fit surface |

### `make_figures.R` — empirical figures

Reads `irw/poly_compare_data/multi_*.rds`:

| Output PDF | Description |
|-----------|-------------|
| `fig_pair_imv.pdf` | Pairwise IMV across datasets |
| `fig_scatter_pairs.pdf` | Scatter of model pairs |
| `fig_gap_by_k.pdf` | IMV gap by number of categories |
| `fig_imv_1pl.pdf` | IMV for 1PL-constrained models |
| `fig_grm_vs_gpcm.pdf` | GRM vs. GPCM comparison |
| `fig_asym.pdf` | Asymmetry analysis |
