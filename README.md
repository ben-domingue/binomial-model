# Reproducing figures — polytomous IRT / binomial approximation paper

This directory contains **all** code needed to reproduce the manuscript figures.
As of 2026-07-12 the Overleaf working directory
(`/home/ben/Dropbox/Apps/Overleaf/binomial_model/`) holds no R source at all —
it's LaTeX + bibliography + figure/table output only. Every script here writes
its output there via an absolute path, so all scripts should be run from this
`src/` directory (or by path — none depend on the caller's working directory).
Three path constants near the top of `make_figures.R` control its input/output
locations:

```r
OUT <- "/home/ben/Dropbox/Apps/Overleaf/binomial_model"   # figure output
SIM <- "/home/ben/Dropbox/projects/binomial_model/src/sim" # simulation data
VIG <- "/home/ben/Dropbox/projects/binomial_model/src/irw" # empirical data
```

`fig_np_boundary_mixture.R`, `fig_np_extreme_categories.R`, and
`make_si_dataset_table.R` hardcode the same Overleaf output path and the
`irw/poly_compare_data` input path directly (no shared `OUT`/`VIG` constants),
since they're run standalone rather than sourced by `make_figures.R`.

## Directory layout

```
src/
  make_figures.R               # generates most simulation + empirical figures
  isoerf_distance2.R           # generates Figure 1 (isoERF distance landscape)
  isoerf_distance2_altcontour.R # alternate-contour variant of Figure 1
  npspline.R                   # generates fig_kl_category.pdf
  fig_single_item_example.R    # generates fig_single_item_example.pdf
  fig_np_boundary_mixture.R    # generates fig_np_boundary_mixture.pdf
  fig_np_extreme_categories.R  # generates fig_np_extreme_categories.pdf (SI)
  make_si_dataset_table.R      # generates si_dataset_table.tex (SI table)
  sim/                         # simulation source code and output data
    gpcm_binomial_report.qmd   # main simulation report (renders to HTML)
    gpcm_binomial_sweep.R      # GPCM DGP simulation runner
    gpcm_binomial_sweep_np.R   # NP DGP simulation runner
    gpcm_binomial_results.rds        # GPCM DGP results (read by make_figures.R)
    gpcm_binomial_results_np.rds     # NP DGP results (read by make_figures.R)
    gpcm_binomial_results_np_oracle.rds  # NP DGP oracle results (read by make_figures.R)
    gpcm_binomial_add_tutz.R
    gpcm_binomial_extend.R
    gpcm_binomial_known_theta.R
    gpcm_binomial_viz.R
  irw/                         # empirical vignette code
    poly_compare.qmd
    poly_compare_multi_compute.R  # runs empirical fits, writes multi_*.rds
    poly_compare_sim.R
    poly_compare_data/            # per-dataset RDS files (read by make_figures.R,
                                   # make_si_dataset_table.R)
```

## Run order

Steps 1–3 are computationally expensive and only need to be rerun if simulation
parameters or the empirical dataset collection changes. The outputs of steps 1–2
are tracked in git (`sim/*.rds`); the outputs of step 3 (`irw/poly_compare_data/`)
are gitignored due to size but are available from the author. Steps 4–10 are fast
and can run in any order or in parallel; all should be run from this `src/`
directory.

| Step | Command | Output |
|------|---------|--------|
| 1 | `Rscript sim/gpcm_binomial_sweep.R` | `sim/gpcm_binomial_results.rds` |
| 2 | `Rscript sim/gpcm_binomial_sweep_np.R` | `sim/gpcm_binomial_results_np.rds`, `sim/gpcm_binomial_results_np_oracle.rds` |
| 3 | `Rscript irw/poly_compare_multi_compute.R` | `irw/poly_compare_data/multi_*.rds` |
| 4 | `Rscript make_figures.R` | simulation + empirical figures → OUT (see inventory below) |
| 5 | `Rscript npspline.R` | `fig_kl_category.pdf` → OUT |
| 6 | `Rscript fig_single_item_example.R` | `fig_single_item_example.pdf` → OUT |
| 7 | `Rscript isoerf_distance2.R` | `fig_isoerf_distance2.pdf` → OUT |
| 8 | `Rscript isoerf_distance2_altcontour.R` | `fig_isoerf_distance2_altcontour.pdf` → OUT |
| 9 | `Rscript fig_np_boundary_mixture.R` | `fig_np_boundary_mixture.pdf` → OUT |
| 10 | `Rscript fig_np_extreme_categories.R` | `fig_np_extreme_categories.pdf` → OUT |
| 11 | `Rscript make_si_dataset_table.R` | `si_dataset_table.tex` → OUT |

## Figure inventory

### Self-contained scripts (no data dependencies)

| Script | Output PDF | Manuscript figure |
|--------|-----------|-------------------|
| `isoerf_distance2.R` | `fig_isoerf_distance2.pdf` | Figure 1 — isoERF distance landscape |
| `isoerf_distance2_altcontour.R` | `fig_isoerf_distance2_altcontour.pdf` | Figure 1 (alt-contour variant) |
| `npspline.R` | `fig_kl_category.pdf` | |
| `fig_single_item_example.R` | `fig_single_item_example.pdf` | |
| `fig_np_boundary_mixture.R` | `fig_np_boundary_mixture.pdf` | |

### `make_si_dataset_table.R` (reads `irw/poly_compare_data/`)

| Output | Description |
|--------|-------------|
| `si_dataset_table.tex` | SI longtable, one row per IRW dataset used empirically |

### `fig_np_extreme_categories.R` (reads nothing; simulates its own DGP draws)

| Output PDF | Description |
|-----------|-------------|
| `fig_np_extreme_categories.pdf` | SI figure — extreme- vs. interior-category share as NP-DGP flex increases |

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
| `fig_np_imv.pdf` | IMV under NP DGP, vs. GPCM |
| `fig_implied_bounds.pdf` | Implied parameter bounds |
| `fig_cross_section.pdf` | Cross-section of fit surface |

Reads `sim/gpcm_binomial_results_np_oracle.rds` (NP DGP, oracle reference):

| Output PDF | Description |
|-----------|-------------|
| `fig_np_imv_oracle.pdf` | IMV vs. the oracle (true NP DGP, true θ) — SI Figure, sec:si-oracle. Note: the original generating script for this figure was lost; the block in `make_figures.R` was reconstructed 2026-07-12 from the `_oracle` results columns (`imv_c_oracle`/`imv_t_oracle`) and verified to reproduce the original PDF pixel-for-pixel (same curves, values, and layout). |

Also generates the SI figure `step1_isoerf.pdf`/`.png` (isoERF injectivity proof, step 1) — self-contained, no data dependency.

### `make_figures.R` — empirical figures

Reads `irw/poly_compare_data/multi_*.rds`:

| Output PDF | Description |
|-----------|-------------|
| `fig_rmse_empirical.pdf` | CV RMSE distribution across models |
| `fig_ranks.pdf` | Mean rank across datasets |
| `fig_imv_empirical.pdf` | Mean IMV across datasets, vs. GPCM |
| `fig_pair_imv.pdf` | Pairwise IMV across datasets |
| `fig_scatter_pairs.pdf` | Scatter of model pairs |
| `fig_gap_by_k.pdf` | IMV gap by number of categories |
| `fig_imv_1pl.pdf` | IMV for 1PL-constrained models |
| `fig_grm_vs_gpcm.pdf` | GRM vs. GPCM comparison |
| `fig_asym.pdf` | Asymmetry analysis |

## `scratch/` — diagnostic investigations (not needed to reproduce any figure)

`scratch/` holds one-off scripts written to answer specific "is this pattern real
or an artifact?" questions during analysis. None of them produce a manuscript
figure or table; their conclusions are already reflected in the manuscript prose.
Kept for the reasoning trail behind those claims, not for reproducibility.

| Script | Question it answers |
|--------|---------------------|
| `diag_grid_resolution.R` | Does the 101-point θ grid artificially collapse IMV differences between binomial models at high flex? |
| `diag_threshold_monotonicity.R` | Does the NP DGP produce non-monotone threshold functions P(X>k\|θ) at high flex? |
| `diag_imvt_decomp.R` | Where does the non-monotone ω_t-vs-flex pattern (see `fig_np_imv_oracle.pdf`) come from, decomposed by threshold k? |
| `diag_flex03_crfs.R` | Single-rep sanity check of true-vs-fitted CRFs at flex = 0.3 (writes `diag_flex03_data.rds`) |

## Not included here (superseded drafts, deleted 2026-07-12)

`isoerf.R`, `isoerf_distance.R` (superseded by `isoerf_distance2.R`),
`fig_np_crf_example.R` (its two candidate figures lost out to
`fig_single_item_example.R`), and `adhoc_np_extreme_category.R` (an earlier,
5-panel draft of `fig_np_extreme_categories.R`) were deleted from the Overleaf
working copy as dead drafts fully superseded by files already in this repo.

## Resolved decisions

Fig 7 (`fig_np_imv.pdf`, main text) stays vs. GPCM rather than being replaced
by the oracle version — decided 2026-07-12. Rationale: the main text already
cross-references the SI's oracle-referenced version (`fig_np_imv_oracle.pdf`,
`sec:si-oracle`) to explain the ω_t non-monotonicity artefact, and vs-GPCM
keeps Fig 7 structurally parallel to Study 1's fitted-competitor framing.
`preview_fig7_oracle.R`, the decision-support script for this question, has
been deleted from the Overleaf working copy.
