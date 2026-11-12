# Old vs. new ω_c and ω_t for the two Section 4.3 datasets

Comparison of the current manuscript IMV code (`poly_compare_multi_compute.R`) with the new `imv_c()` and `imv_t()` on the `polytomous-imv` branch of `ben-domingue/imv`.

**The new formulas barely change Section 4.3 if test responses stay pooled across items as they are now.** Averaging by item, which is what the new `imv()` mirt method does, moves the numbers more. It matters mainly for LSHS, where ω_c drops by about 0.026 and ω_t rises by about 0.037.

"The two tables" was taken to mean the two IRW datasets in Section 4.3, since that section has no LaTeX tables. All numbers are Binom 2PL → GPCM, with old and new computed on identical predictions.

## Average over 20 random 80/20 splits

SD across splits in parentheses.

| Dataset | Metric | Old | New, pooled | New, by-item mean |
|---|---|---|---|---|
| `assessment_time_fournier_2025_bscs` | ω_c | −0.0026 (0.004) | −0.0025 (0.004) | −0.0061 (0.005) |
| | ω_t | −0.0010 (0.001) | −0.0014 (0.002) | −0.0026 (0.003) |
| `LSHS-E_Quidaja_2022` | ω_c | 0.204 (0.008) | 0.204 (0.008) | 0.178 (0.008) |
| | ω_t | 0.105 (0.011) | 0.107 (0.012) | 0.142 (0.014) |

## Using the seed from `fig_empirical_illustrative.R` (20260717)

| Dataset | Metric | Paper | Old, rerun | New, pooled | New, by-item mean |
|---|---|---|---|---|---|
| `assessment_time_fournier_2025_bscs` | ω_c | −0.009 | 0.0004 | 0.0001 | −0.0079 |
| | ω_t | −0.004 | −0.0003 | −0.0005 | −0.0055 |
| `LSHS-E_Quidaja_2022` | ω_c | 0.209 | 0.2086 | 0.2085 | 0.1795 |
| | ω_t | 0.111 | 0.1114 | 0.1134 | 0.1443 |

## Notes

- **Why pooled barely moves:** from the code, the old and new versions score the same binary tasks with the same `imv.binary()`, so only the task weights differ. The largest shift is LSHS ω_t, up about 0.002, and it goes the same direction in all 20 splits.
- **Why by-item moves more:** LSHS gets a lower ω_c and a higher ω_t, and both shifts have the same sign in all 20 splits. So the aggregation choice matters more than the new weights.
- **Paper value not reproduced:** the LSHS values match, but the Fournier values (−0.009 and −0.004) don't. The saved pipeline run has −0.0056 and −0.0022, and the split-to-split SD is about 0.004. The sign of that "ties GPCM" number is therefore within holdout noise, and the 20-split average may be a better thing to report.
- **Bug in the branch:** `.omega_t` in `R/imv_poly.R` overflows R integers at larger N, so `imv_t()` crashed on LSHS with its roughly 8,500 pooled responses. The weight line multiplies integer counts; making `Nk` double fixes it. The fix was applied only inside the diagnostic script, not in the `imv` repo. Single-item calls would hit it too, at a few thousand responses.
- **Fair-coin floor:** with item-level tasks, some binary tasks score worse than a fair coin and get floored at w = 0.5, which sets that task's IMV to 0 when both models are floored. That happens less often once responses are pooled, and it may contribute to the by-item differences.

## Open question

Should the paper use pooled or by-item aggregation? The recommendation is to keep pooled for comparability with the current results, unless the paper should match what the package reports. By-item would change the LSHS numbers and likely the full 436-dataset summaries. Rerunning those requires refitting because the saved results don't include predicted probabilities.

## Code

`diag_imv_new_vs_old.R`, in this folder. Usage: `Rscript diag_imv_new_vs_old.R <imv polytomous-imv checkout> <outdir>`.
