/* Adjacent-category (GPCM-type) model with a free link, matched to a target ERF.

   P(X = k | θ) ∝ exp( Σ_{j ≤ k} (η(θ) − τ_j) ),  k = 0..K−1.

   For fixed step offsets τ, E[X] is strictly increasing in η, so for any target
   expected score m in (0, K−1) there is exactly one η giving E[X] = m. Solving
   that pointwise in θ yields the unique member of this family with step
   structure τ and the target ERF. Used by the fixed-ERF threshold explorer and
   by the CRF editor's Randomize button. Works for any K (= τ.length + 1). */

function acProbs(eta, delta) {
  const K = delta.length + 1, ln = [0];
  let c = 0;
  for (let k = 0; k < K - 1; k++) { c += eta - delta[k]; ln.push(c); }
  const mx = Math.max(...ln), ex = ln.map(v => Math.exp(v - mx)), sm = ex.reduce((s, v) => s + v, 0);
  return ex.map(v => v / sm);
}
function acMean(eta, delta) { return acProbs(eta, delta).reduce((s, p, k) => s + k * p, 0); }

/* Bisection on η for E[X | η] = m. */
function acSolveEta(m, delta) {
  const K = delta.length + 1;
  m = Math.min(Math.max(m, 1e-9), K - 1 - 1e-9);
  let lo = -20, hi = 20;
  while (acMean(lo, delta) > m) lo -= 20;
  while (acMean(hi, delta) < m) hi += 20;
  for (let it = 0; it < 80; it++) {
    const mid = (lo + hi) / 2;
    if (acMean(mid, delta) < m) lo = mid; else hi = mid;
  }
  return (lo + hi) / 2;
}

/* Given a target ERF on a θ grid, return η(θ) and probs[i][k]. */
function acMatchERF(erf, delta) {
  const eta = erf.map(m => acSolveEta(m, delta));
  const probs = eta.map(e => acProbs(e, delta));
  return { eta, probs };
}

/* Step offsets that make this family the binomial: τ_k = log(k / (K − k)). */
function binomDelta(K) { return Array.from({ length: K - 1 }, (_, i) => Math.log((i + 1) / (K - 1 - i))); }

if (typeof module !== 'undefined') module.exports = { acProbs, acMean, acSolveEta, acMatchERF, binomDelta };
