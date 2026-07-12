/* Shared math/chart helpers for the isoERF degeneracy explorer, isoERF manifold
   explorer, and CRF editor pages. GRM/GPCM/Binomial CRF formulas, expected
   response function helpers, and Chart.js plotting utilities used by more than
   one of those pages live here so they aren't duplicated three times. */

const K = 5, TH = Array.from({ length: 120 }, (_, i) => -3.5 + i * 7 / 119);
const CATS = [0, 1, 2, 3, 4];
const COL_A = '#378ADD', COL_B = '#D85A30', COL_REF = '#1D9E75';
const CAT_COLS = ['#378ADD', '#1D9E75', '#BA7517', '#D85A30', '#7F77DD'];
const CAT_NAMES = ['k=0', 'k=1', 'k=2', 'k=3', 'k=4'];
const MF_N = 60, MF_RANGE = [-2.5, 2.5];

function logis(x) { return 1 / (1 + Math.exp(-x)); }
function cumSums(b) { const s = []; let c = 0; b.forEach(v => { c += v; s.push(c); }); return s; }
function allPerms(arr) {
  if (arr.length <= 1) return [arr];
  const res = [];
  arr.forEach((v, i) => {
    const rest = [...arr.slice(0, i), ...arr.slice(i + 1)];
    allPerms(rest).forEach(p => res.push([v, ...p]));
  });
  return res;
}
function rmse(a, b) { return Math.sqrt(a.reduce((t, v, i) => t + (v - b[i]) ** 2, 0) / a.length); }

function gpcmCRFraw(th, a, b) {
  return th.map(t => {
    const ln = [0]; let c = 0;
    for (let k = 0; k < K - 1; k++) { c += a * (t - b[k]); ln.push(c); }
    const mx = Math.max(...ln), ex = ln.map(v => Math.exp(v - mx)), sm = ex.reduce((s, v) => s + v, 0);
    return ex.map(v => v / sm);
  });
}
function grmCRFraw(th, a, b) {
  const bs = [...b].sort((x, y) => x - y);
  return th.map(t => CATS.map(k => {
    const lo = k === 0 ? 1 : logis(a * (t - bs[k - 1])), hi = k === K - 1 ? 0 : logis(a * (t - bs[k]));
    return Math.max(0, lo - hi);
  }));
}
function binomCRFraw(th, a, b) {
  return th.map(t => {
    const p = Math.min(Math.max(logis(a * t - b[0]), 1e-9), 1 - 1e-9);
    return CATS.map(k => {
      let c = 1;
      for (let i = 0; i < k; i++) c *= (K - 1 - i) / (i + 1);
      return c * Math.pow(p, k) * Math.pow(1 - p, K - 1 - k);
    });
  });
}

function computeFull(par, mdl) {
  const { a, b } = par;
  let raw, bound;
  if (mdl === 'grm') {
    raw = grmCRFraw(TH, a, b);
    bound = b.map(bk => TH.map(t => logis(a * (t - bk))));
  } else if (mdl === 'gpcm') {
    raw = gpcmCRFraw(TH, a, b);
    bound = b.map((_, j) => TH.map((_, i) => { let c = 0; for (let k = j + 1; k < K; k++) c += raw[i][k]; return c; }));
  } else {
    raw = binomCRFraw(TH, a, b);
    bound = b.map((_, j) => TH.map((_, i) => { let c = 0; for (let k = j + 1; k < K; k++) c += raw[i][k]; return c; }));
  }
  const crf = CATS.map(k => TH.map((_, i) => raw[i][k]));
  const erf = TH.map((_, i) => CATS.reduce((s, k) => s + k * raw[i][k], 0));
  return { bound, crf, erf };
}
function grmERFfn(a, b) { return TH.map(t => b.reduce((s, bk) => s + logis(a * (t - bk)), 0)); }
function gpcmERFfn(a, b) { const raw = gpcmCRFraw(TH, a, b); return TH.map((_, i) => CATS.reduce((s, k) => s + k * raw[i][k], 0)); }

const charts = {};
function makeChart(id, opts) {
  const c = document.getElementById(id);
  if (charts[id]) charts[id].destroy();
  charts[id] = new Chart(c, { type: 'line', ...opts });
}
const baseOpts = (yLabel, yMax = 1) => ({
  options: {
    animation: false, responsive: true, maintainAspectRatio: false,
    plugins: { legend: { display: false }, tooltip: { enabled: false } },
    scales: {
      x: { type: 'linear', min: -3.5, max: 3.5, ticks: { font: { size: 9 }, color: '#888', callback: v => Number.isInteger(v) ? v : null }, grid: { color: 'rgba(128,128,128,0.1)' }, title: { display: true, text: 'θ', font: { size: 10 }, color: '#888' } },
      y: { min: 0, max: yMax, ticks: { font: { size: 9 }, color: '#888' }, grid: { color: 'rgba(128,128,128,0.1)' }, title: { display: true, text: yLabel, font: { size: 10 }, color: '#888' } }
    }
  }
});
function ds(data, col, dash = []) { return { data: TH.map((t, i) => ({ x: t, y: data[i] })), borderColor: col, borderWidth: 1.5, pointRadius: 0, borderDash: dash }; }

/* Manifold-page-only helpers, but generic enough to share */
function buildGrid(erfFn, refB1, refB2, b3, b4, a) {
  const refERF = erfFn(a, [refB1, refB2, b3, b4]);
  const grid = [];
  for (let iy = 0; iy < MF_N; iy++) {
    const row = [];
    for (let ix = 0; ix < MF_N; ix++) {
      const b1 = MF_RANGE[0] + (MF_RANGE[1] - MF_RANGE[0]) * ix / (MF_N - 1);
      const b2 = MF_RANGE[0] + (MF_RANGE[1] - MF_RANGE[0]) * iy / (MF_N - 1);
      row.push(rmse(refERF, erfFn(a, [b1, b2, b3, b4])));
    }
    grid.push(row);
  }
  return grid;
}

function drawBaseHeatmap(ctx, W, H, grid, admissibleFn) {
  const flat = [].concat(...grid).filter(v => isFinite(v)).sort((a, b) => a - b);
  const maxVal = flat[Math.floor(flat.length * 0.95)] || 1;
  const toVal = (frac) => MF_RANGE[0] + frac * (MF_RANGE[1] - MF_RANGE[0]);
  for (let iy = 0; iy < MF_N; iy++) {
    for (let ix = 0; ix < MF_N; ix++) {
      const b1 = toVal(ix / (MF_N - 1)), b2 = toVal(iy / (MF_N - 1));
      const px = Math.round(ix * W / MF_N), py = Math.round((MF_N - 1 - iy) * H / MF_N);
      const pw = Math.ceil(W / MF_N) + 1, ph = Math.ceil(H / MF_N) + 1;
      if (admissibleFn && !admissibleFn(b1, b2)) { ctx.fillStyle = 'rgba(120,120,120,0.18)'; ctx.fillRect(px, py, pw, ph); continue; }
      const v = Math.min(grid[iy][ix] / maxVal, 1);
      ctx.fillStyle = `rgb(${Math.round(v * 186 + (1 - v) * 8)},${Math.round(v * 117 + (1 - v) * 158)},${Math.round(v * 23 + (1 - v) * 117)})`;
      ctx.fillRect(px, py, pw, ph);
    }
  }
  const thresh = maxVal * 0.04;
  ctx.strokeStyle = 'rgba(255,255,255,0.85)'; ctx.lineWidth = 2;
  for (let iy = 0; iy < MF_N - 1; iy++) {
    for (let ix = 0; ix < MF_N - 1; ix++) {
      const b1 = toVal(ix / (MF_N - 1)), b2 = toVal(iy / (MF_N - 1));
      if (admissibleFn && !admissibleFn(b1, b2)) continue;
      const inC = (v) => v < thresh;
      if ((inC(grid[iy][ix]) || inC(grid[iy][ix + 1]) || inC(grid[iy + 1]?.[ix] ?? 1) || inC(grid[iy + 1]?.[ix + 1] ?? 1)) &&
        !(inC(grid[iy][ix]) && inC(grid[iy][ix + 1]) && inC(grid[iy + 1]?.[ix] ?? 1) && inC(grid[iy + 1]?.[ix + 1] ?? 1))) {
        ctx.beginPath(); ctx.arc(Math.round((ix + 0.5) * W / MF_N), Math.round((MF_N - 1.5 - iy) * H / MF_N), 1.8, 0, Math.PI * 2); ctx.stroke();
      }
    }
  }
  ctx.fillStyle = 'rgba(255,255,255,0.6)'; ctx.font = '10px sans-serif';
  ctx.fillText('b₁ →', W / 2 - 14, H - 3);
  ctx.save(); ctx.translate(8, H / 2 + 12); ctx.rotate(-Math.PI / 2); ctx.fillText('b₂ →', 0, 0); ctx.restore();
  return maxVal;
}

function toPixel(v, W_or_H) { return (v - MF_RANGE[0]) / (MF_RANGE[1] - MF_RANGE[0]) * (W_or_H); }

function renderMfCompareCharts(refPar, selPar, mdl, crfId, erfId) {
  const rRef = computeFull(refPar, mdl), rSel = computeFull(selPar, mdl);
  makeChart(crfId, { data: { datasets: [...rRef.crf.map((c, k) => ds(c, CAT_COLS[k])), ...rSel.crf.map((c, k) => ds(c, CAT_COLS[k], [4, 3]))] }, ...baseOpts('P(X=k|θ)') });
  makeChart(erfId, { data: { datasets: [ds(rRef.erf, COL_REF), ds(rSel.erf, COL_B, [4, 3])] }, ...baseOpts('E[X|θ]', K - 1) });
}
