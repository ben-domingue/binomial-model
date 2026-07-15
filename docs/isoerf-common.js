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

/* Pale-yellow -> orange -> dark-maroon sequential ramp (matches the RMISE
   colorbar used in the manuscript's static isoERF-distance figure). */
const HEAT_STOPS = [
  [0.0, [254, 248, 190]], [0.1, [253, 234, 167]], [0.2, [249, 216, 130]], [0.3, [248, 193, 77]],
  [0.4, [245, 167, 33]], [0.5, [242, 138, 0]], [0.6, [238, 104, 0]], [0.7, [223, 68, 0]],
  [0.8, [196, 37, 0]], [0.9, [162, 12, 14]], [1.0, [134, 2, 33]]
];
function heatColor(v) {
  v = Math.min(Math.max(v, 0), 1);
  for (let i = 0; i < HEAT_STOPS.length - 1; i++) {
    const [f0, c0] = HEAT_STOPS[i], [f1, c1] = HEAT_STOPS[i + 1];
    if (v <= f1) {
      const t = (v - f0) / (f1 - f0);
      return `rgb(${Math.round(c0[0] + t * (c1[0] - c0[0]))},${Math.round(c0[1] + t * (c1[1] - c0[1]))},${Math.round(c0[2] + t * (c1[2] - c0[2]))})`;
    }
  }
  return `rgb(${HEAT_STOPS[HEAT_STOPS.length - 1][1].join(',')})`;
}

/* Simplified marching squares: returns line segments [[x1,y1],[x2,y2]] where
   `grid` crosses `level`. `cellOk(ix,iy)` can veto a cell (e.g. inadmissible
   region) so contour lines don't bleed into the grayed-out area. */
function marchingSquares(grid, W, H, level, cellOk) {
  const N = MF_N;
  const px = ix => ix * W / (N - 1), py = iy => (N - 1 - iy) * H / (N - 1);
  const interp = (x1, y1, v1, x2, y2, v2) => { const t = (level - v1) / (v2 - v1); return [x1 + t * (x2 - x1), y1 + t * (y2 - y1)]; };
  const segs = [];
  for (let iy = 0; iy < N - 1; iy++) {
    for (let ix = 0; ix < N - 1; ix++) {
      if (cellOk && !cellOk(ix, iy)) continue;
      const vbl = grid[iy][ix], vbr = grid[iy][ix + 1], vtr = grid[iy + 1][ix + 1], vtl = grid[iy + 1][ix];
      const xbl = px(ix), ybl = py(iy), xbr = px(ix + 1), ybr = py(iy), xtr = px(ix + 1), ytr = py(iy + 1), xtl = px(ix), ytl = py(iy + 1);
      const bl = vbl >= level, br = vbr >= level, tr = vtr >= level, tl = vtl >= level;
      const c = (bl ? 1 : 0) | (br ? 2 : 0) | (tr ? 4 : 0) | (tl ? 8 : 0);
      if (c === 0 || c === 15) continue;
      const a = () => interp(xbl, ybl, vbl, xbr, ybr, vbr), b = () => interp(xbr, ybr, vbr, xtr, ytr, vtr);
      const cc = () => interp(xtr, ytr, vtr, xtl, ytl, vtl), d = () => interp(xtl, ytl, vtl, xbl, ybl, vbl);
      const avg = (vbl + vbr + vtr + vtl) / 4;
      switch (c) {
        case 1: segs.push([d(), a()]); break;
        case 2: segs.push([a(), b()]); break;
        case 3: segs.push([d(), b()]); break;
        case 4: segs.push([b(), cc()]); break;
        case 5: if (avg >= level) { segs.push([d(), cc()]); segs.push([a(), b()]); } else { segs.push([d(), a()]); segs.push([cc(), b()]); } break;
        case 6: segs.push([a(), cc()]); break;
        case 7: segs.push([d(), cc()]); break;
        case 8: segs.push([cc(), d()]); break;
        case 9: segs.push([cc(), a()]); break;
        case 10: if (avg >= level) { segs.push([a(), d()]); segs.push([b(), cc()]); } else { segs.push([a(), b()]); segs.push([d(), cc()]); } break;
        case 11: segs.push([cc(), b()]); break;
        case 12: segs.push([b(), d()]); break;
        case 13: segs.push([b(), a()]); break;
        case 14: segs.push([a(), d()]); break;
      }
    }
  }
  return segs;
}
function strokeSegs(ctx, segs) { ctx.beginPath(); segs.forEach(([[x1, y1], [x2, y2]]) => { ctx.moveTo(x1, y1); ctx.lineTo(x2, y2); }); ctx.stroke(); }

/* White-halo markers so ✕/● stay legible against every part of the
   pale-yellow -> dark-maroon heatmap ramp, not just the light end. */
function drawHaloCross(ctx, x, y, r, color, lineWidth = 2.5) {
  ctx.lineWidth = lineWidth + 2.5; ctx.strokeStyle = 'rgba(255,255,255,0.95)';
  ctx.beginPath(); ctx.moveTo(x - r, y - r); ctx.lineTo(x + r, y + r); ctx.moveTo(x + r, y - r); ctx.lineTo(x - r, y + r); ctx.stroke();
  ctx.lineWidth = lineWidth; ctx.strokeStyle = color;
  ctx.beginPath(); ctx.moveTo(x - r, y - r); ctx.lineTo(x + r, y + r); ctx.moveTo(x + r, y - r); ctx.lineTo(x - r, y + r); ctx.stroke();
}
function drawHaloDot(ctx, x, y, r, color, haloWidth = 1.5) {
  ctx.fillStyle = 'rgba(255,255,255,0.95)'; ctx.beginPath(); ctx.arc(x, y, r + haloWidth, 0, Math.PI * 2); ctx.fill();
  ctx.fillStyle = color; ctx.beginPath(); ctx.arc(x, y, r, 0, Math.PI * 2); ctx.fill();
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
      ctx.fillStyle = heatColor(grid[iy][ix] / maxVal);
      ctx.fillRect(px, py, pw, ph);
    }
  }
  const cellOk = admissibleFn ? (ix, iy) => {
    const b1a = toVal(ix / (MF_N - 1)), b2a = toVal(iy / (MF_N - 1)), b1b = toVal((ix + 1) / (MF_N - 1)), b2b = toVal((iy + 1) / (MF_N - 1));
    return admissibleFn(b1a, b2a) && admissibleFn(b1b, b2a) && admissibleFn(b1a, b2b) && admissibleFn(b1b, b2b);
  } : null;
  ctx.strokeStyle = 'rgba(90,90,90,0.55)'; ctx.lineWidth = 1; ctx.setLineDash([]);
  for (let k = 1; k <= 8; k++) strokeSegs(ctx, marchingSquares(grid, W, H, maxVal * k / 9, cellOk));
  const thresh = maxVal * 0.04;
  ctx.strokeStyle = COL_B; ctx.lineWidth = 2; ctx.setLineDash([5, 3]);
  strokeSegs(ctx, marchingSquares(grid, W, H, thresh, cellOk));
  ctx.setLineDash([]);
  ctx.font = '10px sans-serif'; ctx.lineWidth = 3; ctx.strokeStyle = 'rgba(255,255,255,0.9)'; ctx.fillStyle = 'rgba(55,55,55,0.9)';
  ctx.strokeText('b₁ →', W / 2 - 14, H - 3); ctx.fillText('b₁ →', W / 2 - 14, H - 3);
  ctx.save(); ctx.translate(8, H / 2 + 12); ctx.rotate(-Math.PI / 2);
  ctx.strokeText('b₂ →', 0, 0); ctx.fillText('b₂ →', 0, 0); ctx.restore();
  return maxVal;
}

function toPixel(v, W_or_H) { return (v - MF_RANGE[0]) / (MF_RANGE[1] - MF_RANGE[0]) * (W_or_H); }

function renderMfCompareCharts(refPar, selPar, mdl, crfId, erfId) {
  const rRef = computeFull(refPar, mdl), rSel = computeFull(selPar, mdl);
  makeChart(crfId, { data: { datasets: [...rRef.crf.map((c, k) => ds(c, CAT_COLS[k])), ...rSel.crf.map((c, k) => ds(c, CAT_COLS[k], [4, 3]))] }, ...baseOpts('P(X=k|θ)') });
  makeChart(erfId, { data: { datasets: [ds(rRef.erf, COL_REF), ds(rSel.erf, COL_B, [4, 3])] }, ...baseOpts('E[X|θ]', K - 1) });
}
