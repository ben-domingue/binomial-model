/* Single source of truth for the widget pages' navigation bar and numbering.
   Widgets are numbered in the order they are first referenced in the manuscript's footnotes;
   index.html lists them in the same order. Each page includes
   <nav class="related-nav" id="site-nav"></nav> and this script. */
(function () {
  const PAGES = [
    ['threshold-structure-explorer.html', 'Threshold independence'],
    ['crf-editor.html', 'CRF editor'],
    ['isoerf-manifold-explorer.html', 'isoERF manifold'],
    ['cross-model-erf-explorer.html', 'Same ERF, different models'],
    ['fixed-erf-threshold-explorer.html', 'Fixed-ERF thresholds'],
    ['category-loss-explorer.html', 'Category loss & KL']
  ];
  const nav = document.getElementById('site-nav');
  if (!nav) return;
  const here = location.pathname.split('/').pop() || 'index.html';
  nav.innerHTML = PAGES.map(([href, label], i) =>
    `<a href="${href}"${href === here ? ' class="current" aria-current="page"' : ''}>${i + 1}. ${label}</a>`).join('');
})();
