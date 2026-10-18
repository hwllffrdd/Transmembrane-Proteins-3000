// Logic for templates/batch_results.html. The results payload and the
// protein_sequence URL template are injected by the server as a JSON data
// island / data-* attribute (see #batch-results-data in that template)
// rather than inlined here, since this is a static file Flask/Jinja
// doesn't process.
(function() {
  const PAGE_SIZE = 100;

  const VIZ_METHODS = [
    { key: 'uniprot',   label: 'UniProt',   cls: 'uniprot-region'   },
    { key: 'phobius',   label: 'Phobius',   cls: 'phobius-region'   },
    { key: 'deeptmhmm', label: 'DeepTMHMM', cls: 'deeptmhmm-region' },
    { key: 'tmbed',     label: 'TMBED',     cls: 'tmbed-region'     },
    { key: 'topcons',   label: 'TOPCONS',   cls: 'topcons-region'   },
  ];

  // Maps a th's data-sort-key to the result field it sorts by. Single
  // source of truth for both the sort comparator and the click handler.
  const SORT_FIELD_MAP = {
    'name':            { field: 'name',                type: 'text' },
    'uniprot-count':   { field: 'uniprot_count',        type: 'num'  },
    'uniprot-pct':     { field: 'uniprot_percentage',   type: 'num'  },
    'phobius-count':   { field: 'phobius_count',        type: 'num'  },
    'phobius-pct':     { field: 'phobius_percentage',   type: 'num'  },
    'deeptmhmm-count': { field: 'deeptmhmm_count',      type: 'num'  },
    'deeptmhmm-pct':   { field: 'deeptmhmm_percentage', type: 'num'  },
    'tmbed-count':     { field: 'tmbed_count',           type: 'num'  },
    'tmbed-pct':       { field: 'tmbed_percentage',      type: 'num'  },
    'topcons-count':   { field: 'topcons_count',         type: 'num'  },
    'topcons-pct':     { field: 'topcons_percentage',    type: 'num'  },
  };

  // Method count fields/colors mirror the "Summary Statistics" table on
  // the full batch-stats page (templates/batch_stats.html) — this is
  // just its first three columns (Method / IMPs / % IMPs), computed
  // client-side from the same result data rather than round-tripping to
  // the stats endpoint.
  const QUICK_STATS_METHODS = {
    UniProt:   { field: 'uniprot_count',   color: '#3498db' },
    Phobius:   { field: 'phobius_count',   color: '#e74c3c' },
    DeepTMHMM: { field: 'deeptmhmm_count', color: '#2ecc71' },
    TMBED:     { field: 'tmbed_count',     color: '#9b59b6' },
    TOPCONS:   { field: 'topcons_count',   color: '#f39c12' },
  };

  function renderQuickStats(results) {
    const tbody = document.getElementById('quick-stats-tbody');
    if (!tbody) return;
    const n = results.length;
    tbody.innerHTML = Object.keys(QUICK_STATS_METHODS).map(method => {
      const { field, color } = QUICK_STATS_METHODS[method];
      const withTM = results.filter(r => parseInt(r[field], 10) > 0).length;
      const pct = n > 0 ? (withTM / n * 100).toFixed(1) : '0.0';
      return `
        <tr>
          <td><span class="legend-dot" style="background:${color}"></span>${method}</td>
          <td>${withTM}</td>
          <td>${pct}%</td>
        </tr>`;
    }).join('');
  }

  // "At least" figures for the summary line above the table: since the 5
  // methods disagree on individual proteins, the only IMP/non-IMP counts
  // that hold no matter which method you trust are the most conservative
  // ones — the method with the fewest IMP calls, and the method with the
  // fewest non-IMP calls (not necessarily the same method).
  function renderQuerySummary(results) {
    const el = document.getElementById('query-summary');
    if (!el) return;
    const total = results.length;
    if (total === 0) {
      el.textContent = 'Your query included 0 valid accessions.';
      return;
    }
    const methodImpCounts = Object.values(QUICK_STATS_METHODS).map(({ field }) =>
      results.filter(r => parseInt(r[field], 10) > 0).length
    );
    const minImps = Math.min(...methodImpCounts);
    const minNonImps = Math.min(...methodImpCounts.map(c => total - c));
    const pct = (minImps / total * 100).toFixed(1);
    el.textContent = `Your query included ${total} accessions, at least ${minImps} (${pct}%) are predicted to be IMPs, at least ${minNonImps} are not IMPs.`;
  }

  // A row counts as "0 TM helices" for the filter below when every one of
  // the 5 methods reports zero (or has no data at all) — i.e. all methods
  // agree there's nothing there, not just the default sort field.
  const TM_COUNT_FIELDS = ['uniprot_count', 'phobius_count', 'deeptmhmm_count', 'tmbed_count', 'topcons_count'];
  function isAllMethodsZero(result) {
    return TM_COUNT_FIELDS.every(field => !(parseInt(result[field], 10) > 0));
  }

  const dataEl = document.getElementById('batch-results-data');
  const allResults = JSON.parse(dataEl.textContent);
  renderQuickStats(allResults);
  renderQuerySummary(allResults);

  // "Quick Summary" floating panel: toggled by its trigger button, and
  // dismissed either by its own close control or a click anywhere outside
  // it (standard popover behavior).
  const quickSummaryBtn = document.getElementById('quick-summary-btn');
  const quickSummaryPanel = document.getElementById('quick-summary-panel');
  const quickSummaryClose = document.getElementById('quick-summary-close');
  if (quickSummaryBtn && quickSummaryPanel) {
    quickSummaryBtn.addEventListener('click', (evt) => {
      evt.stopPropagation();
      quickSummaryPanel.classList.toggle('open');
    });
    if (quickSummaryClose) {
      quickSummaryClose.addEventListener('click', () => {
        quickSummaryPanel.classList.remove('open');
      });
    }
    document.addEventListener('click', (evt) => {
      if (!quickSummaryPanel.classList.contains('open')) return;
      if (quickSummaryPanel.contains(evt.target) || quickSummaryBtn.contains(evt.target)) return;
      quickSummaryPanel.classList.remove('open');
    });
  }
  const hideZeroTmCheckbox = document.getElementById('hide-zero-tm');
  let filteredResults = hideZeroTmCheckbox.checked ? allResults.filter(r => !isAllMethodsZero(r)) : allResults.slice();
  let sortedResults = filteredResults.slice();
  // No default sort — results start in the order the user entered them
  // (as returned by the backend). sortKey stays null until a column
  // header is clicked.
  let sortKey = null;
  let sortDir = 'desc';
  let currentPage = 1;
  let currentPageSlice = [];

  const tbody = document.querySelector('.batch-results-table tbody');
  const tableWrapper = document.querySelector('.batch-results-table-wrapper');
  let hasRenderedOnce = false;
  const sortableThs = document.querySelectorAll('.sortable-th');
  const pagePrevBtns = [document.getElementById('page-prev'), document.getElementById('page-prev-bottom')];
  const pageNextBtns = [document.getElementById('page-next'), document.getElementById('page-next-bottom')];
  const pageRangeStartEls = [document.getElementById('page-range-start'), document.getElementById('page-range-start-bottom')];
  const pageRangeEndEls = [document.getElementById('page-range-end'), document.getElementById('page-range-end-bottom')];
  const pageRangeTotalEls = [document.getElementById('page-range-total'), document.getElementById('page-range-total-bottom')];
  const pageNumbersEls = [document.getElementById('page-numbers'), document.getElementById('page-numbers-bottom')];
  const pageJumpEls = [document.getElementById('page-jump'), document.getElementById('page-jump-bottom')];
  const pageJumpInputs = [document.getElementById('page-jump-input'), document.getElementById('page-jump-input-bottom')];
  const pageJumpGoBtns = [document.getElementById('page-jump-go'), document.getElementById('page-jump-go-bottom')];

  function escapeHtml(value) {
    if (value === null || value === undefined) return '';
    return String(value)
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;')
      .replace(/'/g, '&#39;');
  }

  function compareResults(a, b, field, type, dir) {
    let result;
    if (type === 'text') {
      const av = (a[field] || '').toLowerCase();
      const bv = (b[field] || '').toLowerCase();
      result = av.localeCompare(bv);
    } else {
      const av = parseFloat(a[field]);
      const bv = parseFloat(b[field]);
      const aMissing = isNaN(av);
      const bMissing = isNaN(bv);
      // Rows without a value for this method always sink to the bottom,
      // regardless of sort direction — there's nothing to rank them by.
      if (aMissing && bMissing) return 0;
      if (aMissing) return 1;
      if (bMissing) return -1;
      result = av - bv;
    }
    return dir === 'asc' ? result : -result;
  }

  function sortResults() {
    if (!sortKey) return; // no column selected — keep insertion order
    const { field, type } = SORT_FIELD_MAP[sortKey];
    sortedResults.sort((a, b) => compareResults(a, b, field, type, sortDir));
  }

  // Re-derives sortedResults from filteredResults (e.g. after the filter
  // checkbox or a column sort changes), reapplying the current sort key.
  function rebuildSortedResults() {
    sortedResults = filteredResults.slice();
    sortResults();
  }

  function updateSortIndicators() {
    sortableThs.forEach(th => {
      th.classList.toggle('sorted', th.dataset.sortKey === sortKey);
    });
  }

  function rowHtml(result, i) {
    const vizId = 'viz-' + i;
    const uniprotUrl = 'https://www.uniprot.org/uniprotkb/' + encodeURIComponent(result.accession) + '/entry';
    return `
      <tr class="data-row" data-target="${vizId}" data-rendered="false" data-idx="${i}">
        <td class="protein-identity">
          <a class="prot-name" href="${uniprotUrl}" target="_blank" rel="noopener noreferrer" onclick="event.stopPropagation()">${escapeHtml(result.name)}</a>
          <span class="prot-gene">Gene: ${escapeHtml(result.gene_name)}</span>
          <span class="prot-acc-line"><a class="prot-acc" href="${uniprotUrl}" target="_blank" rel="noopener noreferrer" onclick="event.stopPropagation()">${escapeHtml(result.accession)}</a>${result.entry_name ? ` <span class="prot-entry">(${escapeHtml(result.entry_name)})</span>` : ''}</span>
          <span class="expand-indicator">&#9660; show visualization</span>
        </td>
        <td>${escapeHtml(result.uniprot_count)}</td>
        <td>${escapeHtml(result.uniprot_regions)}</td>
        <td>${escapeHtml(result.uniprot_percentage)}</td>
        <td>${escapeHtml(result.phobius_count)}</td>
        <td>${escapeHtml(result.phobius_regions_str)}</td>
        <td>${escapeHtml(result.phobius_percentage)}</td>
        <td>${escapeHtml(result.deeptmhmm_count)}</td>
        <td>${escapeHtml(result.deeptmhmm_regions_str)}</td>
        <td>${escapeHtml(result.deeptmhmm_percentage)}</td>
        <td>${escapeHtml(result.tmbed_count)}</td>
        <td>${escapeHtml(result.tmbed_regions_str)}</td>
        <td>${escapeHtml(result.tmbed_percentage)}</td>
        <td>${escapeHtml(result.topcons_count)}</td>
        <td>${escapeHtml(result.topcons_regions_str)}</td>
        <td>${escapeHtml(result.topcons_percentage)}</td>
      </tr>
      <tr class="visualizer-row" id="${vizId}">
        <td colspan="16">
          <div class="viz-inner">
            <div class="protein-visualization" id="viz-content-${i}">
              <!-- rendered by JS on first expand -->
            </div>
          </div>
        </td>
      </tr>
    `;
  }

  function totalPages() {
    return Math.max(1, Math.ceil(sortedResults.length / PAGE_SIZE));
  }

  function renderPage() {
    const tp = totalPages();
    currentPage = Math.min(Math.max(1, currentPage), tp);

    const start = (currentPage - 1) * PAGE_SIZE;
    currentPageSlice = sortedResults.slice(start, start + PAGE_SIZE);

    tbody.innerHTML = currentPageSlice.map(rowHtml).join('');

    pageRangeStartEls.forEach(el => { el.textContent = sortedResults.length === 0 ? 0 : start + 1; });
    pageRangeEndEls.forEach(el => { el.textContent = start + currentPageSlice.length; });
    pageRangeTotalEls.forEach(el => { el.textContent = sortedResults.length; });
    pagePrevBtns.forEach(btn => { btn.disabled = currentPage <= 1; });
    pageNextBtns.forEach(btn => { btn.disabled = currentPage >= tp; });
    renderPageNumbers(currentPage, tp);

    // The table now sits directly in the page (no bounded-height scroll
    // frame), so switching pages no longer resets an internal scrollbar —
    // instead scroll the whole page back to the top of the table, so the
    // new page's rows are visible under the sticky header rather than
    // wherever the previous page happened to leave the scroll position.
    // Skipped on the very first render so loading the page doesn't jump
    // past the summary/quick-stats section above the table.
    if (tableWrapper && hasRenderedOnce) {
      tableWrapper.scrollIntoView({ behavior: 'smooth', block: 'start' });
    }
    hasRenderedOnce = true;
  }

  // Builds a compact page list like "1 2 3 4 5 … 15 16 17 18 19 … 46 47
  // 48 49 50" (active page 17 of 50): the first/last BOUNDARY pages, plus
  // a window of DELTA pages on each side of the current page, with a "…"
  // wherever there's a gap between the pieces.
  function getPageNumbers(current, total) {
    const BOUNDARY = 5;
    const DELTA = 2;
    const shown = new Set();
    for (let p = 1; p <= Math.min(BOUNDARY, total); p++) shown.add(p);
    for (let p = Math.max(1, total - BOUNDARY + 1); p <= total; p++) shown.add(p);
    for (let p = Math.max(1, current - DELTA); p <= Math.min(total, current + DELTA); p++) shown.add(p);

    const sorted = Array.from(shown).sort((a, b) => a - b);
    const items = [];
    let prev = null;
    sorted.forEach(p => {
      if (prev !== null && p - prev > 1) items.push(null); // null = ellipsis
      items.push(p);
      prev = p;
    });
    return items;
  }

  function renderPageNumbers(current, total) {
    const items = getPageNumbers(current, total);
    const html = items.map(p => {
      if (p === null) return '<span class="page-ellipsis">&hellip;</span>';
      const active = p === current;
      return `<button type="button" class="page-number${active ? ' active' : ''}" data-page="${p}"${active ? ' aria-current="page"' : ''}>${p}</button>`;
    }).join('');
    pageNumbersEls.forEach(el => { el.innerHTML = html; });

    // "Jump to page" only earns its keep once the strip has a gap (…) —
    // i.e. it isn't already listing every page.
    const hasGap = items.includes(null);
    pageJumpEls.forEach(el => { el.style.display = hasGap ? '' : 'none'; });
    pageJumpInputs.forEach(input => { input.max = total; });
  }

  // Sequences are fetched lazily, one protein at a time, only when a row is
  // expanded — and cached here so re-expanding the same row is instant. This
  // keeps result pages light: a 4000-protein batch ships no sequence data
  // until the user actually asks to see a visualisation.
  const SEQUENCE_URL_TEMPLATE = dataEl.dataset.sequenceUrlTemplate;
  const sequenceCache = {};

  function fetchSequence(accession) {
    if (Object.prototype.hasOwnProperty.call(sequenceCache, accession)) {
      return Promise.resolve(sequenceCache[accession]);
    }
    const url = SEQUENCE_URL_TEMPLATE.replace('__ACC__', encodeURIComponent(accession));
    return fetch(url)
      .then(resp => (resp.ok ? resp.json() : null))
      .then(data => {
        const seq = (data && data.sequence) ? data.sequence : null;
        sequenceCache[accession] = seq;
        return seq;
      })
      .catch(() => {
        sequenceCache[accession] = null;
        return null;
      });
  }

  // Once the sequence arrives, rewrite each region's tooltip to include the
  // actual residues. If the sequence isn't available the position-only
  // tooltip written at render time simply stays as it is.
  function enrichRegionsWithSequence(container, accession) {
    const bars = container.querySelectorAll('.tm-region');
    const note = container.querySelector('.sequence-note');
    if (!bars.length) {
      if (note) note.remove();
      return;
    }

    fetchSequence(accession).then(sequence => {
      if (!sequence) {
        if (note) note.remove();
        return;
      }
      bars.forEach(bar => {
        const start = parseInt(bar.dataset.start, 10);
        const end   = parseInt(bar.dataset.end, 10);
        if (!start || !end) return;
        // Region coordinates are 1-based inclusive; slice() is 0-based, end-exclusive.
        const residues = sequence.slice(start - 1, end);
        bar.dataset.copyText = bar.dataset.posLabel + '\n' + residues;
        bar.title = bar.dataset.copyText + '\n(click to copy)';
      });
      if (note) {
        note.textContent = 'Hover any region to see its amino-acid sequence — click a region to copy it.';
      }
    });
  }

  function renderViz(dataRow) {
    const idx = parseInt(dataRow.dataset.idx, 10);
    const result = currentPageSlice[idx];
    const targetId = dataRow.dataset.target;
    const container = document.getElementById('viz-content-' + targetId.replace('viz-', ''));
    const seqLen = result.sequence_length;

    if (!seqLen) {
      container.innerHTML = '<p class="hint">Sequence length unavailable — cannot render visualisation.</p>';
      return;
    }

    const regions = {
      uniprot:   result.transmembrane_regions,
      phobius:   result.phobius_regions,
      deeptmhmm: result.deeptmhmm_regions,
      tmbed:     result.tmbed_regions,
      topcons:   result.topcons_regions,
    };

    let html = `<div class="sequence-info">Protein length: ${seqLen} amino acids<br>Move cursor to predicted
TM seq to see its aminoacid position and sequence</div>`;

    VIZ_METHODS.forEach(method => {
      const methodRegions = regions[method.key] || [];
      let barsHtml = '';
      methodRegions.forEach(region => {
        const start = region[0], end = region[1];
        const left  = (start / seqLen * 100).toFixed(3);
        const width = ((end - start + 1) / seqLen * 100).toFixed(3);
        const posLabel = `Position ${start}–${end} (${end - start + 1} aa)`;
        barsHtml += `<div class="tm-region ${method.cls}"
          style="left:${left}%;width:${width}%;"
          data-start="${start}" data-end="${end}" data-pos-label="${posLabel}" data-copy-text="${posLabel}"
          title="${posLabel} (click to copy)"></div>`;
      });
      html += `<div class="region-label ${method.key}-label">${method.label}</div>
               <div class="protein-bar">${barsHtml}</div>`;
    });

    html += `<div class="sequence-note">Loading sequence…</div>`;

    container.innerHTML = html;

    enrichRegionsWithSequence(container, result.accession);
  }

  // Clicking a TM-region bar copies its position (and, once loaded, its
  // amino-acid sequence) to the clipboard — same text the hover tooltip
  // shows. Briefly flashes the bar so the click has visible feedback.
  function copyRegionText(bar) {
    const text = bar.dataset.copyText;
    if (!text || !navigator.clipboard) return;
    navigator.clipboard.writeText(text).then(() => {
      bar.classList.add('copied');
      setTimeout(() => bar.classList.remove('copied'), 500);
    }).catch(() => {});
  }

  // Event delegation: rows are re-created on every sort/page change, so a
  // single listener on the (persistent) tbody handles clicks for all of
  // them instead of rebinding per row.
  tbody.addEventListener('click', (evt) => {
    const bar = evt.target.closest('.tm-region');
    if (bar) {
      evt.stopPropagation();
      copyRegionText(bar);
      return;
    }

    const dataRow = evt.target.closest('tr.data-row');
    if (!dataRow) return;

    // If the user just made a text selection (e.g. dragging across a
    // Regions cell to copy it), don't treat that as a row-toggle click.
    const selection = window.getSelection();
    if (selection && selection.toString().length > 0 && dataRow.contains(selection.anchorNode)) {
      return;
    }

    const vizRow = document.getElementById(dataRow.dataset.target);
    const indicator = dataRow.querySelector('.expand-indicator');
    const isOpen = dataRow.classList.contains('open');

    if (isOpen) {
      dataRow.classList.remove('open');
      vizRow.classList.remove('open');
      indicator.innerHTML = '&#9660; show visualisation';
    } else {
      if (dataRow.dataset.rendered === 'false') {
        renderViz(dataRow);
        dataRow.dataset.rendered = 'true';
      }
      dataRow.classList.add('open');
      vizRow.classList.add('open');
      indicator.innerHTML = '&#9650; hide visualisation';
    }
  });

  sortableThs.forEach(th => {
    th.addEventListener('click', () => {
      const key = th.dataset.sortKey;
      const { type } = SORT_FIELD_MAP[key];
      sortDir = sortKey === key
        ? (sortDir === 'asc' ? 'desc' : 'asc')
        : (type === 'text' ? 'asc' : 'desc');
      sortKey = key;
      rebuildSortedResults();
      updateSortIndicators();
      currentPage = 1;
      renderPage();
    });
  });

  hideZeroTmCheckbox.addEventListener('change', () => {
    filteredResults = hideZeroTmCheckbox.checked
      ? allResults.filter(r => !isAllMethodsZero(r))
      : allResults.slice();
    rebuildSortedResults();
    currentPage = 1;
    renderPage();
  });

  pagePrevBtns.forEach(btn => { btn.addEventListener('click', () => { currentPage -= 1; renderPage(); }); });
  pageNextBtns.forEach(btn => { btn.addEventListener('click', () => { currentPage += 1; renderPage(); }); });

  pageNumbersEls.forEach(el => {
    el.addEventListener('click', (evt) => {
      const btn = evt.target.closest('.page-number');
      if (!btn) return;
      currentPage = parseInt(btn.dataset.page, 10);
      renderPage();
    });
  });

  function jumpToPage(input) {
    const tp = totalPages();
    const page = parseInt(input.value, 10);
    if (!page || page < 1 || page > tp) {
      input.value = '';
      return;
    }
    currentPage = page;
    input.value = '';
    renderPage();
  }
  pageJumpGoBtns.forEach((btn, i) => {
    btn.addEventListener('click', () => jumpToPage(pageJumpInputs[i]));
  });
  pageJumpInputs.forEach(input => {
    input.addEventListener('keydown', (evt) => {
      if (evt.key === 'Enter') jumpToPage(input);
    });
  });

  // Initial view: unsorted (insertion order), zero-TM entries hidden, page 1.
  updateSortIndicators();
  renderPage();
})();

function updateStickyHeaderOffsets() {
  // Measure a first-row th without rowspan, so the rowspan'd "Protein"
  // cell (which spans both rows) doesn't skew the row-1 height reading.
  var row1Th = document.querySelector('.batch-results-table thead tr:first-child th:not([rowspan])');
  if (row1Th) {
    document.documentElement.style.setProperty('--thead-row1-height', row1Th.offsetHeight + 'px');
  }
}
updateStickyHeaderOffsets();
window.addEventListener('load', updateStickyHeaderOffsets);
window.addEventListener('resize', updateStickyHeaderOffsets);
