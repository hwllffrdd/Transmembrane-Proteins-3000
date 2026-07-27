// batch_stats.js — chart logic for the batch statistics page
// RAW and PROTEOME_* constants are injected by the inline <script> block in batch_stats.html

// ════════════════════════════════════════════════════════════════════════════
// CONSTANTS
// ════════════════════════════════════════════════════════════════════════════

const METHODS = ['UniProt', 'Phobius', 'DeepTMHMM', 'TMBED', 'TOPCONS'];

const COLORS = {
    'UniProt':   '#3498db',
    'Phobius':   '#e74c3c',
    'DeepTMHMM': '#2ecc71',
    'TMBED':     '#9b59b6',
    'TOPCONS':   '#f39c12',
};

const DOT_CLASS = {
    'UniProt':   'dot-uniprot',
    'Phobius':   'dot-phobius',
    'DeepTMHMM': 'dot-deeptmhmm',
    'TMBED':     'dot-tmbed',
    'TOPCONS':   'dot-topcons',
};

const MAX_CAT  = 12;
const NUM_BINS = MAX_CAT + 1;

const BIN_SIZE = 5;

const PCT_NUM_BINS = Math.max(...Object.values(PROTEOME_PCT_DIST).map(a => a.length));
const PCT_BIN_MAX  = (PCT_NUM_BINS - 1) * BIN_SIZE;

// ════════════════════════════════════════════════════════════════════════════
// HELPERS
// ════════════════════════════════════════════════════════════════════════════

function getCount(r, method) {
    const v = r[`${method.toLowerCase()}_count`];
    return (v === null || v === undefined || v === '') ? 0 : parseInt(v, 10);
}
function getCounts(method) { return RAW.map(r => getCount(r, method)); }

function getPct(r, method) {
    const v = r[`${method.toLowerCase()}_percentage`];
    if (v === null || v === undefined || v === '') return null;
    const f = parseFloat(v);
    return isNaN(f) ? null : f;
}

function makePieLabelsPlugin(id, totalRef) {
    return {
        id,
        afterDatasetsDraw(chart) {
            const { ctx, data } = chart;
            const total = totalRef();
            chart.getDatasetMeta(0).data.forEach((arc, i) => {
                const count = data.datasets[0].data[i];
                if (!count) return;
                const pct = (count / total * 100).toFixed(1);
                const { x, y } = arc.tooltipPosition();
                ctx.save();
                ctx.fillStyle = 'white'; ctx.font = 'bold 12px Arial';
                ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
                ctx.fillText(`${pct}%`, x, y);
                ctx.restore();
            });
        }
    };
}

function buildPctBins(values) {
    const bins = new Array(PCT_NUM_BINS).fill(0);
    values.forEach(v => {
        if (v >= PCT_BIN_MAX) bins[PCT_NUM_BINS - 1]++;
        else bins[Math.min(Math.floor(v / BIN_SIZE), PCT_NUM_BINS - 2)]++;
    });
    return bins;
}

// ── Protein list panel helpers ────────────────────────────────────────────

const _openSlice = { agree: null, cons: null };

function openPanel(panelKey, title, rows) {
    const panel   = document.getElementById(`${panelKey}-protein-panel`);
    const titleEl = document.getElementById(`${panelKey}-panel-title`);
    const tbodyEl = document.getElementById(`${panelKey}-panel-tbody`);

    titleEl.textContent = title;
    tbodyEl.innerHTML = rows.map(r => `
        <tr>
            <td><a href="https://www.uniprot.org/uniprotkb/${r.accession}/entry"
                   target="_blank" rel="noopener noreferrer">${r.accession}</a></td>
            <td>${r.gene_name}</td>
            <td>${r.name}</td>
            <td>${r.label}</td>
        </tr>
    `).join('');

    panel.classList.add('open');
    setTimeout(() => panel.scrollIntoView({ behavior: 'smooth', block: 'nearest' }), 50);
}

function closePanel(panelKey) {
    document.getElementById(`${panelKey}-protein-panel`).classList.remove('open');
    _openSlice[panelKey] = null;
}

window.closePanel = closePanel;

// ════════════════════════════════════════════════════════════════════════════
// SHARED PRE-COMPUTED AVERAGES
// ════════════════════════════════════════════════════════════════════════════

const datasetCountAvg = new Array(NUM_BINS).fill(0);
METHODS.forEach(method => {
    const counts = getCounts(method);
    for (let b = 0; b < MAX_CAT; b++) datasetCountAvg[b] += counts.filter(c => c === b + 1).length;
    datasetCountAvg[NUM_BINS - 1] += counts.filter(c => c > MAX_CAT).length;
});
for (let b = 0; b < NUM_BINS; b++) datasetCountAvg[b] /= METHODS.length;

const proteomeCountAvg = new Array(NUM_BINS).fill(0);
let methodsWithCountData = 0;
METHODS.forEach(method => {
    if (PROTEOME_SIZE[method] === 0) return;
    methodsWithCountData++;
    PROTEOME_DIST[method].forEach((v, b) => { proteomeCountAvg[b] += v; });
});
if (methodsWithCountData > 0)
    for (let b = 0; b < NUM_BINS; b++) proteomeCountAvg[b] /= methodsWithCountData;

const datasetPctAvg = new Array(PCT_NUM_BINS).fill(0);
METHODS.forEach(method => {
    const vals = RAW.map(r => getPct(r, method)).filter(v => v !== null && v > 0);
    buildPctBins(vals).forEach((v, b) => { datasetPctAvg[b] += v; });
});
for (let b = 0; b < PCT_NUM_BINS; b++) datasetPctAvg[b] /= METHODS.length;

const proteomePctAvg = new Array(PCT_NUM_BINS).fill(0);
let methodsWithPctData = 0;
METHODS.forEach(method => {
    const dist = PROTEOME_PCT_DIST[method];
    if (!dist || dist.every(v => v === 0)) return;
    methodsWithPctData++;
    dist.forEach((v, b) => { if (b < PCT_NUM_BINS) proteomePctAvg[b] += v; });
});
if (methodsWithPctData > 0)
    for (let b = 0; b < PCT_NUM_BINS; b++) proteomePctAvg[b] /= methodsWithPctData;

const proteomePctTotal = proteomePctAvg.reduce((a, b) => a + b, 0);
const datasetPctTotal  = datasetPctAvg.reduce((a, b) => a + b, 0);

const datasetFractionAvg = new Array(NUM_BINS).fill(0);
METHODS.forEach(method => {
    const counts = getCounts(method);
    const imps   = counts.filter(c => c > 0);   // denominator: IMPs only, not whole dataset
    const total  = imps.length;
    if (!total) return;
    for (let b = 0; b < MAX_CAT; b++) datasetFractionAvg[b] += counts.filter(c => c === b + 1).length / total;
    datasetFractionAvg[NUM_BINS - 1] += counts.filter(c => c > MAX_CAT).length / total;
});
for (let b = 0; b < NUM_BINS; b++) datasetFractionAvg[b] /= METHODS.length;

const proteomeFractionAvg = new Array(NUM_BINS).fill(0);
let mwpd = 0;
METHODS.forEach(method => {
    const size = PROTEOME_SIZE[method];
    if (!size) return; mwpd++;
    PROTEOME_DIST[method].forEach((v, b) => { proteomeFractionAvg[b] += v / size; });
});
if (mwpd > 0) for (let b = 0; b < NUM_BINS; b++) proteomeFractionAvg[b] /= mwpd;

const datasetPctFrac  = datasetPctAvg.map(v  => datasetPctTotal  > 0 ? +(v / datasetPctTotal  * 100).toFixed(3) : 0);
const proteomePctFrac = proteomePctAvg.map(v => proteomePctTotal > 0 ? +(v / proteomePctTotal * 100).toFixed(3) : 0);

// ════════════════════════════════════════════════════════════════════════════
// 1 · SUMMARY TABLE
// ════════════════════════════════════════════════════════════════════════════

const allSizes = Object.values(PROTEOME_SIZE);
const noteEl   = document.getElementById('proteome-note');
if (allSizes.every(v => v === 0)) {
    noteEl.textContent = 'reference size not yet set';
} else if (allSizes.every(v => v === allSizes[0])) {
    noteEl.textContent = `of ${allSizes[0].toLocaleString()} reference proteins`;
} else {
    noteEl.textContent = 'of per-tool reference sizes (see tooltip)';
}

const tbody = document.getElementById('summary-tbody');
METHODS.forEach(method => {
    const counts = getCounts(method), n = counts.length;
    const withTM = counts.filter(c => c > 0).length;
    const proteomeSize = PROTEOME_SIZE[method];
    const pctProteome  = proteomeSize > 0 ? (withTM / proteomeSize * 100).toFixed(1) : '—';
    const tr = document.createElement('tr');
    tr.innerHTML = `
        <td><span class="method-dot ${DOT_CLASS[method]}"></span>${method}</td>
        <td>${withTM}</td>
        <td>${n > 0 ? (withTM / n * 100).toFixed(1) : '0.0'}%</td>
        <td title="${proteomeSize > 0 ? `Reference: ${proteomeSize.toLocaleString()} proteins` : 'Not set'}">${pctProteome}${pctProteome !== '—' ? '%' : ''}</td>
        <td>${Math.max(0, ...counts)}</td>
    `;
    tbody.appendChild(tr);
});

// ════════════════════════════════════════════════════════════════════════════
// 2 · DISTRIBUTION — RAW COUNTS
// ════════════════════════════════════════════════════════════════════════════

const distLabels = [];
for (let i = 1; i <= MAX_CAT; i++) distLabels.push(String(i));
distLabels.push(`${MAX_CAT}+`);

new Chart(document.getElementById('dist-chart'), {
    type: 'bar',
    data: {
        labels: distLabels,
        datasets: METHODS.map(method => {
            const counts = getCounts(method);
            const data = [];
            for (let i = 1; i <= MAX_CAT; i++) data.push(counts.filter(c => c === i).length);
            data.push(counts.filter(c => c > MAX_CAT).length);
            return { label: method, data, backgroundColor: COLORS[method] + 'cc', borderColor: COLORS[method], borderWidth: 1 };
        })
    },
    options: {
        responsive: true, maintainAspectRatio: false,
        plugins: {
            legend: { position: 'top' },
            tooltip: { callbacks: {
                title: ctx => `${ctx[0].label} TM helix${ctx[0].label === '1' ? '' : 'es'}`,
                label: ctx => ` ${ctx.dataset.label}: ${ctx.parsed.y} protein${ctx.parsed.y !== 1 ? 's' : ''}`,
            }},
        },
        scales: {
            x: { title: { display: true, text: 'Number of Transmembrane Helices', font: { weight: 'bold' } }, grid: { display: false } },
            y: { title: { display: true, text: 'Number of Proteins', font: { weight: 'bold' } }, beginAtZero: true, ticks: { stepSize: 1, precision: 0 }, grid: { color: 'rgba(0,0,0,0.06)' } },
        },
    },
    plugins: [{ id: 'barLabels', afterDatasetsDraw(chart) {
        const { ctx } = chart;
        chart.data.datasets.forEach((ds, i) => {
            chart.getDatasetMeta(i).data.forEach((bar, j) => {
                const val = ds.data[j];
                if (val > 0) { ctx.save(); ctx.fillStyle = '#333'; ctx.font = '10px Arial'; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom'; ctx.fillText(val, bar.x, bar.y - 2); ctx.restore(); }
            });
        });
    }}]
});

// ════════════════════════════════════════════════════════════════════════════
// 3+4 · TM HELIX COUNT — UNIFIED TOGGLED CHART (absolute ↔ fractions)
// ════════════════════════════════════════════════════════════════════════════

(function () {
    // Pre-format data arrays once
    const absDataset   = datasetCountAvg.map(v => +v.toFixed(2));
    const absProteome  = proteomeCountAvg.map(v => +v.toFixed(2));
    const fracDataset  = datasetFractionAvg.map(v => +(v * 100).toFixed(3));
    const fracProteome = proteomeFractionAvg.map(v => +(v * 100).toFixed(3));

    const BLUE = '#3498db', GREY = 'rgba(80,80,80,0.45)';

    const countChart = new Chart(document.getElementById('count-unified-chart'), {
        type: 'bar',
        data: {
            labels: distLabels,
            datasets: [
                { label: 'This dataset',   data: absDataset,  backgroundColor: BLUE + 'bb', borderColor: BLUE, borderWidth: 1 },
                { label: 'Human proteome', data: absProteome, backgroundColor: GREY,         borderColor: '#555', borderWidth: 1 },
            ]
        },
        options: {
            responsive: true, maintainAspectRatio: false,
            animation: { duration: 500, easing: 'easeInOutQuart' },
            plugins: { legend: { display: false }, tooltip: { callbacks: {
                title: ctx => `${ctx[0].label} TM helix${ctx[0].label === '1' ? '' : 'es'}`,
                label: ctx => ` ${ctx.dataset.label}: ${ctx.parsed.y.toFixed(2)} proteins (avg)`,
            }}},
            scales: {
                x: { title: { display: true, text: 'Number of Transmembrane Helices', font: { weight: 'bold' } }, grid: { display: false } },
                y: { title: { display: true, text: 'Number of Proteins (avg across methods)', font: { weight: 'bold' } }, beginAtZero: true, ticks: { precision: 1 }, grid: { color: 'rgba(0,0,0,0.06)' } },
            },
        },
    });

    // Tables
    const absTbody  = document.getElementById('count-abs-tbody');
    const fracTbody = document.getElementById('count-frac-tbody');

    distLabels.forEach((label, b) => {
        const pv = proteomeCountAvg[b], dv = datasetCountAvg[b];
        const cov = pv > 0 ? (dv / pv * 100).toFixed(1) : '—';
        const tr = document.createElement('tr');
        tr.innerHTML = `<td class="col-label">${label} TM helix${label === '1' ? '' : 'es'}</td><td>${pv.toFixed(1)}</td><td>${dv.toFixed(1)}</td><td>${cov}${cov !== '—' ? '%' : ''}</td>`;
        absTbody.appendChild(tr);
    });

    distLabels.forEach((label, b) => {
        const df = datasetFractionAvg[b] * 100, pf = proteomeFractionAvg[b] * 100;
        const ratio = pf > 0 ? (df / pf).toFixed(2) : '—';
        const tr = document.createElement('tr');
        tr.innerHTML = `<td class="col-label">${label} TM helix${label === '1' ? '' : 'es'}</td><td>${df.toFixed(2)}%</td><td>${pf.toFixed(2)}%</td><td>${ratio}</td>`;
        fracTbody.appendChild(tr);
    });

    // Toggle logic
    let countMode = 'abs'; // 'abs' | 'frac'
    const btnAbs  = document.getElementById('count-btn-abs');
    const btnFrac = document.getElementById('count-btn-frac');
    const titleEl = document.getElementById('count-unified-title');
    const noteEl2 = document.getElementById('count-unified-note');
    const absTable  = document.getElementById('count-abs-table');
    const fracTable = document.getElementById('count-frac-table');

    function setCountMode(mode) {
        countMode = mode;
        const isAbs = mode === 'abs';

        countChart.data.datasets[0].data = isAbs ? absDataset  : fracDataset;
        countChart.data.datasets[1].data = isAbs ? absProteome : fracProteome;

        countChart.options.scales.y.title.text = isAbs
            ? 'Number of Proteins (avg across methods)'
            : 'Fraction of IMP population (%)';
        countChart.options.plugins.tooltip.callbacks.label = isAbs
            ? (ctx => ` ${ctx.dataset.label}: ${ctx.parsed.y.toFixed(2)} proteins (avg)`)
            : (ctx => ` ${ctx.dataset.label}: ${ctx.parsed.y.toFixed(2)}%`);
        countChart.options.scales.y.ticks = isAbs
            ? { precision: 1 }
            : { callback: v => v + '%', precision: 2 };

        countChart.update();

        titleEl.textContent = isAbs
            ? 'Proteome Coverage by TM Helix Count'
            : 'TM Helix Count Distribution: Dataset vs. Human Proteome (Fractions)';
        noteEl2.textContent = isAbs
            ? 'Absolute protein counts averaged across all five methods. "Dataset coverage" column shows what fraction of each proteome bin is present in this dataset.'
            : 'Fraction of the IMP population in each TM helix count bin, averaged across all five methods. Ratio > 1 means the dataset over-represents that topology class.';

        absTable.style.display  = isAbs ? '' : 'none';
        fracTable.style.display = isAbs ? 'none' : '';

        btnAbs.classList.toggle('active',  isAbs);
        btnFrac.classList.toggle('active', !isAbs);
    }

    btnAbs.addEventListener('click',  () => setCountMode('abs'));
    btnFrac.addEventListener('click', () => setCountMode('frac'));

    setCountMode('abs'); // initial state
})();

// ════════════════════════════════════════════════════════════════════════════
// 5+6 · % SEQUENCE COVERAGE — UNIFIED TOGGLED CHART (absolute ↔ fractions)
// ════════════════════════════════════════════════════════════════════════════

(function () {
    const pctBinLabels = [];
    for (let lo = 0; lo < PCT_BIN_MAX; lo += BIN_SIZE) pctBinLabels.push(`${lo}–${lo + BIN_SIZE}`);
    pctBinLabels.push(`${PCT_BIN_MAX}+`);

    const absDataset   = datasetPctAvg.map(v => +v.toFixed(2));
    const absProteome  = proteomePctAvg.map(v => +v.toFixed(2));
    const fracDataset  = datasetPctFrac;
    const fracProteome = proteomePctFrac;

    const BLUE = '#3498db', GREY = 'rgba(80,80,80,0.45)';

    const pctChart = new Chart(document.getElementById('pct-unified-chart'), {
        type: 'bar',
        data: {
            labels: pctBinLabels,
            datasets: [
                { label: 'This dataset',   data: absDataset,  backgroundColor: BLUE + 'bb', borderColor: BLUE, borderWidth: 1 },
                { label: 'Human proteome', data: absProteome, backgroundColor: GREY,         borderColor: '#555', borderWidth: 1 },
            ]
        },
        options: {
            responsive: true, maintainAspectRatio: false,
            animation: { duration: 500, easing: 'easeInOutQuart' },
            plugins: { legend: { display: false }, tooltip: { callbacks: {
                title: ctx => `${ctx[0].label}% of sequence`,
                label: ctx => ` ${ctx.dataset.label}: ${ctx.parsed.y.toFixed(2)} proteins (avg)`,
            }}},
            scales: {
                x: { title: { display: true, text: '% of sequence occupied by TM helices', font: { weight: 'bold' } }, ticks: { maxRotation: 45, minRotation: 45 }, grid: { display: false } },
                y: { title: { display: true, text: 'Number of Proteins (avg across methods)', font: { weight: 'bold' } }, beginAtZero: true, ticks: { precision: 1 }, grid: { color: 'rgba(0,0,0,0.06)' } },
            },
        },
    });

    // Tables
    const absTbody  = document.getElementById('pct-abs-tbody');
    const fracTbody = document.getElementById('pct-frac-tbody');

    pctBinLabels.forEach((label, b) => {
        const pv = proteomePctAvg[b], dv = datasetPctAvg[b];
        const cov = pv > 0 ? (dv / pv * 100).toFixed(1) : '—';
        const tr = document.createElement('tr');
        tr.innerHTML = `<td class="col-label">${label}%</td><td>${pv.toFixed(1)}</td><td>${dv.toFixed(1)}</td><td>${cov}${cov !== '—' ? '%' : ''}</td>`;
        absTbody.appendChild(tr);
    });

    pctBinLabels.forEach((label, b) => {
        const df = datasetPctFrac[b], pf = proteomePctFrac[b];
        const ratio = pf > 0 ? (df / pf).toFixed(2) : '—';
        const tr = document.createElement('tr');
        tr.innerHTML = `<td class="col-label">${label}%</td><td>${df.toFixed(2)}%</td><td>${pf.toFixed(2)}%</td><td>${ratio}</td>`;
        fracTbody.appendChild(tr);
    });

    // Toggle logic
    let pctMode = 'abs';
    const btnAbs  = document.getElementById('pct-btn-abs');
    const btnFrac = document.getElementById('pct-btn-frac');
    const titleEl = document.getElementById('pct-unified-title');
    const noteEl2 = document.getElementById('pct-unified-note');
    const absTable  = document.getElementById('pct-abs-table');
    const fracTable = document.getElementById('pct-frac-table');

    function setPctMode(mode) {
        pctMode = mode;
        const isAbs = mode === 'abs';

        pctChart.data.datasets[0].data = isAbs ? absDataset  : fracDataset;
        pctChart.data.datasets[1].data = isAbs ? absProteome : fracProteome;

        pctChart.options.scales.y.title.text = isAbs
            ? 'Number of Proteins (avg across methods)'
            : 'Fraction of IMP population (%)';
        pctChart.options.plugins.tooltip.callbacks.label = isAbs
            ? (ctx => ` ${ctx.dataset.label}: ${ctx.parsed.y.toFixed(2)} proteins (avg)`)
            : (ctx => ` ${ctx.dataset.label}: ${ctx.parsed.y.toFixed(2)}% of IMPs`);
        pctChart.options.scales.y.ticks = isAbs
            ? { precision: 1 }
            : { callback: v => v + '%', precision: 2 };

        pctChart.update();

        titleEl.textContent = isAbs
            ? 'Distribution of TM Region Coverage (% of Sequence)'
            : 'TM Region Coverage Distribution: Dataset vs. Human Proteome (Fractions)';
        noteEl2.textContent = isAbs
            ? 'Absolute protein counts per % sequence coverage bin, averaged across all five methods. Only proteins with at least one predicted TM region are included.'
            : 'Fraction of the IMP population in each % coverage bin, averaged across all five methods. Ratio > 1 means the dataset over-represents that coverage class.';

        absTable.style.display  = isAbs ? '' : 'none';
        fracTable.style.display = isAbs ? 'none' : '';

        btnAbs.classList.toggle('active',  isAbs);
        btnFrac.classList.toggle('active', !isAbs);
    }

    btnAbs.addEventListener('click',  () => setPctMode('abs'));
    btnFrac.addEventListener('click', () => setPctMode('frac'));

    setPctMode('abs'); // initial state
})();

// ════════════════════════════════════════════════════════════════════════════
// 7 · AGREEMENT PIE  +  clickable protein list
// ════════════════════════════════════════════════════════════════════════════

const agreementCounts = new Array(6).fill(0);
const proteinAgreeLevel = RAW.map(r => {
    let n = 0;
    METHODS.forEach(m => { if (getCount(r, m) > 0) n++; });
    agreementCounts[n]++;
    return n;
});

const nTotal       = RAW.length;
const AGREE_LABELS = ['1 method', '2 methods', '3 methods', '4 methods', '5 methods (all)'];
const AGREE_COLORS = ['#e74c3c', '#e67e22', '#f39c12', '#2ecc71', '#27ae60'];
const agreeData    = agreementCounts.slice(1);
const nInPie       = agreeData.reduce((a, b) => a + b, 0);

new Chart(document.getElementById('agreement-pie'), {
    type: 'pie',
    data: { labels: AGREE_LABELS, datasets: [{ data: agreeData, backgroundColor: AGREE_COLORS, borderColor: '#fff', borderWidth: 3 }] },
    options: {
        responsive: true, maintainAspectRatio: false,
        plugins: { legend: { display: false }, tooltip: { callbacks: {
            label: ctx => { const c = ctx.parsed, pct = nTotal > 0 ? (c / nTotal * 100).toFixed(1) : '0.0'; return ` ${c} protein${c !== 1 ? 's' : ''} (${pct}% of all)`; }
        }}},
        onClick(evt, elements) {
            if (!elements.length) return;
            const sliceIdx  = elements[0].index;
            const agreeLevel = sliceIdx + 1;
            if (_openSlice.agree === sliceIdx) { closePanel('agree'); return; }
            _openSlice.agree = sliceIdx;
            const label    = AGREE_LABELS[sliceIdx];
            const matching = RAW.filter((_, i) => proteinAgreeLevel[i] === agreeLevel)
                .map(r => ({ ...r, label: `${agreeLevel} method${agreeLevel !== 1 ? 's' : ''}` }));
            openPanel('agree', `${label} — ${matching.length} protein${matching.length !== 1 ? 's' : ''}`, matching);
        },
    },
    plugins: [makePieLabelsPlugin('agreePieLabels', () => nInPie)]
});

const agreeLegDiv = document.getElementById('agreement-legend');
AGREE_LABELS.forEach((lbl, i) => {
    const cnt = agreeData[i], pct = nTotal > 0 ? (cnt / nTotal * 100).toFixed(1) : '0.0';
    const row = document.createElement('div'); row.className = 'summary-row';
    row.innerHTML = `<span><span class="legend-dot" style="background:${AGREE_COLORS[i]}"></span>${lbl}</span><span>${cnt} (${pct}%)</span>`;
    agreeLegDiv.appendChild(row);
});
const aLeast1 = agreementCounts.slice(1).reduce((a, b) => a + b, 0);
const aLeast3 = agreementCounts.slice(3).reduce((a, b) => a + b, 0);
const aLeast4 = agreementCounts.slice(4).reduce((a, b) => a + b, 0);
const aAll5   = agreementCounts[5];
document.getElementById('agreement-stats').innerHTML = `
    <div class="summary-row"><span>Total proteins</span><span>${nTotal}</span></div>
    <div class="summary-row"><span>Not classified by any method</span><span>${agreementCounts[0]}</span></div>
    <div class="summary-row"><span>≥1 method</span><span>${aLeast1} (${(aLeast1 / nTotal * 100).toFixed(1)}%)</span></div>
    <div class="summary-row"><span>≥3 methods (majority)</span><span>${aLeast3} (${(aLeast3 / nTotal * 100).toFixed(1)}%)</span></div>
    <div class="summary-row"><span>≥4 methods (strong)</span><span>${aLeast4} (${(aLeast4 / nTotal * 100).toFixed(1)}%)</span></div>
    <div class="summary-row"><span>All 5 methods</span><span>${aAll5} (${(aAll5 / nTotal * 100).toFixed(1)}%)</span></div>
`;

// ════════════════════════════════════════════════════════════════════════════
// 8 · CONSENSUS PIE  +  clickable protein list
// ════════════════════════════════════════════════════════════════════════════

// Full classification scheme:
//   modeFreq 5 → Full Consensus
//   modeFreq 4 → High Confidence
//   modeFreq 3 → Moderate Confidence
//   modeFreq 2 → Low Confidence      (two methods agree on one count, rest differ)
//   modeFreq 1 → No Consensus        (all five methods predict different counts)
//
// "No Consensus" is theoretically possible but rare; it only appears in the pie
// and summary box if at least one protein falls into it.

// All five levels defined upfront; No Consensus appended only if present
const CONS_ALL_LABELS = ['Full Consensus', 'High Confidence', 'Moderate Confidence', 'Low Confidence', 'No Consensus'];
const CONS_ALL_COLORS = ['#27ae60', '#2ecc71', '#f39c12', '#e74c3c', '#95a5a6'];
const consCounts = Object.fromEntries(CONS_ALL_LABELS.map(l => [l, 0]));

const proteinConsLevel = RAW.map(r => {
    const preds = METHODS.map(m => getCount(r, m));
    const freq  = {};
    preds.forEach(p => { freq[p] = (freq[p] || 0) + 1; });
    let modeVal = preds[0], modeFreq = 0;
    Object.entries(freq).forEach(([v, c]) => {
        if (c > modeFreq || (c === modeFreq && Number(v) > modeVal)) { modeFreq = c; modeVal = Number(v); }
    });
    if (modeVal === 0) return null; // non-membrane

    let level;
    if      (modeFreq === 5) level = 'Full Consensus';
    else if (modeFreq === 4) level = 'High Confidence';
    else if (modeFreq === 3) level = 'Moderate Confidence';
    else if (modeFreq === 2) level = 'Low Confidence';
    else                     level = 'No Consensus';   // modeFreq === 1: every method differs

    consCounts[level]++;
    return level;
});

// Build the pie arrays — include "No Consensus" only when it has proteins
const hasNoConsensus  = consCounts['No Consensus'] > 0;
const CONS_ORDER  = hasNoConsensus ? CONS_ALL_LABELS : CONS_ALL_LABELS.slice(0, 4);
const CONS_COLORS = hasNoConsensus ? CONS_ALL_COLORS : CONS_ALL_COLORS.slice(0, 4);
const consData    = CONS_ORDER.map(l => consCounts[l]);
const nMembrane   = consData.reduce((a, b) => a + b, 0);

new Chart(document.getElementById('consensus-pie'), {
    type: 'pie',
    data: { labels: CONS_ORDER, datasets: [{ data: consData, backgroundColor: CONS_COLORS, borderColor: '#fff', borderWidth: 3 }] },
    options: {
        responsive: true, maintainAspectRatio: false,
        plugins: { legend: { display: false }, tooltip: { callbacks: {
            label: ctx => { const c = ctx.parsed, pct = nMembrane > 0 ? (c / nMembrane * 100).toFixed(1) : '0.0'; return ` ${c} protein${c !== 1 ? 's' : ''} (${pct}%)`; }
        }}},
        onClick(evt, elements) {
            if (!elements.length) return;
            const sliceIdx = elements[0].index;
            const level    = CONS_ORDER[sliceIdx];
            if (_openSlice.cons === sliceIdx) { closePanel('cons'); return; }
            _openSlice.cons = sliceIdx;
            const matching = RAW.filter((_, i) => proteinConsLevel[i] === level)
                .map(r => ({ ...r, label: level }));
            openPanel('cons', `${level} — ${matching.length} protein${matching.length !== 1 ? 's' : ''}`, matching);
        },
    },
    plugins: [makePieLabelsPlugin('consPieLabels', () => nMembrane)]
});

const consLegDiv = document.getElementById('consensus-legend');
CONS_ORDER.forEach((lbl, i) => {
    const cnt = consCounts[lbl], pct = nMembrane > 0 ? (cnt / nMembrane * 100).toFixed(1) : '0.0';
    const row = document.createElement('div'); row.className = 'summary-row';
    row.innerHTML = `<span><span class="legend-dot" style="background:${CONS_COLORS[i]}"></span>${lbl}</span><span>${cnt} (${pct}%)</span>`;
    consLegDiv.appendChild(row);
});

const hOrP = consCounts['Full Consensus'] + consCounts['High Confidence'];
const noConsLine = hasNoConsensus
    ? `<div class="summary-row"><span>No Consensus</span><span>${consCounts['No Consensus']} (${(consCounts['No Consensus'] / nMembrane * 100).toFixed(1)}%)</span></div>`
    : '';

document.getElementById('consensus-stats').innerHTML = `
    <div class="summary-row"><span>Membrane proteins total</span><span>${nMembrane}</span></div>
    <div class="summary-row"><span>Full + High Confidence</span><span>${hOrP} (${nMembrane > 0 ? (hOrP / nMembrane * 100).toFixed(1) : '0.0'}%)</span></div>
    <div class="summary-row"><span>Moderate Confidence</span><span>${consCounts['Moderate Confidence']} (${nMembrane > 0 ? (consCounts['Moderate Confidence'] / nMembrane * 100).toFixed(1) : '0.0'}%)</span></div>
    <div class="summary-row"><span>Low Confidence</span><span>${consCounts['Low Confidence']} (${nMembrane > 0 ? (consCounts['Low Confidence'] / nMembrane * 100).toFixed(1) : '0.0'}%)</span></div>
    ${noConsLine}
`;

// ════════════════════════════════════════════════════════════════════════════
// 9 · PER-METHOD % COVERAGE SMALL MULTIPLES
// ════════════════════════════════════════════════════════════════════════════
// Uses a fixed 70% upper bound (14 bins at 5 pp) so the x-axis stays
// consistent across all five panels regardless of the PROTEOME_PCT_DIST bin count.

const HIST_BIN_SIZE = 5;
const HIST_BIN_MAX  = 70;
const HIST_NUM_BINS = HIST_BIN_MAX / HIST_BIN_SIZE + 1;

const histBinLabels = [];
for (let lo = 0; lo < HIST_BIN_MAX; lo += HIST_BIN_SIZE) histBinLabels.push(`${lo}–${lo + HIST_BIN_SIZE}`);
histBinLabels.push(`${HIST_BIN_MAX}+`);

function buildHistBins(values) {
    const bins = new Array(HIST_NUM_BINS).fill(0);
    values.forEach(v => {
        if (v >= HIST_BIN_MAX) bins[HIST_NUM_BINS - 1]++;
        else bins[Math.min(Math.floor(v / HIST_BIN_SIZE), HIST_NUM_BINS - 2)]++;
    });
    return bins;
}

const methodPcts = {};
METHODS.forEach(method => { methodPcts[method] = RAW.map(r => getPct(r, method)).filter(v => v !== null && v > 0); });

const allPcts    = Object.values(methodPcts).flat();
const pooledNorm = buildHistBins(allPcts).map(v => v / METHODS.length);
const histYMax   = Math.ceil(Math.max(0, ...METHODS.map(m => Math.max(0, ...buildHistBins(methodPcts[m])))) * 1.15);

const histGrid = document.getElementById('method-hist-grid');
METHODS.forEach(method => {
    const vals = methodPcts[method], bins = buildHistBins(vals), n = vals.length;
    const mean   = n > 0 ? (vals.reduce((a, b) => a + b, 0) / n).toFixed(1) : '—';
    const median = n > 0 ? (() => {
        const sorted = [...vals].sort((a, b) => a - b);
        const mid = Math.floor(sorted.length / 2);
        return (sorted.length % 2 === 0 ? (sorted[mid - 1] + sorted[mid]) / 2 : sorted[mid]).toFixed(1);
    })() : '—';

    const panel = document.createElement('div');
    panel.className = 'method-hist-panel';
    panel.innerHTML = `
        <div class="method-hist-title"><span class="h-dot" style="background:${COLORS[method]}"></span>${method}</div>
        <div class="method-hist-canvas"><canvas id="hist-${method}"></canvas></div>
        <div class="method-hist-stats">n=${n} &nbsp;·&nbsp; mean=${mean}% &nbsp;·&nbsp; median=${median}%</div>
    `;
    histGrid.appendChild(panel);

    new Chart(document.getElementById(`hist-${method}`), {
        type: 'bar',
        data: {
            labels: histBinLabels,
            datasets: [
                { label: 'All methods (avg)', data: pooledNorm, backgroundColor: 'rgba(180,180,180,0.35)', borderColor: 'rgba(150,150,150,0.5)', borderWidth: 1, barPercentage: 1.0, categoryPercentage: 1.0, order: 2 },
                { label: method, data: bins, backgroundColor: COLORS[method] + 'bb', borderColor: COLORS[method], borderWidth: 1, barPercentage: 0.75, categoryPercentage: 1.0, order: 1 },
            ]
        },
        options: {
            responsive: true, maintainAspectRatio: false,
            plugins: {
                legend: { display: false },
                tooltip: { callbacks: {
                    title: ctx => `${ctx[0].label}%`,
                    label: ctx => ctx.datasetIndex === 0
                        ? ` All-method avg: ${ctx.parsed.y.toFixed(1)}`
                        : `${method}: ${ctx.parsed.y} protein${ctx.parsed.y !== 1 ? 's' : ''}`,
                }},
            },
            scales: {
                x: { title: { display: true, text: '% of sequence', font: { size: 10 }, color: '#888' }, ticks: { font: { size: 9 }, maxRotation: 45, minRotation: 45 }, grid: { display: false } },
                y: { title: { display: true, text: 'Proteins', font: { size: 10 }, color: '#888' }, beginAtZero: true, max: histYMax, ticks: { font: { size: 9 }, precision: 0, stepSize: Math.max(1, Math.floor(histYMax / 5)) }, grid: { color: 'rgba(0,0,0,0.05)' } },
            }
        }
    });
});
