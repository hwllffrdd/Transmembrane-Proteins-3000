// Logic for templates/batch_stats.html. RAW is injected by the server
// as a JSON data island (see #batch-stats-raw-data in that template) rather
// than inlined here, since this is a static file Flask/Jinja doesn't process.
const RAW = JSON.parse(document.getElementById('batch-stats-raw-data').textContent);

// Chart.js defaults all of its own text (axis titles, tick labels, native
// legends) to a mid-grey. Applied globally, before any chart is built, so
// every chart on this page uses plain black instead without needing a
// color override repeated in each chart's own options.
Chart.defaults.color = '#000';

// ════════════════════════════════════════════════════════════════════════
// HARDCODED REFERENCE DATA — update after each database refresh
// ════════════════════════════════════════════════════════════════════════

const PROTEOME_SIZE_UniProt   = 5215;
const PROTEOME_SIZE_Phobius   = 5428;
const PROTEOME_SIZE_DeepTMHMM = 4955;
const PROTEOME_SIZE_TMBED     = 4973;
const PROTEOME_SIZE_TOPCONS   = 5633;

const PROTEOME_SIZE = {
  'UniProt':   PROTEOME_SIZE_UniProt,
  'Phobius':   PROTEOME_SIZE_Phobius,
  'DeepTMHMM': PROTEOME_SIZE_DeepTMHMM,
  'TMBED':     PROTEOME_SIZE_TMBED,
  'TOPCONS':   PROTEOME_SIZE_TOPCONS,
};

const PROTEOME_DIST_UniProt   = [2390, 323, 159, 407, 76, 276, 948, 83, 62, 125, 79, 192, 95];
const PROTEOME_DIST_Phobius   = [2501, 455, 223, 354, 138, 207, 858, 137, 88, 94, 100, 168, 105];
const PROTEOME_DIST_DeepTMHMM = [2156, 317, 121, 445, 53, 240, 926, 96, 55, 117, 69, 239, 121];
const PROTEOME_DIST_TMBED     = [2197, 278, 100, 422, 71, 269, 945, 69, 63, 126, 73, 230, 130];
const PROTEOME_DIST_TOPCONS   = [2441, 569, 253, 435, 151, 355, 695, 159, 90, 101, 92, 179, 114];

const PROTEOME_DIST = {
  'UniProt':   PROTEOME_DIST_UniProt,
  'Phobius':   PROTEOME_DIST_Phobius,
  'DeepTMHMM': PROTEOME_DIST_DeepTMHMM,
  'TMBED':     PROTEOME_DIST_TMBED,
  'TOPCONS':   PROTEOME_DIST_TOPCONS,
};

const PROTEOME_PCT_DIST_UniProt   = [1308, 846, 378, 341, 269, 246, 246, 325, 359, 634, 156, 54, 31, 12, 10];
const PROTEOME_PCT_DIST_Phobius   = [1364, 923, 446, 352, 290, 250, 249, 305, 310, 410, 373, 100, 38, 12, 6];
const PROTEOME_PCT_DIST_DeepTMHMM = [1333, 738, 339, 342, 249, 203, 236, 260, 320, 507, 328, 65, 26, 6, 3];
const PROTEOME_PCT_DIST_TMBED     = [1240, 702, 354, 280, 260, 225, 194, 267, 254, 364, 413, 352, 42, 19, 7];
const PROTEOME_PCT_DIST_TOPCONS   = [1428, 1040, 456, 400, 292, 265, 331, 343, 345, 463, 166, 55, 31, 7, 11];

const PROTEOME_PCT_DIST = {
  'UniProt':   PROTEOME_PCT_DIST_UniProt,
  'Phobius':   PROTEOME_PCT_DIST_Phobius,
  'DeepTMHMM': PROTEOME_PCT_DIST_DeepTMHMM,
  'TMBED':     PROTEOME_PCT_DIST_TMBED,
  'TOPCONS':   PROTEOME_PCT_DIST_TOPCONS,
};

// ════════════════════════════════════════════════════════════════════════
// END OF HARDCODED REFERENCE DATA
// ════════════════════════════════════════════════════════════════════════

const METHODS = ['UniProt', 'Phobius', 'DeepTMHMM', 'TMBED', 'TOPCONS'];
const COLORS  = {
  UniProt:   '#3498db',
  Phobius:   '#e74c3c',
  DeepTMHMM: '#2ecc71',
  TMBED:     '#9b59b6',
  TOPCONS:   '#f39c12',
};
const COUNT_KEY = {
  UniProt: 'uniprot_count', Phobius: 'phobius_count', DeepTMHMM: 'deeptmhmm_count',
  TMBED: 'tmbed_count', TOPCONS: 'topcons_count',
};
const PCT_KEY = {
  UniProt: 'uniprot_percentage', Phobius: 'phobius_percentage', DeepTMHMM: 'deeptmhmm_percentage',
  TMBED: 'tmbed_percentage', TOPCONS: 'topcons_percentage',
};

function getCount(r, m) {
  const v = r[COUNT_KEY[m]];
  return (v === null || v === undefined || v === '') ? 0 : parseInt(v, 10);
}
function getPct(r, m) {
  const v = r[PCT_KEY[m]];
  if (v === null || v === undefined || v === '') return null;
  const f = parseFloat(v); return isNaN(f) ? null : f;
}

// A slice's on-chart label needs a minimum wedge of angular (= pixel) room
// to avoid colliding with its neighbors' labels. Below this fraction of the
// full circle, the label is skipped entirely rather than drawn overlapping
// — the exact number is still available via hover tooltip and the
// legend/summary box, so nothing is actually lost.
const MIN_PIE_LABEL_FRACTION = 0.05;

function makePiePlugin(id, nRef, valuesRef) {
  return { id, afterDatasetsDraw(chart) {
    const { ctx } = chart, n = nRef();
    // valuesRef supplies the TRUE (un-inflated) values when the chart's own
    // data array has been adjusted for display — see inflateSmallPieSlices()
    // — so the printed percentage always reflects the real number, even
    // though the slice geometry (and thus arc.tooltipPosition()) is drawn
    // from the inflated array.
    const values = valuesRef ? valuesRef() : chart.data.datasets[0].data;
    // The label's available room is governed by the DRAWN (possibly
    // inflated) arc, not the true value, so collision-avoidance has to be
    // based on that same array.
    const displayData = chart.data.datasets[0].data;
    const displayTotal = displayData.reduce((a, b) => a + b, 0);
    chart.getDatasetMeta(0).data.forEach((arc, i) => {
      const v = values[i]; if (!v) return;
      const arcShare = displayTotal > 0 ? displayData[i] / displayTotal : 0;
      if (arcShare < MIN_PIE_LABEL_FRACTION) return;
      const pct = (v/n*100).toFixed(1), { x, y } = arc.tooltipPosition();
      ctx.save(); ctx.fillStyle='black'; ctx.font='bold 12px Work Sans, sans-serif';
      ctx.textAlign='center'; ctx.textBaseline='middle';
      ctx.fillText(`${pct}%`, x, y); ctx.restore();
    });
  }};
}

// Chart.js pie slices below ~0.6% of the total become too thin to reliably
// hover or click. For DISPLAY ONLY, boost any nonzero slice under that
// threshold up to 0.6% of the total, taking the difference proportionally
// from the slices that are comfortably above it, so the arcs still sum to
// a full circle. The real counts/percentages (tooltips, on-slice labels,
// legend, summary box) are computed from the original values, not this
// adjusted array — only the drawn/clickable geometry uses it. A slice at
// exactly 0 is left at 0, so it's still not drawn at all.
function inflateSmallPieSlices(values, minFraction) {
  minFraction = minFraction || 0.006;
  const total = values.reduce((a, b) => a + b, 0);
  if (total <= 0) return values.slice();
  const minValue = minFraction * total;

  const adjusted = values.slice();
  let boostAmount = 0;
  values.forEach((v, i) => {
    if (v > 0 && v < minValue) {
      boostAmount += minValue - v;
      adjusted[i] = minValue;
    }
  });
  if (boostAmount <= 0) return adjusted;

  const donors = values.map((v, i) => ({ v, i })).filter(({ v }) => v >= minValue);
  const donorTotal = donors.reduce((a, { v }) => a + v, 0);
  if (donorTotal > 0) {
    donors.forEach(({ v, i }) => {
      adjusted[i] = Math.max(0, v - (v / donorTotal) * boostAmount);
    });
  }
  return adjusted;
}

// ── Panel helpers ─────────────────────────────────────────────────────
const _openSlice = { agree: null, cons: null, dist: null, count: null, pct: null };
const PANEL_ROWS = { agree: [], cons: [], dist: [], count: [], pct: [] };

const VIZ_METHODS = [
  { key: 'uniprot',   label: 'UniProt',   cls: 'uniprot-region'   },
  { key: 'phobius',   label: 'Phobius',   cls: 'phobius-region'   },
  { key: 'deeptmhmm', label: 'DeepTMHMM', cls: 'deeptmhmm-region' },
  { key: 'tmbed',     label: 'TMBED',     cls: 'tmbed-region'     },
  { key: 'topcons',   label: 'TOPCONS',   cls: 'topcons-region'   },
];

function escapeHtml(value) {
  if (value === null || value === undefined) return '';
  return String(value)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#39;');
}

function panelRowHtml(row, i, panelKey) {
  const vizId = `${panelKey}-viz-${i}`;
  const uniprotUrl = 'https://www.uniprot.org/uniprotkb/' + encodeURIComponent(row.accession) + '/entry';
  return `
    <tr class="data-row" data-target="${vizId}" data-rendered="false" data-idx="${i}" data-panel="${panelKey}">
      <td><a href="${uniprotUrl}" target="_blank" rel="noopener noreferrer" onclick="event.stopPropagation()">${escapeHtml(row.accession)}</a></td>
      <td>${escapeHtml(row.gene_name)}</td>
      <td>${escapeHtml(row.name)}<span class="expand-indicator">&#9660; show visualisation</span></td>
      <td>${escapeHtml(row.label)}</td>
    </tr>
    <tr class="visualizer-row" id="${vizId}">
      <td colspan="4">
        <div class="viz-inner">
          <div class="protein-visualization" id="${vizId}-content"></div>
        </div>
      </td>
    </tr>`;
}

function renderPanelViz(dataRow) {
  const panelKey = dataRow.dataset.panel;
  const idx = parseInt(dataRow.dataset.idx, 10);
  const row = PANEL_ROWS[panelKey][idx];
  const container = document.getElementById(dataRow.dataset.target + '-content');
  const seqLen = row.sequence_length;

  if (!seqLen) {
    container.innerHTML = '<p class="hint">Sequence length unavailable — cannot render visualisation.</p>';
    return;
  }

  const regions = {
    uniprot:   row.transmembrane_regions,
    phobius:   row.phobius_regions,
    deeptmhmm: row.deeptmhmm_regions,
    tmbed:     row.tmbed_regions,
    topcons:   row.topcons_regions,
  };

  let html = `<div class="sequence-info">Protein length: ${seqLen} amino acids</div>`;
  VIZ_METHODS.forEach(method => {
    const methodRegions = regions[method.key] || [];
    let barsHtml = '';
    methodRegions.forEach(region => {
      const start = region[0], end = region[1];
      const left  = (start / seqLen * 100).toFixed(3);
      const width = ((end - start + 1) / seqLen * 100).toFixed(3);
      barsHtml += `<div class="tm-region ${method.cls}"
        style="left:${left}%;width:${width}%;"
        title="Position ${start}–${end} (${end - start + 1} aa)"></div>`;
    });
    html += `<div class="region-label ${method.key}-label">${method.label}</div>
             <div class="protein-bar">${barsHtml}</div>`;
  });
  container.innerHTML = html;
}

// Event delegation: rows are re-created every time a slice is clicked, so
// a single listener on each (persistent) panel tbody handles all of them.
function handlePanelTbodyClick(evt) {
  const dataRow = evt.target.closest('tr.data-row');
  if (!dataRow) return;

  const selection = window.getSelection();
  if (selection && selection.toString().length > 0 && dataRow.contains(selection.anchorNode)) return;

  const vizRow = document.getElementById(dataRow.dataset.target);
  const indicator = dataRow.querySelector('.expand-indicator');
  const isOpen = dataRow.classList.contains('open');

  if (isOpen) {
    dataRow.classList.remove('open');
    vizRow.classList.remove('open');
    indicator.innerHTML = '&#9660; show visualisation';
  } else {
    if (dataRow.dataset.rendered === 'false') {
      renderPanelViz(dataRow);
      dataRow.dataset.rendered = 'true';
    }
    dataRow.classList.add('open');
    vizRow.classList.add('open');
    indicator.innerHTML = '&#9650; hide visualisation';
  }
}
document.getElementById('agree-panel-tbody').addEventListener('click', handlePanelTbodyClick);
document.getElementById('cons-panel-tbody').addEventListener('click', handlePanelTbodyClick);
document.getElementById('dist-panel-tbody').addEventListener('click', handlePanelTbodyClick);
document.getElementById('count-panel-tbody').addEventListener('click', handlePanelTbodyClick);
document.getElementById('pct-panel-tbody').addEventListener('click', handlePanelTbodyClick);

function openPanel(panelKey, title, rows) {
  const panel=document.getElementById(`${panelKey}-protein-panel`);
  const titleEl=document.getElementById(`${panelKey}-panel-title`);
  const tbodyEl=document.getElementById(`${panelKey}-panel-tbody`);
  titleEl.textContent=title;
  PANEL_ROWS[panelKey] = rows;
  tbodyEl.innerHTML = rows.map((r, i) => panelRowHtml(r, i, panelKey)).join('');
  panel.classList.add('open');
  setTimeout(()=>panel.scrollIntoView({behavior:'smooth',block:'nearest'}),50);
}
function closePanel(k) { document.getElementById(`${k}-protein-panel`).classList.remove('open'); _openSlice[k]=null; }
window.closePanel = closePanel;

// ── Pre-computed averages ─────────────────────────────────────────────
const MAX_CAT=12, NUM_BINS=MAX_CAT+1;

const datasetCountAvg = new Array(NUM_BINS).fill(0);
METHODS.forEach(m => {
  const c=RAW.map(r=>getCount(r,m));
  for(let b=0;b<MAX_CAT;b++) datasetCountAvg[b]+=c.filter(x=>x===b+1).length;
  datasetCountAvg[NUM_BINS-1]+=c.filter(x=>x>MAX_CAT).length;
});
for(let b=0;b<NUM_BINS;b++) datasetCountAvg[b]/=METHODS.length;

const proteomeCountAvg = new Array(NUM_BINS).fill(0);
let mwcd=0;
METHODS.forEach(m => {
  if(!PROTEOME_SIZE[m]) return; mwcd++;
  PROTEOME_DIST[m].forEach((v,b)=>{ proteomeCountAvg[b]+=v; });
});
if(mwcd>0) for(let b=0;b<NUM_BINS;b++) proteomeCountAvg[b]/=mwcd;

const PCT_NUM_BINS = Math.max(...Object.values(PROTEOME_PCT_DIST).map(a=>a.length));
const PCT_BIN_SIZE=5, PCT_BIN_MAX=(PCT_NUM_BINS-1)*PCT_BIN_SIZE;

function pctBinIndexFor(v) {
  if (v >= PCT_BIN_MAX) return PCT_NUM_BINS - 1;
  return Math.min(Math.floor(v / PCT_BIN_SIZE), PCT_NUM_BINS - 2);
}

function buildPctBins(vals) {
  const bins=new Array(PCT_NUM_BINS).fill(0);
  vals.forEach(v => { bins[pctBinIndexFor(v)]++; });
  return bins;
}

const datasetPctAvg = new Array(PCT_NUM_BINS).fill(0);
METHODS.forEach(m => {
  const vals=RAW.map(r=>getPct(r,m)).filter(v=>v!==null&&v>0);
  buildPctBins(vals).forEach((v,b)=>{ datasetPctAvg[b]+=v; });
});
for(let b=0;b<PCT_NUM_BINS;b++) datasetPctAvg[b]/=METHODS.length;

const proteomePctAvg = new Array(PCT_NUM_BINS).fill(0);
let mwpd=0;
METHODS.forEach(m => {
  const dist=PROTEOME_PCT_DIST[m];
  if(!dist||dist.every(v=>v===0)) return; mwpd++;
  dist.forEach((v,b)=>{ if(b<PCT_NUM_BINS) proteomePctAvg[b]+=v; });
});
if(mwpd>0) for(let b=0;b<PCT_NUM_BINS;b++) proteomePctAvg[b]/=mwpd;

const proteomePctTotal=proteomePctAvg.reduce((a,b)=>a+b,0);
const datasetPctTotal =datasetPctAvg.reduce((a,b)=>a+b,0);
const datasetPctFrac  =datasetPctAvg.map(v=>datasetPctTotal>0?+(v/datasetPctTotal*100).toFixed(3):0);
const proteomePctFrac =proteomePctAvg.map(v=>proteomePctTotal>0?+(v/proteomePctTotal*100).toFixed(3):0);

const datasetFractionAvg=new Array(NUM_BINS).fill(0);
METHODS.forEach(m => {
  const c=RAW.map(r=>getCount(r,m));
  const tot=c.filter(x=>x>0).length;  // denominator: IMPs only, not whole dataset
  if(!tot) return;
  for(let b=0;b<MAX_CAT;b++) datasetFractionAvg[b]+=c.filter(x=>x===b+1).length/tot;
  datasetFractionAvg[NUM_BINS-1]+=c.filter(x=>x>MAX_CAT).length/tot;
});
for(let b=0;b<NUM_BINS;b++) datasetFractionAvg[b]/=METHODS.length;

const proteomeFractionAvg=new Array(NUM_BINS).fill(0);
let mwfd=0;
METHODS.forEach(m => {
  const sz=PROTEOME_SIZE[m]; if(!sz) return; mwfd++;
  PROTEOME_DIST[m].forEach((v,b)=>{ proteomeFractionAvg[b]+=v/sz; });
});
if(mwfd>0) for(let b=0;b<NUM_BINS;b++) proteomeFractionAvg[b]/=mwfd;

// Standard, vivid blue (Bootstrap's own "primary") rather than the muted
// navy this used to be — kept distinct from COLORS.UniProt (#3498db) since
// these bars mean "this dataset, averaged across all 5 methods", not the
// UniProt method specifically, and reusing that exact hue would blur that
// distinction. The explanatory <p class="chart-note"> next to each chart
// colors its "blue bar" text with this same DIST_BLUE, so they stay in sync.
const DIST_BLUE  = '#0d6efd';
const DIST_GREY  = 'rgba(80,80,80,0.45)';
const BORDER_BLUE = DIST_BLUE;
const BORDER_GREY = '#6a5959';

// 1 · Summary table
const allSizes=Object.values(PROTEOME_SIZE);
const noteEl=document.getElementById('proteome-note');
if(allSizes.every(v=>v===0)) noteEl.textContent='reference size not yet set';
else if(allSizes.every(v=>v===allSizes[0])) noteEl.textContent=`of ${allSizes[0].toLocaleString()} reference proteins`;
else noteEl.textContent='of reference proteomes';

const tbody=document.getElementById('summary-tbody');
METHODS.forEach(m => {
  const counts=RAW.map(r=>getCount(r,m)), n=counts.length, withTM=counts.filter(c=>c>0).length;
  const sz=PROTEOME_SIZE[m], pct=sz>0?(withTM/sz*100).toFixed(1):'—';
  const tr=document.createElement('tr');
  tr.innerHTML=`
    <td><span class="legend-dot" style="background:${COLORS[m]}"></span>${m}</td>
    <td>${withTM}</td>
    <td>${n>0?(withTM/n*100).toFixed(1):'0.0'}%</td>
    <td title="${sz>0?`Reference: ${sz.toLocaleString()} proteins`:'Not set'}">${pct}${pct!=='—'?'%':''}</td>
    <td>${Math.max(0,...counts)}</td>`;
  tbody.appendChild(tr);
});

// 2 · Distribution bar chart
const distLabels=[];
for(let i=1;i<=MAX_CAT;i++) distLabels.push(String(i));
distLabels.push(`${MAX_CAT}+`);

const distChart=new Chart(document.getElementById('dist-chart'), {
  type:'bar',
  data:{ labels:distLabels, datasets:METHODS.map(m=>{
    const c=RAW.map(r=>getCount(r,m)), d=[];
    for(let i=1;i<=MAX_CAT;i++) d.push(c.filter(x=>x===i).length);
    d.push(c.filter(x=>x>MAX_CAT).length);
    return{label:m,data:d,backgroundColor:COLORS[m]+'cc',borderColor:COLORS[m],borderWidth:1};
  })},
  options:{
    responsive:true,maintainAspectRatio:false,
    // On screen the legend is the HTML checkbox list (#dist-legend) built
    // below; this canvas legend is switched on only while exporting a PNG
    // (see data-export-legend), listing just the methods currently shown.
    plugins:{legend:{display:false,position:'top',labels:{filter:item=>!item.hidden}},tooltip:{callbacks:{
      title:ctx=>`${ctx[0].label} TM helix${ctx[0].label==='1'?'':'es'}`,
      label:ctx=>` ${ctx.dataset.label}: ${ctx.parsed.y} protein${ctx.parsed.y!==1?'s':''}`,
    }}},
    scales:{
      x:{title:{display:true,text:'Number of Transmembrane Helices',font:{weight:'bold'}},grid:{display:false}},
      y:{title:{display:true,text:'Number of Proteins',font:{weight:'bold'}},beginAtZero:true,ticks:{precision:0},grid:{color:'rgba(0,0,0,0.06)'}},
    },
    onClick(evt, elements) {
      if (!elements.length) return;
      const { datasetIndex, index: bin } = elements[0];
      const method = METHODS[datasetIndex];
      const isOverflow = bin === MAX_CAT;
      const exactCount = bin + 1;
      const key = `${datasetIndex}-${bin}`;
      if (_openSlice.dist === key) { closePanel('dist'); return; }
      _openSlice.dist = key;
      const rows = RAW.filter(r => {
        const v = getCount(r, method);
        return isOverflow ? v > MAX_CAT : v === exactCount;
      }).map(r => ({ ...r, label: String(getCount(r, method)) }));
      const binLabel = distLabels[bin];
      openPanel('dist',
        `${method} — ${binLabel} TM helix${binLabel === '1' ? '' : 'es'} — ${rows.length} protein${rows.length !== 1 ? 's' : ''}`,
        rows);
    },
  },
  plugins:[{id:'barLabels',afterDatasetsDraw(chart){
    const{ctx}=chart;
    chart.data.datasets.forEach((ds,i)=>{
      if(!chart.isDatasetVisible(i)) return; // skip methods hidden via the legend
      chart.getDatasetMeta(i).data.forEach((bar,j)=>{
        const v=ds.data[j]; if(v>0){ctx.save();ctx.fillStyle='#333';ctx.font='10px sans-serif';ctx.textAlign='center';ctx.textBaseline='bottom';ctx.fillText(v,bar.x,bar.y-2);ctx.restore();}
      });
    });
  }}]
});

// Checkbox legend: one entry per method, ticked = shown in the chart
(function(){
  const legend=document.getElementById('dist-legend');
  distChart.data.datasets.forEach((ds,i)=>{
    const item=document.createElement('label');
    item.className='chart-legend-item chart-legend-toggle';
    item.innerHTML=`<input type="checkbox" checked>`+
      `<span class="chart-legend-swatch" style="border:1px solid ${ds.borderColor}; background:${ds.backgroundColor};"></span>`+
      `<span>${escapeHtml(ds.label)}</span>`;
    item.querySelector('input').addEventListener('change',e=>{
      distChart.setDatasetVisibility(i,e.target.checked);
      distChart.update();
    });
    legend.appendChild(item);
  });
})();

// The toggled charts below switch the y axis between absolute counts and
// percentages, whose tick labels differ in width ("1,250" vs "12.5%"). Chart.js
// sizes the axis to its labels, so every toggle nudged the plot area (and the
// bars) sideways. This afterFit hook pins the axis to the widest width seen so
// far; each chart's setMode() is run through both modes at creation, so the
// width is settled before the user ever clicks a toggle.
function stableYAxisWidth() {
  let maxWidth = 0;
  return scale => {
    maxWidth = Math.max(maxWidth, scale.width);
    scale.width = maxWidth;
  };
}

// 3+4 · TM helix count — unified toggled chart
(function(){
  const absDs  = datasetCountAvg.map(v=>+v.toFixed(2));
  const absPr  = proteomeCountAvg.map(v=>+v.toFixed(2));
  const fracDs = datasetFractionAvg.map(v=>+(v*100).toFixed(3));
  const fracPr = proteomeFractionAvg.map(v=>+(v*100).toFixed(3));
  const BLUE=DIST_BLUE+'cc', GREY=DIST_GREY;

  const ch=new Chart(document.getElementById('count-unified-chart'),{
    type:'bar',
    data:{labels:distLabels,datasets:[
      {label:'This dataset',  data:absDs, backgroundColor:BLUE,borderColor:BORDER_BLUE,borderWidth:1},
      {label:'Human proteome',data:absPr, backgroundColor:GREY,borderColor:BORDER_GREY,borderWidth:1},
    ]},
    options:{
      responsive:true,maintainAspectRatio:false,
      animation:{duration:500,easing:'easeInOutQuart'},
      plugins:{legend:{display:false},tooltip:{callbacks:{
        title:ctx=>`${ctx[0].label} TM helix${ctx[0].label==='1'?'':'es'}`,
        label:ctx=>` ${ctx.dataset.label}: ${ctx.parsed.y.toFixed(2)} proteins (avg)`,
      }}},
      scales:{
        x:{title:{display:true,text:'Number of Transmembrane Helices',font:{weight:'bold'}},grid:{display:false}},
        y:{title:{display:true,text:'Number of Proteins (avg across methods)',font:{weight:'bold'}},beginAtZero:true,ticks:{precision:1},grid:{color:'rgba(0,0,0,0.06)'},afterFit:stableYAxisWidth()},
      },
      onHover(evt, elements) {
        evt.native.target.style.cursor = (elements.length && elements[0].datasetIndex === 0) ? 'pointer' : 'default';
      },
      onClick(evt, elements) {
        if (!elements.length) return;
        const { datasetIndex, index: bin } = elements[0];
        if (datasetIndex !== 0) return; // "Human proteome" — reference data, no per-protein rows
        const key = `${bin}`;
        if (_openSlice.count === key) { closePanel('count'); return; }
        _openSlice.count = key;
        const isOverflow = bin === MAX_CAT;
        const exactCount = bin + 1;
        const rows = [];
        METHODS.forEach(m => {
          RAW.forEach(r => {
            const v = getCount(r, m);
            const matches = isOverflow ? v > MAX_CAT : v === exactCount;
            if (matches) rows.push({ ...r, label: m });
          });
        });
        const binLabel = distLabels[bin];
        openPanel('count',
          `Your dataset — proteins predicted to contain ${binLabel} TM helix${binLabel === '1' ? '' : 'es'}`,
          rows);
      },
    }
  });

  const absTb=document.getElementById('count-abs-tbody');
  const frTb =document.getElementById('count-frac-tbody');
  distLabels.forEach((lbl,b)=>{
    const pv=proteomeCountAvg[b],dv=datasetCountAvg[b],cov=pv>0?(dv/pv*100).toFixed(1):'—';
    const tr=document.createElement('tr');
    tr.innerHTML=`<td class="col-label">${lbl} TM helix${lbl==='1'?'':'es'}</td><td>${pv.toFixed(1)}</td><td>${dv.toFixed(1)}</td><td>${cov}${cov!=='—'?'%':''}</td>`;
    absTb.appendChild(tr);
  });
  distLabels.forEach((lbl,b)=>{
    const df=datasetFractionAvg[b]*100,pf=proteomeFractionAvg[b]*100,ratio=pf>0?(df/pf).toFixed(2):'—';
    const tr=document.createElement('tr');
    tr.innerHTML=`<td class="col-label">${lbl} TM helix${lbl==='1'?'':'es'}</td><td>${df.toFixed(2)}%</td><td>${pf.toFixed(2)}%</td><td>${ratio}</td>`;
    frTb.appendChild(tr);
  });

  const btnA=document.getElementById('count-btn-abs'),btnF=document.getElementById('count-btn-frac');
  const descEl=document.getElementById('count-unified-desc'),noteEl=document.getElementById('count-unified-note');
  const tabA=document.getElementById('count-abs-table'),tabF=document.getElementById('count-frac-table');

  function setMode(isAbs){
    ch.data.datasets[0].data=isAbs?absDs:fracDs;
    ch.data.datasets[1].data=isAbs?absPr:fracPr;
    ch.options.scales.y.title.text=isAbs?'Number of Proteins (avg across methods)':'Fraction of IMP population (%)';
    ch.options.plugins.tooltip.callbacks.label=isAbs
      ?(ctx=>` ${ctx.dataset.label}: ${ctx.parsed.y.toFixed(2)} proteins (avg)`)
      :(ctx=>` ${ctx.dataset.label}: ${ctx.parsed.y.toFixed(2)}%`);
    ch.options.scales.y.ticks=isAbs?{precision:1}:{callback:v=>v+'%',precision:2};
    ch.update();
    descEl.textContent=isAbs
      ?'Total counts of proteins predicted to contain specific number of TM segments (avg numbers for all prediction tools are shown)'
      :'Fractions of IMPs containing specific numbers of TM segments (expressed as % of all IMPs in your dataset or whole human proteome)';
    noteEl.textContent=isAbs
      ?'Absolute protein counts averaged across all five methods. "Dataset coverage" column shows what fraction of each proteome bin is present in this dataset.'
      :'Fraction of the IMP population in each TM helix count bin, averaged across all five methods. Ratio > 1 means the dataset over-represents that topology class.';
    tabA.style.display=isAbs?'':'none'; tabF.style.display=isAbs?'none':'';
    btnA.classList.toggle('active',isAbs); btnF.classList.toggle('active',!isAbs);
  }
  btnA.addEventListener('click',()=>setMode(true));
  btnF.addEventListener('click',()=>setMode(false));
  setMode(false); setMode(true);   // lay out both modes once so stableYAxisWidth() knows the widest axis
})();

// 5+6 · % sequence coverage — unified toggled chart
(function(){
  const pctBinLabels=[];
  for(let lo=0;lo<PCT_BIN_MAX;lo+=PCT_BIN_SIZE) pctBinLabels.push(`${lo}–${lo+PCT_BIN_SIZE}`);
  pctBinLabels.push(`${PCT_BIN_MAX}+`);

  const absDs  = datasetPctAvg.map(v=>+v.toFixed(2));
  const absPr  = proteomePctAvg.map(v=>+v.toFixed(2));
  const fracDs = datasetPctFrac;
  const fracPr = proteomePctFrac;
  const BLUE=DIST_BLUE+'cc', GREY=DIST_GREY;

  const ch=new Chart(document.getElementById('pct-unified-chart'),{
    type:'bar',
    data:{labels:pctBinLabels,datasets:[
      {label:'This dataset',  data:absDs, backgroundColor:BLUE,borderColor:BORDER_BLUE,borderWidth:1},
      {label:'Human proteome',data:absPr, backgroundColor:GREY,borderColor:BORDER_GREY,borderWidth:1},
    ]},
    options:{
      responsive:true,maintainAspectRatio:false,
      animation:{duration:500,easing:'easeInOutQuart'},
      plugins:{legend:{display:false},tooltip:{callbacks:{
        title:ctx=>`${ctx[0].label}% of sequence`,
        label:ctx=>` ${ctx.dataset.label}: ${ctx.parsed.y.toFixed(2)} proteins (avg)`,
      }}},
      scales:{
        x:{title:{display:true,text:'% of sequence occupied by TM helices',font:{weight:'bold'}},ticks:{maxRotation:45,minRotation:45},grid:{display:false}},
        y:{title:{display:true,text:'Number of Proteins (avg across methods)',font:{weight:'bold'}},beginAtZero:true,ticks:{precision:1},grid:{color:'rgba(0,0,0,0.06)'},afterFit:stableYAxisWidth()},
      },
      onHover(evt, elements) {
        evt.native.target.style.cursor = (elements.length && elements[0].datasetIndex === 0) ? 'pointer' : 'default';
      },
      onClick(evt, elements) {
        if (!elements.length) return;
        const { datasetIndex, index: bin } = elements[0];
        if (datasetIndex !== 0) return; // "Human proteome" — reference data, no per-protein rows
        const key = `${bin}`;
        if (_openSlice.pct === key) { closePanel('pct'); return; }
        _openSlice.pct = key;
        const rows = [];
        METHODS.forEach(m => {
          RAW.forEach(r => {
            const v = getPct(r, m);
            if (v === null || v <= 0) return;
            if (pctBinIndexFor(v) === bin) rows.push({ ...r, label: m });
          });
        });
        const binLabel = pctBinLabels[bin];
        openPanel('pct',
          `Your dataset — proteins predicted to contain ${binLabel}% of sequence inside the membrane`,
          rows);
      },
    }
  });

  const absTb=document.getElementById('pct-abs-tbody');
  const frTb =document.getElementById('pct-frac-tbody');
  pctBinLabels.forEach((lbl,b)=>{
    const pv=proteomePctAvg[b],dv=datasetPctAvg[b],cov=pv>0?(dv/pv*100).toFixed(1):'—';
    const tr=document.createElement('tr');
    tr.innerHTML=`<td class="col-label">${lbl}%</td><td>${pv.toFixed(1)}</td><td>${dv.toFixed(1)}</td><td>${cov}${cov!=='—'?'%':''}</td>`;
    absTb.appendChild(tr);
  });
  pctBinLabels.forEach((lbl,b)=>{
    const df=datasetPctFrac[b],pf=proteomePctFrac[b],ratio=pf>0?(df/pf).toFixed(2):'—';
    const tr=document.createElement('tr');
    tr.innerHTML=`<td class="col-label">${lbl}%</td><td>${df.toFixed(2)}%</td><td>${pf.toFixed(2)}%</td><td>${ratio}</td>`;
    frTb.appendChild(tr);
  });

  const btnA=document.getElementById('pct-btn-abs'),btnF=document.getElementById('pct-btn-frac');
  const descEl=document.getElementById('pct-unified-desc'),noteEl=document.getElementById('pct-unified-note');
  const tabA=document.getElementById('pct-abs-table'),tabF=document.getElementById('pct-frac-table');

  function setMode(isAbs){
    ch.data.datasets[0].data=isAbs?absDs:fracDs;
    ch.data.datasets[1].data=isAbs?absPr:fracPr;
    ch.options.scales.y.title.text=isAbs?'Number of Proteins (avg across methods)':'Fraction of IMP population (%)';
    ch.options.plugins.tooltip.callbacks.label=isAbs
      ?(ctx=>` ${ctx.dataset.label}: ${ctx.parsed.y.toFixed(2)} proteins (avg)`)
      :(ctx=>` ${ctx.dataset.label}: ${ctx.parsed.y.toFixed(2)}% of IMPs`);
    ch.options.scales.y.ticks=isAbs?{precision:1}:{callback:v=>v+'%',precision:2};
    ch.update();
    descEl.textContent=isAbs
      ?'Percentage of IMP protein aminoacid sequence predicted to be embedded in membrane in your dataset and whole human proteome (amionoacids in TM segments/total aminoacids). (Absolute counts of proteins are shown)'
      :'Percentage of IMP protein aminoacid sequence predicted to be embedded in membrane. Here shown as a fraction of all predicted IMPs in your dataset and whole human proteome (amionoacids in TM segments/total aminoacids)';
    noteEl.textContent=isAbs
      ?'Absolute protein counts per % sequence coverage bin, averaged across all five methods. Only proteins with at least one predicted TM region are included.'
      :'Fraction of the IMP population in each % coverage bin, averaged across all five methods. Ratio > 1 means the dataset over-represents that coverage class.';
    tabA.style.display=isAbs?'':'none'; tabF.style.display=isAbs?'none':'';
    btnA.classList.toggle('active',isAbs); btnF.classList.toggle('active',!isAbs);
  }
  btnA.addEventListener('click',()=>setMode(true));
  btnF.addEventListener('click',()=>setMode(false));
  setMode(false); setMode(true);   // lay out both modes once so stableYAxisWidth() knows the widest axis
})();

// 7 · Agreement pie + panel
const agreeCounts=new Array(6).fill(0);
const proteinAgreeLevel=RAW.map(r=>{
  let n=0; METHODS.forEach(m=>{if(getCount(r,m)>0)n++;}); agreeCounts[n]++; return n;
});
const nTotal=RAW.length;
const AGREE_LABELS=['1 method','2 methods','3 methods','4 methods','5 methods (all)'];
const AGREE_COLORS=['#e74c3c','#e67e22','#f39c12','#2ecc71','#27ae60'];
// The pie (and so every percentage in this section) covers only proteins that
// at least one method predicts as membrane; proteins no method classifies are
// listed as a plain count in the summary box.
const agreeData=agreeCounts.slice(1), nInPie=agreeData.reduce((a,b)=>a+b,0);
const agreeDisplayData=inflateSmallPieSlices(agreeData);
const pctOfInPie=v=>nInPie>0?(v/nInPie*100).toFixed(1):'0.0';

new Chart(document.getElementById('agreement-pie'),{
  type:'pie',
  data:{labels:AGREE_LABELS,datasets:[{data:agreeDisplayData,backgroundColor:AGREE_COLORS,borderColor:'#fff',borderWidth:3}]},
  options:{responsive:true,maintainAspectRatio:false,
    plugins:{legend:{display:false},tooltip:{callbacks:{label:ctx=>{
      const v=agreeData[ctx.dataIndex];
      return ` ${v} protein${v!==1?'s':''} (${pctOfInPie(v)}%)`;
    }}}},
    onClick(evt,elements){
      if(!elements.length) return;
      const si=elements[0].index, al=si+1;
      if(_openSlice.agree===si){closePanel('agree');return;} _openSlice.agree=si;
      const m=RAW.filter((_,i)=>proteinAgreeLevel[i]===al).map(r=>({...r,label:`${al} method${al!==1?'s':''}`}));
      openPanel('agree',`${AGREE_LABELS[si]} — ${m.length} protein${m.length!==1?'s':''}`,m);
    },
  },
  plugins:[makePiePlugin('agreePie',()=>nInPie,()=>agreeData)]
});

const agreeLegDiv=document.getElementById('agreement-legend');
AGREE_LABELS.forEach((lbl,i)=>{
  const cnt=agreeData[i],pct=pctOfInPie(cnt);
  const row=document.createElement('div');row.className='stat-row';
  row.innerHTML=`<span><span class="legend-dot" style="background:${AGREE_COLORS[i]}"></span>${lbl}</span><span>${cnt} (${pct}%)</span>`;
  agreeLegDiv.appendChild(row);
});
const aL3=agreeCounts.slice(3).reduce((a,b)=>a+b,0),
      aL4=agreeCounts.slice(4).reduce((a,b)=>a+b,0),aA5=agreeCounts[5];
document.getElementById('agreement-stats').innerHTML=`
  <div class="stat-row"><span>Total proteins</span><span>${nTotal}</span></div>
  <div class="stat-row"><span>Not classified by any method</span><span>${agreeCounts[0]}</span></div>
  <div class="stat-row"><span>Membrane proteins (≥1 method)</span><span>${nInPie}</span></div>
  <div class="stat-row"><span>≥3 methods (majority)</span><span>${aL3} (${pctOfInPie(aL3)}%)</span></div>
  <div class="stat-row"><span>≥4 methods (strong)</span><span>${aL4} (${pctOfInPie(aL4)}%)</span></div>
  <div class="stat-row"><span>All 5 methods</span><span>${aA5} (${pctOfInPie(aA5)}%)</span></div>`;

// 8 · Consensus pie + panel
const CONS_ALL_LABELS=['Full Consensus','High Confidence','Moderate Confidence','Low Confidence','No Consensus'];
const CONS_ALL_COLORS=['#27ae60','#2ecc71','#f39c12','#e74c3c','#95a5a6'];
const consCounts=Object.fromEntries(CONS_ALL_LABELS.map(l=>[l,0]));
const proteinConsLevel=RAW.map(r=>{
  const preds=METHODS.map(m=>getCount(r,m)), freq={};
  preds.forEach(p=>{freq[p]=(freq[p]||0)+1;});
  let mv=preds[0],mf=0;
  Object.entries(freq).forEach(([v,c])=>{if(c>mf||(c===mf&&Number(v)>mv)){mf=c;mv=Number(v);}});
  if(mv===0) return null;
  let level;
  if(mf===5) level='Full Consensus'; else if(mf===4) level='High Confidence';
  else if(mf===3) level='Moderate Confidence'; else if(mf===2) level='Low Confidence';
  else level='No Consensus';
  consCounts[level]++; return level;
});
const hasNC=consCounts['No Consensus']>0;
const CONS_ORDER=hasNC?CONS_ALL_LABELS:CONS_ALL_LABELS.slice(0,4);
const CONS_COLORS=hasNC?CONS_ALL_COLORS:CONS_ALL_COLORS.slice(0,4);
const consData=CONS_ORDER.map(l=>consCounts[l]), nMembrane=consData.reduce((a,b)=>a+b,0);
const consDisplayData=inflateSmallPieSlices(consData);

new Chart(document.getElementById('consensus-pie'),{
  type:'pie',
  data:{labels:CONS_ORDER,datasets:[{data:consDisplayData,backgroundColor:CONS_COLORS,borderColor:'#fff',borderWidth:3}]},
  options:{responsive:true,maintainAspectRatio:false,
    plugins:{legend:{display:false},tooltip:{callbacks:{label:ctx=>{
      const v=consData[ctx.dataIndex];
      return ` ${v} protein${v!==1?'s':''} (${nMembrane>0?(v/nMembrane*100).toFixed(1):'0.0'}%)`;
    }}}},
    onClick(evt,elements){
      if(!elements.length) return;
      const si=elements[0].index, level=CONS_ORDER[si];
      if(_openSlice.cons===si){closePanel('cons');return;} _openSlice.cons=si;
      const m=RAW.filter((_,i)=>proteinConsLevel[i]===level).map(r=>({...r,label:level}));
      openPanel('cons',`${level} — ${m.length} protein${m.length!==1?'s':''}`,m);
    },
  },
  plugins:[makePiePlugin('consPie',()=>nMembrane,()=>consData)]
});

const consLegDiv=document.getElementById('consensus-legend');
CONS_ORDER.forEach((lbl,i)=>{
  const cnt=consCounts[lbl],pct=nMembrane>0?(cnt/nMembrane*100).toFixed(1):'0.0';
  const row=document.createElement('div');row.className='stat-row';
  row.innerHTML=`<span><span class="legend-dot" style="background:${CONS_COLORS[i]}"></span>${lbl}</span><span>${cnt} (${pct}%)</span>`;
  consLegDiv.appendChild(row);
});
const hOrP=consCounts['Full Consensus']+consCounts['High Confidence'];
const noConsLine=hasNC?`<div class="stat-row"><span>No Consensus</span><span>${consCounts['No Consensus']} (${(consCounts['No Consensus']/nMembrane*100).toFixed(1)}%)</span></div>`:'';
document.getElementById('consensus-stats').innerHTML=`
  <div class="stat-row"><span>Membrane proteins total</span><span>${nMembrane}</span></div>
  <div class="stat-row"><span>Full + High Confidence</span><span>${hOrP} (${nMembrane>0?(hOrP/nMembrane*100).toFixed(1):'0.0'}%)</span></div>
  <div class="stat-row"><span>Moderate Confidence</span><span>${consCounts['Moderate Confidence']} (${nMembrane>0?(consCounts['Moderate Confidence']/nMembrane*100).toFixed(1):'0.0'}%)</span></div>
  <div class="stat-row"><span>Low Confidence</span><span>${consCounts['Low Confidence']} (${nMembrane>0?(consCounts['Low Confidence']/nMembrane*100).toFixed(1):'0.0'}%)</span></div>
  ${noConsLine}`;

// ── Section tabs ──────────────────────────────────────────────────────
// Inactive panels are hidden via the position:absolute + visibility:hidden
// .tab-panel-hidden class (see its CSS comment in batch_stats.html) rather
// than display:none, so every chart already has correct dimensions from
// the moment it's created — no resize-on-activate step needed here.
(function () {
  const tabs = document.querySelectorAll('.stats-tab');
  const panels = document.querySelectorAll('[data-tab-panel]');

  function activateTab(name) {
    tabs.forEach(tab => {
      const isActive = tab.dataset.tab === name;
      tab.classList.toggle('active', isActive);
      tab.setAttribute('aria-selected', isActive ? 'true' : 'false');
    });
    panels.forEach(panel => {
      panel.classList.toggle('tab-panel-hidden', panel.dataset.tabPanel !== name);
    });
  }

  tabs.forEach(tab => {
    tab.addEventListener('click', () => activateTab(tab.dataset.tab));
  });
})();

// ── Chart PNG export ─────────────────────────────────────────────────
// Chart.js keeps the canvas's current pixel content addressable as a PNG
// data URL via toBase64Image() — no extra library needed. Works for
// whichever mode (absolute/fractions) a toggled chart is currently in,
// since it just reads the canvas as rendered.
document.querySelectorAll('.chart-download-btn').forEach(btn => {
  btn.addEventListener('click', () => {
    const canvas = document.getElementById(btn.dataset.chart);
    const chart = canvas && Chart.getChart(canvas);
    if (!chart) return;
    // Charts whose on-screen legend is HTML (outside the canvas) get their
    // canvas legend drawn just for the export, so the PNG still says which
    // colour is which. update('none') redraws synchronously — provided no
    // animation is still running (e.g. right after toggling a checkbox),
    // hence stop() first — so the page never visibly shows the temporary legend.
    const exportLegend = 'exportLegend' in btn.dataset;
    if (exportLegend) { chart.stop(); chart.options.plugins.legend.display = true; chart.update('none'); }
    const link = document.createElement('a');
    link.href = chart.toBase64Image('image/png', 1);
    if (exportLegend) { chart.options.plugins.legend.display = false; chart.update('none'); }
    link.download = btn.dataset.filename || 'chart.png';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  });
});
