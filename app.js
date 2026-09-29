'use strict';
/* LIVIA Atlas — search AlphaFold-predicted protein interactions by protein.
 * One page per protein per species, gathering its predictions from every screen of that species: each prediction
 * keeps its screen, chain order and rank. The page runs LIVIA cLIP on them (contact residue frequency, clustered
 * interaction fingerprint, cluster info, interaction residues, 3D structure, interaction scatter plot) and adds a
 * partner overview, a network and a partner table.
 * Static: each screen is a folder (manifest.json, proteins.json, edges.tsv, b/<id>.zip, s/<id>.fa) — or one
 * uncompressed zip read by byte range (Zenodo) — and each species an index over its screens (proteins.json keyed by
 * UniProt accession, or by FlyBase gene for fly, whose bundles gather every construct of a gene), all in datasets.json. */

const DEV = location.hostname === 'localhost' || location.hostname === '127.0.0.1';   // local preview: LIVIA on :8000, screens from ../data/
const LIVIA = DEV ? 'http://localhost:8000/' : 'https://flyark.github.io/LIVIA/';
const CUT = { 10: 0.223, 5: 0.339, 1: 0.551 };
// Common names and abbreviations of viruses the taxonomy names otherwise: alias → taxonomy names (search, virus list filter)
const VALIAS = (() => { const m = {}; for (const [name, xs] of Object.entries({
  'monkeypox virus': ['mpox', 'mpox virus', 'mpxv'], 'variola virus': ['smallpox'], 'epstein-barr virus': ['ebv', 'hhv-4', 'hhv4'],
  'human herpesvirus 1': ['hsv-1', 'hsv1', 'herpes simplex', 'herpes simplex virus', 'herpes simplex virus 1'], 'human herpesvirus 2': ['hsv-2', 'hsv2', 'herpes simplex virus 2'],
  'human herpesvirus 8 type p': ['kshv', 'hhv-8', 'hhv8'], 'human cytomegalovirus': ['hcmv', 'cmv', 'hhv-5'], 'varicella-zoster virus': ['vzv', 'chickenpox', 'hhv-3'],
  'human papillomavirus type 16': ['hpv16', 'hpv-16'], 'hepatitis b virus': ['hbv'], 'zaire ebolavirus': ['ebola', 'ebola virus'] })) for (const x of xs) (m[x] ||= []).push(name);
  return m; })();
const cutNote = (f) => `<i class="kc">iLIS ≥ ${CUT[f].toFixed(3)}</i>`;   // the cutoff under a "past … FPR" count
const BAND = { 1: '#6B21A8', 5: '#0C735C', 10: '#875F00', 0: '#A7B2BF' };   // LIVIA's band colors (js/livia-core.js ilisColor)
const bandOf = (v) => (v >= CUT[1] ? 1 : v >= CUT[5] ? 5 : v >= CUT[10] ? 10 : 0);
const bandLabel = { 1: '1% FPR', 5: '5% FPR', 10: '10% FPR', 0: 'below' };
// Benchmarked cutoffs at 10 / 5 / 1% FPR for single models and for the average over a pair's models: AFM-LIS
// thresholds_data_yfh_lipdockq.xlsx ("total group"; Y2H reference sets in yeast, fly and human — Kim et al. 2026, FlyPredictome).
const FPR = { iLIS: [0.223, 0.339, 0.551], ipTM: [0.48, 0.59, 0.72], iLIA: [620.3, 1247.4, 3078.8], iLISA: [143.9, 360.6, 1241.0],
  LIS: [0.168, 0.257, 0.439], cLIS: [0.298, 0.449, 0.716], ipSAE: [0.165, 0.363, 0.615], actifpTM: [0.745, 0.880, 0.963],
  pDockQ: [0.293, 0.385, 0.567], LIpDockQ: [0.108, 0.178, 0.332], pDockQ2: [0.022, 0.036, 0.158], LIpDockQ2: [0.058, 0.132, 0.377] };   // best model, same benchmark as iLIS
const FPR_AVG = { iLIS: [0.072, 0.120, 0.268], ipTM: [0.292, 0.336, 0.442], iLISA: [40.0, 101.7, 332.0], LIS: [0.057, 0.096, 0.201], cLIS: [0.084, 0.158, 0.358],
  ipSAE: [0.039, 0.085, 0.241], actifpTM: [0.419, 0.504, 0.738], pDockQ: [0.288, 0.343, 0.438], LIpDockQ: [0.041, 0.064, 0.128], pDockQ2: [0.014, 0.019, 0.052], LIpDockQ2: [0.025, 0.044, 0.131] };   // mean of five models, same benchmark
// A UniProt domain box under the pointer: its name and span (the plots draw D1, D2 … or a clipped name)
const domTip = (d) => { const a = d.s != null ? d.s : d.start, z = d.e != null ? d.e : d.end; return `<b>${d.idx ? `D${d.idx} ` : ''}${esc(d.name)}</b><br>residues ${a}–${z} (${z - a + 1} aa) · UniProt domain`; };
// Each screen's run settings (registry: models, recycles). The cutoffs were calibrated on five models with five recycles
// per pair; a screen run otherwise says so where its settings are shown.
const CALIB = 'Cutoffs and average-iLIS cutoffs were calibrated on AlphaFold-Multimer runs with five models and five recycles per pair. For screens run with fewer models or other settings, they are a guide rather than a measured error rate.';
const runSettings = (d) => { if (!d || !d.models) return '';
  const off = !(d.models === 5 && d.recycles === 5), t = `${d.models} model${d.models === 1 ? '' : 's'}${d.recycles ? ` × ${d.recycles} recycles` : ''} per pair${d.modelsNote ? ` (${d.modelsNote})` : ''}`;
  return `<span class="runset${off ? ' off' : ''}"${off ? ` title="${esc(CALIB)}"` : ''}>${esc(t)}</span>`; };
// A pair folded in several runs counts once, by its run with the highest iLIS: which run a view shows, and what else holds it.
function overlapNote(sp, B, preds, P) {
  const by = new Map();
  for (const p of preds) { if (!by.has(p.run)) by.set(p.run, { n: 0, best: 0, rep: p.rep }); const r = by.get(p.run); r.n++; r.best = Math.max(r.best, p.iLIS || 0); }
  const shown = [...by].filter(([, r]) => !r.rep), also = [...by].filter(([, r]) => r.rep);
  if (!shown.length || !also.length) return '';
  const lab = ([rid, r]) => `${runLabel(sp, B, rid, P)} (${r.n} model${r.n === 1 ? '' : 's'}; highest iLIS ${r.best.toFixed(3)})`;
  return `Shown: ${shown.map(lab).join(', ')}. Also in: ${also.map(lab).join(', ')}, listed as ${also.length === 1 ? 'a repeat' : 'repeats'}. A pair folded more than once counts once, by its run with the highest iLIS.`;
}
const bandIn = (cuts, v) => (v >= cuts[2] ? 1 : v >= cuts[1] ? 5 : v >= cuts[0] ? 10 : 0);
// Text shades of the band colors, at least 4.5:1 on white and on the light band chips (BAND itself stays for plots and swatches)
const BAND_TXT = { 1: '#6B21A8', 5: '#0C735C', 10: '#875F00', 0: '#5F6771' };
const bandCol = (cuts, v) => BAND_TXT[bandIn(cuts, v)];   // a value's color = its FPR band under its own metric's cutoff
const BAND_W = { 1: 700, 5: 600, 10: 500, 0: 400 };   // weight grows with the band: 1% boldest, below 10% plain
const bandSty = (cuts, v) => `color:${bandCol(cuts, v)};font-weight:${BAND_W[bandIn(cuts, v)]}`;
const ARCHIVE = { doi: '10.5281/zenodo.22964479', url: 'https://doi.org/10.5281/zenodo.22964479' };   // the atlas's data record: the concept DOI, always the latest version
// "Cite this view": the page, its link and the date, with the Atlas data record, copied for a methods section or a legend
const citeBtn = (title, ds = []) => `<button class="cite-link" type="button" data-cite="${esc(title)}" data-ds="${esc(ds.join(','))}" title="copy a citation of this page: its title, link and date, the Atlas data version it reads, and the screens its predictions come from">Cite</button>`;
// A page's citation: the page (title, link, date), the Atlas data version each of its screens is read from (the Zenodo
// version DOI, not the concept DOI that always opens the newest), the method, and the source of each screen shown, which
// the screens' licenses (CC BY) require.
const REC_VERSION = { 22964480: '1.0', 22967610: '1.1', 22968056: '1.2', 22984781: '0.1.3' };   // the Atlas record's versions on Zenodo
const archiveOf = (recs) => (recs.length ? recs.map((r) => `LIVIA Atlas version ${REC_VERSION[r] || '?'}, Zenodo, https://doi.org/10.5281/zenodo.${r}`).join('; ') : `LIVIA Atlas, Zenodo, https://doi.org/${ARCHIVE.doi}`);
const archiveLine = () => { const recs = [...new Set(((REG && REG.datasets) || []).filter((d) => d.status === 'live').map(recOf).filter(Boolean))].sort();   // the versions this site reads
  return recs.map((r) => { const ids = [...new Set(REG.datasets.filter((d) => recOf(d) === r).map((d) => d.short))]; return `LIVIA Atlas version ${REC_VERSION[r] || '?'}, <i>Zenodo</i>, <a href="https://doi.org/10.5281/zenodo.${r}" target="_blank" rel="noopener">doi:10.5281/zenodo.${r}</a> (${ids.map(esc).join(', ')})`; }).join('; ') || `<a href="${ARCHIVE.url}" target="_blank" rel="noopener">doi:${ARCHIVE.doi}</a>`; };
const recOf = (d) => { const m = /records\/(\d+)\//.exec((d && d.zip && d.zip.url) || ''); return m ? m[1] : null; };
const citeText = (title, ids = []) => {
  const ds = ids.map((id) => ((REG && REG.datasets) || []).find((d) => d.id === id)).filter(Boolean);
  const recs = [...new Set(ds.map(recOf).filter(Boolean))];
  const doiUrl = (u) => { const m = /(10\.\d{4,9}\/[^\s?#]+?)(v\d+)?(\.full(\.pdf)?)?$/.exec(u || ''); return m && !/^http.*nvidia/.test(u) ? `https://doi.org/${m[1]}` : u; };   // a preprint cited by its DOI
  const src = [...new Set(ds.filter((d) => d.paper && !d.paper.includes(ARCHIVE.doi)).map((d) => `${d.source.replace(' · ', ', ')}, ${doiUrl(d.paper)}`))];
  return `${title}. LIVIA Atlas, ${location.origin}${location.pathname.replace(/index\.html$/, '')}${location.hash} (accessed ${new Date().toISOString().slice(0, 10)}). `
    + `Kim, A.-R. & Perrimon, N. (2026). LIVIA: a browser-based tool for assessing and visualizing predicted protein interactions. bioRxiv. https://doi.org/${REF.livia[1]}. `
    + `Data archive: ${archiveOf(recs)}.` + (src.length ? ` Screens: ${src.join('; ')}.` : '');
};
document.addEventListener('click', async (e) => {
  const b = e.target.closest && e.target.closest('[data-cite]'); if (!b) return;
  const t = citeText(b.dataset.cite, (b.dataset.ds || '').split(',').filter(Boolean)), box = b.parentElement.parentElement.querySelector('.cite-box');
  try { await navigator.clipboard.writeText(t); b.textContent = 'Copied'; setTimeout(() => { b.textContent = 'Cite this view'; }, 2000); if (box) box.remove(); }
  catch (err) {   // no clipboard (an insecure page, a locked-down browser): show the text to copy by hand
    const ta = box || Object.assign(document.createElement('textarea'), { className: 'cite-box', readOnly: true, rows: 3 });
    ta.value = t; if (!box) b.parentElement.after(ta); ta.focus(); ta.select(); }
});
const REF = {
  livia: ['Kim & Perrimon (2026) LIVIA, bioRxiv', '10.64898/2026.05.01.721633'],
  flypredictome: ['Kim et al. (2026) FlyPredictome, bioRxiv', '10.64898/2026.04.14.718529'],
  afmlis: ['Kim et al. (2024) AFM-LIS, bioRxiv', '10.1101/2024.02.19.580970'],
  biogrid: ['Oughtred et al. 2021', '10.1002/pro.3978'],
  uniprot: ['UniProt Consortium 2025', '10.1093/nar/gkae1010'], flybase: ['Öztürk-Çolak et al. 2024', '10.1093/genetics/iyad211'],
  afdb: ['Varadi et al. 2024', '10.1093/nar/gkad1011'], alphafold: ['Jumper et al. 2021', '10.1038/s41586-021-03819-2'],
  alphamissense: ['Cheng et al. 2023', '10.1126/science.adg7492'], afm: ['Evans et al. 2021', '10.1101/2021.10.04.463034'],
  colabfold: ['Mirdita et al. 2022', '10.1038/s41592-022-01488-1'], silhouette: ['Rousseeuw 1987', '10.1016/0377-0427(87)90125-7'],
  molstar: ['Sehnal et al. 2021', '10.1093/nar/gkab314'], d3: ['Bostock et al. 2011', '10.1109/TVCG.2011.185'],
  leiden: ['Traag et al. 2019', '10.1038/s41598-019-41695-z'], zenodo: ['Zenodo', '10.25495/7GXK-RD71'],
  schmid2025: ['Schmid et al. 2025', '10.1101/2025.11.10.687652'], kim2025: ['Kim et al. 2025', '10.1101/2025.10.10.681672'],
};
const cite = (k) => `<a href="https://doi.org/${REF[k][1]}" target="_blank" rel="noopener">${REF[k][0]}</a>`;

const $ = (sel, el = document) => el.querySelector(sel);
const short = (name) => String(name || '').replace(/\s*\((?:EC [^)]*|[^)]*)\).*$/, '').trim() || String(name || '');
const esc = (s) => String(s ?? '').replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
const fmtInt = (n) => Number(n).toLocaleString('en-US');
const fmtNum = (v, d) => (Number.isFinite(v) ? v.toFixed(d) : '–');
const mean = (a) => (a.length ? a.reduce((s, v) => s + v, 0) / a.length : 0);
const el = (html) => { const t = document.createElement('template'); t.innerHTML = html.trim(); return t.content.firstElementChild; };
const app = $('#app');
const TIP = document.body.appendChild(Object.assign(document.createElement('div'), { className: 'gtip', hidden: true }));
function showTip(html, x, y) {
  TIP.innerHTML = html; TIP.hidden = false;
  const w = TIP.offsetWidth, h = TIP.offsetHeight, m = 10;
  let left = x - w / 2, top = y - h - 14;
  if (top < m) top = y + 18;
  left = Math.max(m, Math.min(window.innerWidth - w - m, left));
  TIP.style.left = left + 'px'; TIP.style.top = top + 'px';
}
const hideTip = () => { TIP.hidden = true; };
const uniprotLink = (acc, label) => (acc ? `<a href="https://www.uniprot.org/uniprotkb/${esc(acc)}" target="_blank" rel="noopener">${label || 'UniProt ' + esc(acc)}</a>` : '');

/* ── data: registry, screens, species, merged bundles, sequences ──────────────────────────────────────── */
let REG = null;
const DSC = {}, SPC = {};
let ROUTE = 0;   // bumped by every navigation: a view still loading for an earlier page stops at its next await
const stale = (gen) => gen !== ROUTE;
// The site's and the data repos' files. GitHub Pages answers 5xx now and then (right after a deploy): one retry, then a
// message a reader can act on, never a parse error.
async function getFile(url, as) {
  for (let i = 0; ; i++) {
    let res = null; try { res = await fetch(url); } catch (e) { /* offline, or blocked */ }
    if (res && res.ok) return as === 'text' ? res.text() : res.json();
    if (i || (res && res.status < 500)) throw new Error(`The atlas data could not be loaded${res ? ` (HTTP ${res.status})` : ''}. Try reloading the page.`);
    await new Promise((r) => setTimeout(r, 1500));
  }
}
const getJSON = (url) => getFile(url, 'json'), getText = (url) => getFile(url, 'text');
async function registry() { if (!REG) REG = await getJSON('datasets.json'); return REG; }
const regDataset = async (id) => (await registry()).datasets.find((d) => d.id === id);
const regSpecies = async (id) => ((await registry()).species || []).find((s) => s.id === id);

async function dataset(id) {   // one screen: its manifest; proteins.json only when a screen page needs it
  if (DSC[id]) return DSC[id];
  const reg = await regDataset(id);
  if (!reg || !reg.base) throw new Error(`Unknown dataset “${id}”.`);
  const base = new URL(DEV && reg.dev ? reg.dev : reg.base, location.href).href;   // each screen is its own Pages repo
  DSC[id] = { id, reg, base, manifest: await getJSON(base + 'manifest.json'), raw: new Map(), rows: null };
  return DSC[id];
}
// Thematic sets of a dataset (sets.json): views over the one dataset, each a list of its prediction runs — a screen the
// source names (FlyPredictome's kinase–TF screen, with its own paper) or a source category. A run can be in several.
function setsOf(ds) {
  if (!ds.setsJob) ds.setsJob = (async () => {
    const f = ds.manifest.files && ds.manifest.files.sets; if (!f) return null;
    try {
      const j = await getJSON(new URL(f, ds.base).href), runSets = new Map();
      for (const s of j.sets) { SET_COL[s.short] = s.color; for (const r of s.runs) { if (!runSets.has(r)) runSets.set(r, []); runSets.get(r).push(s); } }
      return { list: j.sets, byId: new Map(j.sets.map((s) => [s.id, s])), runSets, ds };
    } catch (e) { return null; }
  })();
  return ds.setsJob;
}
const primarySet = (tags) => (tags && tags.length ? tags.find((s) => s.type === 'screen') || tags[0] : null);   // a named screen before its category
// A screen's per-protein files: in its folder, or — for a screen kept as one uncompressed zip (Zenodo) — one HTTP Range
// read each, at the byte range the offsets map kept with the site gives. Never the whole archive.
const OFFS = new Map();
// Zenodo can answer slowly or stall, and a large interactome takes seconds to cluster. Each archive read and each clustering
// shows in a status pill (what, how long, a note when slow). A read times out after 45 s and is tried twice; then the pill
// says Zenodo did not answer, with a button to try again (a clustering failure is reported on its card).
const LOADS = new Map(); let LOADN = 0, LOADT = null, LOADFAIL = '';
function loadBar() { let el = document.getElementById('loadbar'); if (!el) { el = document.createElement('div'); el.id = 'loadbar'; el.setAttribute('role', 'status'); el.setAttribute('aria-live', 'polite'); el.hidden = true; document.body.appendChild(el); } return el; }
function paintLoads() {
  const el = loadBar(), list = [...LOADS.values()];
  if (!list.length) { clearInterval(LOADT); LOADT = null;
    if (LOADFAIL) { el.hidden = false; el.className = 'fail'; el.innerHTML = `${esc(LOADFAIL)} <button type="button" id="loadbar-retry">Try again</button>`; $('#loadbar-retry').onclick = () => location.reload(); }
    else el.hidden = true;
    return; }
  const now = Date.now(), sec = Math.round((now - Math.min(...list.map((x) => x.t0))) / 1000);
  const reads = [...new Set(list.filter((x) => x.kind === 'zenodo').map((x) => x.what))], tasks = [...new Set(list.filter((x) => x.kind !== 'zenodo').map((x) => x.what))];
  const slowRead = list.some((x) => x.kind === 'zenodo' && now - x.t0 >= 12000), slowTask = list.some((x) => x.kind !== 'zenodo' && now - x.t0 >= 8000);
  const text = [...(reads.length ? [`Reading ${reads.join(', ')} from Zenodo`] : []), ...tasks].join(' · ');
  el.hidden = false; el.className = slowRead ? 'slow' : '';
  el.innerHTML = `<span class="spin" aria-hidden="true"></span>${esc(text)} · ${sec} s${slowRead ? '<br><small>Zenodo is answering slowly. The page fills in as soon as the data arrive.</small>'
    : slowTask ? '<br><small>A large interactome takes a little longer to cluster; the page fills in when it is done.</small>' : ''}`;
}
function trackLoad(what, job, kind = 'zenodo') {   // kind 'zenodo': "Reading <what> from Zenodo"; any other kind: <what> as written
  const id = ++LOADN; if (kind === 'zenodo') LOADFAIL = ''; LOADS.set(id, { what, kind, t0: Date.now() }); if (!LOADT) LOADT = setInterval(paintLoads, 1000); paintLoads();
  return job.then((v) => { LOADS.delete(id); paintLoads(); return v; }, (e) => { LOADS.delete(id); if (kind === 'zenodo') LOADFAIL = e.message; paintLoads(); throw e; });
}
async function rangeRead(url, range) {   // → the bytes of one Range read, or null when the host did not answer 206
  for (let i = 1; ; i++) {
    const ac = new AbortController(), t = setTimeout(() => ac.abort(), 45000);
    try {
      const res = await fetch(url, { headers: { Range: range }, signal: ac.signal });
      if (res.status !== 206) { try { if (res.body) res.body.cancel(); } catch (e) { /* nothing to cancel */ } return null; }
      return await res.arrayBuffer();   // the body under the same clock: a read can stall after the headers
    } catch (e) { if (i >= 2) throw new Error(`Zenodo did not answer in time (tried twice). It is often slow for a few minutes.`); }
    finally { clearTimeout(t); }
  }
}
async function screenFile(ds, rel) {
  const z = ds.reg.zip;
  if (!z) { const res = await fetch(new URL(rel.split('/').map(encodeURIComponent).join('/'), ds.manifest.bundleBase || ds.base).href); return res.ok ? res : null; }
  if (!OFFS.has(ds.id)) OFFS.set(ds.id, fetch(new URL(z.offsets, location.href).href).then((r) => (r.ok ? r.json() : {})).catch(() => ({})));
  const at = (await OFFS.get(ds.id))[rel]; if (!at) return null;
  const t = ds.reg.title || ds.id, what = /^(Human|Zebrafish|Yeast|Fly|Worm)\b/.test(t) ? t.charAt(0).toLowerCase() + t.slice(1) : t;   // “the human kinase–kinase screen”, “C. elegans …”, “FlyPredictome”
  const buf = await trackLoad(what, rangeRead(DEV && z.dev ? new URL(z.dev, location.href).href : z.url, `bytes=${at[0]}-${at[0] + at[1] - 1}`));
  return buf ? new Response(buf) : null;
}
async function datasetRows(ds) {
  if (ds.rows) return ds.rows;
  const prot = await getJSON(ds.base + 'proteins.json'), c = Object.fromEntries(prot.columns.map((k, i) => [k, i]));
  ds.rows = prot.rows.map((r) => ({ id: r[c.id], gene: r[c.gene], acc: r[c.acc], partners: r[c.partners], pos10: r[c.pos10], pos5: r[c.pos5], pos1: r[c.pos1] }));
  return ds.rows;
}
const SPM = {};
function speciesManifest(id) {   // a species' counts and screens only (the home page needs no more)
  if (!SPM[id]) SPM[id] = (async () => { const reg = await regSpecies(id); if (!reg) throw new Error(`Unknown species “${id}”.`); return getJSON(new URL(reg.base, location.href).href + 'manifest.json'); })()
    .catch((e) => { delete SPM[id]; throw e; });
  return SPM[id];
}
function species(id) {   // one load per species however many callers ask at once
  if (!SPC[id]) SPC[id] = speciesIndex(id).catch((e) => { delete SPC[id]; throw e; });
  return SPC[id];
}
async function speciesIndex(id) {   // the species index: one row per protein over every screen, keyed by UniProt accession (fly: FlyBase gene)
  const reg = await regSpecies(id);
  if (!reg) throw new Error(`Unknown species “${id}”.`);
  const base = new URL(reg.base, location.href).href;
  const [manifest, prot] = await Promise.all([speciesManifest(id), getJSON(base + 'proteins.json')]);
  const c = Object.fromEntries(prot.columns.map((k, i) => [k, i]));
  const rows = prot.rows.map((r, i) => {
    const occ = String(r[c.occ] || '').split(' ').filter(Boolean).map((t) => { const k = t.indexOf(':'); return { di: +t.slice(0, k), name: t.slice(k + 1) }; });
    return { i, key: r[c.key], gene: r[c.gene] || r[c.name] || r[c.key], acc: r[c.acc], name: r[c.name], syn: r[c.syn], len: r[c.len], clen: r[c.clen], status: r[c.status], occ,
      id: occ.length ? occ[0].name : r[c.key], partners: r[c.partners], pos10: r[c.pos10], pos5: r[c.pos5], pos1: r[c.pos1], best: r[c.bestIlis],
      src: occ.reduce((m, o) => m | (1 << o.di), 0) };
  });
  const byKey = new Map(rows.map((r) => [r.key, r])), byName = new Map(), byGene = new Map();
  for (const r of rows) { for (const o of r.occ) byName.set(o.name, r); if (!byGene.has(r.gene)) byGene.set(r.gene, r); }
  const keys = rows.map((r) => ({ gene: r.gene.toLowerCase(), key: r.key.toLowerCase(), acc: (r.acc || '').toLowerCase(), ids: r.occ.map((o) => o.name.toLowerCase()),
    syn: (r.syn || '').toLowerCase().split(/\s+/).filter(Boolean), name: (r.name || '').toLowerCase() }));
  const dsIds = manifest.datasets.map((d) => d.id), regs = await Promise.all(dsIds.map(regDataset));
  let viruses = null;   // the viral species: which virus (taxon) each protein belongs to, for the virus pages and search
  if (manifest.files && manifest.files.viruses) {
    const V = await getJSON(base + manifest.files.viruses);
    const c = Object.fromEntries(V.columns.map((k, i) => [k, i])), both = (r, f) => (c[`hetero_pos${f}`] == null ? null : r[c[`hetero_pos${f}`]] + r[c[`homo_pos${f}`]]);
    viruses = V.viruses.map((r) => ({ taxid: r[c.taxid], name: r[c.name], n: r[c.proteins], het: r[c.hetero_pairs], hom: r[c.homo_pairs], hpos: r[c.hetero_pos10], mpos: r[c.homo_pos10],
      pos5: both(r, 5), pos1: both(r, 1), pwp: r[c.proteins_with_partner], members: r[c.rows],
      family: c.family != null ? r[c.family] || '' : '', genus: c.genus != null ? r[c.genus] || '' : '', host: c.host != null ? r[c.host] || '' : '',   // ICTV VMR MSL40
      species: c.species != null ? r[c.species] || '' : '', spTax: c.species_taxid != null ? r[c.species_taxid] : null, lineage: c.lineage != null ? r[c.lineage] || [] : [] }));
    for (const v of viruses) for (const i of v.members) rows[i].virus = v;
  }
  return { id, reg, base, manifest, rows, byKey, byName, byGene, keys, viruses, dsIds, dsShort: manifest.datasets.map((d) => d.short),
    dsColor: regs.map((r) => (r && r.color) || '#5B6B7F'), one: regs.every((r) => r && r.models === 1), edges: null, cache: new Map() };   // one: one model per pair, so best = average
}
// Pairs reported in BioGRID among a species' proteins, from data/species/<sp>/biogrid.tsv on the site: row numbers of the
// species index and the number of publications reporting a physical and a genetic interaction. The file names the index it
// was built for; a file for another index is not used, so a pair is never marked by a stale match. null: no file, or not this index.
// One parser for the whole file and for a shard (the same header: release, index size and hash of the species index)
async function parseReported(sp, t) {
  const head = t.slice(0, 600).split('\n'), idx = (head.find((l) => l.startsWith('# index ')) || '').split(' '), N = sp.rows.length;
  if (+idx[2] !== N) return null;
  const dig = await crypto.subtle.digest('SHA-1', new TextEncoder().encode(sp.rows.map((r) => r.key).join('\n')));
  if ([...new Uint8Array(dig)].map((b) => b.toString(16).padStart(2, '0')).join('').slice(0, 12) !== idx[3]) return null;
  const m = new Map(), gm = new Map(), rm = new Map();   // rm: the pair's publications, where the file carries them (per-protein shards)
  for (let p = 0; p < t.length;) { let e = t.indexOf('\n', p); if (e < 0) e = t.length;
    if (t[p] !== '#' && t[p] !== 'i') { const [a, b, ph, ge, pr, gr] = t.slice(p, e).split('\t'), k = +a * N + +b;
      if (+ph) m.set(k, +ph); if (+ge) gm.set(k, +ge); if (pr || gr) rm.set(k, { p: pr ? pr.split(',') : [], g: gr ? gr.split(',') : [] }); }
    p = e + 1; }
  const release = ((head[0] || '').match(/BioGRID ([\d.]+)/) || [])[1] || '', at = (a, b) => Math.min(a, b) * N + Math.max(a, b);
  return { release, pubs: (a, b) => m.get(at(a, b)) || 0, gen: (a, b) => gm.get(at(a, b)) || 0, refs: (a, b) => rm.get(at(a, b)) || null };   // publications: physical, genetic; their ids
}
function reported(sp) {   // every reported pair of the species: the network views
  if (!(sp.manifest.files || {}).biogridShards) return Promise.resolve(null);   // a species without BioGRID records (zebrafish, viruses) has no file to ask for
  if (!sp.known) sp.known = (async () => { let t; try { t = await getText(sp.base + 'biogrid.tsv'); } catch (e) { return null; } return parseReported(sp, t); })();
  return sp.known;
}
// The reported pairs of one protein (row i): its shard of biogrid.tsv (biogrid/<i % S>.tsv holds every pair of the proteins
// in it), so a protein page reads ~100 kB instead of the whole file; species without shards read the whole file.
function reportedOf(sp, i) {
  const S = (sp.manifest.files || {}).biogridShards; if (!S || i == null) return reported(sp);
  // always the protein's shard, even when the whole file is loaded: the shard carries each pair's publications
  sp.knownShard = sp.knownShard || new Map(); const k = i % S;
  if (!sp.knownShard.has(k)) sp.knownShard.set(k, (async () => { let t; try { t = await getText(`${sp.base}biogrid/${k}.tsv`); } catch (e) { return reported(sp); } return parseReported(sp, t); })());
  return sp.knownShard.get(k);
}
// Every pair folded in the species' screens (data/species/<sp>/tested.bin, built by tools/tested_index.py): whether a pair
// was tested at all, below the cutoff or not. The file names the index it was built for (rows and key hash); a file for
// another index is not used. → { has(a, b), P } or null (no file, or not this index)
function tested(sp) {
  const f = (sp.manifest.files || {}).tested; if (!f) return Promise.resolve(null);
  if (!sp.testedP) sp.testedP = (async () => {
    const res = await fetch(sp.base + f); if (!res.ok) return null; const buf = await res.arrayBuffer(), dv = new DataView(buf), td = new TextDecoder();
    const ver = td.decode(new Uint8Array(buf, 0, 5)); if (ver !== 'LVTP1' && ver !== 'LVTP2') return null;   // LVTP2 adds each pair's best iLIS
    const N = dv.getUint32(6, true), P = dv.getUint32(10, true), sha = td.decode(new Uint8Array(buf, 14, 12));
    if (N !== sp.rows.length) return null;
    const dig = await crypto.subtle.digest('SHA-1', new TextEncoder().encode(sp.rows.map((r) => r.key).join('\n')));
    if ([...new Uint8Array(dig)].map((b) => b.toString(16).padStart(2, '0')).join('').slice(0, 12) !== sha) return null;
    const off = new Uint32Array(buf, 32, N + 1), nb = new Uint16Array(buf, 32 + 4 * (N + 1), P), sc = ver === 'LVTP2' ? new Uint8Array(buf, 32 + 4 * (N + 1) + 2 * P, P) : null;
    const find = (a, b) => { const i = Math.min(a, b), j = Math.max(a, b); if (!(i >= 0 && j < N)) return -1; let lo = off[i], hi = off[i + 1];
      while (lo < hi) { const m = (lo + hi) >> 1; if (nb[m] < j) lo = m + 1; else hi = m; } return lo < off[i + 1] && nb[lo] === j ? lo : -1; };
    return { has: (a, b) => find(a, b) >= 0, score: (a, b) => { const x = sc ? find(a, b) : -1; return x >= 0 ? sc[x] / 255 : NaN; }, P };   // score: best iLIS to about ±0.002
  })().catch(() => null);
  return sp.testedP;
}
// Paralog families (data/species/<sp>/paralogs.json, tools/paralog_index.py): row → family, for counting two paralogs of one
// complex as one connection. Checked against the index it was built for. → { of(i), note } or null
function paralogs(sp) {
  const f = (sp.manifest.files || {}).paralogs; if (!f) return Promise.resolve(null);
  if (!sp.paraP) sp.paraP = (async () => { const J = await getJSON(sp.base + f);
    const dig = await crypto.subtle.digest('SHA-1', new TextEncoder().encode(sp.rows.map((r) => r.key).join('\n')));
    if ([...new Uint8Array(dig)].map((b) => b.toString(16).padStart(2, '0')).join('').slice(0, 12) !== J.index) return null;
    const m = new Map(); J.families.forEach((fam, n) => fam.forEach((i) => m.set(i, n)));
    return { of: (i) => (m.has(i) ? 'f' + m.get(i) : 'r' + i), note: J.identity, n: J.families.length }; })().catch(() => null);
  return sp.paraP;
}
const srcBadges = (sp, mask) => sp.dsIds.map((_, di) => (mask & (1 << di) ? `<span class="src" style="--c:${sp.dsColor[di]}">${esc(sp.dsShort[di])}</span>` : '')).join('');

// edges.tsv of a species: every pair past 10% FPR in any screen — best iLIS over every model, mean iLIS, screens
async function edges(sp, setId = '') {   // setId: that thematic set's edges (same row numbers as the species index)
  sp.edgesBy = sp.edgesBy || new Map();
  if (!sp.edgesBy.has(setId)) sp.edgesBy.set(setId, (async () => {
    let url = sp.base + 'edges.tsv', map = null;
    if (setId && sp.dsIds.includes(setId)) {   // one screen: its own best / average over its models, rows mapped onto the species index
      const ds = await dataset(setId), rows = await datasetRows(ds);
      map = rows.map((r) => { const R = sp.byName.get(r.id) || sp.byKey.get(r.id); return R ? R.i : -1; }); url = ds.base + 'edges.tsv';
    } else if (setId) for (const id of sp.dsIds) { const TS = await setsOf(await dataset(id)).catch(() => null), S = TS && TS.byId.get(setId); if (S) { url = new URL(S.files.edges, TS.ds.base).href; break; } }
    // Every protein's neighbours side by side in typed arrays (scores in thousandths, exact to the file's 3 decimals):
    // a few MB where a Map per protein took tens. adj.get(i) builds one protein's Map when a network asks for it.
    const text = await getText(url), nl = text.indexOf('\n'), head = text.slice(0, nl).split('\t');
    const cb = head.indexOf('iLIS_best'), ca = head.indexOf('iLIS_avg'), cs = head.indexOf('src'), cp = head.indexOf('ipTM_best'), N = sp.rows.length;
    const I = [], J = [], BE = [], AV = [], SR = [], IP = [];
    for (let p = nl + 1; p < text.length;) {
      let e = text.indexOf('\n', p); if (e < 0) e = text.length;
      if (e > p) { const t = text.slice(p, e).split('\t'), i = map ? map[+t[0]] : +t[0], j = map ? map[+t[1]] : +t[1];
        if (i >= 0 && j >= 0 && i < N && j < N) { I.push(i); J.push(j); BE.push(Math.round(+t[cb] * 1000)); AV.push(Math.round(+t[ca] * 1000)); SR.push(cs >= 0 ? +t[cs] : 1); IP.push(cp >= 0 ? Math.round(+t[cp] * 100) : -1); } }
      p = e + 1;
    }
    const off = new Uint32Array(N + 1);
    for (let k = 0; k < I.length; k++) { off[I[k] + 1]++; off[J[k] + 1]++; }
    for (let v = 0; v < N; v++) off[v + 1] += off[v];
    const at = off.slice(0, N), nb = new Int32Array(2 * I.length), best = new Int16Array(2 * I.length), avg = new Int16Array(2 * I.length), src = new Uint8Array(2 * I.length), ipt = new Int16Array(2 * I.length);
    const put = (a, b, k) => { const x = at[a]++; nb[x] = b; best[x] = BE[k]; avg[x] = AV[k]; src[x] = SR[k]; ipt[x] = IP[k]; };
    for (let k = 0; k < I.length; k++) { put(I[k], J[k], k); put(J[k], I[k], k); }
    const memo = new Map();
    const adj = { get(v) {
      if (!(v >= 0 && v < N) || off[v] === off[v + 1]) return undefined;
      if (!memo.has(v)) { if (memo.size > 500) memo.clear(); const m = new Map();
        for (let x = off[v]; x < off[v + 1]; x++) m.set(nb[x], { best: best[x] / 1000, avg: avg[x] / 1000, src: src[x], iptm: ipt[x] >= 0 ? ipt[x] / 100 : NaN }); memo.set(v, m); }
      return memo.get(v); } };
    return { adj };
  })().catch((e) => { sp.edgesBy.delete(setId); throw e; }));   // a failed read is not kept: the next draw tries again
  return sp.edgesBy.get(setId);
}
function parseCSV(text) {   // lis.py CSV: quoted only where a residue list holds commas. Fields are sliced, never built a character at a time
  const out = [], n = text.length;
  for (let i = 0; i < n;) {
    let e = text.indexOf('\n', i); if (e < 0) e = n;
    const end = e > i && text.charCodeAt(e - 1) === 13 ? e - 1 : e, row = [];
    for (let k = i; k <= end;) {
      if (text.charCodeAt(k) === 34 && k < end) {   // "…", with "" for a quote
        let f = '', m = k + 1;
        for (;;) { const q = text.indexOf('"', m); if (q < 0 || q >= end) { f += text.slice(m, end); k = end; break; }
          f += text.slice(m, q); if (text.charCodeAt(q + 1) === 34) { f += '"'; m = q + 2; } else { k = q + 1; break; } }
        row.push(f); if (k < end && text.charCodeAt(k) === 44) k++; else break;
      } else { let c = text.indexOf(',', k); if (c < 0 || c > end) c = end; row.push(text.slice(k, c)); k = c + 1; }
    }
    out.push(row); i = e + 1;
  }
  return out;
}
// A string cut from a large text keeps that whole text alive (V8 slices point into their parent): what the predictions
// keep is copied out, so a bundle's CSV can be freed once it is read.
const own = (s) => (typeof s === 'string' && s.length > 12 ? (' ' + s).slice(1) : s);
const rangesOf = (a) => { if (!a.length) return '[]'; const t = []; let x = a[0], y = a[0];
  for (const r of a.slice(1)) { if (r === y + 1) { y = r; continue; } t.push(x === y ? String(x) : `${x}-${y}`); x = y = r; } t.push(x === y ? String(x) : `${x}-${y}`); return '[' + t.join(',') + ']'; };
const expand = (s) => { const out = []; for (const t of String(s || '').replace(/[\[\]\s]/g, '').split(',')) { if (!t) continue;
  const [a, b] = t.split('-').map(Number); if (b >= a) for (let x = a; x <= b; x++) out.push(x); else if (a) out.push(a); } return out; };
function parseFasta(text) { const m = new Map(); let id = null, buf = [];
  for (const line of (text || '').split('\n')) { if (line.startsWith('>')) { if (id) m.set(id, buf.join('')); id = line.slice(1).trim(); buf = []; } else if (line) buf.push(line.trim()); }
  if (id) m.set(id, buf.join('')); return m; }

const bundleRel = (ds, name) => (ds.manifest.files.bundle || 'b/{id}.zip').replace('{id}', name);
const bundleUrl = (ds, name) => new URL(bundleRel(ds, name).split('/').map(encodeURIComponent).join('/'), ds.manifest.bundleBase || ds.base).href;
// A gene-keyed bundle names its constructs: name → gene key, kind, label, length, offset on the gene's reference sequence
function parseCons(text) {
  const m = new Map();
  for (const line of text.split('\n').slice(1)) { if (!line) continue; const [name, key, kind, label, len, off, exact, mut, seg] = line.split('\t');
    m.set(name, { key, kind, label, len: +len, off: off === '' || off == null ? null : +off, exact: exact === 'y', mut: mut || '',
      seg: seg ? seg.split(',').map((t) => t.split(':').map(Number)) : null }); }
  return m;
}
function bundleRaw(ds, name) {   // one screen's cLIP bundle for one protein: lis.py rows + FASTA (+ its construct table)
  const rel = bundleRel(ds, name);
  if (!ds.raw.has(rel)) {
    const job = (async () => {
      const res = await screenFile(ds, rel);
      if (!res) throw new Error(`No interaction data for ${name}.`);
      const bytes = await res.arrayBuffer(), zip = await JSZip.loadAsync(bytes), files = Object.keys(zip.files);
      const csvName = files.find((f) => /\.csv$/i.test(f) && !/identity[_-]?map/i.test(f)), faName = files.find((f) => /\.(fa|fasta)$/i.test(f));
      const conName = files.find((f) => /(^|\/)constructs\.tsv$/.test(f));
      return { rel, bytes, csvName, faName, ...(await csvRows(zip, csvName)),
        seqs: parseFasta(faName ? await zip.file(faName).async('string') : ''), cons: conName ? parseCons(await zip.file(conName).async('string')) : null,
        isoforms: files.includes('isoforms.json') ? JSON.parse(await zip.file('isoforms.json').async('string')) : null };   // v1.2: its other isoforms are files of their own
    })();
    ds.raw.set(rel, job); job.catch(() => ds.raw.delete(rel));
    while (ds.raw.size > 16) ds.raw.delete(ds.raw.keys().next().value);   // the most recent bundles only
  } else { const job = ds.raw.get(rel); ds.raw.delete(rel); ds.raw.set(rel, job); }
  return ds.raw.get(rel);
}
async function csvRows(zip, csvName) {   // a bundle's lis.py rows; merged() reads them once and lets them go (a big bundle's rows run to hundreds of MB)
  const t = parseCSV(await zip.file(csvName).async('string'));
  return { header: t[0], H: Object.fromEntries(t[0].map((k, i) => [k, i])), rows: t.slice(1).filter((r) => r.length > 10) };
}
// Every prediction of one protein across the screens of its species. A pair predicted in two screens, or both ways
// round in one, keeps every model; each keeps its screen and chain order (its "run") and rank. A gene-keyed bundle
// (FlyPredictome) holds every construct of the gene: its construct table says which construct belongs to which gene,
// and each prediction also keeps its screen category (set) and prediction run (batch).
function merged(sp, P, scope = '', whole = false) {   // scope: one screen of the species (its dataset id) or one thematic set; whole: every isoform file
  const ck = P.key + (scope ? '?' + scope : '') + (whole ? '#whole' : ''), onlyDi = scope ? sp.dsIds.indexOf(scope) : -1, setId = onlyDi >= 0 ? '' : scope;
  if (!sp.cache.has(ck)) {
    const job = (async () => {
      const errs = [], parts = (await Promise.all(P.occ.filter((o) => onlyDi < 0 || o.di === onlyDi).map(async (o) => {   // a screen that cannot be reached is left out
        try { const ds = await dataset(sp.dsIds[o.di]); return { di: o.di, name: o.name, ds, raw: await bundleRaw(ds, o.name), TS: await setsOf(ds) }; } catch (e) { errs.push(e); return null; }
      }))).filter(Boolean);
      if (!parts.length) throw (errs.find((e) => /Zenodo/.test(e.message)) || new Error(`No interaction data for ${P.gene}.`));
      const split = parts.find((x) => x.raw.isoforms) || null;   // atlas v1.2: this gene's other isoforms are files of their own
      if (split && whole) parts.push(...await Promise.all(split.raw.isoforms.choices.map(async (c) => ({ ...split, name: c.file, raw: await bundleRaw(split.ds, c.file) }))));
      const all = await assemble(sp, P, parts, scope, onlyDi, setId);
      if (split && !whole) isoformFiles(all, split, sp, P, scope, onlyDi, setId);
      if (scope && !all.preds.length && !(all.choices || []).length) throw new Error(`${P.gene} has no predictions in this ${onlyDi >= 0 ? 'screen' : 'set'}.`);
      return all;
    })();
    sp.cache.set(ck, job); job.catch(() => sp.cache.delete(ck));
    while (sp.cache.size > 6) sp.cache.delete(sp.cache.keys().next().value);   // the proteins read most recently
  } else { const job = sp.cache.get(ck); sp.cache.delete(ck); sp.cache.set(ck, job); }
  return sp.cache.get(ck);
}
// The models a view of just these predictions counts: one run per sequence pair, its highest-iLIS run.
function countOnce(list) {
  const best = new Map(), n = new Map();
  for (const p of list) { n.set(p.run, (n.get(p.run) || 0) + 1); const b = best.get(p.sp); if (!b || (p.iLIS || 0) > b.v) best.set(p.sp, { v: p.iLIS || 0, run: p.run }); }
  let s = 0; for (const b of best.values()) s += n.get(b.run); return s;
}
// A pair of the viral release that the AlphaFold Database does not display: its model and PAE read by byte range from the
// release archive at EBI (zstd members of chunk_N.tar; addresses from the Atlas index), unpacked in the browser (fzstd,
// fflate) and handed to LIVIA's prediction page as a bundle. The tab opens at the click so no popup blocker stops it.
const ARCH_COLS = ['chunk', 'cif_off', 'cif_len', 'pae_off', 'pae_len'];
const ARCH_URL = (n) => `https://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/pandemic_prep/chunk_${n}.tar`;
// A virus pair's structure: its model in the AlphaFold Database (pairs the database displays) or in the release archive
// (the rest), from its virus's pair table; null for a pair below the 10% FPR cutoff (not in the table) or outside viruses.
// A virus's pairs file (data/species/virus/pairs/<taxid>.tsv), read once: each pair's model, its archive address and every
// score lis.py wrote, for the virus protein page's partner table (the same columns as the virus page's pairs table).
const VMCOL = ['iLIS', 'iLISA', 'ipSAE', 'actifpTM', 'ipTM', 'LIS', 'cLIS', 'pDockQ', 'LIpDockQ', 'pDockQ2', 'LIpDockQ2'];
const vmfmt = (k, x) => (!Number.isFinite(x) ? '–' : k === 'iLISA' ? x.toFixed(1) : k === 'ipTM' ? x.toFixed(2) : x.toFixed(3));
function virusRows(sp, taxid) {
  sp.vpairs = sp.vpairs || new Map();
  if (!sp.vpairs.has(taxid)) sp.vpairs.set(taxid, getText(sp.base + `pairs/${taxid}.tsv`).then((t) => {
    const L = t.trim().split('\n'), ix = Object.fromEntries(L[0].split('\t').map((k, i) => [k, i]));
    return L.slice(1).map((l) => { const r = l.split('\t'), ad = {}; for (const k of ARCH_COLS) if (r[ix[k]] !== '' && r[ix[k]] != null) ad[k] = +r[ix[k]];
      const m = {}; for (const k of VMCOL) if (ix[k] != null && r[ix[k]] !== '') m[k] = +r[ix[k]];
      return { a: +r[0], b: +r[1], model: r[ix.model] || '', shown: r[ix.afdb] === '1', addr: ARCH_COLS.every((k) => Number.isFinite(ad[k])) ? ad : null, m }; });
  }).catch(() => []));
  return sp.vpairs.get(taxid);
}
async function virusStruct(sp, ia, ib) {
  if (!sp.viruses || ia == null || ib == null) return null;
  const v = sp.viruses.find((x) => (x.members || []).includes(ia) && (x.members || []).includes(ib)); if (!v) return null;
  const rows = await virusRows(sp, v.taxid);
  return rows.find((x) => (x.a === ia && x.b === ib) || (x.a === ib && x.b === ia)) || null;
}
// The link to a virus pair's model in LIVIA, wired: the database's copy, or a byte-range read of the release archive.
function structLink(host, x, cls) {
  if (!host) return; host.innerHTML = '';
  if (!x || !x.model || !(x.shown || x.addr)) return;
  const a = document.createElement('a'); a.className = cls; a.textContent = 'Structure in LIVIA ↗';
  if (x.shown) { a.href = `${LIVIA}dimer.html?id=${encodeURIComponent(x.model)}`; a.target = '_blank'; a.rel = 'noopener'; a.title = `${x.model} in LIVIA, from the AlphaFold Database`; }
  else { a.href = '#'; a.title = `${x.model} in LIVIA, read from the release archive at EBI (${((x.addr.cif_len + x.addr.pae_len) / 1048576).toFixed(1)} MB; not displayed by the AlphaFold Database)`;
    a.onclick = (e) => { e.preventDefault(); openFromArchive(x.model, x.addr, a); }; }
  host.appendChild(a);
}
async function openFromArchive(model, ad, link) {
  const w = window.open(`${LIVIA}universal.html?post=1`, '_blank'); if (!w) return;
  const txt = link.textContent; link.textContent = 'reading…';
  const get = async (o, n) => { for (let t = 0; ; t++) { try { const r = await fetch(ARCH_URL(ad.chunk), { headers: { Range: `bytes=${o}-${o + n - 1}` } });
    if (r.status !== 206) throw new Error(`the archive answered ${r.status}`); const b = new Uint8Array(await r.arrayBuffer()); if (b.length !== n) throw new Error('short read'); return b; }
    catch (e) { if (t >= 2) throw e; await new Promise((res) => setTimeout(res, 700 * (t + 1))); } } };   // EBI sometimes refuses a connection: retry
  try {
    const [[{ decompress }, { zipSync }], cif, pae] = await Promise.all([Promise.all([import('https://cdn.jsdelivr.net/npm/fzstd@0.1.1/+esm'), import('https://cdn.jsdelivr.net/npm/fflate@0.8.2/+esm')]),
      get(ad.cif_off, ad.cif_len), get(ad.pae_off, ad.pae_len)]);
    const zip = zipSync({ [`${model}-model_v1.cif`]: decompress(cif), [`${model}-predicted_aligned_error_v1.json`]: decompress(pae) });
    handTo(w, { type: 'livia-load', name: `${model}_viral.zip`, data: zip.buffer }); link.textContent = txt;
  } catch (e) { link.textContent = 'not read'; link.title = `The archive could not be read (${e.message || e}); try again.`; try { w.close(); } catch (_) {} }
}
// Hand data to a LIVIA tab opened with ?post=1 (cLIP, network): ping until it says it is ready, then post (its handshake).
function handTo(w, msg) {
  const origin = new URL(LIVIA, location.href).origin; let done = false, n = 0;
  const onMsg = (ev) => { if (ev.source !== w || !ev.data || ev.data.type !== 'livia-ready' || done) return; done = true; window.removeEventListener('message', onMsg); clearInterval(t); w.postMessage(msg, origin); };
  window.addEventListener('message', onMsg);
  const t = setInterval(() => { if (done || w.closed || ++n > 120) { clearInterval(t); window.removeEventListener('message', onMsg); return; } try { w.postMessage({ type: 'livia-ping' }, origin); } catch (e) { /* not loaded yet */ } }, 250);
}
// A split gene's bundle as one file again (the download): its rows from every isoform file, its constructs and sequences.
async function wholeBundle(ds, raw) {
  const subs = await Promise.all(raw.isoforms.choices.map((c) => bundleRaw(ds, c.file)));
  const texts = await Promise.all([raw, ...subs].map(async (r) => (await JSZip.loadAsync(r.bytes)).file(r.csvName).async('string')));
  const main = await JSZip.loadAsync(raw.bytes), zip = new JSZip(), nl = texts[0].indexOf('\n') + 1;
  zip.file(raw.csvName, texts[0].slice(0, nl) + texts.map((t) => t.slice(t.indexOf('\n') + 1)).join(''));
  for (const f of ['constructs.tsv', raw.faName].filter(Boolean)) zip.file(f, await main.file(f).async('string'));
  return zip.generateAsync({ type: 'blob', compression: 'DEFLATE', compressionOptions: { level: 6 } });
}
// Clustering choices in page order: the reference first when it holds at least 5% of the models (the page opens on it),
// then the others by models, ties by name; a nearly empty reference goes by its count like the rest.
function orderChoices(choices) {
  const total = choices.reduce((a, c) => a + c.n, 0), ref = choices.find((c) => c.id === '');
  choices.sort((x, y) => y.n - x.n || (x.id < y.id ? -1 : x.id > y.id ? 1 : 0));
  if (ref && ref.n >= 0.05 * total) { choices.splice(choices.indexOf(ref), 1); choices.unshift(ref); }
  return choices;
}
// The predictions of the given bundles merged into one view of the protein: its partners, its clustering choices, cLIP rows.
async function assemble(sp, P, parts, scope, onlyDi, setId) {
      const preds = [], runs = new Map(), seqs = new Map(), cons = new Map(), sets = [], TS = (parts.find((x) => x.TS) || {}).TS || null;
      const seqNo = new Map(), sno = (s) => { let i = seqNo.get(s); if (i == null) seqNo.set(s, i = seqNo.size); return i; };
      // Without a FASTA, a construct placed exactly on its gene's reference (no mutation) is that stretch of the reference,
      // whatever it is called: 'Cdk1' and 'FBgn0004106' folded as the same 297 residues are one sequence, not two isoforms.
      const seqId = (C, nm, len) => { const c = C && C.get(nm);
        return c && c.exact && c.off != null && !c.mut && c.kind !== 'mutant' && c.kind !== 'variant' ? `${c.key}@${c.off}+${c.len || len}` : `${nm}:${len}`; };
      const pairKey = (S, q, o, ql, ol, C) => { const x = S && S.get(q), y = S && S.get(o);   // the two sequences, the query's first
        return x && y ? sno(x) + '|' + sno(y) : `${seqId(C, q, ql)}|${seqId(C, o, ol)}`; };   // (their placement, or names and lengths, without a FASTA)
      for (const part of parts) {
        if (!part.raw.rows) Object.assign(part.raw, await csvRows(await JSZip.loadAsync(part.raw.bytes), part.raw.csvName));   // read again: an earlier view let them go
        const C = part.raw.cons;
        if (C) for (const [nm, c] of C) cons.set(nm, c);
        const keyOf = (nm) => { const c = C && C.get(nm); if (c) return c.key; const r = sp.byName.get(nm); return r ? r.key : nm; };
        for (const [nm, s] of part.raw.seqs) seqs.set(C ? nm : keyOf(nm), s);   // gene-keyed: the reference under the gene key, constructs by name
        const H = part.raw.H, hdr = part.raw.header.map(own), num = (r, c) => (H[c] == null || r[H[c]] === '' ? NaN : +r[H[c]]);
        for (const r of part.raw.rows) {
          const nm = r[H.name], at = nm.indexOf('___'), a = nm.slice(0, at), b = nm.slice(at + 3), qi = C ? keyOf(a) === P.key : a === part.name;
          if (!qi && (C ? keyOf(b) !== P.key : b !== part.name)) continue;
          const tags = part.TS && H.batch != null ? part.TS.runSets.get(+r[H.batch]) || [] : null;   // the run's thematic sets
          if (setId && !(tags && tags.some((t) => t.id === setId))) continue;
          const ps = primarySet(tags), qc = own(qi ? a : b), pc = own(qi ? b : a), key = keyOf(pc), set = ps ? ps.short : H.set != null ? own(r[H.set]) : '';
          const rid = part.di + '|' + nm + (H.batch != null ? '|' + r[H.batch] : '');
          if (!runs.has(rid)) runs.set(rid, { id: rid, di: part.di, qi, key, qc, pc, set, tags: tags ? tags.map((t) => t.id) : [] });
          if (set && !sets.includes(set)) sets.push(set);
          const s = (x, y) => (qi ? x : y);
          const p = { partner: key, run: rid, di: part.di, qi, qc, pc, set, tags: runs.get(rid).tags, rank: num(r, 'rank'), iLIS: num(r, 'iLIS'), iLIA: num(r, 'iLIA'), iLISA: num(r, 'iLISA'), ipTM: num(r, 'ipTM'),
            pTM: num(r, 'pTM'), LIS: num(r, 'LIS'), cLIS: num(r, 'cLIS'), LIA: num(r, 'LIA'), cLIA: num(r, 'cLIA'), ipSAE: num(r, 'ipSAE'), actifpTM: num(r, 'actifpTM'),
            qPl: num(r, s('pLDDT_i', 'pLDDT_j')), pPl: num(r, s('pLDDT_j', 'pLDDT_i')), qLIR: num(r, s('LIR_i', 'LIR_j')), pLIR: num(r, s('LIR_j', 'LIR_i')),
            qcLIR: num(r, s('cLIR_i', 'cLIR_j')), pcLIR: num(r, s('cLIR_j', 'cLIR_i')), qLen: num(r, s('len_i', 'len_j')), pLen: num(r, s('len_j', 'len_i')),
            qL: own(r[H[s('LIR_indices_i', 'LIR_indices_j')]]), pL: own(r[H[s('LIR_indices_j', 'LIR_indices_i')]]),
            qC: own(r[H[s('cLIR_indices_i', 'cLIR_indices_j')]]), pC: own(r[H[s('cLIR_indices_j', 'cLIR_indices_i')]]), row: null, hdr: null };
          if (p.iLIS >= CUT[10]) { p.row = r.map(own); p.hdr = hdr; }   // only rows past the lowest cutoff are ever clustered
          if (!Number.isFinite(p.iLISA)) p.iLISA = (p.iLIS || 0) * (p.iLIA || 0);
          p.sp = pairKey(part.raw.seqs, qi ? a : b, qi ? b : a, p.qLen, p.pLen, C);
          preds.push(p);
        }
        part.raw.rows = null;   // read: the predictions keep what they need (and the rows past the cutoff, below)
      }
      // A sequence pair counts once in a view (the atlas rule, 2026-09-27): the same two sequences folded in several runs
      // (two screens or sets, both chain orders, a re-run) keep the run with the highest iLIS for every residue view and
      // count; the other runs stay listed as repeats (the pair page shows their models) but add no contacts or models.
      const top = new Map();   // run → its best iLIS and sequence pair
      for (const p of preds) { const t = top.get(p.run); if (!t) top.set(p.run, { best: p.iLIS || 0, sp: p.sp }); else if ((p.iLIS || 0) > t.best) t.best = p.iLIS || 0; }
      const better = (x, y) => { const X = top.get(x), Y = top.get(y), rx = runs.get(x), ry = runs.get(y);   // ties: the earlier screen, the query first, the id
        return X.best !== Y.best ? X.best > Y.best : rx.di !== ry.di ? rx.di < ry.di : rx.qi !== ry.qi ? rx.qi : x < y; };
      const kept = new Map();
      for (const [rid, t] of top) { const k = kept.get(t.sp); if (k == null || better(rid, k)) kept.set(t.sp, rid); }
      for (const [rid, t] of top) { const k = kept.get(t.sp); runs.get(rid).repeatOf = k === rid ? null : k; }
      for (const p of preds) p.rep = !!runs.get(p.run).repeatOf;
      const aggregate = (list) => {   // predictions → one row per partner: its runs, its best and average scores
        const byP = new Map();
        for (const p of list) { if (!byP.has(p.partner)) byP.set(p.partner, []); byP.get(p.partner).push(p); }
        return [...byP].map(([key, ps]) => {
          const ids = [...new Set(ps.map((p) => p.run))].sort((x, y) => runs.get(x).di - runs.get(y).di || (runs.get(y).qi ? 1 : 0) - (runs.get(x).qi ? 1 : 0) || (x < y ? -1 : x > y ? 1 : 0));   // the same order however the rows were packed
          ps.sort((x, y) => ids.indexOf(x.run) - ids.indexOf(y.run) || x.rank - y.rank);
          const cs = ps.some((p) => !p.rep) ? ps.filter((p) => !p.rep) : ps;   // scores from the runs this view counts
          const il = cs.map((p) => p.iLIS || 0), ip = cs.map((p) => p.ipTM || 0);
          const rep = !ps.some((p) => !p.rep), of = rep ? runs.get(runs.get(ps[0].run).repeatOf) : null;   // every run a repeat: same sequences as another partner
          return { id: key, row: sp.byKey.get(key) || null, preds: ps, counted: cs, rep, repOf: of ? of.key : null, runs: ids, src: ps.reduce((m, p) => m | (1 << p.di), 0), best: Math.max(...il), avg: mean(il),
            ilisaBest: Math.max(...cs.map((p) => p.iLISA || 0)), iptmBest: Math.max(...ip), iptmAvg: mean(ip), contacts: Math.max(...cs.map((p) => p.qcLIR || 0)),
            sets: [...new Set(ps.map((p) => p.set).filter(Boolean))] };
        });
      };
      const partners = aggregate(preds);
      for (const pt of partners) for (const rid of pt.runs) { const ru = runs.get(rid); ru.twin = pt.runs.some((o) => o !== rid && runs.get(o).di === ru.di); }
      if (cons.size) for (const pt of partners) {   // a run of a gene-keyed screen: its category, and the constructs when they are not the genes themselves
        const og = (sp.byKey.get(pt.id) || {}).gene || pt.id, seen = new Map(), k = new Map();
        for (const rid of pt.runs) { const ru = runs.get(rid), qc = cons.get(ru.qc), pc = cons.get(ru.pc), ql = qc ? qc.label : P.gene, pl = pc ? pc.label : og;
          ru.base = (ru.set || sp.dsShort[ru.di]) + (ql !== P.gene || pl !== og ? ` · ${ql} – ${pl}` : ''); seen.set(ru.base, (seen.get(ru.base) || 0) + 1); }
        for (const rid of pt.runs) { const ru = runs.get(rid); const i = (k.get(ru.base) || 0) + 1; k.set(ru.base, i); ru.label = seen.get(ru.base) > 1 ? `${ru.base} · run ${i}` : ru.base; }
      }
      // What cLIP clusters. A gene-keyed screen (FlyPredictome) offers choices: the gene's reference sequence, with
      // every construct placed on it exactly (the full length, a trimmed or tiled part, a phosphosite window) or mapped
      // residue by residue (another isoform sharing 95% of its residues) — or any other construct, in its own numbering.
      // The default is the choice with the most models; point mutants and variants are never clustered. Other screens:
      // one construct, the reference-length one (UniProt) when a screen used it, else the one most predictions used.
      const R0 = P.len || 0;
      const onRef = (p) => { const c = cons.get(p.qc);
        return !!(c && c.len === p.qLen && c.kind !== 'mutant' && c.kind !== 'variant' && ((c.exact && c.off != null && c.off + p.qLen <= R0) || c.seg)); };
      const choices = [];
      if (cons.size && R0) {
        const own = new Map(); let refN = 0;
        for (const p of preds) { if (onRef(p)) { refN++; continue; } const c = cons.get(p.qc); if (c && c.kind !== 'mutant' && c.kind !== 'variant') own.set(p.qc, (own.get(p.qc) || 0) + 1); }
        if (refN) choices.push({ id: '', label: `${P.gene}, its reference (${fmtInt(R0)} aa)`, n: refN, len: R0 });
        for (const [name, n] of own) { const c = cons.get(name); choices.push({ id: name, label: /\(\d/.test(c.label) ? c.label : `${c.label} (${fmtInt(c.len)} aa)`, n, len: c.len }); }
        orderChoices(choices);
      }
      const toRef = (c, r) => { if (c.exact && c.off != null) return r + c.off; if (c.seg) for (const [a, b, n] of c.seg) if (r >= a && r < a + n) return b + r - a; return null; };
      const onRefRanges = (v, c) => rangesOf([...new Set(expand(v).map((r) => toRef(c, r)).filter((x) => x != null))].sort((x, y) => x - y));
      let fallback = null;
      if (!choices.length) {
        const full = (p) => { const c = cons.get(p.qc); return !c || c.kind === 'gene' || c.kind === 'isoform'; };
        const pool = preds.some(full) ? preds.filter(full) : preds;
        const lenN = new Map(); for (const p of pool) lenN.set(p.qLen, (lenN.get(p.qLen) || 0) + 1);
        const qL = [...lenN].sort((x, y) => (y[0] === P.len) - (x[0] === P.len) || y[1] - x[1])[0][0];
        fallback = { qLen: qL, inClip: (p) => p.qLen === qL && (pool === preds || full(p)) };
      }
      const labels = new Map(), labelOf = new Map();
      for (const pt of partners) pt.runs.forEach((rid, i) => { const lab = pt.runs.length > 1 || pt.id === P.key ? `${pt.id}~${i + 1}` : pt.id; labelOf.set(rid, lab); labels.set(lab, { key: pt.id, run: rid }); });
      for (const p of preds) p.label = labelOf.get(p.run);
      const clipFor = (choice) => {   // → the rows cLIP clusters for one choice, its axis, and what is left out
        const ref = choices.length > 0 && choice === '', inClip = !choices.length ? fallback.inClip : ref ? onRef : (p) => p.qc === choice;
        const qLen = !choices.length ? fallback.qLen : ref ? R0 : cons.get(choice).len;
        const rows = [], aside = new Map(), nameN = new Map(), bg = new Uint32Array(qLen + 2); let bgN = 0;
        for (const p of preds) {
          if (p.rep) continue;   // a repeat of a sequence pair adds no contacts
          if (!inClip(p)) { const c = cons.get(p.qc), k = c ? p.qc : 'screen ' + p.di; if (!aside.has(k)) aside.set(k, { di: p.di, len: p.qLen, con: c || null, n: 0 }); aside.get(k).n++; continue; }
          nameN.set(p.qc, (nameN.get(p.qc) || 0) + 1);
          if (Number.isFinite(p.iLIS) && p.iLIS < CUT[10]) {   // the background: every model below the 10% cutoff, its contacts on this axis
            bgN++; const c = ref ? cons.get(p.qc) : null, map = c && !(c.exact && c.off === 0);
            for (const r of expand(p.qC)) { const x = map ? toRef(c, r) : r; if (x != null && x >= 1 && x <= qLen) bg[x]++; }
            continue; }
          if (!(p.iLIS >= CUT[10]) || !p.row) continue;
          const o = {}; p.hdr.forEach((c, i) => { o[c] = p.row[i]; });
          if (ref) {   // the query side in the reference's numbering
            const q = p.qi ? 'i' : 'j', c = cons.get(p.qc);
            if (!(c.exact && c.off === 0)) { o['LIR_indices_' + q] = onRefRanges(o['LIR_indices_' + q], c); o['cLIR_indices_' + q] = onRefRanges(o['cLIR_indices_' + q], c); }
            o['len_' + q] = String(R0);
          }
          o.name = p.qi ? `${P.key}___${p.label}` : `${p.label}___${P.key}`;
          rows.push(o);
        }
        const qName = choices.length ? choice : [...nameN].sort((x, y) => y[1] - x[1]).map(([n]) => n)[0] || '';
        return { choice, rows, qLen, qName, aside, bg, bgN };
      };
      const C0 = clipFor(choices.length ? choices[0].id : '');
      const all = { parts, preds, counted: preds.filter((p) => !p.rep), partners, runs, seqs, cons, sets, TS, setId: scope, choices, clipFor, C0, clipRows: C0.rows, qLabel: P.key, labels, qLen: C0.qLen, qName: C0.qName, aside: C0.aside, iso: null };
      // One choice's predictions only (the reference with everything placed on it, or one other construct): a gene whose
      // isoforms were folded separately is read one isoform at a time, every card on the same predictions.
      const views = new Map(), others = new Set(choices.filter((c) => c.id).map((c) => c.id));
      all.only = (id) => {
        if (choices.length < 2) return all;
        if (!views.has(id)) {
          const ps = preds.filter((p) => (id ? p.qc === id : !others.has(p.qc))), v = { ...all, preds: ps, counted: ps.filter((p) => !p.rep), partners: aggregate(ps), iso: id };
          let C = null; const clip = () => (C ||= clipFor(id));   // built when this isoform is the one shown, not for the comparison table
          Object.defineProperties(v, { C0: { get: clip }, clipRows: { get: () => clip().rows }, qLen: { get: () => clip().qLen }, qName: { get: () => clip().qName }, aside: { get: () => clip().aside } });
          views.set(id, v);
        }
        return views.get(id);
      };
      return all;
}
// A gene split into isoform files (v1.2): the reference view comes from its own file; the other choices are listed from
// the file's isoforms.json (models in scope, for the switch and the Isoforms card) and each is read when it is chosen.
function isoformFiles(all, split, sp, P, scope, onlyDi, setId) {
  const S = split.raw.isoforms, ref = (all.choices || []).find((c) => c.id === '');
  const scoped = (t) => (setId ? (t.bySet || {})[setId] || { models: 0, partners: 0, p10: 0, p5: 0, p1: 0, top: [] } : t);   // a set-scoped page counts that set only
  const others = S.choices.filter((c) => scoped(c).models > 0).map((c) => ({ id: c.id, label: /\(\d/.test(c.label) ? c.label : `${c.label} (${fmtInt(c.len)} aa)`, n: scoped(c).models, len: c.len, file: c.file, stats: scoped(c) }));
  all.choices = orderChoices([...(ref ? [ref] : []), ...others]);
  all.split = { summary: S, part: split, others, allStats: scoped(S.all) };
  const views = new Map();
  all.load = (id) => {
    if (!id) return Promise.resolve(all);
    if (!views.has(id)) {
      const c = S.choices.find((x) => x.id === id); if (!c) return Promise.reject(new Error(`No isoform “${id}”.`));
      views.set(id, (async () => { const raw = await bundleRaw(split.ds, c.file), v = await assemble(sp, P, [{ ...split, name: c.file, raw }], scope, onlyDi, setId);
        return Object.assign(v, { iso: id, choices: all.choices, split: all.split }); })().catch((e) => { views.delete(id); throw e; }));
    }
    return views.get(id);
  };
  all.only = (id) => (id ? null : all);
}
const runLabel = (sp, B, rid, P) => { const ru = B.runs.get(rid); if (!ru) return ''; if (ru.label) return ru.label;
  const O = sp.byKey.get(ru.key), og = O ? O.gene : ru.key;
  return sp.dsShort[ru.di] + (ru.twin ? ` · ${ru.qi ? P.gene + '–' + og : og + '–' + P.gene}` : ''); };
// FlyPredictome's source categories (the paper's), as badges
const SET_COL = { 'Ligand–receptor': '#2B7A78', Kinase: '#B04A73', 'Kinase–TF': '#9A4A8F', 'Literature-derived': '#8C6D31', 'Large-scale proteomics': '#2B5F8E',
  'Subcellular organelle': '#5E7D2B', 'Inferred from orthologs': '#B5543C', Other: '#5B6B7F' };
const srcPapers = (sp, mask) => { const seen = new Set();   // the paper behind each screen in `mask`, linked
  return sp.dsIds.filter((_, di) => mask & (1 << di)).map((id) => ((REG && REG.datasets) || []).find((d) => d.id === id)).filter((d) => d && d.paper && !seen.has(d.source) && seen.add(d.source))
    .map((d) => `<a href="${esc(d.paper)}" target="_blank" rel="noopener">${esc(d.source.split(' · ')[0])} ↗</a>`).join(', '); };
const setBadges = (sets) => sets.map((s) => `<span class="src" style="--c:${SET_COL[s] || '#5B6B7F'}">${esc(s)}</span>`).join('');
const runColor = (sp, p) => (p.set ? SET_COL[p.set] || '#5B6B7F' : sp.dsColor[p.di]);
// A protein's predicted sequence: from a bundle FASTA when one carries it, else a screen's per-protein s/<name>.fa.
const SEQS = new Map();
function seqOf(sp, R, B) {
  if (!R) return Promise.resolve('');
  const s = B && B.seqs.get(R.key); if (s) return Promise.resolve(s);
  return (async () => {
    for (const o of R.occ || []) {
      let ds; try { ds = await dataset(sp.dsIds[o.di]); } catch (e) { continue; }
      const f = ds.manifest.files && ds.manifest.files.sequence; if (!f) continue;
      const rel = f.replace('{id}', o.name), ck = ds.id + '|' + rel;
      if (!SEQS.has(ck)) SEQS.set(ck, screenFile(ds, rel).then((r) => (r ? r.text() : '')).then((t) => [...parseFasta(t).values()][0] || '').catch(() => ''));
      const seq = await SEQS.get(ck); if (seq) return seq;
    }
    return '';
  })();
}

/* ── LIVIA's shared modules, loaded from the LIVIA site so the atlas and clip.html resolve, draw and export alike ── */
const LIVIA_JS = new Map();
function liviaScript(path) {
  if (!LIVIA_JS.has(path)) LIVIA_JS.set(path, new Promise((resolve, reject) => {
    const s = document.createElement('script'); s.src = LIVIA + path; s.onload = resolve;
    s.onerror = () => { LIVIA_JS.delete(path); reject(new Error(`Could not load ${path} from LIVIA.`)); };
    document.head.appendChild(s);
  }));
  return LIVIA_JS.get(path);
}
const liviaReady = () => Promise.all([liviaScript('js/clip-resolver.js'), liviaScript('js/livia-viewer.js')]);
const exportsReady = () => liviaScript('js/canvas2svg.js').then(() => liviaScript('js/livia-maps.js'));
const DOMS = new Map();   // accession → Promise<[{start, end, name}]> (UniProt Domain / DNA binding / Zinc finger, as clip.html)
const domainsOf = (acc) => { if (!DOMS.has(acc)) DOMS.set(acc, liviaReady().then(() => CLIPResolver.fetchDomains(acc)).catch(() => [])); return DOMS.get(acc); };
const AFDB = new Map();   // accession → Promise<{cifUrl, seq, amUrl} | null>
function afdbEntry(acc) {
  if (!AFDB.has(acc)) AFDB.set(acc, fetch(`https://alphafold.ebi.ac.uk/api/prediction/${encodeURIComponent(acc)}`).then((r) => {
    if (r.status === 404) return null; if (!r.ok) throw new Error(`AlphaFold DB answered ${r.status}`); return r.json(); }).then((d) => {
    const list = Array.isArray(d) ? d : d ? [d] : [];
    const e = list.find((x) => (x.uniprotAccession === acc) && (x.uniprotStart || x.sequenceStart || 1) === 1) || list[0];
    return e && e.cifUrl ? { cifUrl: e.cifUrl, seq: e.uniprotSequence || e.sequence || '', amUrl: e.amAnnotationsUrl || null } : null;
  }).catch(() => { AFDB.delete(acc); return { failed: true }; }));   // null: AFDB has no model; failed: the request did not get through (not remembered, so it can be retried)
  return AFDB.get(acc);
}
async function alphaMissense(url) {   // mean pathogenicity over all substitutions at each residue (as LIVIA's resolver)
  try {
    const lines = (await (await fetch(url)).text()).split('\n'), sum = [], cnt = [];
    for (let i = 1; i < lines.length; i++) { const c = lines[i].split(','), m = /^[A-Z](\d+)[A-Z]$/.exec(c[0]), p = parseFloat(c[1]); if (!m || isNaN(p)) continue; const k = +m[1]; sum[k] = (sum[k] || 0) + p; cnt[k] = (cnt[k] || 0) + 1; }
    const out = []; sum.forEach((s, k) => { out[k] = s / cnt[k]; }); return out.length ? out : null;
  } catch (e) { return null; }
}
// LIVIA's figure export bar (↓ SVG · ↓ PNG · W · H · font ×) under a canvas figure; the SVG re-runs the figure's own draw.
function attachExport(canvasId, name, redraw) {
  exportsReady().then(() => { if (document.getElementById(canvasId)) LiviaMaps.attachExportBar(canvasId, { name, svg: redraw, png: true }); }).catch(() => {});
}
// The same bar for figures that are SVG already (overview scatter, network): serialize for SVG, rasterize for PNG.
function svgExport(host, name, getSvg) {
  const old = host.querySelector(':scope > .xbar'); if (old) old.remove();
  const bar = el(`<div class="xbar"><button type="button" data-k="svg">↓ SVG</button><button type="button" data-k="png">↓ PNG</button>
    <span>W</span><input list="lm-exp-dim" placeholder="auto" aria-label="Export width in pixels"><span>H</span><input list="lm-exp-dim" placeholder="auto" aria-label="Export height in pixels"><span>font ×</span><input list="lm-exp-font" value="1" aria-label="Export font scale"></div>`);
  host.appendChild(bar);
  exportsReady().then(() => { if (!document.getElementById('lm-exp-dim')) { const mk = (id, vals) => document.body.appendChild(el(`<datalist id="${id}">${vals.map((v) => `<option value="${v}"></option>`).join('')}</datalist>`)); mk('lm-exp-dim', ['auto', '600', '900', '1200', '1600']); mk('lm-exp-font', ['1', '1.25', '1.5', '2', '0.8']); } }).catch(() => {});
  const [wI, hI, fI] = bar.querySelectorAll('input');
  bar.onclick = async (e) => {
    const b = e.target.closest('button'); if (!b) return;
    const svgEl = getSvg(); if (!svgEl) return;
    await exportsReady();
    const c = svgEl.cloneNode(true), W = svgEl.clientWidth || +svgEl.getAttribute('width'), H = svgEl.clientHeight || +svgEl.getAttribute('height'), f = +fI.value || 1;
    c.setAttribute('xmlns', 'http://www.w3.org/2000/svg'); c.setAttribute('width', W); c.setAttribute('height', H); c.setAttribute('font-family', 'IBM Plex Sans, Helvetica, Arial, sans-serif');
    const bg = document.createElementNS('http://www.w3.org/2000/svg', 'rect'); bg.setAttribute('width', W); bg.setAttribute('height', H); bg.setAttribute('fill', '#fff'); c.insertBefore(bg, c.firstChild);
    if (f !== 1) c.querySelectorAll('[font-size]').forEach((n) => n.setAttribute('font-size', (parseFloat(n.getAttribute('font-size')) * f).toFixed(2)));
    LiviaMaps.setExportOpts({ font: 1, width: +wI.value || 0, height: +hI.value || 0 });
    const svg = new XMLSerializer().serializeToString(c);
    if (b.dataset.k === 'svg') { LiviaMaps.downloadSVGFile(name, svg); return; }
    const out = LiviaMaps.applyExportOpts(svg), m = out.match(/<svg\b[^>]*\bwidth="([\d.]+)"[^>]*\bheight="([\d.]+)"/);
    const tw = m ? +m[1] : W, th = m ? +m[2] : H, scale = (+wI.value || +hI.value) ? 1 : 2;
    const img = new Image();
    img.onload = () => { const cv = document.createElement('canvas'); cv.width = Math.round(tw * scale); cv.height = Math.round(th * scale);
      const g = cv.getContext('2d'); g.drawImage(img, 0, 0, cv.width, cv.height);
      const a = document.createElement('a'); a.href = cv.toDataURL('image/png'); a.download = name.replace(/[^\w.-]+/g, '_') + '.png'; a.click(); };
    img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(out);
  };
}

/* ── search ──────────────────────────────────────────────────────────────────────────────────────────── */
function scoreProteins(sp, raw, limit = 10) {   // → [{ row, score }], best first; an exact-case symbol first (fly: tor is torso, Tor is TOR)
  const q = raw.trim().toLowerCase(); if (!q) return [];
  const scored = [];
  for (let i = 0; i < sp.rows.length; i++) {
    const k = sp.keys[i], r = sp.rows[i]; let s = 0;
    if (r.gene === raw.trim() || r.key === raw.trim()) s = 130;                        // exact case first: fly's Tor (mTor) is not tor (torso);
    else if ((' ' + (r.syn || '') + ' ').includes(' ' + raw.trim() + ' ')) s = 115;   // tiers stay further apart than the partner bonus (≤ 10)
    else if (k.gene === q || k.acc === q || k.key === q || k.ids.includes(q)) s = 100;
    else if (k.syn.includes(q)) s = 80;
    else if (k.gene.startsWith(q)) s = 60 - Math.min(20, k.gene.length - q.length);
    else if (k.ids.some((x) => x.startsWith(q))) s = 45;
    else if (q.length >= 3 && k.name.includes(q)) s = 25;
    else if (q.length >= 3 && k.syn.some((x) => x.startsWith(q))) s = 20;
    if (s) { const bonus = Math.min(10, Math.log10(1 + sp.rows[i].pos10) * 3); scored.push([s + bonus, i, s, bonus]); }
  }
  scored.sort((a, b) => b[0] - a[0]);
  return scored.slice(0, limit).map(([score, i, base, bonus]) => ({ row: sp.rows[i], score, base, bonus }));
}
function mountSearch(host, { big = false, spId = null, autofocus = false, only = null, set = '' } = {}) {   // spId null: every species; only: a set's keys
  host.innerHTML = `<div class="search ${big ? 'big' : ''}">
      <svg class="glass" width="${big ? 20 : 17}" height="${big ? 20 : 17}" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round"><circle cx="11" cy="11" r="7"/><path d="m20 20-3.6-3.6"/></svg>
      <input type="search" placeholder="${big ? (spId ? (spId === 'virus' ? 'Virus, protein name or accession' : 'Gene, accession or protein name') : 'Gene, UniProt accession, FlyBase ID, protein or virus') : 'Search a protein or virus'}" aria-label="Search proteins" autocomplete="off" spellcheck="false">
      <div class="suggest" hidden></div></div>`;
  const input = $('input', host), box = $('.suggest', host);
  let items = [], on = -1, timer = null;
  const go = (it) => { box.hidden = true; input.value = ''; location.hash = it.v ? `#/${it.sp.id}/taxon/${it.v.taxid}` : `#/${it.sp.id}/${it.row.key}${set ? '?set=' + encodeURIComponent(set) : ''}`; };
  const paint = () => { [...box.children].forEach((c, k) => c.classList.toggle('on', k === on)); };
  async function update() {
    const q = input.value, sps = spId ? [await species(spId)] : await Promise.all(((await registry()).species || []).map((x) => species(x.id).catch(() => null)));
    const many = sps.filter(Boolean).length > 1;
    const t = q.trim().toLowerCase(), al = VALIAS[t] || [], vir = (t.length < 3 && !al.length) || only ? [] : sps.filter((sp) => sp && sp.viruses).flatMap((sp) => sp.viruses.map((v) => {
      const nm = v.name.toLowerCase(), spn = (v.species || '').toLowerCase(), score = nm === t || al.includes(nm) || spn === t ? 3 : nm.startsWith(t) || spn.startsWith(t) ? 2 : nm.includes(t) ? 1 : 0; return { v, sp, score }; }))
      .filter((x) => x.score).sort((a, b) => b.score - a.score || b.v.hpos - a.v.hpos).slice(0, 4);   // a named virus first
    const hits = sps.filter(Boolean).flatMap((sp) => scoreProteins(sp, q, only ? 400 : 10).filter((h) => !only || only.has(h.row.key)).map((h) => ({ ...h, sp })));
    let prot;
    if (many) {   // across species an exact match counts whatever its case (the case rule separates genes within a species: fly's Tor
      // is not tor), so exact matches are ordered by how many partners pass ("p53": human TP53 before fly p53); each species' own
      // hits keep their order
      const x = (h) => (h.base >= 80 ? 100 : h.base) + h.bonus, bySp = new Map();
      for (const h of [...hits].sort((a, b) => b.score - a.score)) { if (!bySp.has(h.sp)) bySp.set(h.sp, []); bySp.get(h.sp).push(h); }
      prot = [...hits].sort((a, b) => x(b) - x(a)).map((h) => bySp.get(h.sp).shift());
    } else prot = hits.sort((a, b) => b.score - a.score);
    // a virus named exactly (or by a common name) first, then exact gene matches, then viruses named in part ("Tor" is a gene first)
    const exact = prot.filter((h) => h.base >= 80), rest = prot.filter((h) => h.base < 80);
    items = [...vir.filter((x) => x.score === 3), ...exact, ...vir.filter((x) => x.score < 3), ...rest].slice(0, 10);
    on = items.length ? 0 : -1;
    box.innerHTML = items.map(({ row: r, sp, v }) => (v ? `<div class="sg"><b>${esc(v.name)}</b><span class="nm">virus · ${fmtInt(v.n)} proteins</span><span class="ct">${fmtInt(v.hpos + v.mpos)} pairs</span>
        <span class="sub">${many ? `<span class="sp-tag">${esc(sp.reg.label)}</span>` : ''}its proteins' network · taxon ${v.taxid}</span></div>`
      : `<div class="sg"><b>${esc(r.gene)}</b><span class="nm">${esc(short(r.name))}</span><span class="ct">${fmtInt(r.pos10)} / ${fmtInt(r.partners)}</span>
        <span class="sub">${many ? `<span class="sp-tag">${esc(sp.reg.label)}</span>` : ''}${r.virus ? `${esc(r.virus.name)} · ` : ''}${esc(r.acc || '—')} · ${esc(r.id)}</span></div>`)).join('')
      || (q.trim() ? '<div class="sg-note">No match. Try a gene symbol, UniProt accession, FlyBase ID or protein name.</div>' : '');
    box.hidden = !box.innerHTML;
    [...box.querySelectorAll('.sg')].forEach((d, k) => d.onclick = () => go(items[k]));
    paint();
  }
  input.addEventListener('input', () => { clearTimeout(timer); timer = setTimeout(update, 60); });
  input.addEventListener('keydown', (e) => {
    if (e.key === 'ArrowDown') { on = Math.min(items.length - 1, on + 1); paint(); e.preventDefault(); }
    else if (e.key === 'ArrowUp') { on = Math.max(0, on - 1); paint(); e.preventDefault(); }
    else if (e.key === 'Enter' && items[on >= 0 ? on : 0]) { go(items[on >= 0 ? on : 0]); }
    else if (e.key === 'Escape') { box.hidden = true; }
  });
  const outside = (e) => { if (!host.isConnected) { document.removeEventListener('click', outside); return; } if (!host.contains(e.target)) box.hidden = true; };
  document.addEventListener('click', outside);   // gone with the page
  if (autofocus) setTimeout(() => input.focus(), 50);
  return input;
}

/* ── shared drawing: cluster palette (LIVIA cLIP's 'auto'), residue axes, interface tracks, sequences ──── */
const TAB10 = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'];
const TAB20 = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5'];
const TAB20B = ['#393b79', '#5254a3', '#6b6ecf', '#9c9ede', '#637939', '#8ca252', '#b5cf6b', '#cedb9c', '#8c6d31', '#bd9e39', '#e7ba52', '#e7cb94', '#843c39', '#ad494a', '#d6616b', '#e7969c', '#7b4173', '#a55194', '#ce6dbd', '#de9ed6'];
const TAB20C = ['#3182bd', '#6baed6', '#9ecae1', '#c6dbef', '#e6550d', '#fd8d3c', '#fdae6b', '#fdd0a2', '#31a354', '#74c476', '#a1d99b', '#c7e9c0', '#756bb1', '#9e9ac8', '#bcbddc', '#dadaeb', '#636363', '#969696', '#bdbdbd', '#d9d9d9'];
const TAB60 = TAB20.concat(TAB20B, TAB20C);
function hslDistinct(i) {   // hex, so Mol* (MVS) can take it too
  const h = (i * 137.508) % 360, s = (62 + (i % 2) * 16) / 100, l = (42 + (i % 3) * 13) / 100, a = s * Math.min(l, 1 - l);
  const f = (n) => { const k = (n + h / 30) % 12, c = l - a * Math.max(-1, Math.min(k - 3, 9 - k, 1)); return Math.round(255 * c).toString(16).padStart(2, '0'); };
  return `#${f(0)}${f(8)}${f(4)}`;
}
const clusterColor = (id, k) => { const pal = k <= 10 ? TAB10 : k <= 20 ? TAB20 : k <= 60 ? TAB60 : null; return pal ? pal[(id - 1) % pal.length] : hslDistinct(id - 1); };
const clusterLabel = (c, brief) => (brief ? 'C' : 'Cluster ') + c;
function amCol(v) {   // AlphaMissense: benign (blue) → ambiguous (gray) → pathogenic (red), as clip.html
  const t = Math.max(0, Math.min(1, v)), lerp = (a, b, u) => [0, 1, 2].map((i) => Math.round(a[i] + (b[i] - a[i]) * u));
  const c = t < 0.5 ? lerp([44, 123, 182], [224, 224, 224], t / 0.5) : lerp([224, 224, 224], [215, 25, 28], (t - 0.5) / 0.5);
  return '#' + c.map((x) => x.toString(16).padStart(2, '0')).join('');
}
const plddtCol = (b) => (b > 90 ? '#0053D6' : b > 70 ? '#65CBF3' : b > 50 ? '#FFDB13' : '#FF7D45');
const PAIR_COL = { q: { base: '#E0E0E0', lir: '#80CBC4', clir: '#00897B' }, p: { base: '#E0E0E0', lir: '#FFAB91', clir: '#E64A19' } };

// A 2D context sized in CSS px. On screen it is scaled by the device pixel ratio; while LIVIA's exporter re-runs a draw
// against canvas2svg (which ignores transforms), it is drawn at 1× so the SVG comes out at the on-screen size.
function canvasCtx(cv, W, H) {
  let g = cv.getContext('2d');
  const svg = typeof CanvasRenderingContext2D !== 'undefined' && !(g instanceof CanvasRenderingContext2D), dpr = svg ? 1 : (window.devicePixelRatio || 1);
  cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); cv.style.width = W + 'px'; cv.style.height = H + 'px';
  if (svg) g = cv.getContext('2d');
  g.setTransform(dpr, 0, 0, dpr, 0, 0); g.clearRect(0, 0, W, H); g.__dpr = dpr;
  return g;
}
function runs(fp) { const s = [...fp].sort((a, b) => a - b), out = []; if (!s.length) return out; let a = s[0], p = s[0];
  for (let i = 1; i < s.length; i++) { if (s[i] === p + 1) p = s[i]; else { out.push([a, p]); a = p = s[i]; } } out.push([a, p]); return out; }
// Residue ticks: 1, round steps, and the last residue — one set for every residue axis on a page, so plots line up.
// `want` (the page's x-ticks box) asks for about that many ticks; blank = fit the width.
function resTicks(L, plotW, want) {
  want = want || Math.max(3, Math.min(14, Math.floor(plotW / 64)));
  const step = [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 2500, 5000].find((s) => L / s <= want) || 10000;
  const t = [1]; for (let r = step; r < L; r += step) if (r - 1 > step * 0.35 && L - r > step * 0.35) t.push(r);
  if (L > 1) t.push(L); return t;
}
function drawTicks(g, ticks, y, xc, W) {   // xc(r): tick x for residue r; labels stay inside the canvas
  g.save(); g.fillStyle = '#5A697C'; g.strokeStyle = '#B9C4CF'; g.lineWidth = 1; g.font = '10.5px "IBM Plex Mono", ui-monospace, monospace'; g.textBaseline = 'top'; g.textAlign = 'left';
  for (const r of ticks) { const x = Math.round(xc(r)) + 0.5, s = String(r), w = g.measureText(s).width;
    g.beginPath(); g.moveTo(x, y); g.lineTo(x, y + 4); g.stroke(); g.fillText(s, Math.max(1, Math.min(W - w - 1, x - w / 2)), y + 6); }
  g.restore();
}
function lanes(items, x0Of, x1Of) {   // stack overlapping boxes into lanes; returns lane count
  const end = []; for (const d of items) { const x0 = x0Of(d); let l = 0; while (l < end.length && end[l] > x0 + 0.5) l++; d.lane = l; end[l] = x1Of(d); } return end.length;
}
// Interaction and contact residues of one prediction on both proteins: gray = neither, light = LIR (interaction), dark = cLIR (contact); UniProt domains above.
function drawIfaceTracks(cv, tracks) {
  const W = cv.parentElement.clientWidth, pad = 8, gap = 18;
  for (const t of tracks) { t.x = (r) => pad + (r - 1) / t.len * (W - 2 * pad); t.nl = t.doms.length ? lanes(t.doms, (d) => t.x(d.start), (d) => t.x(d.end + 1)) : 0; }
  const bh = (t) => 22 + t.nl * 15 + (t.nl ? 4 : 0) + 22 + 20;
  const H = tracks.reduce((s, t) => s + bh(t), 0) + gap * (tracks.length - 1);
  const g = canvasCtx(cv, W, H);
  let y = 0; cv._blocks = []; cv._doms = [];
  for (const t of tracks) {
    const L = t.len, bw = Math.max(1.3, (W - 2 * pad) / L), name = t.label || t.gene;
    g.font = '600 13.5px "IBM Plex Sans", system-ui, sans-serif'; g.fillStyle = t.col.clir; g.textBaseline = 'alphabetic'; g.textAlign = 'left'; g.fillText(name, pad, y + 15);
    const gw = g.measureText(name).width; g.font = '12px "IBM Plex Mono", ui-monospace, monospace'; g.fillStyle = '#5A697C';
    const extent = t.span ? `residues ${fmtInt(t.span[0])}–${fmtInt(t.span[1])} of ${fmtInt(L)}` : t.own ? `${fmtInt(L)} aa construct, its own numbering` : `${fmtInt(L)} aa`;
    g.fillText(`${extent} · ${t.lir.length} interaction · ${t.clir.length} contact`, pad + gw + 10, y + 15);
    y += 22;
    for (const d of t.doms) { const x0 = t.x(Math.max(1, d.start)), x1 = Math.max(x0 + 2, t.x(Math.min(L, d.end) + 1)), yy = y + d.lane * 15;
      g.fillStyle = '#E3E9F1'; g.fillRect(x0, yy, x1 - x0, 13); g.strokeStyle = '#9FB0C4'; g.lineWidth = 0.6; g.strokeRect(x0 + 0.3, yy + 0.3, x1 - x0 - 0.6, 12.4);
      g.font = '10.5px "IBM Plex Sans", system-ui, sans-serif'; g.fillStyle = '#34445A'; let s = d.name; while (s.length > 2 && g.measureText(s).width > x1 - x0 - 6) s = s.slice(0, -2) + '…';
      if (s.length > 2 && g.measureText(s).width <= x1 - x0 - 6) g.fillText(s, x0 + 3, yy + 10); cv._doms.push({ d, x0, x1, y0: yy, y1: yy + 13 }); }
    y += t.nl * 15 + (t.nl ? 4 : 0);
    const by = y;
    if (t.span) {   // a fragment or window on its whole gene: the gene faint, the folded part in the usual gray
      g.fillStyle = '#F2F4F7'; g.fillRect(pad, y + 3, W - 2 * pad, 14);
      g.fillStyle = t.col.base; g.fillRect(t.x(t.span[0]), y + 3, Math.max(bw, t.x(t.span[1] + 1) - t.x(t.span[0])), 14);
    } else { g.fillStyle = t.col.base; g.fillRect(pad, y + 3, W - 2 * pad, 14); }
    g.fillStyle = t.col.lir; for (const r of t.lir) g.fillRect(t.x(r), y + 3, bw, 14);
    g.fillStyle = t.col.clir; for (const r of t.clir) g.fillRect(t.x(r), y, bw, 20);
    cv._blocks.push({ t, y0: by, y1: by + 20 });
    y += 22;
    drawTicks(g, resTicks(L, W - 2 * pad), y, (r) => t.x(r) + bw / 2, W);
    y += 20 + gap;
  }
  cv.onmousemove = (e) => { const b = cv.getBoundingClientRect(), x = e.clientX - b.left, yy = e.clientY - b.top;
    const hd = cv._doms.find((k) => x >= k.x0 && x <= k.x1 && yy >= k.y0 && yy <= k.y1); if (hd) return showTip(domTip(hd.d), e.clientX, e.clientY);   // a domain box
    const blk = cv._blocks.find((k) => yy >= k.y0 - 4 && yy <= k.y1 + 4); if (!blk) return hideTip();
    const t = blk.t, r = Math.floor((x - pad) / (W - 2 * pad) * t.len) + 1; if (r < 1 || r > t.len) return hideTip();
    const st = t.span && (r < t.span[0] || r > t.span[1]) ? 'not in the folded construct' : t.clirSet.has(r) ? 'contact (cLIR)' : t.lirSet.has(r) ? 'interaction (LIR)' : 'not an interaction residue';
    const dom = t.doms.filter((d) => r >= d.start && r <= d.end).map((d) => d.name).join(', ');
    showTip(`<b>${esc(t.gene)}</b> ${t.seq && t.seq[r - 1] ? t.seq[r - 1] : ''}${r} · ${st}${dom ? `<br>${esc(dom)}` : ''}`, e.clientX, e.clientY); };
  cv.onmouseleave = hideTip;
}
// A sequence fills its line: the font (11–14.5 px) is chosen so a whole number of 10-residue groups spans the width,
// on a phone as on a wide screen, instead of leaving a ragged gap at the right.
const SEQFIT = typeof ResizeObserver === 'function' ? new ResizeObserver((ents) => ents.forEach((e) => fitSeq(e.target))) : null;
function fitSeq(f) {
  const W = f.clientWidth; if (!W || !f.isConnected) return;
  const probe = document.createElement('span'); probe.textContent = 'MMMMMMMMMM'; probe.style.cssText = 'position:absolute;visibility:hidden;white-space:nowrap;font-size:100px';
  f.appendChild(probe); const unit = probe.getBoundingClientRect().width / 1000; probe.remove();   // a character's width per px of font size
  const per = 10 * unit + 0.7;                                   // one group and its 0.7em gap, per px of font size
  let pick = null;
  for (let n = 1; n <= 20; n++) { const px = Math.floor((W - 2) / n / per * 10) / 10; if (px >= 11 && px <= 14.5 && (!pick || Math.abs(px - 12.5) < Math.abs(pick - 12.5))) pick = px; }
  if (pick && f.style.fontSize !== pick + 'px') f.style.fontSize = pick + 'px';
}
function fitSeqs(root) { root.querySelectorAll('.seq-flow').forEach((f) => { fitSeq(f); if (SEQFIT) SEQFIT.observe(f); }); }
// One sequence as a continuous, searchable flow in 10-residue groups (position numbers are CSS, not text), as clip.html.
function seqPanel(label, seq, lir, clir, col, len, span) {   // span: the folded part of a longer gene; the rest is dimmed
  const ext = span ? `residues ${fmtInt(span[0])}–${fmtInt(span[1])} of ${fmtInt(len)}` : `${fmtInt(len || seq.length)} residues`;
  const TXT = { '#00897B': '#00776B', '#E64A19': '#C13E15' };   // the chain colors as text (at least 4.5:1 on white)
  const head = `<div class="seqh"><b style="color:${TXT[String(col.clir).toUpperCase()] || col.clir}">${esc(label)}</b> <span class="muted">${ext} · ${lir.length} interaction · ${clir.length} contact</span></div>`;
  if (!seq || (len && seq.length !== len)) return head + '<p class="muted" style="margin:6px 0 0">The predicted sequence is not available for this construct.</p>';
  const L = new Set(lir), C = new Set(clir); let body = '';
  for (let i = 0; i < seq.length; i++) {
    if (i % 10 === 0) { if (i) body += '</span><wbr>'; const n = Math.min(i + 10, seq.length);   // a 1–2 residue last group: its number runs right, clear of the one before
      body += `<span class="g10${n - i < 3 ? ' tail' : ''}" data-n="${n}">`; }
    const r = i + 1, ch = seq[i];
    body += C.has(r) ? `<span class="c" style="background:${col.clir}">${ch}</span>` : L.has(r) ? `<span class="l" style="background:${col.lir}">${ch}</span>`
      : span && (r < span[0] || r > span[1]) ? `<span class="out">${ch}</span>` : ch;
  }
  return head + `<div class="seq-flow">${body}</span></div>`;
}
// Interaction Residues for one prediction (protein page card and pair page), with LIVIA's export bar under the map.
// A construct placed on its gene's reference sequence (a fragment, a phosphosite window) is drawn on the whole gene in
// reference numbering, its extent shaded; a construct that is not placed keeps its own numbering, and says so.
async function ifaceView(host, { sp, P, O, pred, B, canvasId }) {
  const tok = (host._tok = (host._tok || 0) + 1);
  const side = (R, name, len) => {
    const c = B.cons && B.cons.get(name), label = c ? c.label : R.gene;
    const placed = c && c.off != null && c.len === len && R.len && c.off + len <= R.len;
    if (!c || (placed && c.off === 0 && len === R.len)) return { label, len: len || R.clen || 1, shift: 0, span: null, own: false, name };
    if (placed) return { label, len: R.len, shift: c.off, span: [c.off + 1, c.off + len], own: false, name };
    return { label, len: len || 1, shift: 0, span: null, own: true, name };
  };
  const q = side(P, pred.qc, pred.qLen), o = side(O, pred.pc, pred.pLen), sh = (s, k) => expand(s).map((r) => r + k);
  const qL = sh(pred.qL, q.shift), qC = sh(pred.qC, q.shift), pL = sh(pred.pL, o.shift), pC = sh(pred.pC, o.shift);
  const none = !qL.length && !pL.length && !qC.length && !pC.length;
  host.innerHTML = `${none ? `<p class="note">This model has no residue pair with PAE ≤ 12 Å between the two proteins, so it has no interaction residues.</p>` : ''}
    <div class="ifmap"><canvas id="${canvasId}"></canvas></div><div class="seqpanels"><div class="seqp"></div><div class="seqp"></div></div>`;
  const cv = $('canvas', host);
  const T = [{ gene: P.gene, label: q.label, len: q.len, span: q.span, own: q.own, lir: qL, clir: qC, lirSet: new Set(qL), clirSet: new Set(qC), col: PAIR_COL.q, doms: [] },
    { gene: O.gene, label: o.label, len: o.len, span: o.span, own: o.own, lir: pL, clir: pC, lirSet: new Set(pL), clirSet: new Set(pC), col: PAIR_COL.p, doms: [] }];
  host._redraw = () => drawIfaceTracks(cv, T);
  host._redraw();
  attachExport(canvasId, `atlas_${P.gene}_vs_${O.gene}_residues`, host._redraw);
  const seqFor = (R, s) => (s.own ? Promise.resolve(B.seqs.get(s.name) || '') : seqOf(sp, R, B));   // the reference, or the construct's own
  const [qs, os] = await Promise.all([seqFor(P, q), seqFor(O, o)]);
  if (host._tok !== tok) return;
  T[0].seq = qs; T[1].seq = os;
  const panels = host.querySelectorAll('.seqp');
  panels[0].innerHTML = seqPanel(q.label, qs, qL, qC, PAIR_COL.q, T[0].len, q.span);
  panels[1].innerHTML = seqPanel(o.label, os, pL, pC, PAIR_COL.p, T[1].len, o.span);
  fitSeqs(panels[0]); fitSeqs(panels[1]);
  // UniProt domain coordinates fit only the UniProt sequence: a full-length human construct; for fly, a FlyBase
  // reference that is UniProt's sequence (checked against the AlphaFold DB entry)
  const domOK = async (R, s, seq) => { if (!R || !R.acc || s.own) return false;
    if (!sp.manifest.keyedBy) return !!(R.len && R.len === s.len);
    const e = await afdbEntry(R.acc); return !!(e && !e.failed && seq && e.seq === seq); };
  const [qd, od] = await Promise.all([domOK(P, q, qs).then((ok) => (ok ? domainsOf(P.acc) : [])), domOK(O, o, os).then((ok) => (ok ? domainsOf(O.acc) : []))]);
  if (host._tok !== tok) return;
  T[0].doms = (qd || []).map((d) => ({ ...d })); T[1].doms = (od || []).map((d) => ({ ...d }));
  if (T[0].doms.length || T[1].doms.length) host._redraw();
}
function pearson(xs, ys) { const n = xs.length; if (n < 3) return NaN; let mx = 0, my = 0; for (let i = 0; i < n; i++) { mx += xs[i]; my += ys[i]; } mx /= n; my /= n;
  let sxy = 0, sxx = 0, syy = 0; for (let i = 0; i < n; i++) { const dx = xs[i] - mx, dy = ys[i] - my; sxy += dx * dy; sxx += dx * dx; syy += dy * dy; } return sxy / Math.sqrt(sxx * syy); }
function spearman(xs, ys) { const n = xs.length; if (n < 3) return NaN;
  const rank = (a) => { const idx = a.map((v, i) => [v, i]).sort((p, q) => p[0] - q[0]), r = new Array(n); let i = 0;
    while (i < n) { let j = i; while (j + 1 < n && idx[j + 1][0] === idx[i][0]) j++; const av = (i + j) / 2 + 1; for (let k = i; k <= j; k++) r[idx[k][1]] = av; i = j + 1; } return r; };
  return pearson(rank(xs), rank(ys)); }
// Axes on a canvas (so the whole plot exports through canvas2svg): frame, ticks, labels, titles.
function canvasAxes(g, xs, ys, m, W, H, xTitle, yTitle) {
  g.save(); g.strokeStyle = '#D5DDE6'; g.lineWidth = 1;
  g.beginPath(); g.moveTo(m.l + 0.5, m.t); g.lineTo(m.l + 0.5, H - m.b + 0.5); g.lineTo(W - m.r, H - m.b + 0.5); g.stroke();
  g.fillStyle = '#5A697C'; g.font = '11px "IBM Plex Mono", ui-monospace, monospace';
  const xt = xs.ticks(Math.max(3, Math.floor((W - m.l - m.r) / 110))), xf = xs.tickFormat(xt.length);
  g.textAlign = 'center'; g.textBaseline = 'top';
  for (const t of xt) { const x = Math.round(xs(t)) + 0.5; g.beginPath(); g.moveTo(x, H - m.b); g.lineTo(x, H - m.b + 5); g.stroke(); g.fillText(xf(t), x, H - m.b + 8); }
  const yt = ys.ticks(6), yf = ys.tickFormat(yt.length);
  g.textAlign = 'right'; g.textBaseline = 'middle';
  for (const t of yt) { const y = Math.round(ys(t)) + 0.5; g.beginPath(); g.moveTo(m.l - 5, y); g.lineTo(m.l, y); g.stroke(); g.fillText(yf(t), m.l - 8, y); }
  g.fillStyle = '#34445A'; g.font = '12.5px "IBM Plex Sans", system-ui, sans-serif'; g.textAlign = 'center'; g.textBaseline = 'alphabetic';
  g.fillText(xTitle, (m.l + W - m.r) / 2, H - 8);
  g.translate(15, (m.t + H - m.b) / 2); g.rotate(-Math.PI / 2); g.textBaseline = 'middle'; g.fillText(yTitle, 0, 0);
  g.restore();
}

/* ── views ───────────────────────────────────────────────────────────────────────────────────────────── */
const TRY = { human: ['TP53', 'MDM2', 'CTNNB1', 'MAPK3', 'SMAD4', 'KRAS'], fly: ['arm', 'dsh', 'Ras85D', 'N', 'yki', 'Akt'],   // example proteins on the home page
  zebrafish: ['ksr1a', 'map3k7'], yeast: ['PRR2', 'YAK1'], worm: ['ksr-1', 'par-1'],
  virus: [['Herpes simplex virus 1', 'taxon/10299'], ['Epstein-Barr virus', 'taxon/10377'], ['Poliovirus', 'taxon/12081'], ['Mpox virus', 'taxon/619591']] };   // viruses open their network page; each shows reported complexes at the default 1% FPR
async function viewHome() {
  // Totals come from the species manifests alone; the species indexes (search) load once the page is idle.
  const gen = ROUTE, reg = await registry(), all = await Promise.all((reg.species || []).map((x) => speciesManifest(x.id)));
  if (stale(gen)) return;
  const tot = (k) => all.reduce((s, m) => s + (m.counts[k] || 0), 0), nScreens = all.reduce((s, m) => s + (m.datasets || []).length, 0);
  app.innerHTML = `
    <section class="hero hero-center">
      <h1>Where does <span class="ini">each&nbsp;partner</span> bind?</h1>
      <p class="lede">AlphaFold-Multimer predicts not only whether two proteins bind, but through which residues. Pooled across large-scale
        screens, those interfaces map a protein's many partners onto its sequence. On every protein page, <a href="${LIVIA}clip.html" target="_blank" rel="noopener">cLIP</a>
        (<span class="q">c</span>lustered <span class="q">L</span>ocal <span class="q">I</span>nteraction <span class="q">P</span>rofiler) groups the partners by the residues they contact, so partners that share a binding site fall into one cluster.<br>Search a protein
        to see who is predicted to bind it, how confidently, and where.</p>
      <p class="lede-src">cLIP is part of <a href="${LIVIA}" target="_blank" rel="noopener">LIVIA</a>, a browser-based tool for assessing and visualizing predicted protein interactions
        (Kim &amp; Perrimon, 2026, <a href="https://doi.org/${REF.livia[1]}" target="_blank" rel="noopener">bioRxiv</a>).</p>
      <div id="home-search" class="hero-search"></div>
      <div class="totals"><span><b>${fmtInt(tot('runs'))}</b> predictions</span><span><b>${fmtInt(tot('predictions'))}</b> models</span><span><b>${fmtInt(tot('pairs'))}</b> protein pairs</span>
        <span><b>${fmtInt(tot('proteins'))}</b> proteins</span><span><b>${fmtInt(nScreens)}</b> screen${nScreens === 1 ? '' : 's'}</span></div>
      <div class="chips">${(reg.species || []).map((x, i) => `<span class="chip-group">${i ? '' : '<span class="lbl">Try</span>'}${reg.species.length > 1 ? `<a class="lbl sp-link" href="#/${x.id}" title="every ${esc(x.label.toLowerCase())} protein, screen and network">${esc(x.label)}</a>` : ''}`
        + (TRY[x.id] || []).map((g) => (Array.isArray(g) ? `<a class="chip" href="#/${x.id}/${g[1]}">${esc(g[0])}</a>` : `<a class="chip" href="#/${x.id}/${encodeURIComponent(g)}">${esc(g)}</a>`)).join('') + '</span>').join('')}</div>
      <div class="showcase" id="showcase" aria-roledescription="carousel" aria-label="Example proteins"></div>
    </section>`;
  mountSearch($('#home-search'), { big: true, autofocus: true });
  showcase();
  const idle = window.requestIdleCallback || ((f) => setTimeout(f, 1500));
  idle(() => { if (!stale(gen)) for (const x of reg.species || []) species(x.id).catch(() => {}); }, { timeout: 4000 });   // search is instant by the first keystroke
}
async function fillThemes() {   // home: each theme's species and totals, from its members' counts
  const reg = await registry(), box = $('#themes'); if (!box) return;
  const cards = await Promise.all((reg.themes || []).map(async (T) => { const rows = await themeMembers(T), sum = (k) => rows.reduce((a, r) => a + (r.counts[k] || 0), 0);
    return `<div class="ds live"><span class="badge on">Theme</span><h3><a href="#/themes/${T.id}">${esc(T.title)}</a></h3>
      <div class="sp">${rows.map((r) => esc(r.S.label)).join(' · ')}</div>
      <div class="stats"><div><b>${fmtInt(sum('proteins'))}</b><span>proteins</span></div><div><b>${fmtInt(sum('pairs'))}</b><span>pairs</span></div><div><b>${fmtInt(sum('pairsFpr10'))}</b><span>past 10% FPR${cutNote(10)}</span></div></div>
      <div class="src-line">${esc(T.about || '')}</div></div>`; }));
  box.innerHTML = cards.join('');
}
function dsCard(d) {
  const live = d.status === 'live', spx = ((REG && REG.species) || []).find((s) => s.id === d.species);
  const head = `<span class="badge ${live ? 'on' : 'soon'}">${live ? 'Searchable' : d.status === 'external' ? 'Separate site' : 'Planned'}</span>
    <h3>${live ? `<a href="#/datasets/${d.id}">${esc(d.title)}</a>` : d.url ? `<a href="${esc(d.url)}" target="_blank" rel="noopener">${esc(d.title)}</a>` : esc(d.title)}</h3>
    <div class="sp">${esc(spx ? spx.name : d.species)}${live && d.short ? ` · <span class="src" style="--c:${d.color}">${esc(d.short)}</span>` : ''}</div>${live && d.models ? `<div class="sp">${runSettings(d)}</div>` : ''}`;
  const src = d.paper ? `<a href="${esc(d.paper)}" target="_blank" rel="noopener">${esc(d.source)} ↗</a>` : esc(d.source);   // every paper reference links to the paper
  return `<div class="ds ${live ? 'live' : ''}" data-ds="${d.id}">${head}${live ? '<div class="stats"></div>' : ''}
    <div class="src-line">${src}${d.note ? ` — ${esc(d.note)}` : ''}</div></div>`;
}
async function fillDsStats() {
  for (const d of (await registry()).datasets.filter((x) => x.status === 'live')) {
    let s; try { s = await dataset(d.id); } catch (e) { continue; }
    const k = s.manifest.counts, box = app.querySelector(`[data-ds="${d.id}"] .stats`);
    if (box) box.innerHTML = `<div><b>${fmtInt(k.proteins)}</b><span>proteins</span></div><div><b>${fmtInt(k.pairs)}</b><span>pairs</span></div><div><b>${fmtInt(k.pairsFpr10)}</b><span>past 10% FPR${cutNote(10)}</span></div>`;
    const TS = await setsOf(s).catch(() => null), scr = TS ? TS.list.filter((x) => x.type === 'screen') : [];   // its named screens, as thematic sets
    if (box && scr.length && !box.parentElement.querySelector('.setchips')) box.insertAdjacentHTML('afterend', `<div class="setchips"><span>Sets</span>${scr.map((x) =>
      `<a class="src" style="--c:${x.color}" href="#/datasets/${d.id}/${x.id}">${esc(x.short)}</a>`).join('')}<a href="#/datasets/${d.id}">all ${TS.list.length} ›</a></div>`);
  }
}
// The home banner: example proteins from every species, precomputed by LIVIA cLIP (atlas build/make_showcase.js), one
// every 7 s, each with its contact residue frequency over its clustered fingerprints. Hover or focus pauses it; with
// reduced motion it only moves when asked.
async function showcase() {
  const gen = ROUTE, box = $('#showcase'); if (!box) return;
  let data; try { data = await getJSON('data/showcase.json'); } catch (e) { box.remove(); return; }
  const S = (data && data.slides) || []; if (stale(gen) || !S.length) { if (!S.length) box.remove(); return; }
  box.innerHTML = `<div class="sc-top"><div class="sc-name"><a id="sc-gene"></a><span class="sc-sp" id="sc-sp"></span></div><span id="sc-note"></span></div>
    <div class="sc-stage"><canvas id="sc-cv"></canvas></div>
    <div class="sc-foot"><span id="sc-cap"></span><a id="sc-open"></a></div>
    <div class="sc-own">Have your own screen? <a href="${LIVIA}clip.html" target="_blank" rel="noopener">Use cLIP in LIVIA →</a></div>
    <div class="sc-nav"><button type="button" class="sc-arrow" id="sc-prev" aria-label="Previous protein"><svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"><path d="m15 18-6-6 6-6"/></svg></button>
      <div class="sc-dots">${S.map((s, k) => `<button type="button" data-k="${k}" aria-label="${esc(s.gene)}, ${esc(s.spLabel)}"></button>`).join('')}</div>
      <button type="button" class="sc-arrow" id="sc-next" aria-label="Next protein"><svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"><path d="m9 18 6-6-6-6"/></svg></button></div>`;
  const reduce = window.matchMedia && matchMedia('(prefers-reduced-motion: reduce)').matches;
  let at = 0, timer = null, paused = false;
  const draw = (s) => {   // as on a protein page: frequency (bars in their most frequent cluster's color) over the fingerprint (contacts in navy, clusters in the strip), one residue axis
    const cv = $('#sc-cv'), W = cv.clientWidth, FH = 78, GAP = 6, PH = 112, H = FH + GAP + PH + 16, g = canvasCtx(cv, W, H), L = s.L, X0 = 12, bw = (W - X0) / L, x = (r) => X0 + (r - 1) * bw;
    const max = Math.max(1, ...s.tot);
    for (let r = 1; r <= L; r++) { const n = s.tot[r - 1]; if (!n) continue; const h = n / max * (FH - 4); g.fillStyle = clusterColor(s.dom[r - 1] || 1, s.k); g.fillRect(x(r), FH - h, Math.max(1, bw), h); }
    g.fillStyle = '#D5DDE6'; g.fillRect(X0, FH, W - X0, 1);
    const n = s.fp.length, y0 = FH + GAP, rh = PH / n;
    g.fillStyle = '#F7FBFF'; g.fillRect(X0, y0, W - X0, PH);
    s.fp.forEach(([c, ...res], i) => { const y = y0 + i * rh, hh = Math.max(1, rh);
      g.fillStyle = clusterColor(c, s.k); g.fillRect(0, y, X0 - 4, hh);
      g.fillStyle = '#08306B'; for (const r of res) g.fillRect(x(r), y, Math.max(1, bw), hh); });
    g.strokeStyle = '#D5DDE6'; g.lineWidth = 1; g.strokeRect(X0 + 0.5, y0 + 0.5, W - X0 - 1, PH - 1);
    g.fillStyle = '#627085'; g.font = '10.5px "IBM Plex Mono", monospace'; g.textBaseline = 'alphabetic';
    g.textAlign = 'left'; g.fillText('1', X0, H - 3); g.textAlign = 'right'; g.fillText(fmtInt(L), W, H - 3);
  };
  const drawNet = (s) => {   // a virus slide: its network at the cutoff, as on the virus page (gray edges by iLIS, nodes colored by community, black rings for homodimers)
    const cv = $('#sc-cv'), W = cv.clientWidth, H = 206, g = canvasCtx(cv, W, H);
    if (!s._pos) {   // communities as on the virus page (Louvain, iLIS-weighted; groups of three or more colored), each group laid out apart
      const cm = louvain(s.nodes.length, s.edges), csize = new Map(); cm.forEach((c) => csize.set(c, (csize.get(c) || 0) + 1));
      const big = [...csize].filter(([, z]) => z >= 3).map(([c]) => c).sort((a, b) => a - b), center = new Map();
      big.forEach((c, n) => { if (n === 0) center.set(c, [0, 0]); else { const t = (n - 1) / Math.max(1, big.length - 1) * 2 * Math.PI; center.set(c, [150 * Math.cos(t), 70 * Math.sin(t)]); } });
      const linked = new Set(s.edges.flatMap(([a, b]) => [a, b]));
      s._col = cm.map((c, i) => (!linked.has(i) ? '#C3CCD6' : big.includes(c) && big.indexOf(c) < TAB10.length ? TAB10[big.indexOf(c)] : '#1A5276'));
      const home = (d) => center.get(cm[d.i]) || [0, 0], nodes = s.nodes.map((l, i) => ({ i, x: (center.get(cm[i]) || [0, 0])[0] + 20 * Math.cos(i), y: (center.get(cm[i]) || [0, 0])[1] + 20 * Math.sin(i) }));
      const links = s.edges.map(([a, b, w]) => ({ source: a, target: b, w, same: cm[a] === cm[b] }));
      const sim = d3.forceSimulation(nodes).force('link', d3.forceLink(links).distance(26).strength((l) => (l.same ? 0.6 : 0.02))).force('charge', d3.forceManyBody().strength(-38))
        .force('x', d3.forceX((d) => home(d)[0]).strength((d) => (center.has(cm[d.i]) ? 0.3 : 0.05))).force('y', d3.forceY((d) => home(d)[1]).strength((d) => (center.has(cm[d.i]) ? 0.3 : 0.09))).stop();
      for (let t = 0; t < 320; t++) sim.tick();
      const xs = nodes.map((n) => n.x), ys = nodes.map((n) => n.y), x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys), y1 = Math.max(...ys);
      s._pos = nodes.map((n) => [(n.x - x0) / Math.max(1, x1 - x0), (n.y - y0) / Math.max(1, y1 - y0)]); }
    const P = s._pos.map(([u, v]) => [24 + u * (W - 48), 14 + v * (H - 28)]), deg = new Map(); s.edges.forEach(([a, b]) => { deg.set(a, (deg.get(a) || 0) + 1); deg.set(b, (deg.get(b) || 0) + 1); });
    for (const [a, b, w] of s.edges) { g.strokeStyle = EGRAY(w); g.lineWidth = EWID(w); g.beginPath(); g.moveTo(...P[a]); g.lineTo(...P[b]); g.stroke(); }
    const homo = new Set(s.homo);
    P.forEach(([x, y], i) => { g.beginPath(); g.arc(x, y, 4.5, 0, 2 * Math.PI); g.fillStyle = s._col[i]; g.fill(); g.lineWidth = homo.has(i) ? 1.8 : 1; g.strokeStyle = homo.has(i) ? HOMO_RING : '#fff'; g.stroke(); });
    g.font = '600 10px "IBM Plex Sans", system-ui, sans-serif'; g.textAlign = 'center'; g.textBaseline = 'alphabetic';
    P.forEach(([x, y], i) => { if ((deg.get(i) || 0) < 3) return; g.lineWidth = 3; g.strokeStyle = 'rgba(255,255,255,0.9)'; g.strokeText(s.nodes[i], x, y - 7); g.fillStyle = '#17263A'; g.fillText(s.nodes[i], x, y - 7); });
  };
  const show = (k, user) => {
    if (stale(gen) || !box.isConnected) { clearInterval(timer); return; }
    at = (k + S.length) % S.length; const s = S[at], href = `#/${s.sp}/${encodeURIComponent(s.key)}`;
    const paint = () => {
      if (s.kind === 'virus') {   // a virus: its predicted network
        const vh = `#/virus/taxon/${s.taxid}`;
        $('#sc-gene').textContent = s.gene; $('#sc-gene').href = vh; $('#sc-sp').textContent = s.spLabel;
        $('#sc-note').innerHTML = `virus network · ${fmtInt(s.proteins)} proteins · ${fmtInt(s.edges.length)} pairs and ${fmtInt(s.homo.length)} homodimers past 5% FPR (iLIS ≥ ${s.cut})`;
        $('#sc-cap').innerHTML = 'every protein pair of the virus folded and scored: an edge per pair past the cutoff, darker for higher iLIS; ringed proteins bind themselves';
        $('#sc-open').textContent = `Open ${s.gene} →`; $('#sc-open').href = vh;
        drawNet(s); box.querySelectorAll('.sc-dots button').forEach((b, i) => b.setAttribute('aria-current', i === at ? 'true' : 'false'));
        box.classList.remove('sc-out'); return;
      }
      $('#sc-gene').textContent = s.gene; $('#sc-gene').href = href; $('#sc-sp').textContent = s.spLabel;
      $('#sc-note').innerHTML = `<span class="q">c</span>lustered <span class="q">L</span>ocal <span class="q">I</span>nteraction <span class="q">P</span>rofiler (<span class="q">cLIP</span>) · ${s.k} clusters · ${fmtInt(s.partners)} partners past 10% FPR`;
      $('#sc-cap').innerHTML = '<span class="q">cLIP</span> groups partners by the residues they contact: '   // the full name is on the line above
        + ['contacts per residue above,', 'one row per partner below'].map((t) => t.replace(/ /g, '&nbsp;')).join(' ');   // each legend phrase stays on one line
      $('#sc-open').textContent = `Open ${s.gene} →`; $('#sc-open').href = href;
      draw(s); box.querySelectorAll('.sc-dots button').forEach((b, i) => b.setAttribute('aria-current', i === at ? 'true' : 'false'));
      box.classList.remove('sc-out');
    };
    if (reduce || !user && at === 0 && !box.dataset.shown) { paint(); box.dataset.shown = '1'; return; }
    box.classList.add('sc-out'); setTimeout(paint, 180);
  };
  const run = () => { clearInterval(timer); if (!reduce) timer = setInterval(() => { if (!paused) show(at + 1); }, 7000); };
  box.querySelectorAll('.sc-dots button').forEach((b) => b.onclick = () => { show(+b.dataset.k, true); run(); });
  const step = (d) => { show(at + d, true); run(); };   // arrows, the arrow keys and a sideways swipe move one protein
  $('#sc-prev').onclick = () => step(-1); $('#sc-next').onclick = () => step(1);
  box.addEventListener('keydown', (e) => { if (e.key === 'ArrowLeft' || e.key === 'ArrowRight') { e.preventDefault(); step(e.key === 'ArrowLeft' ? -1 : 1); } });
  let t0 = null;
  box.addEventListener('touchstart', (e) => { const t = e.touches[0]; t0 = [t.clientX, t.clientY]; }, { passive: true });
  box.addEventListener('touchend', (e) => { if (!t0) return; const t = e.changedTouches[0], dx = t.clientX - t0[0], dy = t.clientY - t0[1]; t0 = null;
    if (Math.abs(dx) > 40 && Math.abs(dx) > 1.5 * Math.abs(dy)) step(dx < 0 ? 1 : -1); }, { passive: true });
  box.addEventListener('mouseenter', () => { paused = true; }); box.addEventListener('mouseleave', () => { paused = false; });
  box.addEventListener('focusin', () => { paused = true; }); box.addEventListener('focusout', () => { paused = false; });
  window.onresize = () => { if (box.isConnected) (S[at].kind === 'virus' ? drawNet : draw)(S[at]); };
  show(0); run();
}

async function viewDatasets() {
  const gen = ROUTE, reg = await registry();
  if (stale(gen)) return;
  app.innerHTML = `<div class="crumbs"><a href="#/">Atlas</a> / <a href="#/datasets">Datasets</a></div>
    ${(reg.themes || []).length ? '<h2 class="section-h" style="margin-top:4px">Themes</h2><div class="datasets live-row" id="themes"></div>' : ''}
    <h2 class="section-h"${(reg.themes || []).length ? '' : ' style="margin-top:4px"'}>Datasets</h2>
    <div class="datasets live-row" id="ds-cards">${reg.datasets.filter((d) => d.status !== 'planned').map(dsCard).join('')}</div>`;
  fillDsStats(); fillThemes(); setCards(gen, reg);
}
// A screen kept as a set inside a larger dataset (the fly kinase–kinase and kinase–TF screens in FlyPredictome) gets its
// own card, placed after the other screens of its kind; it cites what they cite.
async function setCards(gen, reg) {
  for (const d of reg.datasets.filter((x) => x.status === 'live')) {
    let TS; try { TS = await setsOf(await dataset(d.id)); } catch (e) { continue; }
    if (stale(gen) || !TS) continue;
    const spx = (reg.species || []).find((s) => s.id === d.species);
    for (const s of TS.list.filter((x) => x.type === 'screen')) {
      const box = $('#ds-cards'); if (!box || box.querySelector(`[data-set="${d.id}/${s.id}"]`)) continue;
      const kin = reg.datasets.filter((x) => x.short === s.short && x.id !== d.id), cite = kin[0], k = s.counts || {};
      const src = cite ? (cite.paper ? `<a href="${esc(cite.paper)}" target="_blank" rel="noopener">${esc(cite.source)} ↗</a>` : esc(cite.source)) : '';
      const card = `<div class="ds live" data-set="${d.id}/${s.id}"><span class="badge on">Searchable</span>
        <h3><a href="#/datasets/${d.id}/${s.id}">${esc(spx ? spx.label : d.species)} ${esc(s.title.charAt(0).toLowerCase() + s.title.slice(1))}</a></h3>
        <div class="sp">${esc(spx ? spx.name : d.species)} · <span class="src" style="--c:${s.color}">${esc(s.short)}</span></div>${d.models ? `<div class="sp">${runSettings(d)}</div>` : ''}
        <div class="stats"><div><b>${fmtInt(k.proteins)}</b><span>proteins</span></div><div><b>${fmtInt(k.pairs)}</b><span>pairs</span></div><div><b>${fmtInt(k.pairsFpr10)}</b><span>past 10% FPR${cutNote(10)}</span></div></div>
        <div class="src-line">${src}${src ? ' · ' : ''}a set within <a href="#/datasets/${d.id}">${esc(d.title)}</a></div></div>`;
      const after = kin.map((x) => box.querySelector(`[data-ds="${x.id}"]`)).filter(Boolean)[0];
      if (after) after.insertAdjacentHTML('afterend', card); else box.insertAdjacentHTML('beforeend', card);
    }
  }
}

async function viewAbout() {
  const gen = ROUTE, reg = await registry(); if (stale(gen)) return;
  const ref = (k, text) => `<li>${text} <a href="https://doi.org/${REF[k][1]}" target="_blank" rel="noopener">doi.org/${REF[k][1]}</a></li>`;
  const AFM = '<a href="https://github.com/flyark/AFM-LIS" target="_blank" rel="noopener">AFM-LIS</a>';
  // each live screen's source, screens that share one grouped on one line; a screen not listed here shows its registry source
  const SRC = { 'human-predictomes': ['schmid2025', 'Schmid, E. W. et al. (2025). Proteome-wide in silico screening for human protein-protein interactions. <i>bioRxiv</i>.'],
    'human-kinase-tf': ['kim2025', 'Kim, A.-R. et al. (2025). A structure-guided kinase–transcription factor interactome atlas reveals docking landscapes of the kinome. <i>bioRxiv</i>.'],
    flypredictome: ['flypredictome', 'Kim, A.-R. et al. (2026). FlyPredictome: a structural atlas of predicted protein-protein interactions in <i>Drosophila</i>. <i>bioRxiv</i>.'] };
  const kk = 'Kim, A.-R. &amp; Perrimon, N. (2026). LIVIA: a browser-based tool for assessing and visualizing predicted protein interactions. <i>bioRxiv</i>. The predictions are in the LIVIA Atlas record on Zenodo, cited above.';
  const groups = new Map();
  for (const d of (reg.datasets || []).filter((x) => x.status === 'live')) {
    const [k, text] = SRC[d.id] || (/kinase-kinase$/.test(d.id) ? ['livia', kk] : [null, `${esc(d.source || '')}${d.paper ? ` <a href="${esc(d.paper)}" target="_blank" rel="noopener">${esc(d.paper.replace(/^https?:\/\//, ''))}</a>` : ''}`]);
    const spx = (reg.species || []).find((x) => x.id === d.species), g = groups.get(k || d.id) || { k, text, short: d.short, titles: [], sps: [] };
    g.titles.push(d.title); g.sps.push(spx ? (/^[A-Z]\. /.test(spx.label) ? spx.label : spx.label.toLowerCase()) : d.species); groups.set(k || d.id, g); }
  const andJoin = (xs) => (xs.length < 3 ? xs.join(' and ') : `${xs.slice(0, -1).join(', ')} and ${xs[xs.length - 1]}`);
  const screenList = [...groups.values()].map((g) => { const t = `<b>${esc(g.titles.length > 1 ? `${g.short} screens (${andJoin(g.sps)})` : g.titles[0])}</b>: ${g.text}`;
    return g.k ? ref(g.k, t) : `<li>${t}</li>`; }).join('');
  app.innerHTML = `<div class="reading about"><div class="crumbs"><a href="#/">Atlas</a> / <a href="#/about">About</a></div>
    <div class="card" style="margin-top:6px"><h2>What this is</h2>
      <p>LIVIA Atlas makes large AlphaFold-Multimer interaction screens searchable at the level of residues. Every prediction is scored with
      <b>iLIS</b>, the integrated local interaction score, computed by lis.py (${AFM} on GitHub) over residue pairs with predicted aligned error of at most 12 Å
      (LIS), and over those that are also in contact, Cβ–Cβ distance of at most 8 Å (cLIS): iLIS = √(LIS × cLIS). <b>iLISA</b> weighs iLIS by
      the size of the interface: iLISA = iLIS × iLIA, where iLIA = √(LIA × cLIA) and LIA and cLIA count the residue pairs that enter LIS and cLIS.</p>
      <table class="cuts"><caption>Cutoffs at a 10%, 5% and 1% false-positive rate, from binary reference sets of the kind used to benchmark proteome-scale Y2H screens (positive reference sets of literature-curated interactions; random reference sets as negatives), in yeast, fly and human (357 positive and 1,255 negative pairs; ${cite('flypredictome')}), for the top-ranked model of a pair and for the average of its five models. The Atlas applies the top-ranked cutoffs to a pair's best model; on the same set that gives 10.4%, 5.3% and 1.1%. ipTM as reported in the run logs (two decimals). Cutoff tables: <a href="https://github.com/flyark/AFM-LIS/blob/main/thresholds_data_yfh_lipdockq.xlsx" target="_blank" rel="noopener">AFM-LIS repository ↗</a></caption>
        <thead><tr><th>Metric</th><th>Model</th><th>10% FPR</th><th>5% FPR</th><th>1% FPR</th></tr></thead>
        <tbody>${[['iLIS', 3], ['iLISA', 1], ['LIS', 3], ['cLIS', 3], ['ipTM', 2], ['ipSAE', 3], ['actifpTM', 3], ['pDockQ', 3], ['LIpDockQ', 3], ['pDockQ2', 3], ['LIpDockQ2', 3]]
          .map(([m, d]) => [['top-ranked', FPR[m]], ['average of 5', FPR_AVG[m]]].map(([l, c], k) => `<tr${k ? ' class="avg"' : ''}>${k ? '' : `<th rowspan="2">${m}</th>`}<td class="mdl">${l}</td>${c ? c.map((v, j) => `<td style="color:${BAND_TXT[[10, 5, 1][j]]}">≥ ${v.toFixed(d)}</td>`).join('') : '<td>—</td><td>—</td><td>—</td>'}</tr>`).join('')).join('')}</tbody></table>
      <p class="muted">${CALIB} Where two human screens both called a pair, they placed the interface on nearly the same residues; which pairs they call agrees less often. Each screen's models and recycles are listed on its <a href="#/datasets">dataset page</a>. The screens in the Atlas were scored with lis.py before 26 September 2026, which left out values exactly at a cutoff (a PAE of 12 Å, a Cβ distance of 8 Å); LIVIA includes them, and iLIS differs by at most about 0.003.</p>
      <p>Each protein has one page per species that gathers its predictions from every screen. A pair predicted in two screens, or both ways round,
      keeps every model with its source. The interaction residues on both proteins are kept for every prediction, and each protein page runs
      <a href="${LIVIA}clip.html" target="_blank" rel="noopener">LIVIA cLIP</a> in the browser: partners are clustered by their interaction fingerprints, and the clusters
      are mapped onto the AlphaFold DB structure with pLDDT and, for human proteins, AlphaMissense.</p>
      <p>Fly pages are organized by FlyBase gene. FlyPredictome folded many constructs of a gene: isoforms, fragments, phosphosite windows,
      point mutants. Every name is resolved to its gene, and each construct is placed on the gene's reference sequence by its folded sequence, so its
      contacts are drawn in the gene's residue numbering. An isoform too unlike the reference to be placed is shown on its own: its page has an
      Isoform row, and every card follows the isoform chosen. The table of construct names and their genes is on the FlyPredictome page.</p>
      <p>Interactions reported in BioGRID (release 5.0.261, MIT license; ${cite('biogrid')}) are marked, matched to each species' proteins by UniProt
      accession, official symbol or systematic name. In networks, edges are shaded in gray by the pair's best iLIS, and the pairs BioGRID reports
      (physical, genetic or either, as the reader chooses) can be colored by what was reported (physical red, genetic green, both purple); the width is the average iLIS. On a protein page, a partner
      with a reported physical interaction is ringed in the Overview, and in the Clusters and Partners lists a reported partner's name is marked in the same colors, light.
      The matched pairs are a file on this site for each species, so nothing is looked up elsewhere.</p></div>
    <div class="card"><h2>Cite</h2>
      <ul class="refs">
        ${ref('livia', 'LIVIA: Kim, A.-R. &amp; Perrimon, N. (2026). LIVIA: a browser-based tool for assessing and visualizing predicted protein interactions. <i>bioRxiv</i>.')}
        <li>Data archive: ${archiveLine()}</li>
        ${ref('flypredictome', 'iLIS and its cutoffs: Kim, A.-R. et al. (2026). FlyPredictome: a structural atlas of predicted protein-protein interactions in <i>Drosophila</i>. <i>bioRxiv</i>.')}
        ${ref('afmlis', 'LIS and AFM-LIS: Kim, A.-R. et al. (2024). Enhanced protein-protein interaction discovery via AlphaFold-Multimer. <i>bioRxiv</i>.')}
      </ul>
      <p class="muted" style="font-size:14px;margin-bottom:0">Please also cite the source of each screen you use, and the resources below that your work draws on.</p></div>
    <div class="card"><h2>Screens</h2>
      <ul class="refs">${screenList}</ul></div>
    <div class="card"><h2>Data and software</h2>
      <h3 class="refs-h">Structures and annotations</h3>
      <ul class="refs">
        ${ref('uniprot', 'Protein names, sequences and domains (Domain, DNA-binding and zinc-finger features): The UniProt Consortium (2025). UniProt: the Universal Protein Knowledgebase in 2025. <i>Nucleic Acids Res.</i> 53, D609–D617.')}
        ${ref('flybase', 'Fly genes, symbols and synonyms (release FB2026_03): Öztürk-Çolak, A. et al. (2024). FlyBase: updates to the <i>Drosophila</i> genes and genomes database. <i>Genetics</i> 227, iyad211.')}
        ${ref('afdb', 'Single-chain models in the 3D view: Varadi, M. et al. (2024). AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. <i>Nucleic Acids Res.</i> 52, D368–D375.')}
        ${ref('alphafold', 'Those models and their pLDDT: Jumper, J. et al. (2021). Highly accurate protein structure prediction with AlphaFold. <i>Nature</i> 596, 583–589.')}
        ${ref('alphamissense', 'Missense pathogenicity on human proteins: Cheng, J. et al. (2023). Accurate proteome-wide missense variant effect prediction with AlphaMissense. <i>Science</i> 381, eadg7492.')}
        ${ref('biogrid', 'Reported physical and genetic interactions (release 5.0.261): Oughtred, R. et al. (2021). The BioGRID database: a comprehensive biomedical resource of curated protein, genetic, and chemical interactions. <i>Protein Sci.</i> 30, 187–200.')}
      </ul>
      <h3 class="refs-h">Prediction, scoring and display</h3>
      <ul class="refs">
        ${ref('afm', 'Every screen was folded with AlphaFold-Multimer: Evans, R. et al. (2021). Protein complex prediction with AlphaFold-Multimer. <i>bioRxiv</i>.')}
        ${ref('colabfold', 'The kinase–kinase screens were run through ColabFold: Mirdita, M. et al. (2022). ColabFold: making protein folding accessible to all. <i>Nat. Methods</i> 19, 679–682.')}
        <li>Scores: lis.py from ${AFM} (LIS, cLIS, iLIS; cited above), on every model.</li>
        ${ref('silhouette', 'Binding sites: LIVIA cLIP, run in the browser, picks the number of clusters by silhouette: Rousseeuw, P. J. (1987). Silhouettes: a graphical aid to the interpretation and validation of cluster analysis. <i>J. Comput. Appl. Math.</i> 20, 53–65.')}
        ${ref('molstar', '3D view, in LIVIA\'s Mol* page: Sehnal, D. et al. (2021). Mol* Viewer: modern web app for 3D visualization and analysis of large biomolecular structures. <i>Nucleic Acids Res.</i> 49, W431–W437.')}
        ${ref('d3', 'Plots and networks, D3 7.9: Bostock, M., Ogievetsky, V. &amp; Heer, J. (2011). D³ data-driven documents. <i>IEEE Trans. Vis. Comput. Graph.</i> 17, 2301–2309.')}
        ${ref('leiden', 'Communities in LIVIA\'s network page: Traag, V. A., Waltman, L. &amp; van Eck, N. J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. <i>Sci. Rep.</i> 9, 5233.')}
        <li>Reading prediction bundles in the browser: JSZip 3.10. <a href="https://stuk.github.io/jszip/" target="_blank" rel="noopener">stuk.github.io/jszip</a></li>
        ${ref('zenodo', 'Hosting: the site and some screens on GitHub Pages; the others in the LIVIA Atlas record on Zenodo, read by byte range, so a page loads only its own bundles. European Organization for Nuclear Research &amp; OpenAIRE (2013). <i>Zenodo</i>. CERN.')}
        <li>Visits are counted with GoatCounter, without cookies.</li>
      </ul></div></div>`;
}

// The counts of a screen, a species or a set, in one row: a prediction is a pair folded once (its ranked models are
// counted as models); protein pairs are unique pairs. One model per pair (the same three numbers): one tile, so marked
const kpiRow = (k) => { const one = k.predictions === k.pairs && (!k.runs || k.runs === k.pairs);
  return `<div class="kpirow"><div class="kpi"><b>${fmtInt(k.proteins)}</b><span>proteins</span></div><div class="kpi"><b>${fmtInt(k.pairs)}</b><span>protein pairs${one ? '<i class="kc t">one model each</i>' : ''}</span></div>
  ${one ? '' : `${k.runs ? `<div class="kpi"><b>${fmtInt(k.runs)}</b><span>predictions</span></div>` : ''}<div class="kpi"><b>${fmtInt(k.predictions)}</b><span>models</span></div>`}
  <div class="kpi f10"><b>${fmtInt(k.pairsFpr10)}</b><span>pairs past 10% FPR${cutNote(10)}</span></div><div class="kpi f5"><b>${fmtInt(k.pairsFpr5)}</b><span>past 5% FPR${cutNote(5)}</span></div>
  <div class="kpi f1"><b>${fmtInt(k.pairsFpr1)}</b><span>past 1% FPR${cutNote(1)}</span></div></div>`; };
const shortCite = (src) => `${src.citation.split(' ')[0]} et al. ${(src.citation.match(/\((\d{4})\)/) || [])[1] || ''}`.trim();
// Thematic sets of a dataset, for its dataset and species pages: one table, screens first, then source categories;
// the bar is each set's share of the dataset's predictions (sets overlap, so shares need not add up)
function setsCard(TS) {
  if (!TS) return '';
  const total = TS.ds.manifest.counts.predictions || Math.max(...TS.list.map((x) => x.counts.predictions));
  const row = (S) => `<tr><td class="set-name"><i style="background:${S.color}"></i><a href="#/datasets/${TS.ds.id}/${S.id}">${esc(S.title)}</a>${S.source
    ? `<a class="set-cite" href="${esc(S.source.url)}" target="_blank" rel="noopener">${esc(shortCite(S.source))} ↗</a>` : ''}</td>
    <td class="n">${fmtInt(S.counts.proteins)}</td><td class="n">${fmtInt(S.counts.pairs)}</td><td class="n">${fmtInt(S.counts.predictions)}</td>
    <td class="set-bar" title="${(100 * S.counts.predictions / total).toFixed(1)}% of all models"><div><span style="width:${Math.max(1, 100 * S.counts.predictions / total).toFixed(1)}%;background:${S.color}"></span></div></td></tr>`;
  const grp = (label, list) => (list.length ? `<tr class="set-grp"><th colspan="5">${label}</th></tr>${list.map(row).join('')}` : '');
  const scr = TS.list.filter((x) => x.type === 'screen'), cat = TS.list.filter((x) => x.type !== 'screen');
  return `<div class="card"><div class="card-head"><h2>Thematic sets</h2><span class="muted">subsets of ${esc(TS.ds.reg.short)} · open one, or choose it at the top of a protein page</span></div>
    <div class="tbl-wrap"><table class="sets"><thead><tr><th>Set</th><th class="n">Proteins</th><th class="n">Pairs</th><th class="n">Models</th><th>Share of models</th></tr></thead>
      <tbody>${grp('Screens', scr)}${grp('Source categories', cat)}</tbody></table></div></div>`;
}
/* themes: one biological question across species — each member a screen of a species, or a thematic set of a screen */
async function themeMembers(T) {
  return (await Promise.all(T.members.map(async (m) => {
    const S = await regSpecies(m.species); let ds; try { ds = await dataset(m.dataset); } catch (e) { return null; }
    if (m.set) { const TS = await setsOf(ds).catch(() => null), st = TS && TS.byId.get(m.set); if (!st) return null;
      return { S, title: st.title, within: ds.reg.title, counts: st.counts, href: `#/datasets/${ds.id}/${st.id}`, source: st.source || null }; }
    return { S, title: ds.reg.title, within: '', counts: ds.manifest.counts, href: `#/datasets/${ds.id}`, source: ds.manifest.source || null };
  }))).filter(Boolean);
}
async function viewTheme(id) {
  const gen = ROUTE, reg = await registry(), T = (reg.themes || []).find((t) => t.id === id);
  if (stale(gen)) return;
  if (!T) { app.innerHTML = `<div class="empty">No theme “${esc(id)}”. <a href="#/">Go to the atlas home</a></div>`; return; }
  document.title = `${T.title} · LIVIA Atlas`;
  const rows = await themeMembers(T), sum = (k) => rows.reduce((a, r) => a + (r.counts[k] || 0), 0);
  if (stale(gen)) return;
  const cite = (src) => (src && src.url ? `<a href="${esc(src.url)}" target="_blank" rel="noopener">${esc(shortCite(src))} ↗</a>` : '');
  app.innerHTML = `<div class="crumbs"><a href="#/">Atlas</a> / <a href="#/datasets">Themes</a> / <a href="#/themes/${T.id}">${esc(T.title)}</a></div>
    <div class="dshead"><h1>${esc(T.title)}</h1><div class="pname">${esc(T.about || '')}</div></div>
    ${kpiRow({ proteins: sum('proteins'), pairs: sum('pairs'), runs: sum('runs'), predictions: sum('predictions'), pairsFpr10: sum('pairsFpr10'), pairsFpr5: sum('pairsFpr5'), pairsFpr1: sum('pairsFpr1') })}
    <div class="card"><div class="card-head"><h2>By species</h2><span class="muted">open one to search it; its protein pages show only this theme</span></div>
      <div class="tbl-wrap"><table class="sets"><thead><tr><th>Species</th><th class="n">Proteins</th><th class="n">Protein pairs</th><th class="n">Predictions</th><th class="n">Past 10% FPR</th></tr></thead><tbody>
      ${rows.map((r) => `<tr><td class="set-name"><i style="background:${T.color}"></i><a href="${r.href}">${esc(r.S.label)}</a> <span class="muted" style="font-style:italic">${esc(r.S.name)}</span>${r.within ? ` <span class="muted">· in ${esc(r.within)}</span>` : ''}${r.source ? `<span class="set-cite">${cite(r.source)}</span>` : ''}</td>
        <td class="n">${fmtInt(r.counts.proteins)}</td><td class="n">${fmtInt(r.counts.pairs)}</td><td class="n">${fmtInt(r.counts.runs || r.counts.predictions)}</td><td class="n">${fmtInt(r.counts.pairsFpr10)}</td></tr>`).join('')}
      </tbody></table></div></div>`;
}
/* a thematic set: a view over one dataset (its prediction runs), with its own counts, citation, hubs and search */
async function viewSet(dsId, setId) {
  const gen = ROUTE, ds = await dataset(dsId), TS = await setsOf(ds), S = TS && TS.byId.get(setId), sp = await species(ds.reg.species), m = ds.manifest;
  if (stale(gen)) return;
  if (!S) { app.innerHTML = `<div class="empty">No set “${esc(setId)}” in ${esc(ds.reg.title)}.</div>`; return; }
  document.title = `${S.short} · ${ds.reg.short} · LIVIA Atlas`;
  const prot = await getJSON(new URL(S.files.proteins, ds.base).href), c = Object.fromEntries(prot.columns.map((k, i) => [k, i]));
  if (stale(gen)) return;
  const rows = prot.rows.map((r) => ({ key: r[c.key], pos10: r[c.pos10] })), keys = new Set(rows.map((r) => r.key));
  const hubs = [...rows].sort((a, b) => b.pos10 - a.pos10).slice(0, 24), k = S.counts;
  const link = (u, t) => `<a href="${esc(u)}" target="_blank" rel="noopener">${esc(t)} ↗</a>`;
  app.innerHTML = `<div class="crumbs"><a href="#/">Atlas</a> / <a href="#/datasets">Datasets</a> / <a href="#/datasets/${ds.id}">${esc(ds.reg.title)}</a> / <a href="#/datasets/${ds.id}/${S.id}">${esc(S.short)}</a></div>
    <div class="dshead"><h1>${esc(S.title)}</h1>
      <div class="pname"><span class="src" style="--c:${S.color}">${esc(S.short)}</span> ${S.type === 'screen' ? 'A screen' : 'A source category'} within ${esc(ds.reg.title)} · <i>${esc(m.species.name)}</i></div>${ds.reg.models ? `<div class="pname">${runSettings(ds.reg)}</div>` : ''}
      ${S.source ? `<div class="cite">${link(S.source.url, `${S.source.citation} doi:${S.source.doi}`)}</div>` : ''}
      <div class="cite">Part of ${m.source.url ? link(m.source.url, m.source.citation) : esc(m.source.citation)}</div></div>
    ${kpiRow(k)}
    <div class="card"><h2>Search this set</h2><p class="muted" style="margin:2px 0 10px">Protein pages opened from here show only this set's predictions, with a switch to all of ${esc(ds.reg.short)}.</p><div id="set-search"></div></div>
    <div class="card"><h2>Most connected proteins in this set <span class="muted">partners past the 10% FPR cutoff</span></h2>
      <div class="chips">${hubs.map((r) => { const R = sp.byKey.get(r.key); return `<a class="chip" href="#/${sp.id}/${r.key}?set=${S.id}">${esc(R ? R.gene : r.key)} <span class="num" style="color:var(--ink-3)">${fmtInt(r.pos10)}</span></a>`; }).join('')}</div></div>
    <div class="card"><h2>Data</h2><p class="muted" style="font-size:14px;margin:4px 0 0">This set is part of the ${esc(ds.reg.title)} download in the
      <a href="${ARCHIVE.url}" target="_blank" rel="noopener">LIVIA Atlas record on Zenodo ↗</a>;
      <a href="https://github.com/flyark/livia-atlas/blob/main/tools/extract_set.py" target="_blank" rel="noopener">extract_set.py ↗</a> pulls out just this set as a table.</p></div>`;
  mountSearch($('#set-search'), { spId: sp.id, only: keys, set: S.id });
}
async function viewDataset(dsId) {   // one screen: what it is, its counts and files; proteins link to their species pages
  const gen = ROUTE, ds = await dataset(dsId), m = ds.manifest, k = m.counts, sp = await species(ds.reg.species), TS = await setsOf(ds);
  const rows = await datasetRows(ds), hubs = [...rows].sort((a, b) => b.pos10 - a.pos10).slice(0, 24);
  if (stale(gen)) return;
  const scopeQ = sp.dsIds.length > 1 ? `?set=${encodeURIComponent(ds.id)}` : '';   // one of several screens: its protein pages open in its scope
  app.innerHTML = `<div class="crumbs"><a href="#/">Atlas</a> / <a href="#/datasets">Datasets</a> / <a href="#/datasets/${ds.id}">${esc(ds.reg.title)}</a></div>
    <div class="dshead"><h1>${esc(ds.reg.title)}</h1><div class="pname"><i>${esc(m.species.name)}</i> · ${esc(m.source.method)} · ${esc(m.analysis.tool)}, <span title="${m.analysis.inclusive ? 'cutoffs as LIVIA applies them' : 'as scored: lis.py before 26 Sep 2026 left out values exactly at a cutoff; LIVIA includes them, and iLIS differs by at most about 0.003'}">${m.analysis.inclusive ? `PAE ≤ ${m.analysis.paeCutoff} Å, Cβ ≤ ${m.analysis.cbCutoff} Å` : `PAE ${m.analysis.paeCutoff} Å, Cβ ${m.analysis.cbCutoff} Å, values exactly at a cutoff excluded`}</span></div>${ds.reg.models ? `<div class="pname">${runSettings(ds.reg)}</div>` : ''}
      <div class="cite">${m.source.url ? `<a href="${esc(m.source.url)}" target="_blank" rel="noopener">${esc(m.source.citation)}${m.source.doi ? ` doi:${esc(m.source.doi)}` : ''} ↗</a>` : esc(m.source.citation)}</div></div>
    ${kpiRow(k)}
    <div class="card"><h2>Search</h2><p class="muted" style="margin:2px 0 10px">${sp.dsIds.length > 1 ? `Protein pages opened from here show only this screen, <span class="src" style="--c:${ds.reg.color}">${esc(ds.reg.short)}</span>, with a switch to every ${esc(sp.reg.label.toLowerCase())} screen.`
      : `Protein pages show every prediction of this screen${TS ? ', with a switch to each of its thematic sets' : ''}.`}</p><div id="ds-search"></div></div>
    <div class="card"><h2>Most connected proteins in this screen <span class="muted">partners past the 10% FPR cutoff</span></h2>
      <div class="chips">${hubs.map((r) => { const R = sp.byName.get(r.id); return `<a class="chip" href="#/${sp.id}/${R ? R.key : r.id}${scopeQ}">${esc(r.gene)} <span class="num" style="color:var(--ink-3)">${fmtInt(r.pos10)}</span></a>`; }).join('')}</div></div>
    ${setsCard(TS)}
    <div class="card"><h2>Data</h2><p class="muted" style="font-size:14px;margin:4px 0 0">Every file of this screen is in the
      <a href="${ARCHIVE.url}" target="_blank" rel="noopener">LIVIA Atlas record on Zenodo (doi:${ARCHIVE.doi}) ↗</a>. A protein page's
      <b>Data</b> menu downloads that protein's predictions.${m.files.identity ? ` <a href="${ds.base}${m.files.identity}" download>Construct names and their FlyBase genes</a> (table).` : ''}</p></div>`;
  mountSearch($('#ds-search'), sp.dsIds.length > 1 ? { spId: sp.id, set: ds.id, only: new Set(rows.map((r) => { const R = sp.byName.get(r.id); return R ? R.key : r.id; })) } : { spId: sp.id });
}

/* ── protein page: LIVIA cLIP, natively, over every screen, with a partner overview, a network and a partner table ── */
let CLIPW = null, clipSeq = 0; const clipWait = new Map();
function runClip(rows, gene, cut) {
  if (!CLIPW) { CLIPW = new Worker('clipworker.js?v=20260925p'); CLIPW.onmessage = (e) => { const w = clipWait.get(e.data.id); if (w) { clipWait.delete(e.data.id); e.data.ok ? w.resolve(e.data) : w.reject(new Error(e.data.message)); } }; }
  const id = ++clipSeq;
  return new Promise((resolve, reject) => { clipWait.set(id, { resolve, reject }); CLIPW.postMessage({ id, livia: LIVIA, rows: rows.filter((r) => +r.iLIS >= cut), gene, cut }); });
}
function stopClip() {   // leaving a page mid-clustering: drop its job, so the next page does not wait behind it
  if (!CLIPW || !clipWait.size) return;
  CLIPW.terminate(); CLIPW = null;
  for (const w of clipWait.values()) w.reject(new Error('The page changed.'));
  clipWait.clear();
}
const AXL = 64, AXR = 18;   // shared residue axis of the frequency plot and the fingerprint: dendrogram + cluster strip / y axis live in AXL
const METRICS = { iLIS: 'iLIS', iLISA: 'iLISA', iLIA: 'iLIA', ipTM: 'ipTM', pTM: 'pTM', LIS: 'LIS', cLIS: 'cLIS', LIA: 'LIA', cLIA: 'cLIA', ipSAE: 'ipSAE', actifpTM: 'actifpTM', qPl: 'pLDDT (query)', pPl: 'pLDDT (partner)', _rank: 'global rank' };
// Edge colors of every network, two independent layers the reader switches on or off: the pair's best iLIS on a color
// scale the reader picks (light at the 10% FPR cutoff, dark at 0.85 and above; or one flat gray) and, when chosen, the
// pairs reported in BioGRID colored by what was reported: physical red, genetic green, both purple. Edge width is the
// average iLIS.
const ESCALE = { gray: ['#C5CCD4', '#1E2A38'], blue: ['#C6DBEF', '#08306B'], brown: ['#E8D9C4', '#5B3A1A'] };
const ESCALE_LBL = { gray: 'gray', blue: 'blue', brown: 'brown', flat: 'off (one gray)' };
const EFLAT = '#9AA5B1', KB_COL = { p: '#C62828', g: '#2E7D32', pg: '#6A1B9A' };
const escale = (k) => d3.scaleLinear().domain([CUT[10], 0.85]).range(ESCALE[k] || ESCALE.gray).clamp(true);
const EGRAY = escale('gray');   // the fixed gray scale of the virus networks
const kbPubs = (ph, ge) => [ph ? `physical, ${ph} publication${ph === 1 ? '' : 's'}` : '', ge ? `genetic, ${ge} publication${ge === 1 ? '' : 's'}` : ''].filter(Boolean).join('; ');
const kbHit = (d, ev) => (ev === 'p' ? d.pubs > 0 : ev === 'g' ? d.gen > 0 : d.pubs > 0 || d.gen > 0);
const kbCol = (d) => (d.pubs > 0 && d.gen > 0 ? KB_COL.pg : d.pubs > 0 ? KB_COL.p : KB_COL.g);   // what was reported for the pair
const KB_EV = { pg: 'physical or genetic', p: 'physical', g: 'genetic' };
const edgeCtl = (id, kbOn = false) => `<label class="ctl" title="shade each edge by the best iLIS of the pair: light at the 10% FPR cutoff, dark at 0.85 and above">iLIS scale<select id="${id}-shade" aria-label="iLIS color scale">${Object.entries(ESCALE_LBL).map(([k, l]) => `<option value="${k}">${l}</option>`).join('')}</select></label>
  <label class="ctl" title="color the pairs reported in BioGRID: physical red, genetic green, both purple"><input type="checkbox" id="${id}-kb"${kbOn ? ' checked' : ''}> BioGRID</label><select id="${id}-ev" aria-label="Which BioGRID evidence"${kbOn ? '' : ' disabled'}>${Object.entries(KB_EV).map(([k, l]) => `<option value="${k}">${l}</option>`).join('')}</select>`;
function edgeStyle(id, K) {   // the reader's choice for one network → the color of an edge and whether it is a reported one
  const sc = $(`#${id}-shade`).value, shade = sc !== 'flat', kb = !!K && $(`#${id}-kb`).checked, ev = $(`#${id}-ev`).value, hit = (d) => kb && kbHit(d, ev);
  const ramp = escale(sc);
  return { shade, sc, kb, ev, hit, color: (d) => (hit(d) ? kbCol(d) : shade ? ramp(d.best) : EFLAT) };
}
function edgeKey(st, K, links, extra = []) {   // the key under a network: the iLIS scale (or flat gray) and, when on, the reported pairs by type
  const [lo, hi] = ESCALE[st.sc] || ESCALE.gray;
  const base = st.shade ? `<span><i class="kb-grad" style="background:linear-gradient(90deg, ${lo}, ${hi})"></i>best iLIS, 0.223 to 0.85+</span>` : `<span><i style="background:${EFLAT}"></i>predicted pair</span>`;
  const on = links.filter((d) => kbHit(d, st.ev)), n = (f) => fmtInt(on.filter(f).length);
  const red = !K ? '<span class="muted">no BioGRID records for this species</span>'
    : st.kb ? `<span class="muted">reported in BioGRID ${esc(K.release)} (${fmtInt(on.length)} of ${fmtInt(links.length)} pairs):</span>`
      + (st.ev !== 'g' ? `<span><i style="background:${KB_COL.p}"></i>physical (${n((d) => d.pubs > 0 && !(d.gen > 0))})</span>` : '')
      + (st.ev !== 'p' ? `<span><i style="background:${KB_COL.g}"></i>genetic (${n((d) => d.gen > 0 && !(d.pubs > 0))})</span>` : '')
      + `<span><i style="background:${KB_COL.pg}"></i>both (${n((d) => d.pubs > 0 && d.gen > 0)})</span>`
      + (extra.length ? `<span><i class="kb-dash"></i>reported, not predicted, dashed in the same colors (${fmtInt(extra.length)})</span>` : '') : '';
  return `<span class="kbhead">Edge color</span><div class="kbrow">${base}${red}</div>`;
}
const HOMO_RING = '#111111';   // a black ring (community colors include red) around a protein predicted to form a homodimer (its pair with itself past the cutoff)
const homoKey = (n) => `<span><i style="background:#fff;border:3px solid ${HOMO_RING};border-radius:50%;box-sizing:border-box"></i>predicted homodimer${n == null ? '' : ` (${fmtInt(n)})`}</span>`;
// Communities of a small network by Louvain (modularity, resolution 1): nodes 0..n-1, edges [u, v, weight]. Deterministic
// (nodes visited in index order), so the same network always gives the same grouping. Returns a community index per node,
// communities numbered by size, largest first; isolated nodes get their own.
function louvain(n, edges) {
  let comm = [...Array(n).keys()], N = n, E = edges.map(([u, v, w]) => [u, v, w]);
  const member = [...Array(n).keys()].map((i) => [i]);   // original nodes in each current node
  for (let level = 0; level < 10; level++) {
    const adj = [...Array(N)].map(() => new Map()), k = new Float64Array(N); let m2 = 0;
    for (const [u, v, w] of E) { adj[u].set(v, (adj[u].get(v) || 0) + w); adj[v].set(u, (adj[v].get(u) || 0) + w); k[u] += w; k[v] += w; m2 += 2 * w; }
    if (!m2) break;
    const c = [...Array(N).keys()], tot = Float64Array.from(k);
    let moved = true, any = false, rounds = 0;
    while (moved && rounds++ < 50) { moved = false;
      for (let i = 0; i < N; i++) {
        const ci = c[i], wt = new Map(); for (const [j, w] of adj[i]) if (j !== i) wt.set(c[j], (wt.get(c[j]) || 0) + w);
        tot[ci] -= k[i]; let best = ci, gain = (wt.get(ci) || 0) - tot[ci] * k[i] / m2;
        for (const [cj, w] of wt) { const g = w - tot[cj] * k[i] / m2; if (g > gain + 1e-12) { gain = g; best = cj; } }
        tot[best] += k[i]; if (best !== ci) { c[i] = best; moved = true; any = true; } }
    }
    if (!any) break;
    const ids = new Map(); c.forEach((x) => { if (!ids.has(x)) ids.set(x, ids.size); });
    const next = [...Array(ids.size)].map(() => []); for (let i = 0; i < N; i++) next[ids.get(c[i])].push(...member[i]);
    const agg = new Map(); for (const [u, v, w] of E) { const a = ids.get(c[u]), b = ids.get(c[v]), key = a < b ? a * ids.size + b : b * ids.size + a; agg.set(key, (agg.get(key) || 0) + w); }
    E = [...agg].map(([key, w]) => [Math.floor(key / ids.size), key % ids.size, w]); N = ids.size; member.length = 0; member.push(...next);
  }
  const order = member.map((m, i) => [i, m.length]).sort((a, b) => b[1] - a[1]), out = new Array(n);
  order.forEach(([i], rank) => { for (const x of member[i]) out[x] = rank; });
  return out;
}
// A network's labels in one layer drawn after every node, so no node covers a label; the query's label last. Returns
// the function each tick calls to keep the labels on their nodes.
function liftLabels(g, node) {
  const layer = g.append('g').attr('class', 'net-labels').style('pointer-events', 'none'), wraps = [];
  node.each(function (d) { const t = this.querySelector('text'); if (!t) return; const w = layer.append('g').datum(d).node(); w.appendChild(t); wraps.push(w); });
  wraps.sort((a, b) => (d3.select(a).datum().q ? 1 : 0) - (d3.select(b).datum().q ? 1 : 0)).forEach((w) => layer.node().appendChild(w));
  const ws = d3.selectAll(wraps);
  return () => ws.attr('transform', (d) => `translate(${d.x},${d.y})`);
}
// A Find partner box that suggests partners with their iLIS, in place of the browser's datalist (a second column looks
// different in every browser, and a native popup cannot be checked headless). items() returns the rows to list, each
// {name, best, dot, sty, tip}, in the order to show; the list filters by what is typed (substring, any case), opens on focus,
// click and typing, and picking a partner sets the box and fires input and change, so the box's own handlers run as if it
// had been typed. Arrow keys move, Enter picks, Escape closes. The list sits on the page body, so no card clips it.
const PSUG_ABORT = new Map();   // input id -> the controller of the window listeners its list added on the last draw
function partnerSuggest(input, items) {
  const id = `psug-${input.id}`, old = document.getElementById(id); if (old) old.remove();   // the page is drawn again on each route
  if (PSUG_ABORT.has(id)) PSUG_ABORT.get(id).abort();
  const ac = new AbortController(); PSUG_ABORT.set(id, ac);
  const list = document.createElement('div'); list.className = 'psug'; list.id = id; list.hidden = true; list.setAttribute('role', 'listbox'); document.body.appendChild(list);
  input.setAttribute('role', 'combobox'); input.setAttribute('aria-autocomplete', 'list'); input.setAttribute('aria-expanded', 'false'); input.setAttribute('aria-controls', id); input.autocomplete = 'off';
  let rows = [], at = -1, picking = false, shown = 0;
  const close = () => { list.hidden = true; input.setAttribute('aria-expanded', 'false'); input.removeAttribute('aria-activedescendant'); at = -1; };
  const mark = (n) => { at = n; [...list.querySelectorAll('.psug-row')].forEach((r, i) => { r.classList.toggle('on', i === n); r.setAttribute('aria-selected', i === n ? 'true' : 'false'); });
    if (n >= 0) { const r = list.querySelectorAll('.psug-row')[n]; r.scrollIntoView({ block: 'nearest' }); input.setAttribute('aria-activedescendant', r.id); } else input.removeAttribute('aria-activedescendant'); };
  const pick = (name) => { picking = true; close(); input.value = name; input.dispatchEvent(new Event('input', { bubbles: true })); input.dispatchEvent(new Event('change', { bubbles: true })); picking = false; };
  const place = () => {   // under the box; kept open while the page scrolls or the phone's keyboard resizes the window, closed when the box leaves the screen
    const r = input.getBoundingClientRect(); if (r.bottom < 0 || r.top > innerHeight) { close(); return; }
    const w = Math.max(230, r.width);
    list.style.minWidth = `${w}px`; list.style.left = `${Math.max(8, Math.min(r.left, innerWidth - w - 8))}px`; list.style.top = `${r.bottom + 4}px`;
    list.style.maxHeight = `${Math.max(140, Math.min(340, innerHeight - r.bottom - 16))}px`; };
  const open = () => {
    if (picking) return;
    if (!input.isConnected) { list.remove(); return; }
    const q = input.value.trim().toLowerCase(); rows = items().filter((it) => !q || it.name.toLowerCase().includes(q));
    if (!rows.length) { close(); return; }
    shown = Math.min(rows.length, 300);
    list.innerHTML = rows.slice(0, shown).map((it, i) => `<div class="psug-row" role="option" aria-selected="false" id="${id}-${i}" data-i="${i}" title="${esc(it.tip || '')}"><span class="mdot" style="background:${it.dot}"></span><span class="psug-n">${esc(it.name)}</span><span class="psug-v num" style="${it.sty}">${it.best.toFixed(3)}</span></div>`).join('')
      + (rows.length > shown ? `<div class="psug-more">${fmtInt(rows.length - shown)} more: keep typing</div>` : '');
    list.hidden = false; place(); input.setAttribute('aria-expanded', 'true'); mark(-1);
  };
  input.addEventListener('focus', open); input.addEventListener('click', open); input.addEventListener('input', open);
  input.addEventListener('keydown', (e) => {   // capture: before the box's own key handler, so Enter on a marked row picks it instead of searching
    if (e.key === 'ArrowDown' || e.key === 'ArrowUp') { e.preventDefault(); if (list.hidden) { open(); return; } mark(Math.max(0, Math.min(shown - 1, at + (e.key === 'ArrowDown' ? 1 : -1)))); }
    else if (e.key === 'Enter' && !list.hidden && at >= 0) { e.preventDefault(); e.stopImmediatePropagation(); pick(rows[at].name); }
    else if (e.key === 'Escape' && !list.hidden) { e.preventDefault(); e.stopPropagation(); close(); }   // preventDefault: Chrome's search box would also clear its text
    else if (e.key === 'Tab') close();
  }, true);
  list.addEventListener('mousedown', (e) => { const row = e.target.closest('.psug-row'); if (!row) return; e.preventDefault(); pick(rows[+row.dataset.i].name); });
  input.addEventListener('blur', () => setTimeout(close, 120));
  const again = (e) => { if (!list.hidden && !(e.target instanceof Node && list.contains(e.target))) place(); };   // a resize's target is the window, not a node
  addEventListener('resize', again, { signal: ac.signal }); addEventListener('scroll', again, { capture: true, signal: ac.signal });
}
const EWID = (a) => 0.5 + 4.3 * Math.max(0, Math.min(1, a / 0.8));
function resolveRow(sp, q) {   // a key, any screen's name, an accession, a gene symbol (exact case first), a CG number or an older name
  const h = resolveHow(sp, q); return h ? h.row : null;
}
// The same lookup, with the rule that matched (for the network builder's report). Maps built once per species, so a list of
// thousands of names does not scan every row per name; the first row in index order wins, as a scan would.
// iso: also take a UniProt isoform (P04637-2) as its canonical accession
const ACC_ISO = /^((?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}))-\d+$/i;
function resolveHow(sp, q, iso = false) {
  const keyRule = sp.manifest.keyedBy ? 'FlyBase ID' : 'UniProt accession';
  if (sp.byKey.has(q)) return { row: sp.byKey.get(q), how: keyRule };
  if (sp.byName.has(q)) return { row: sp.byName.get(q), how: 'screen name' };
  if (sp.byGene.has(q)) return { row: sp.byGene.get(q), how: 'gene symbol' };
  if (!sp.lk) { const first = (m, k, i) => { if (k && !m.has(k)) m.set(k, i); };
    const L = { syn: new Map(), gene: new Map(), acc: new Map(), ids: new Map(), lsyn: new Map() };
    sp.rows.forEach((r, i) => { for (const t of String(r.syn || '').split(' ')) first(L.syn, t, i); const k = sp.keys[i];
      first(L.gene, k.gene, i); first(L.acc, k.acc, i); k.ids.forEach((t) => first(L.ids, t, i)); k.syn.forEach((t) => first(L.lsyn, t, i)); });
    sp.lk = L; }
  const synRule = (t) => (/^CG\d+$/i.test(t) ? 'CG number' : 'older name');
  if (sp.lk.syn.has(q)) return { row: sp.rows[sp.lk.syn.get(q)], how: synRule(q) };   // an older name, exact case
  const l = q.toLowerCase(), hits = [['gene', 'gene symbol'], ['acc', 'UniProt accession'], ['ids', 'screen name']].filter(([m]) => sp.lk[m].has(l)).map(([m, how]) => ({ i: sp.lk[m].get(l), how }));
  if (hits.length) { const h = hits.reduce((a, b) => (b.i < a.i ? b : a)); return { row: sp.rows[h.i], how: h.how }; }
  if (sp.lk.lsyn.has(l)) return { row: sp.rows[sp.lk.lsyn.get(l)], how: synRule(q) };
  const m = iso && q.match(ACC_ISO);
  if (m) { const h = resolveHow(sp, m[1].toUpperCase()); if (h) return { row: h.row, how: 'isoform accession' }; }
  return null;
}
// A pasted or loaded list: plain names (commas, spaces, new lines) or a table (tab or comma separated, two or more columns).
// For a table, every column is scored by how many of its cells resolve; a column of numbers is never read as names.
// column titles that are also gene names somewhere (prey, bait, target...): in the first row they are a header, not a protein
const HEADWORD = /^(bait|baits|prey|preys|gene|genes|gene[ _]?(name|symbol|id)|protein|proteins|protein[ _]?[ab12]|gene[ _]?[ab12]|symbol|symbols|name|names|id|ids|uniprot|accession|fbgn|flybase|target|targets|source|query|partner|partners|hit|hits|interactor|interactors)$/i;
const NUMCELL = /^[-+]?(\d+\.?\d*|\.\d+)(e[-+]?\d+)?%?$/i, BLANK = /^(na|nan|n\/a|null|none|-|—|–|\.)$/i;   // blanks a number column may hold
function readIdTable(sp, text) {
  const lines = text.split(/\r?\n/).filter((l) => l.trim());
  const delim = lines.length >= 2 && lines.every((l) => l.includes('\t')) ? '\t' : lines.length >= 2 && lines.filter((l) => l.includes(',')).length >= 0.8 * lines.length && lines.some((l) => /,\s*\S/.test(l) && /\S\s*,/.test(l)) && new Set(lines.map((l) => l.split(',').length)).size <= 2 ? ',' : null;
  if (!delim) return null;
  const splitRow = (l) => { if (delim === '\t' || !l.includes('"')) return l.split(delim);   // CSV cells may quote a comma
    const out = []; let cur = '', q = false; for (let k = 0; k < l.length; k++) { const ch = l[k];
      if (q) { if (ch === '"' && l[k + 1] === '"') { cur += '"'; k++; } else if (ch === '"') q = false; else cur += ch; }
      else if (ch === '"') q = true; else if (ch === ',') { out.push(cur); cur = ''; } else cur += ch; }
    out.push(cur); return out; };
  const cells = lines.map((l) => splitRow(l).map((c) => c.trim().replace(/^"(.*)"$/, '$1').trim()));
  const nc = Math.max(...cells.map((r) => r.length));
  if (nc < 2) return null;
  const sample = cells.slice(0, 400), cols = [];
  for (let c = 0; c < nc; c++) {
    const v = sample.slice(1).map((r) => r[c] || '').filter((x) => x && !BLANK.test(x)), num = v.length && v.every((x) => NUMCELL.test(x));
    const res = num ? [] : v.map((x) => { const h = resolveHow(sp, x, true); return h ? h.row.i : null; }).filter((i) => i != null);
    const hits = new Set(res).size, raw = res.length, headHit = !HEADWORD.test(sample[0][c] || '') && !!(sample[0][c] && resolveHow(sp, sample[0][c], true));
    const few = v.length >= 6 && new Set(v).size <= Math.max(3, 0.1 * v.length);   // a handful of values repeated (bait/prey, up/down): a group column, not names
    cols.push({ c, name: sample[0][c] || `column ${c + 1}`, num, few, hits: few ? 0 : hits, raw: few ? 0 : raw, n: v.length, headHit: few ? false : headHit });
  }
  const header = cols.every((x) => !x.headHit || x.num);   // a first row whose cells are not names is a header
  if (!header) cols.forEach((x) => { if (!x.num) { x.hits += x.headHit ? 1 : 0; x.raw += x.headHit ? 1 : 0; x.n += sample[0][x.c] ? 1 : 0; } });   // no header: the first row is data too
  // hits count different proteins, so a column that repeats one word (a "prey" that is also an old gene name) cannot outscore the names
  // no header, no number column and every column mostly names: a list wrapped over lines, so every name is read;
  // and when no column holds a single name we know, a list too, so each unknown name is reported rather than taken for a header
  if (!header && !cols.some((x) => x.num) && cols.every((x) => !x.n || x.raw >= 0.5 * x.n)) return null;   // per cell: a pair list repeats names
  if (cols.every((x) => x.num || !x.hits) && !cols.some((x) => x.headHit)) return null;
  if (!header) cols.forEach((x) => { x.name = `column ${x.c + 1}`; });
  return { cols, header, rows: header ? cells.slice(1) : cells };
}

// The names a reader gives (network builder, nested network): a plain list, or the chosen column of a table; lines starting
// with # are notes. → { toks, T, col, named, data, note }: the names in order, the table (or null), its name column, how many
// columns hold names, its other columns (numbers or a few categories, for coloring) and a note on how the text was read
function readIdInput(sp, raw, colPref = null) {
  const text = String(raw || '').split(/\r?\n/).filter((l) => !/^\s*#/.test(l)).join('\n'), T = readIdTable(sp, text);
  if (!T) { let body = text, note = ''; const ls = text.split(/\r?\n/).map((l) => l.trim()).filter(Boolean);
    if (ls.length >= 3 && !/[\s,;]/.test(ls[0]) && !resolveHow(sp, ls[0], true)) {   // one column with a header line: the header is not a missing name
      const rest = ls.slice(1, 200); if (rest.filter((x) => resolveHow(sp, x, true)).length >= 0.5 * rest.length) { body = ls.slice(1).join('\n'); note = `first line “${ls[0]}” read as a header`; } }
    return { toks: [...new Set(body.split(/[\s,;]+/).map((t) => t.trim()).filter(Boolean))], T: null, col: null, named: 0, data: [], note }; }
  const named = T.cols.filter((x) => !x.num), best = named.reduce((a, b) => (b.hits > a.hits ? b : a), named[0] || T.cols[0]);
  const strong = named.filter((x) => x.n && x.raw >= 0.5 * x.n);   // columns that are mostly names we know: a pair list (bait, prey), or a symbol beside an accession
  const all = strong.length >= 2 && (colPref == null || colPref === 'all');
  const col = all ? best.c : colPref == null || colPref === 'all' || !T.cols[colPref] || T.cols[colPref].num ? (best ? best.c : 0) : colPref;
  const nums = T.cols.filter((x) => x.num).map((x) => x.name);
  const note = `a table of ${fmtInt(T.rows.length)} row${T.rows.length === 1 ? '' : 's'}, names read from ${all ? strong.map((x) => `“${x.name}”`).join(' and ') : `“${T.cols[col].name}”`}${nums.length ? `; number column${nums.length === 1 ? '' : 's'} not read as names: ${nums.map((n) => `“${n}”`).join(', ')}` : ''}`;
  // a row whose name is unknown is tried by its other ID columns (a UniProt accession beside an old symbol)
  const alts = T.cols.filter((x) => !x.num && x.c !== col && x.hits >= 3);
  const rowName = T.rows.map((r) => { const nm = (r[col] || '').trim(); if (!nm || resolveHow(sp, nm, true)) return nm;
    for (const x of alts) { const a = (r[x.c] || '').trim(); if (a && resolveHow(sp, a, true)) return a; } return nm; });
  // the other columns, to color the proteins by: numbers, or a few categories (a column of names is not one)
  const data = []; for (const x of T.cols) { if (x.c === col) continue; const vals = new Map();
    T.rows.forEach((r, k) => { const nm = rowName[k], v = (r[x.c] || '').trim(), h = nm && v && !BLANK.test(v) && resolveHow(sp, nm, true); if (h && !vals.has(h.row.i)) vals.set(h.row.i, x.num ? parseFloat(v) : v); });
    const kinds = new Set(vals.values());
    if (x.num) data.push({ name: x.name, kind: 'num', vals }); else if (kinds.size >= 2 && kinds.size <= 12 && x.hits <= 0.2 * Math.max(1, x.n)) data.push({ name: x.name, kind: 'cat', vals }); }
  const toks = all ? [...new Set(T.rows.flatMap((r) => strong.map((x) => (r[x.c] || '').trim())).filter(Boolean))] : [...new Set(rowName.filter(Boolean))];
  return { toks, T, col: all ? 'all' : col, strongN: strong.length, named: named.length, data, note };
}
// Residue lookup input: "983", "T983", "T983A", "p.T983A" or "Thr983Ala" → { n, wt, mut } (one-letter codes)
const AA3 = { ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G', HIS: 'H', ILE: 'I', LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S', THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V', TER: '*' };
function parseRes(v) {
  let t = String(v || '').trim().toUpperCase().replace(/^P\./, '');
  t = t.replace(/^([A-Z]{3})(\d+)([A-Z]{3})?$/, (m, a, n, b) => (AA3[a] || '?') + n + (b ? AA3[b] || '?' : ''));
  const m = t.match(/^([A-Z])?(\d+)([A-Z*])?$/);
  return m && +m[2] > 0 ? { n: +m[2], wt: m[1] || '', mut: m[3] || '' } : { n: 0, wt: '', mut: '' };
}
const resLabel = (r) => `${r.wt}${r.n}${r.mut}`;
async function viewProtein(spId, q, setId = '', iso = null) {   // setId: only that thematic set's predictions; iso: one isoform ('reference' or a construct)
  const gen = ROUTE, gone = () => stale(gen);   // after every await: stop if the reader has moved to another page
  /* BioGRID: partners with a reported interaction, from the species file (reported); marked once it arrives */
  let RESV = parseRes(hashPath().q.get('res')), RESQ = RESV.n, resScrolled = false;   // residue lookup: ?res=983 or ?res=T983A in the link
  let KB = null, scKnown = 'pass';   // rings in the Overview: 'pass' (partners past the 10% FPR cutoff), 'all' or 'off'
  const kbOf = (key) => { if (!KB) return null; const r = sp.byKey.get(key); if (!r || r.i == null || r.i === P.i) return null;
    const ph = KB.pubs(P.i, r.i), ge = KB.gen(P.i, r.i); return ph || ge ? { ph, ge } : null; };
  const kbTip = (kb) => `reported in BioGRID ${KB.release}: ${[kb.ph ? `physical, ${kb.ph} publication${kb.ph === 1 ? '' : 's'}` : '', kb.ge ? `genetic, ${kb.ge} publication${kb.ge === 1 ? '' : 's'}` : ''].filter(Boolean).join('; ')}`;
  const qp = new URLSearchParams(); if (setId) qp.set('set', setId); if (iso) qp.set('iso', iso); if (RESQ) qp.set('res', resLabel(RESV));   // kept through the redirect to the canonical key
  const sp = await species(spId), P = resolveRow(sp, q), qs = qp.toString() ? `?${qp}` : '';
  if (gone()) return;
  if (!P) { app.innerHTML = `<div class="empty">No protein “${esc(q)}” in the ${esc(sp.reg.label.toLowerCase())} screens. <a href="#/${sp.id}">Search ${esc(sp.reg.label.toLowerCase())} proteins</a></div>`; return; }
  if (P.key !== q) { location.replace(`#/${sp.id}/${P.key}${qs}`); return; }
  document.title = `${P.gene} · LIVIA Atlas`;
  const flags = [];
  if (P.status === 'renamed') flags.push(`<span class="flag">named ${esc(P.occ.map((o) => o.name).filter((n, i, a) => a.indexOf(n) === i).join(' / '))} in the screens; UniProt renamed it</span>`);
  if (P.status === 'unreviewed') flags.push('<span class="flag">unreviewed UniProt entry</span>');
  if (P.status === 'other species') flags.push('<span class="flag">not a fly protein: folded as a partner of fly proteins</span>');
  if (P.status === 'construct') flags.push('<span class="flag">an engineered construct or a retired gene, kept under its screen name</span>');
  if (P.status === 'obsolete') flags.push('<span class="flag">UniProt has since retired this entry; the sequence is the one the screen folded</span>');
  const fbLink = /^FBgn\d{7}$/.test(P.key) ? `<a href="https://flybase.org/reports/${P.key}" target="_blank" rel="noopener">FlyBase ${P.key}</a>` : '';
  const nav = [['c-partners', 'Overview'], ['c-sites', 'Binding sites'], ['c-info', 'Clusters'], ['c-3d', '3D structure'], ['c-freq', 'Frequency'], ['c-fp', 'Fingerprint'], ['c-res', 'Residues'], ['c-net', 'Network'], ['c-pt', 'Partners']];   // answers first (who, where, which share), then evidence, then tools
  const chips = '<div class="chips cl-chips" data-chips></div>';
  const xticks = '<label class="xt">x-ticks <input type="number" class="xticks" min="2" max="40" placeholder="auto"></label>';
  const occ = (await Promise.all(P.occ.map(async (o) => { try { return { ...o, ds: await dataset(sp.dsIds[o.di]) }; } catch (e) { return null; } }))).filter(Boolean);
  if (gone()) return;
  const dataMenu = occ.map((o, i) => { const u = bundleUrl(o.ds, o.name), label = `${sp.dsShort[o.di]}${occ.filter((x) => x.di === o.di).length > 1 ? ' · ' + o.name : ''}`;
    if (o.ds.reg.zip) return `<div class="dm-row"><span class="src" style="--c:${sp.dsColor[o.di]}">${esc(label)}</span><a href="#" data-clip="${i}">Open in LIVIA cLIP ↗</a><a href="#" data-dl="${i}">Download .zip</a></div>`;   // inside the archive: read here, then handed to cLIP or saved
    return `<div class="dm-row"><span class="src" style="--c:${sp.dsColor[o.di]}">${esc(label)}</span><a href="${LIVIA}clip.html?data=${encodeURIComponent(u)}&gene=${encodeURIComponent(P.gene)}" target="_blank" rel="noopener">Open in LIVIA cLIP ↗</a><a href="${u}" download>Download .zip</a></div>`; }).join('');
  app.innerHTML = `<div class="crumbs"><a href="#/">Atlas</a> / <a href="#/${sp.id}">${esc(sp.reg.label)}</a> / <a href="#/${sp.id}/${P.key}">${esc(P.gene)}</a></div>
    <div class="phead"><div><h1>${esc(P.gene)}</h1><div class="pname" title="${esc(P.name)}">${esc(short(P.name) || P.id)}</div>
      <div class="ids">${fbLink}${uniprotLink(P.acc)}${fbLink ? '' : `<span>${esc(P.id)}</span>`}${P.clen ? `<span>${fmtInt(P.clen)} aa</span>` : ''}</div>
      ${P.virus ? `<div class="srcs">Virus <a href="#/${sp.id}/taxon/${P.virus.taxid}">${esc(P.virus.name)}</a> <span class="muted">· ${fmtInt(P.virus.n)} proteins, every pair folded${P.key.includes('_p') ? ' · a mature peptide of a polyprotein, numbered from its own first residue' : ''}</span></div>` : ''}
      <div class="srcs scope" id="scope">In ${srcBadges(sp, P.src)}</div>
      <div class="flags" id="flags">${flags.join('')}</div>
      <div class="actions"><details class="dmenu"><summary class="btn">Data &amp; LIVIA cLIP ▾</summary><div class="dm-pop">${dataMenu}</div></details>${citeBtn(`${P.gene} (${sp.reg.label})`, setId && sp.dsIds.includes(setId) ? [setId] : sp.dsIds.filter((_, i) => P.src & (1 << i)))}</div></div>
      <div class="kpis"><div class="kpi"><b id="kp-all">${fmtInt(P.partners)}</b><span>partners predicted</span></div><div class="kpi f10"><b id="kp-10">${fmtInt(P.pos10)}</b><span>past 10% FPR${cutNote(10)}</span></div>
        <div class="kpi f5"><b id="kp-5">${fmtInt(P.pos5)}</b><span>past 5% FPR${cutNote(5)}</span></div><div class="kpi f1"><b id="kp-1">${fmtInt(P.pos1)}</b><span>past 1% FPR${cutNote(1)}</span></div></div></div>
    <div class="srcs scope isorow" id="isorow" hidden></div>
    <div class="setbar" id="setbar" hidden></div>
    <nav class="subnav" aria-label="Sections">${nav.map(([t, l]) => `<button data-t="${t}">${l}</button>`).join('')}</nav>
    <div class="card" id="c-iso" hidden></div>
    <div class="card" id="c-partners"><div class="card-head"><div><h2>Overview</h2><div class="muted">each partner by its best model</div></div>
        <div class="controls" style="margin:0"><label>Y <select id="sc-y"></select></label><label>X <select id="sc-x"></select></label>
          <label>show <select id="sc-pts"><option value="partner">one per partner (best model)</option><option value="all">every prediction</option><option value="rank1">rank-1 per pair</option></select></label>
          <input type="search" id="sc-find" placeholder="Find partner" aria-label="Find a partner in the plot" style="width:140px"></div></div>
      <div class="overview"><div><div class="muted" id="sc-rho" style="margin:2px 0 8px"></div><div class="muted" id="sc-def" style="margin:-4px 0 8px" hidden></div><div class="plot" id="scat2"><canvas id="scatter-canvas"></canvas></div><div class="legend" id="sc-legend"></div></div>
        <div><h3>Top partners <span class="muted">by the share of models past the 10% cutoff, then average iLIS</span></h3>
          <div class="tl-head"><span></span><span>Partner</span><span>Cluster</span><span>iLIS<br>best</span><span>iLIS<br>avg</span><span>ipTM<br>best</span><span>ipTM<br>avg</span><span title="models past the 10% FPR cutoff, of the pair's models">models<br>past</span></div>
          <ol class="toplist" id="toplist"></ol>
          <div class="legend tl-key"><span>FPR band, each value by its own benchmarked cutoff</span><span class="tl-keys"><span><i style="background:${BAND[1]}"></i>1%</span><span><i style="background:${BAND[5]}"></i>5%</span><span><i style="background:${BAND[10]}"></i>10%</span><span><i style="background:#A7B2BF"></i>below</span></span></div></div></div></div>
    <div class="card" id="c-sites"><div class="card-head"><div><h2>Predicted binding sites</h2><div class="muted" id="clip-sub">Loading the predictions…</div></div>
        <div class="clip-ctl"><label class="ctl muted" title="a residue number or a variant (983, T983A): which predictions, partners and sites contact it; the link keeps it">Residue<input type="text" id="res-q" placeholder="983 or T983A" spellcheck="false" autocomplete="off" value="${esc(RESQ ? resLabel(RESV) : '')}"></label><div class="ctl"><span class="muted">iLIS cutoff</span><div class="seg" id="cut-seg">${[10, 5, 1].map((f) => `<button data-f="${f}" class="${f === 10 ? 'on' : ''}">${f}% FPR · ${CUT[f].toFixed(3)}</button>`).join('')}</div></div></div></div>
      <p class="sites-answer" id="sites-answer">Finding the binding sites…</p><p class="res-look" id="res-look" hidden></p>
      <div class="plot" id="sites-map"></div><p class="sites-bg muted" id="sites-bg"></p><p class="sites-more" id="sites-more" hidden></p><div class="domlegend" id="sites-domains"></div>
      <p class="muted sites-note">How the sites are found: LIVIA cLIP, run in your browser on the predictions shown, takes the residues of ${esc(P.gene)} that
        each prediction past the cutoff contacts (cLIR: PAE ≤ 12 Å and Cβ ≤ 8 Å) as its fingerprint, compares fingerprints by cosine distance, joins them by
        average linkage and chooses the number of clusters by silhouette. Each cluster is a site, numbered by size; the map has one lane per major site, one
        with at least 3 partners and 2% of all partners, and shows the minor sites when asked. Solid marks
        a site's footprint, the residues at least 30% of its predictions contact; lighter shades, how often the other residues are contacted. Click a lane to
        show only that site on the page; the partners of each site are listed under Clusters. Type a residue number or a variant (for example <span id="res-eg">T983A</span>) to see which
        predictions, partners and sites contact it; the link keeps it.</p><p class="note" id="clip-aside" hidden></p></div>
    <div class="card" id="c-info"><div class="card-head"><h2>Clusters</h2><span class="muted">the partners of each binding site · Cluster n (proteins / predictions) · largest first</span></div><div class="clinfo" id="cluster-info"></div><div class="legend" id="info-kb" hidden></div></div>
    <div class="card" id="c-3d"><div class="card-head"><h2>3D structure</h2><span class="muted" id="struct-badge"></span></div>
      <p class="muted" style="margin:2px 0 6px">${esc(P.gene)} as predicted alone in the AlphaFold Database, residues colored by the cluster that consensus-contacts them · click clusters to isolate.
        <b>C<i>n</i> (N)</b>: N = predictions (AlphaFold ranks) in that cluster.</p>
      <div class="controls"><span>Highlight residues contacted by ≥</span><select id="commonality" aria-label="Share of a cluster's predictions that must contact a residue" title="Fraction of a cluster's members (predictions) that must contact a residue for it to take the cluster's color">${[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9].map((v) => `<option value="${v}"${v === 0.5 ? ' selected' : ''}>${Math.round(v * 100)}%</option>`).join('')}<option value="custom">custom…</option></select>
        <span id="commonality-custom" hidden><input type="number" id="commonality-pct" min="1" max="100" step="1" style="width:64px" title="Any value from 1 to 100%; applied on Enter or when the field loses focus"> %</span>
        <span>of a cluster's members · color only clusters with ≥</span><input type="number" id="min-cluster" min="1" value="5" style="width:60px" aria-label="Smallest cluster shown"><span>predictions</span></div>
      ${chips}
      <div class="controls"><div class="ctl"><span>Color by</span><div class="seg" id="cmode"><button data-m="cluster" class="on">cluster</button><button data-m="plddt">pLDDT</button><button data-m="am" id="cm-am" hidden>AlphaMissense</button></div></div>
        <span id="am-cut-wrap" hidden>AM pathogenicity average ≥ <select id="am-cutoff">${[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9].map((v) => `<option value="${v}">${Math.round(v * 100)}%</option>`).join('')}</select></span></div>
      <div class="viewer3d"><iframe id="viewer3d-frame" title="3D structure viewer"></iframe><div class="v3d-msg" id="v3d-msg">Loading the AlphaFold DB model…</div></div>
      <div class="legend" id="legend-3d"></div></div>
    <div class="card" id="c-freq"><div class="card-head"><h2>Contact residue frequency</h2><div class="card-tools"><span class="muted">predictions contacting each residue, colored by their most frequent cluster</span>${xticks}</div></div>
      ${chips}<div class="plot" id="freq-wrap"></div><div class="domlegend" id="freq-domains"></div><div class="hot" id="hot"></div></div>
    <div class="card" id="c-fp"><div class="card-head"><h2>Clustered interaction fingerprint</h2><div class="card-tools"><span class="muted">one row per prediction, in dendrogram order · hover for the partner, click to open the pair</span>
        <input type="search" id="fp-find" placeholder="Find partner" style="width:140px" aria-label="Find a partner in the fingerprint">${xticks}</div></div>
      ${chips}<div class="plot" id="fp-wrap"></div><p class="sites-more" id="fp-found" hidden></p>
      <div class="legend"><span><i style="background:#08306B"></i>contact residue (cLIR)</span><span><i style="background:#F7FBFF;box-shadow:inset 0 0 0 1px #C9D6E3"></i>no contact</span><span>left: dendrogram and cluster of each prediction</span></div></div>
    <div class="card" id="c-res"><div class="card-head"><h2>Interaction Residues</h2>
        <div class="controls" style="margin:0"><select id="res-partner" aria-label="Partner" style="max-width:300px"></select><input type="search" id="res-find" placeholder="Find partner" aria-label="Find a partner for the residue view" style="width:130px"><select id="res-rank" aria-label="Model" style="max-width:280px"></select><span id="res-struct"></span></div></div>
      <div class="legend" style="margin:2px 0 12px"><span><i style="background:#E0E0E0"></i>not an interaction residue</span><span><i style="background:#80CBC4"></i><i style="background:#FFAB91;margin-left:-2px"></i>interaction residue (LIR: PAE ≤ 12 Å)</span><span><i style="background:#00897B"></i><i style="background:#E64A19;margin-left:-2px"></i>contact (cLIR: also Cβ ≤ 8 Å)</span></div>
      <div id="res-body"></div></div>
    <div class="card" id="c-net"><div class="card-head"><h2>Network</h2>
      <div class="controls" style="margin:0"><label class="ctl">Partners<select id="net-n"><option>30</option><option>60</option><option selected>100</option><option>200</option></select></label>
        <label class="ctl">Cutoff<select id="net-cut"><option value="10">10% FPR · iLIS ${CUT[10]}</option><option value="5">5% FPR · iLIS ${CUT[5]}</option><option value="1">1% FPR · iLIS ${CUT[1]}</option></select></label>${edgeCtl('net', true)}</div></div>
      <p class="muted" style="margin:2px 0 12px">${esc(P.gene)} at the center; partners sit closer the higher their iLIS and are filled with their cluster color. Edges between partners join partners
        predicted to bind each other, in any screen. Drag to move, scroll to zoom, click to open. <a href="#/${sp.id}/network?ids=${encodeURIComponent(P.gene)}&add=top&k=10">Build a network with other proteins →</a></p>
      <div class="net" id="net"><div class="loading" style="padding:20px">Loading the network…</div></div>
      <div class="netkey"><div class="kbkey" id="net-kbkey"></div>
        <div><span>Edge width · ${sp.one ? 'iLIS' : 'average iLIS'}</span><svg id="net-w" width="260" height="30" aria-hidden="true"></svg></div></div>
      <div class="legend" id="net-legend"></div><div id="net-x"></div></div>
    <div class="card" id="c-pt"><div class="card-head"><h2>Partners <span class="muted" id="pt-note"></span></h2>
      <div class="controls" style="margin:0"><select id="pt-src" aria-label="Which screens"><option value="0">all screens</option>${sp.dsIds.map((_, di) => (P.src & (1 << di) ? `<option value="${1 << di}">${esc(sp.dsShort[di])}</option>` : '')).join('')}</select>
        <select id="pt-band" aria-label="Which partners, by FPR band"><option value="10">past 10% FPR (iLIS ${CUT[10]})</option><option value="5">past 5% FPR (iLIS ${CUT[5]})</option><option value="1">past 1% FPR (iLIS ${CUT[1]})</option><option value="0">all predicted</option></select>
        <input type="search" id="pt-filter" placeholder="Filter partners" style="width:180px"></div></div>
      <div class="tbl-wrap"><table class="pt" id="pt"></table></div><div class="pager" id="pager"></div></div>`;
  { const bar = $('.subnav'); let cur = null, raf = 0;
    const spy = () => { raf = 0; if (!bar || !bar.isConnected) { window.removeEventListener('scroll', onScroll); return; }
      const lim = bar.getBoundingClientRect().bottom + 24; let on = null;
      for (const b of bar.querySelectorAll('button')) { const t = document.getElementById(b.dataset.t); if (t && !t.hidden && t.getBoundingClientRect().top <= lim) on = b; }
      if (on === cur) return; cur = on; bar.querySelectorAll('button').forEach((b) => b.classList.toggle('on', b === on));
      if (on && bar.scrollWidth > bar.clientWidth + 1) bar.scrollTo({ left: Math.max(0, on.offsetLeft - (bar.clientWidth - on.offsetWidth) / 2), behavior: 'auto' }); };
    const onScroll = () => { if (!raf) raf = requestAnimationFrame(spy); };
    window.addEventListener('scroll', onScroll, { passive: true }); requestAnimationFrame(spy); }
  app.querySelectorAll('.subnav button').forEach((b) => b.onclick = () => { const t = document.getElementById(b.dataset.t); if (t) window.scrollTo({ top: t.getBoundingClientRect().top + window.scrollY - 112, behavior: 'auto' }); });
  { const dm = $('.dmenu'), shut = (e) => { if (!dm || !dm.isConnected) { document.removeEventListener('click', shut); return; } if (dm.open && !dm.contains(e.target)) dm.open = false; }; document.addEventListener('click', shut); }
  app.querySelectorAll('[data-clip]').forEach((a) => a.onclick = async (e) => {   // a bundle inside a screen archive: read it here and hand its bytes to cLIP
    e.preventDefault(); const o = occ[+a.dataset.clip], w = window.open(`${LIVIA}clip.html?post=1`, '_blank'); if (!w) return;
    try { const raw = await bundleRaw(o.ds, o.name), data = raw.isoforms ? await (await wholeBundle(o.ds, raw)).arrayBuffer() : raw.bytes.slice(0);
      handTo(w, { type: 'livia-load', name: `${o.name}.zip`, data, gene: P.gene }); }
    catch (err) { a.textContent = 'Not available'; try { w.close(); } catch (x) { /* already closed */ } } });
  app.querySelectorAll('[data-dl]').forEach((a) => a.onclick = async (e) => {   // a bundle inside a screen archive: read it, save it
    e.preventDefault(); const o = occ[+a.dataset.dl];
    try { const raw = await bundleRaw(o.ds, o.name), u = URL.createObjectURL(raw.isoforms ? await wholeBundle(o.ds, raw) : new Blob([raw.bytes], { type: 'application/zip' }));
      const d = document.createElement('a'); d.href = u; d.download = `${o.name}.zip`; d.click(); setTimeout(() => URL.revokeObjectURL(u), 5000); }
    catch (err) { a.textContent = 'Not available'; } });

  // Scope: one screen of the species (human kinase–TF …) or one thematic set of a screen (FlyPredictome's kinase–TF
  // screen, a source category …). The "In" row switches it; a banner says what is shown; links keep it.
  const TS0 = (await Promise.all(occ.map((o) => setsOf(o.ds).catch(() => null)))).find(Boolean) || null;
  if (gone()) return;
  const dsScope = setId ? occ.find((o) => o.ds.id === setId) : null;
  const SET = dsScope ? { id: setId, type: 'dataset', title: dsScope.ds.reg.title, short: dsScope.ds.reg.short, color: dsScope.ds.reg.color, source: dsScope.ds.manifest.source }
    : setId && TS0 ? TS0.byId.get(setId) || null : null;
  const scopeQ = SET ? `?set=${encodeURIComponent(SET.id)}` : '';
  let B, B0;
  try { B = await merged(sp, P, SET ? SET.id : ''); B0 = SET ? await merged(sp, P) : B; } catch (e) { if (!gone()) $('#clip-sub').textContent = e.message; return; }
  if (gone()) return;
  // Isoforms folded separately (a construct too unlike the reference to be drawn on it): the page shows one at a time,
  // every card on its predictions; the "Isoform" row switches, the Isoforms card compares them. BA keeps them all.
  const BA = B, pick = BA.choices && (BA.choices.length > 1 || (BA.split && BA.choices[0] && BA.choices[0].id));   // a split gene in a set with isoform models only opens on one
  const ISO = pick ? BA.choices.find((c) => c.id === (iso === 'reference' ? '' : iso)) || BA.choices[0] : null;
  const WORD = ISO && BA.choices.every((c) => !c.id || (BA.cons.get(c.id) || {}).kind === 'isoform') ? 'Isoform' : 'Construct';
  const isoName = (c) => (c.id ? String((BA.cons.get(c.id) || {}).label || c.id).replace(P.gene + ' ', '').replace(/ \(([\d,]+) aa\)$/, ' · $1 aa') : `reference · ${fmtInt(c.len)} aa`);
  const isoHref = (c) => `#/${sp.id}/${P.key}?${SET ? `set=${encodeURIComponent(SET.id)}&` : ''}iso=${c.id ? encodeURIComponent(c.id) : 'reference'}`;
  if (ISO) {
    if (BA.split && ISO.id) {   // a split gene (v1.2): this isoform's predictions are in a file of their own, read now
      $('#clip-sub').textContent = 'Loading this isoform…';
      try { B = await BA.load(ISO.id); } catch (e) { if (!gone()) $('#clip-sub').textContent = e.message; return; }
      if (gone()) return;
    } else B = BA.only(ISO.id);
    const row = $('#isorow'); row.hidden = false;
    row.innerHTML = `<span class="lbl">${WORD}</span>${BA.choices.map((c) => `<a class="src scope-chip${c === ISO ? ' on' : ''}" style="--c:#1A5276" href="${isoHref(c)}"
      title="${esc(c.label)}: ${fmtInt(c.n)} models">${esc(isoName(c))} <span class="n">${fmtInt(c.n)}</span></a>`).join('')}`;
    if (ISO !== BA.choices[0]) document.title = `${P.gene} · ${isoName(ISO)} · LIVIA Atlas`;
  }
  {
    const n = new Map(), by = new Map();   // each screen and set counts its own models, a sequence pair once
    for (const p of B0.preds) for (const k of [sp.dsIds[p.di], ...(p.tags || [])]) { if (!by.has(k)) by.set(k, []); by.get(k).push(p); }
    for (const [k, ps] of by) n.set(k, countOnce(ps));
    let total = B0.counted.length;
    if (B0.split) for (const c of B0.split.summary.choices) {   // the isoform files, counted from the gene's summary
      const d = sp.dsIds[B0.split.part.di]; n.set(d, (n.get(d) || 0) + c.models); total += c.models;
      for (const [t, k] of Object.entries(c.sets || {})) n.set(t, (n.get(t) || 0) + k); }
    const screens = occ.length > 1 ? occ.map((o) => ({ id: o.ds.id, short: o.ds.reg.short, color: o.ds.reg.color, title: o.ds.reg.title })) : [];
    const sets = TS0 ? TS0.list.filter((x) => n.get(x.id)).sort((x, y) => (x.type === 'screen' ? 0 : 1) - (y.type === 'screen' ? 0 : 1)) : [];
    const chip = (x) => `<a class="src scope-chip${SET && SET.id === x.id ? ' on' : ''}" style="--c:${x.color}" href="#/${sp.id}/${P.key}?set=${encodeURIComponent(x.id)}"
      title="${esc(x.title)}: ${fmtInt(n.get(x.id) || 0)} models">${esc(x.short)} <span class="n">${fmtInt(n.get(x.id) || 0)}</span></a>`;
    if (screens.length || sets.length) $('#scope').innerHTML = `<span class="lbl">In</span><a class="src scope-chip all${SET ? '' : ' on'}" href="#/${sp.id}/${P.key}"
      title="every screen and set">All <span class="n">${fmtInt(total)}</span></a>${screens.map(chip).join('')}${screens.length && sets.length ? '<span class="sep"></span>' : ''}${sets.map(chip).join('')}`;
  }
  if (SET) {
    const who = SET.source && SET.source.citation ? `${SET.source.citation.split(' ')[0]} et al. ${(SET.source.citation.match(/\((\d{4})\)/) || [])[1] || ''}` : '';
    const what = SET.type === 'dataset' ? `the ${SET.title.charAt(0).toLowerCase() + SET.title.slice(1)}` : SET.type === 'screen'
      ? `the ${SET.title.charAt(0).toLowerCase() + SET.title.slice(1)} of ${TS0.ds.reg.short}` : `the ${SET.title} category of ${TS0.ds.reg.short}`;
    const about = SET.type === 'dataset' ? `#/datasets/${SET.id}` : `#/datasets/${TS0.ds.id}/${SET.id}`;
    const bar = $('#setbar'); bar.hidden = false;
    bar.innerHTML = `<span class="src" style="--c:${SET.color}">${esc(SET.short)}</span><span>Only ${esc(what)}${who ? ` (<a href="${esc(SET.source.url)}" target="_blank" rel="noopener">${esc(who)}</a>)` : ''}:
      ${fmtInt(BA.counted.length + (BA.split ? BA.split.others.reduce((a, c) => a + c.n, 0) : 0))} of ${fmtInt(B0.counted.length + (B0.split ? B0.split.summary.choices.reduce((a, c) => a + c.models, 0) : 0))} models.</span><span><a href="#/${sp.id}/${P.key}">Show every prediction</a> · <a href="${about}">about this ${SET.type === 'dataset' ? 'screen' : 'set'}</a></span>`;
    document.title = `${P.gene} · ${SET.short} · LIVIA Atlas`;
  }
  if (SET || ISO) {   // the tiles count what the page shows
    const others = B.partners.filter((x) => x.id !== P.key && !x.rep), cnt = (c) => others.filter((x) => x.best >= c).length;
    $('#kp-all').textContent = fmtInt(others.length); $('#kp-10').textContent = fmtInt(cnt(CUT[10])); $('#kp-5').textContent = fmtInt(cnt(CUT[5])); $('#kp-1').textContent = fmtInt(cnt(CUT[1]));
  }
  // One model per pair (the viral release, or a page whose screens folded one model): best = average, so one value is shown
  const ONE = sp.one || B.partners.every((x) => x.counted.length <= 1);
  if (ONE) {
    $('#c-partners .card-head .muted').textContent = 'each partner, one model per pair'; $('#sc-pts').closest('label').style.display = 'none';
    $('#c-partners h3 .muted').textContent = 'by iLIS';
    $('.tl-head').innerHTML = '<span></span><span>Partner</span><span>Cluster</span><span>iLIS</span><span>ipTM</span>';
    $('.tl-head').classList.add('one'); $('#toplist').classList.add('one');
  }
  const refSeq = await seqOf(sp, P, B);   // the reference sequence (UniProt; FlyBase for fly)
  if (gone()) return;
  let CQ = B.C0, qSeq = '';   // the clustered construct (a choice of B.choices): its rows, its axis, and the letters along it
  const setQSeq = () => { qSeq = B.cons.size && CQ.qName ? B.seqs.get(CQ.qName) || '' : refSeq;
    if (qSeq && CQ.qLen && qSeq.length !== CQ.qLen) qSeq = (CQ.qName && B.seqs.get(CQ.qName)) || ''; };
  setQSeq();
  const SETS = B.sets.length ? B.sets : null;   // a screen with categories: the Source column shows them
  const gname = (key) => { const r = sp.byKey.get(key); return r ? r.gene : key; };
  const range = (k) => Array.from({ length: k }, (_, i) => i + 1);
  let cut = 10, M = null, ACTIVE = new Set(), NET = null, infoOpen = new Set();
  const predCluster = new Map(), partnerCluster = new Map();
  // A partner folded as several constructs (a receptor's isoforms, its fragments): each one's scores and cluster, so a
  // partner that binds through one isoform and not another shows it. One row per gene; the constructs open under it.
  const isoCache = new Map();
  const partnerIsos = (pt) => {
    if (!B.cons.size) return null;
    if (!isoCache.has(pt.id)) {
      const by = new Map(); for (const p of pt.counted) { if (!by.has(p.pc)) by.set(p.pc, []); by.get(p.pc).push(p); }   // counted runs only: the same sequence under another name is a repeat, not an isoform
      isoCache.set(pt.id, by.size < 2 ? null : [...by].map(([pc, ps]) => { const c = B.cons.get(pc), il = ps.map((p) => p.iLIS || 0), ip = ps.map((p) => p.ipTM || 0);
        return { pc, label: c ? c.label : pc, kind: c ? c.kind : '', preds: [...ps].sort((a, b) => b.iLIS - a.iLIS), best: Math.max(...il), avg: mean(il), iptmBest: Math.max(...ip), iptmAvg: mean(ip),
          contacts: Math.max(...ps.map((p) => p.qcLIR || 0)), pass: ps.filter((p) => p.iLIS >= CUT[10]).length, n: ps.length }; }).sort((a, b) => b.best - a.best));
    }
    return isoCache.get(pt.id);
  };
  const isoCluster = (x) => { for (const p of x.preds) { const c = predCluster.get(p.label + '|' + p.rank); if (c != null) return c; } return 0; };
  const isoWord = (xs) => (xs.every((x) => x.kind === 'isoform' || x.kind === 'gene') ? 'isoforms' : 'constructs');
  const isoTip = (p, xs) => `${xs.length} ${isoWord(xs)} of ${gname(p.id)} were folded; best iLIS: ${xs.map((x) => `${x.label} ${x.best.toFixed(2)}`).join(' · ')}`;
  const isoOpen = new Set();
  const clustered = () => !!(M && M.fingerprints.length >= 2);
  const allOn = () => !M || ACTIVE.size === M.k;
  const S = { state: 'loading', map: null, mapOK: false };              // AlphaFold DB model of this protein
  const V = { mode: 'cluster', thr: 0.5, minC: 5, amCut: 0, shown: false, want: false };
  const toStruct = (r) => (S.map ? S.map[r - 1] : r);
  const who = (label) => B.labels.get(label) || { key: label, run: null };    // a cLIP partner label → partner key + screen run
  const xtWant = () => { const v = +((app.querySelector('.xticks') || {}).value); return v >= 2 ? Math.min(40, v) : 0; };
  { const byDs = new Map(); for (const p of B.preds) { if (!byDs.has(p.di)) byDs.set(p.di, new Set()); byDs.get(p.di).add(p.qLen); }   // the construct each screen folded
    const lens = new Set(B.preds.map((p) => p.qLen)), ref = sp.manifest.keyedBy ? 'FlyBase reference' : 'UniProt';
    const f = !B.cons.size && lens.size > 1 ? `constructs differ between screens: ${[...byDs].map(([di, ls]) => `${sp.dsShort[di]} ${[...ls].map(fmtInt).join(' / ')} aa`).join(' · ')}${P.len ? `; UniProt ${fmtInt(P.len)} aa` : ''}`
      : !B.cons.size && P.len && CQ.qLen !== P.len ? `clustered construct ${fmtInt(CQ.qLen)} aa, ${ref} ${fmtInt(P.len)} aa; residue numbers follow the construct` : '';
    if (f) $('#flags').insertAdjacentHTML('beforeend', `<span class="flag">${esc(f)}</span>`); }
  const constructNote = () => {   // what the clustering uses, and what it leaves out
    const a = $('#clip-aside');
    if (B.cons.size) {   // a gene-keyed screen: its other constructs by kind
      const KIND = { phosphosite: ['phosphosite window', 'phosphosite windows'], fragment: ['fragment', 'fragments'], mutant: ['point mutant', 'point mutants'],
        variant: ['variant', 'variants'], peptide: ['peptide', 'peptides'], isoform: ['isoform', 'isoforms'], gene: ['construct', 'constructs'] };
      const pickable = new Set((B.choices || []).map((c) => c.id)), say = (list) => { const byKind = new Map();
        for (const x of list) { const k = x.con ? x.con.kind : 'construct'; if (!byKind.has(k)) byKind.set(k, []); byKind.get(k).push(x); }
        return [...byKind].map(([k, xs]) => { const [one, many] = KIND[k] || ['construct', 'constructs'], ex = xs.slice(0, 3).map((x) => (x.con ? x.con.label : '')).filter(Boolean);
          return `${fmtInt(xs.length)} ${xs.length === 1 ? one : many}${ex.length ? ` (${ex.join(', ')}${xs.length > 3 ? ', …' : ''})` : ''}`; }).join(', '); };
      const own = (BA.split ? BA.split.others.map((c) => ({ con: BA.cons.get(c.id) || { kind: 'isoform', label: c.label } }))
        : [...CQ.aside.entries()].filter(([k]) => pickable.has(k)).map(([, x]) => x)).sort((x, y) => ((x.con || {}).label || '').localeCompare((y.con || {}).label || ''));
      const rest = [...CQ.aside.entries()].filter(([k]) => !pickable.has(k)).map(([, x]) => x);
      const row = `the ${WORD} row at the top`;
      a.textContent = CQ.qName ? `Clustering uses ${(B.cons.get(CQ.qName) || {}).label || CQ.qName}, in its own numbering.`
          + (ISO ? ` Every card on this page shows its predictions; the ${P.gene} reference and the other ${WORD.toLowerCase()}s are in ${row}.` : '')
        : `Clustering uses the ${fmtInt(CQ.qLen)} aa reference, with every construct placed or mapped on it.`
          + (own.length ? ` ${say(own)} ${own.length === 1 ? 'shares' : 'share'} too little sequence with the reference to be drawn on it; ${row} shows each on its own.` : '')
          + (rest.length ? ` ${say(rest)} ${rest.length === 1 ? 'is' : 'are'} listed with ${rest.length === 1 ? 'its' : 'their'} partners and on the pair pages, not clustered.` : '');
      a.hidden = !(own.length || rest.length || CQ.qName);
    } else if (CQ.aside.size) { a.hidden = false;
      a.textContent = [...CQ.aside.values()].map((x) => `${sp.dsShort[x.di]} predicted ${P.gene} as a ${fmtInt(x.len)} aa construct`).join('; ') + `, so those models are listed but left out of the clustering, which uses the ${fmtInt(CQ.qLen)} aa construct.`;
    } else a.hidden = true;
  };
  constructNote();
  if (ISO) {   // the isoforms side by side: what each was folded with and what binds it best
    const tally = (V) => { const ot = V.partners.filter((x) => x.id !== P.key && !x.rep), n = (c) => ot.filter((x) => x.best >= c).length;
      return { n: ot.length, p10: n(CUT[10]), p5: n(CUT[5]), p1: n(CUT[1]), top: [...ot].sort((x, y) => y.best - x.best).slice(0, 3) }; };
    const fromSummary = (t) => ({ n: t.partners, p10: t.p10, p5: t.p5, p1: t.p1, top: (t.top || []).map(([id, best]) => ({ id, best })) });
    const all = BA.split ? fromSummary(BA.split.allStats) : tally(BA), card = $('#c-iso');
    const stat = new Map(BA.choices.map((c) => [c, c.stats ? fromSummary(c.stats) : tally(BA.only(c.id))]));
    const listed = BA.choices.filter((c) => c === ISO || stat.get(c).n > 1), few = BA.choices.length - listed.length;   // a construct folded with one partner is left out of the table (its pair page still has it)
    card.hidden = listed.length < 2;
    card.innerHTML = `<div class="card-head"><h2>${WORD}s</h2><span class="muted">each folded separately · open one to see it on this page</span></div>
      <div class="tbl-wrap"><table class="sets isotbl"><thead><tr><th>${WORD}</th><th class="n">Models</th><th class="n">Partners</th><th class="n">Past 10% FPR</th><th class="n">5%</th><th class="n">1%</th><th>Top partners (iLIS)</th></tr></thead><tbody>
      ${listed.map((c) => { const t = stat.get(c);
        return `<tr class="${c === ISO ? 'on' : ''}"><td class="set-name"><a href="${isoHref(c)}">${esc(isoName(c))}</a>${c === ISO ? ' <span class="muted">· shown</span>' : ''}</td>
          <td class="n">${fmtInt(c.n)}</td><td class="n">${fmtInt(t.n)}</td><td class="n">${fmtInt(t.p10)}</td><td class="n">${fmtInt(t.p5)}</td><td class="n">${fmtInt(t.p1)}</td>
          <td class="iso-top">${t.top.map((x) => `<a href="#/${sp.id}/${x.id}${scopeQ}">${esc(gname(x.id))}</a> <span class="num">${x.best.toFixed(2)}</span>`).join(' · ')}</td></tr>`; }).join('')}
      </tbody></table></div>
      <p class="muted" style="margin:10px 0 0;font-size:13.5px">All ${fmtInt(BA.choices.length)} together: ${fmtInt(all.n)} partners, ${fmtInt(all.p10)} past the 10% FPR cutoff.${few ? ` ${fmtInt(few)} ${few === 1 ? 'construct' : 'constructs'} folded with a single partner ${few === 1 ? 'is' : 'are'} not listed.` : ''}</p>`;
    if (!card.hidden) { const nav = $('.subnav'), btn = document.createElement('button'); btn.dataset.t = 'c-iso'; btn.textContent = `${WORD}s`; nav.prepend(btn);
      btn.onclick = () => window.scrollTo({ top: card.getBoundingClientRect().top + window.scrollY - 112, behavior: 'auto' }); }
  }
  app.querySelectorAll('.xticks').forEach((i) => { i.oninput = () => { app.querySelectorAll('.xticks').forEach((o) => { if (o !== i) o.value = i.value; }); drawFreq(); drawHeatmap(); }; });

  /* cLIP ─ clustering + everything drawn from it */
  async function cluster() {
    $('#clip-sub').textContent = 'Clustering…';
    let m = null, err = null;
    try { m = await trackLoad(`Clustering ${P.gene}'s partners (${cut}% FPR)`, runClip(CQ.rows, B.qLabel, CUT[cut]), 'task'); } catch (e) { err = e; }
    if (gone()) return;
    M = m; if (err) $('#clip-sub').textContent = `Clustering failed: ${err.message}`;
    predCluster.clear(); partnerCluster.clear();
    if (M) {
      const best = new Map();
      M.preds.forEach((p, i) => { predCluster.set(p.partner + '|' + p.rank, M.labels[i]); const k = who(p.partner).key, b = best.get(k); if (!b || p.iLIS > b.iLIS) best.set(k, { iLIS: p.iLIS, c: M.labels[i] }); });
      for (const [k, v] of best) partnerCluster.set(k, v.c);
      ACTIVE = new Set(range(M.k));
    }
    const n = M ? M.fingerprints.length : 0, k = clustered() ? M.k : 0;
    if (M) $('#clip-sub').textContent = clustered() ? `${fmtInt(n)} predictions of ${fmtInt(partnerCluster.size)} partners past iLIS ${CUT[cut]} (${cut}% FPR) · query ${fmtInt(M.plen)} aa · cosine distance, average linkage, silhouette`
      : `${n ? 'Only one prediction' : 'No predictions'} past the ${cut}% FPR cutoff, so there is nothing to cluster.`;
    const want = !clustered() && V.mode === 'cluster' ? 'plddt' : clustered() && V.auto && S.mapOK ? 'cluster' : null;   // no clusters: show pLDDT until there are
    if (want) { V.auto = want === 'plddt'; V.mode = want; app.querySelectorAll('#cmode button').forEach((b) => b.classList.toggle('on', b.dataset.m === want)); }
    renderChips(); renderSites(); drawFreq(); drawHeatmap(); renderClusterInfo(); recolor3D(); drawScatter(); drawTopList(); drawTable(); fillPartners(); if (NET) NET.recolor();
  }
  function renderChips() {
    app.querySelectorAll('[data-chips]').forEach((box) => {
      if (!clustered()) { box.innerHTML = ''; return; }
      const sizes = {}; for (const l of M.labels) sizes[l] = (sizes[l] || 0) + 1;
      box.innerHTML = `<button class="cchip all" data-c="all">All clusters (${fmtInt(M.labels.length)})</button>` + range(M.k).map((c) =>
        `<button class="cchip" data-c="${c}"><i style="background:${clusterColor(c, M.k)}"></i>${clusterLabel(c, true)} (${sizes[c] || 0})</button>`).join('');
      box.onclick = (e) => { const b = e.target.closest('.cchip'); if (b) toggleCluster(b.dataset.c); };
    });
    paintChips();
  }
  function toggleCluster(c) {   // as clip.html: from "all", a click isolates that cluster; further clicks add or remove
    if (c === 'all') ACTIVE = new Set(range(M.k));
    else { const id = +c; if (allOn()) ACTIVE = new Set([id]); else if (ACTIVE.has(id)) { ACTIVE.delete(id); if (!ACTIVE.size) ACTIVE = new Set(range(M.k)); } else ACTIVE.add(id); }
    paintChips(); renderSites(); drawFreq(); drawHeatmap(); recolor3D();
  }
  function paintChips() {
    const all = allOn();
    app.querySelectorAll('[data-chips] .cchip').forEach((b) => { const c = b.dataset.c;
      if (c === 'all') b.classList.toggle('sel', all); else { b.classList.toggle('sel', !all && ACTIVE.has(+c)); b.classList.toggle('off', !all && !ACTIVE.has(+c)); } });
  }
  function emptyPlot(host, msg) { host.innerHTML = `<div class="plot-empty">${msg}</div>`; }
  function qDomains(L) {   // UniProt domains sit on the current UniProt sequence: drawn only when the predicted construct is that sequence's length
    if (!S.domains || !S.domains.length || !(P.len && P.len === L)) return [];
    if (sp.manifest.keyedBy && !(qSeq && S.uniSeq === qSeq)) return [];   // fly: only when the FlyBase reference is UniProt's sequence
    return S.domains.map((d) => ({ name: d.name, start: d.start, end: d.end, s: d.start, e: d.end })).filter((d) => d.e >= 1 && d.s <= L).sort((a, b) => a.s - b.s);
  }

  function drawFreq() {
    const host = $('#freq-wrap'); if (!host) return;
    if (!clustered()) { emptyPlot(host, M ? 'Nothing to cluster at this cutoff.' : 'Clustering…'); $('#hot').innerHTML = ''; $('#freq-domains').innerHTML = ''; return; }
    if (!$('#freq', host)) host.innerHTML = '<canvas id="freq"></canvas>';
    const cv = $('#freq', host), L = M.plen, n = M.fingerprints.length, k = M.k;
    const tot = new Uint16Array(L + 2), byC = new Array(L + 2);
    for (let i = 0; i < n; i++) { const lab = M.labels[i]; if (!ACTIVE.has(lab)) continue;
      for (const r of M.fingerprints[i]) { if (r < 1 || r > L) continue; tot[r]++; (byC[r] ||= {})[lab] = (byC[r][lab] || 0) + 1; } }
    let max = 1; for (let r = 1; r <= L; r++) if (tot[r] > max) max = tot[r];
    const W = host.clientWidth, bw = (W - AXL - AXR) / L, xOf = (r) => AXL + (r - 1) * bw, xc = (r) => xOf(r) + bw / 2;
    const doms = qDomains(L); doms.forEach((d, i) => { d.idx = i + 1; });
    const nL = doms.length ? lanes(doms, (d) => xOf(Math.max(1, d.s)), (d) => Math.max(xOf(Math.max(1, d.s)) + 2, xOf(Math.min(L, d.e) + 1))) : 0;
    const DRH = 15, padT = 4 + nL * DRH + (nL ? 8 : 10), fH = 124, fBase = padT + fH;
    const pl = !!(S.plddt && S.mapOK), am = !!(S.am && S.mapOK);
    const pTop = fBase + 18, pH = 46, pBase = pl ? pTop + pH : fBase, aTop = pBase + 18, aH = 42, aBase = am ? aTop + aH : pBase, H = aBase + 26;
    const g = canvasCtx(cv, W, H);
    const yLabel = (text, y0, y1, color) => { g.save(); g.translate(13, (y0 + y1) / 2); g.rotate(-Math.PI / 2); g.textAlign = 'center'; g.textBaseline = 'middle'; g.fillStyle = color; g.font = '11px "IBM Plex Sans", system-ui, sans-serif'; g.fillText(text, 0, 0); g.restore(); };
    const yTicks = (vals, fmt, yOf) => { g.fillStyle = '#5A697C'; g.font = '10.5px "IBM Plex Mono", ui-monospace, monospace'; g.textAlign = 'right'; g.textBaseline = 'middle'; for (const v of vals) g.fillText(fmt(v), AXL - 6, yOf(v)); };
    const frame = (y0, y1) => { g.fillStyle = '#F6F8FB'; g.fillRect(AXL, y0, W - AXL - AXR, y1 - y0); g.strokeStyle = '#D5DDE6'; g.lineWidth = 1; g.beginPath(); g.moveTo(AXL + 0.5, y0); g.lineTo(AXL + 0.5, y1 + 0.5); g.lineTo(W - AXR, y1 + 0.5); g.stroke(); };
    for (const d of doms) { const x0 = xOf(Math.max(1, d.s)), x1 = Math.max(x0 + 2, xOf(Math.min(L, d.e) + 1)), y = 4 + d.lane * DRH;   // domains D1, D2 … (names in the legend), as LIVIA cLIP
      g.fillStyle = '#E3E9F1'; g.fillRect(x0, y, x1 - x0, 12); g.strokeStyle = '#9FB0C4'; g.lineWidth = 0.6; g.strokeRect(x0 + 0.3, y + 0.3, x1 - x0 - 0.6, 11.4);
      g.fillStyle = '#34445A'; g.font = '10px "IBM Plex Sans", system-ui, sans-serif'; g.textAlign = 'center'; g.textBaseline = 'middle'; const lb = 'D' + d.idx; if (g.measureText(lb).width < x1 - x0 - 2) g.fillText(lb, (x0 + x1) / 2, y + 6.5); }
    frame(padT, fBase);   // # predictions contacting the residue; bar color = the residue's most frequent cluster
    g.strokeStyle = '#E3E8EE'; g.lineWidth = 1; for (const f of [0.5, 1]) { const y = Math.round(fBase - f * fH) + 0.5; g.beginPath(); g.moveTo(AXL, y); g.lineTo(W - AXR, y); g.stroke(); }
    for (let r = 1; r <= L; r++) { if (!tot[r]) continue; let dom = 1, b = -1; const cc = byC[r]; for (const c in cc) if (cc[c] > b) { b = cc[c]; dom = +c; }
      const h = tot[r] / max * fH; g.fillStyle = clusterColor(dom, k); g.fillRect(xOf(r), fBase - h, Math.max(1, bw), h); }
    yTicks([...new Set([0, Math.round(max / 2), max])], String, (v) => fBase - v / max * fH); yLabel('# predictions', padT, fBase, '#34445A');
    const line = (vals, yOf, col) => { let prev = null; g.lineWidth = 1.5;
      for (let r = 1; r <= L; r++) { const sr = toStruct(r), v = sr ? vals(sr) : null; if (v == null) { prev = null; continue; } const x = xc(r), y = yOf(v);
        if (prev) { const mx = (prev.x + x) / 2, my = (prev.y + y) / 2; g.strokeStyle = col(prev.v); g.beginPath(); g.moveTo(prev.x, prev.y); g.lineTo(mx, my); g.stroke(); g.strokeStyle = col(v); g.beginPath(); g.moveTo(mx, my); g.lineTo(x, y); g.stroke(); }
        prev = { x, y, v }; } };
    if (pl) { frame(pTop, pBase); const yOf = (v) => pBase - v / 100 * pH; line((sr) => S.plddt.get('A:' + sr), yOf, plddtCol); yTicks([0, 50, 100], String, yOf); yLabel('pLDDT', pTop, pBase, '#2F6FA8'); }
    if (am) { frame(aTop, aBase); const yOf = (v) => aBase - v * aH; line((sr) => S.am[sr], yOf, amCol); yTicks([0, 0.5, 1], (v) => v.toFixed(1), yOf); yLabel('AM path.', aTop, aBase, '#A33'); }
    drawTicks(g, resTicks(L, W - AXL - AXR, xtWant()), aBase + 2, xc, W);
    cv.onmousemove = (e) => { const b = cv.getBoundingClientRect(), r = Math.floor((e.clientX - b.left - AXL) / bw) + 1; if (r < 1 || r > L) return hideTip();
      const my = e.clientY - b.top, hd = my < padT - 4 && doms.find((d) => r >= d.s && r <= d.e && my >= 4 + d.lane * DRH && my <= 16 + d.lane * DRH);
      if (hd) return showTip(domTip(hd), e.clientX, e.clientY);   // a domain box: its name and span
      const cc = byC[r] || {}, parts = Object.keys(cc).map(Number).sort((a, z) => cc[z] - cc[a]).map((c) => `<span style="color:${clusterColor(c, k)}">●</span> ${clusterLabel(c, true)} ${cc[c]}`).join(' · ');
      const sr = toStruct(r), plv = pl && sr ? S.plddt.get('A:' + sr) : null, amv = am && sr ? S.am[sr] : null, dom = doms.filter((d) => r >= d.s && r <= d.e).map((d) => d.name).join(', ');
      showTip(`<b>${qSeq[r - 1] || ''}${r}</b> · ${tot[r]} prediction${tot[r] === 1 ? '' : 's'}${parts ? '<br>' + parts : ''}${dom ? `<br>${esc(dom)}` : ''}${plv != null ? `<br>pLDDT ${plv.toFixed(0)}` : ''}${amv != null ? `${plv != null ? ' · ' : '<br>'}AM ${amv.toFixed(2)}` : ''}`, e.clientX, e.clientY); };
    cv.onmouseleave = hideTip;
    $('#freq-domains').innerHTML = doms.length ? `<b>Domains</b> ${doms.map((d) => `<span><b>D${d.idx}</b> ${esc(d.name)} <span class="muted">(${d.start}–${d.end})</span></span>`).join('')}${P.acc ? uniprotLink(P.acc, `UniProt ${esc(P.acc)} ↗`) : ''}` : '';
    const hot = range(L).filter((r) => tot[r]).sort((a, b) => tot[b] - tot[a]).slice(0, 12);
    $('#hot').innerHTML = hot.length ? `<span class="hot-lbl">Most contacted</span>${hot.map((r) => `<span title="${tot[r]} predictions">${qSeq[r - 1] || ''}${r} · ${tot[r]}</span>`).join('')}` : '';
    attachExport('freq', `atlas_${P.gene}_frequency`, drawFreq);
  }

  function fpHit() {   // the partner typed in the fingerprint's Find partner box: exact name, then prefix, then substring
    const f = ($('#fp-find') || {}).value; const q = (f || '').trim().toLowerCase(); if (!q) return null;
    const gl = (p) => gname(p.id).toLowerCase();
    return B.partners.find((p) => gl(p) === q) || B.partners.find((p) => gl(p).startsWith(q)) || B.partners.find((p) => gl(p).includes(q)) || false;
  }
  function drawHeatmap() {
    const host = $('#fp-wrap'); if (!host) return;
    if (!clustered()) { emptyPlot(host, M ? 'Nothing to cluster at this cutoff.' : 'Clustering…'); return; }
    if (!$('#heatmap', host)) host.innerHTML = '<canvas id="heatmap"></canvas>';
    const cv = $('#heatmap', host), L = M.plen, k = M.k, all = allOn();
    const order = all ? M.order : M.order.filter((i) => ACTIVE.has(M.labels[i])), n = order.length || 1;
    const W = host.clientWidth, rowsH = Math.min(300, Math.max(110, n * 6)), H = rowsH + 26;   // the y axis is capped: many partners share the same height
    const g = canvasCtx(cv, W, H), dpr = g.__dpr, bw = (W - AXL - AXR) / L, xOf = (r) => AXL + (r - 1) * bw;
    g.setTransform(1, 0, 0, 1, 0, 0);                                           // rows in device pixels: integer edges, no seams
    const X = (x) => Math.round(x * dpr), rowsHd = Math.round(rowsH * dpr), rowTop = (row) => Math.round(row * rowsHd / n);
    const sx0 = X(34), sx1 = X(AXL - 8), rx0 = X(AXL), rx1 = X(W - AXR);
    for (let row = 0; row < n; row++) { const i = order[row]; if (i == null) continue; const y0 = rowTop(row), hh = Math.max(1, rowTop(row + 1) - y0);
      g.fillStyle = clusterColor(M.labels[i], k); g.fillRect(sx0, y0, sx1 - sx0, hh);
      g.fillStyle = '#F7FBFF'; g.fillRect(rx0, y0, rx1 - rx0, hh);
      g.fillStyle = '#08306B'; for (const [s, e] of runs(M.fingerprints[i])) { const x0 = X(xOf(s)), x1 = Math.max(x0 + 1, X(xOf(e + 1))); g.fillRect(x0, y0, x1 - x0, hh); } }
    const hit = fpHit(), found = $('#fp-found'), mine = hit ? M.order.filter((i) => who(M.preds[i].partner).key === hit.id) : [];
    $('#fp-find').classList.toggle('nf', hit === false); found.hidden = !hit;
    if (hit) {   // its rows: contacts in orange on a light band, at least 3 px tall however many rows share the height, and a marker at the right
      const on = new Set(order), cl = [...new Set(mine.map((i) => M.labels[i]))].sort((a, b) => a - b), mh = Math.round(3 * dpr), col = '#E67E22';
      order.forEach((i, row) => { if (!mine.includes(i)) return; const yc = (rowTop(row) + rowTop(row + 1)) / 2, y0 = Math.round(yc - mh / 2);
        g.fillStyle = 'rgba(230,126,34,0.2)'; g.fillRect(rx0, y0, rx1 - rx0, mh);
        g.fillStyle = col; for (const [s, e] of runs(M.fingerprints[i])) { const x0 = X(xOf(s)), x1 = Math.max(x0 + 1, X(xOf(e + 1))); g.fillRect(x0, y0, x1 - x0, mh); }
        g.beginPath(); g.moveTo(rx1 + X(2), yc); g.lineTo(rx1 + X(11), yc - X(5)); g.lineTo(rx1 + X(11), yc + X(5)); g.closePath(); g.fill(); });
      const shownN = mine.filter((i) => on.has(i)).length;
      found.innerHTML = !mine.length ? `<b>${esc(gname(hit.id))}</b> has no prediction past the ${cut}% FPR cutoff.`
        : `<b style="color:${col}">${esc(gname(hit.id))}</b>: ${mine.length} prediction${mine.length === 1 ? '' : 's'}, in ${cl.map((c) => `<span style="color:${clusterColor(c, k)}">●</span> ${clusterLabel(c)}`).join(', ')}`
          + (shownN < mine.length ? ` · ${mine.length - shownN} hidden by the cluster choice above` : '') + ' · marked in orange';
    }
    const Z = M.Z || [];
    if (all && Z.length === n - 1 && n >= 2) {                                  // row dendrogram: leaves at the strip, root at the left; x ∝ merge distance
      const yv = new Float64Array(2 * n), xv = new Float64Array(2 * n), d0 = X(3), d1 = X(31), maxD = Z[Z.length - 1][2] || 1;
      for (let row = 0; row < n; row++) { yv[order[row]] = (rowTop(row) + rowTop(row + 1)) / 2; xv[order[row]] = d1; }
      for (let m = 0; m < Z.length; m++) { const a = Z[m][0], b = Z[m][1]; yv[n + m] = (yv[a] + yv[b]) / 2; xv[n + m] = d0 + (d1 - d0) * (1 - Z[m][2] / maxD); }
      g.strokeStyle = '#5B6B7F'; g.lineWidth = Math.max(0.6, dpr * 0.6); g.beginPath();
      for (let m = 0; m < Z.length; m++) { const a = Z[m][0], b = Z[m][1], xm = xv[n + m];
        g.moveTo(xv[a], yv[a]); g.lineTo(xm, yv[a]); g.moveTo(xv[b], yv[b]); g.lineTo(xm, yv[b]); g.moveTo(xm, yv[a]); g.lineTo(xm, yv[b]); }
      g.stroke();
    }
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.strokeStyle = '#D5DDE6'; g.lineWidth = 1; g.strokeRect(AXL + 0.5, 0.5, W - AXL - AXR - 1, rowsH - 1);
    drawTicks(g, resTicks(L, W - AXL - AXR, xtWant()), rowsH + 2, (r) => xOf(r) + bw / 2, W);
    const at = (e) => { const b = cv.getBoundingClientRect(), x = e.clientX - b.left, y = e.clientY - b.top; if (y < 0 || y > rowsH || x < 30) return null;
      const i = order[Math.min(n - 1, Math.floor(y / rowsH * n))]; return i == null ? null : { i, r: Math.floor((x - AXL) / bw) + 1 }; };
    cv.onmousemove = (e) => { const h = at(e); if (!h) return hideTip(); const p = M.preds[h.i], c = M.labels[h.i], w = who(p.partner);
      const res = h.r >= 1 && h.r <= L ? `<br>residue ${qSeq[h.r - 1] || ''}${h.r}${M.fingerprints[h.i].includes(h.r) ? ' · contact' : ''}` : '';
      showTip(`<b>${esc(gname(w.key))}</b> · ${esc(runLabel(sp, B, w.run, P))} · rank ${p.rank} · iLIS ${p.iLIS.toFixed(3)} · <span style="color:${clusterColor(c, k)}">●</span> ${clusterLabel(c)}${res}`, e.clientX, e.clientY); };
    cv.onmouseleave = hideTip;
    cv.onclick = (e) => { const h = at(e); if (!h) return; hideTip(); location.hash = `#/${sp.id}/${P.key}/${who(M.preds[h.i].partner).key}${scopeQ}`; };
    attachExport('heatmap', `atlas_${P.gene}_fingerprint`, drawHeatmap);
  }

  /* Binding sites: each cLIP cluster's footprint on the query, the page's answer to "where do the partners bind" */
  // R4: how often this protein's own predictions below the 10% cutoff (all models of all partners) contact each residue:
  // surfaces the model places any partner on. Descriptive only; hidden with fewer than 50 such predictions.
  const BG_MIN = 50;
  const bgOf = () => { if (!CQ || !CQ.bg || !M) return null; if (CQ._bg && CQ._bg.L === M.plen) return CQ._bg;
    const L = M.plen, n = CQ.bgN, f = new Float32Array(L + 1); let sum = 0, max = 0;
    for (let r = 1; r <= L; r++) { f[r] = n ? (CQ.bg[r] || 0) / n : 0; sum += f[r]; if (f[r] > max) max = f[r]; }
    return (CQ._bg = { L, n, f, mean: sum / L, max, ok: n >= BG_MIN }); };
  const bgRatio = (st) => { const b = bgOf(); if (!b || !b.ok || !st.foot.length || !b.mean) return null;
    return st.foot.reduce((a, r) => a + (b.f[r] || 0), 0) / st.foot.length / b.mean; };
  // The two contact rates over a site's footprint: its own predictions, and the protein's predictions below the cutoff
  // (both as the mean over the footprint's residues). Shown side by side, never as a ratio: the footprint is chosen by
  // the first rate, and a ratio over a near-empty background only grows.
  const bgRates = (st) => { const b = bgOf(); if (!b || !b.ok || !st.foot.length) return null;
    return { own: st.foot.reduce((a, r) => a + (st.hits.get(r) || 0) / st.n, 0) / st.foot.length, bg: st.foot.reduce((a, r) => a + (b.f[r] || 0), 0) / st.foot.length }; };
  const SITE_FRAC = 0.3;   // a site's footprint: residues contacted by at least 30% of its predictions (the cluster footprint of the paper)
  function sitesOf() {     // → sites, most partners first: cluster, predictions, partners by best iLIS, contacts per residue, footprint
    const by = new Map();
    M.preds.forEach((p, i) => { const c = M.labels[i]; let t = by.get(c); if (!t) by.set(c, t = { c, n: 0, hits: new Map(), best: new Map() });
      t.n++; for (const r of M.fingerprints[i]) t.hits.set(r, (t.hits.get(r) || 0) + 1);
      const k = who(p.partner).key; if (k !== P.key) t.best.set(k, Math.max(t.best.get(k) || 0, p.iLIS || 0)); });
    return [...by.values()].map((t) => { const foot = [...t.hits].filter(([, h]) => h / t.n >= SITE_FRAC).map(([r]) => r).sort((a, b) => a - b);
      return { ...t, partners: [...t.best].map(([key, best]) => ({ key, best })).sort((a, b) => b.best - a.best), foot, ranges: runs(foot) }; })
      .sort((a, b) => b.partners.length - a.partners.length || b.n - a.n || a.c - b.c);
  }
  const rangeText = (rs) => rs.map(([a, b]) => (a === b ? `${a}` : `${a}–${b}`)).join(', ');
  const stretches = (rs) => { const out = []; for (const [a, b] of rs) { const l = out[out.length - 1]; if (l && a - l[1] <= 4) l[1] = b; else out.push([a, b]); } return out; };   // gaps of up to 3 residues joined
  const major = (st, np) => st.partners.length >= Math.max(3, Math.ceil(0.02 * np));   // a major site: at least 3 partners and 2% of them
  let sitesAll = false;   // the map shows the major sites until asked for all
  const siteDomains = (st, doms) => doms.filter((d) => st.foot.some((r) => r >= d.s && r <= d.e)).map((d) => d.name).slice(0, 2);
  // Partners past the cutoff (the tiles) that no site holds, and why: cLIP takes a pair only when its rank-1 model passes,
  // and drops a model without contact residues and constructs named phospho- or mutant (its default filter)
  function outOfSites() {
    const miss = B.partners.filter((x) => x.id !== P.key && !x.rep && x.best >= CUT[cut] && !partnerCluster.has(x.id));
    if (!miss.length) return '';
    const r1 = (x) => x.counted.filter((p) => p.rank === 1 && p.iLIS >= CUT[cut]);
    const low = miss.filter((x) => !r1(x).length).length, bare = miss.filter((x) => r1(x).length && r1(x).every((p) => !p.qcLIR)).length, other = miss.length - low - bare;
    const n1 = miss.length === 1, why = [[low, `pass${n1 ? 'es' : ''} only in lower-ranked models (a site takes a pair when its rank-1 model passes)`],
      [bare, `ha${n1 ? 's' : 've'} no contact residues`], [other, `${n1 ? 'is a construct' : 'are constructs'} that cLIP leaves out by name (phospho-, mutant)`]].filter(([k]) => k);
    const head = `${fmtInt(miss.length)} more partner${n1 ? '' : 's'} past the cutoff ${n1 ? 'is' : 'are'} in no site`;
    return ` <span class="muted">${why.length === 1 ? `${head}: ${n1 ? 'it' : 'they'} ${why[0][1]}` : `${head}: ${why.map(([k, t]) => `${fmtInt(k)} ${t}`).join('; ')}`}.</span>`;
  }
  function renderSites() {
    const ans = $('#sites-answer'), map = $('#sites-map'); if (!ans) return;
    if (!clustered()) { ans.innerHTML = M ? `No binding site to show: ${M.fingerprints.length ? 'only one prediction is' : 'no prediction is'} past the ${cut}% FPR cutoff.` : 'Finding the binding sites…';
      map.innerHTML = ''; $('#sites-domains').innerHTML = ''; $('#sites-more').hidden = true; return; }
    const sites = sitesOf(), doms = qDomains(M.plen), K = sites.length, np = partnerCluster.size, big = sites[0], lanesOrder = [...sites].sort((a, b) => a.c - b.c);
    const names = (ps) => ps.map((x) => `<a href="#/${sp.id}/${P.key}/${x.key}${scopeQ}">${esc(gname(x.key))}</a>`);
    const andList = (xs) => (xs.length < 3 ? xs.join(' and ') : `${xs.slice(0, -1).join(', ')} and ${xs[xs.length - 1]}`);
    // the partners named as examples: the largest share of models past the cutoff first, then the average iLIS, so a pair
    // past it in one model of three is not named before one past it in all (3 of 3 counts as 5 of 5)
    const PB = new Map(B.partners.map((x) => [x.id, x])), passN = (t) => (t ? t.counted.filter((p) => p.iLIS >= CUT[cut]).length / Math.max(1, t.counted.length) : 0);
    const steady = (ps) => ps.map((x) => ({ ...x, pass: passN(PB.get(x.key)), avg: (PB.get(x.key) || {}).avg || 0 })).sort((a, b) => b.pass - a.pass || b.avg - a.avg || b.best - a.best);
    const where = (st) => { if (!st.foot.length) return 'no residue shared by 30% of its predictions';
      const r = stretches(st.ranges), dn = siteDomains(st, doms);
      return `${r.length <= 3 ? `residues ${rangeText(r)}` : `${r.length} stretches between residues ${r[0][0]} and ${r[r.length - 1][1]}`}${dn.length ? ` (${esc(dn.join(', '))})` : ''}`; };
    ans.innerHTML = `<b>${esc(P.gene)}</b>: ${fmtInt(np)} partner${np === 1 ? ' is' : 's are'} predicted, past the ${cut}% FPR cutoff, to contact ${K === 1 ? 'one site' : `${K} sites`}.`
      + ` The largest, at ${where(big)}, is contacted by ${fmtInt(big.partners.length)} of them${big.partners.length ? `, including ${andList(names(steady(big.partners).slice(0, 3)))}` : ''}.`
      + outOfSites();
    { const [pos] = [...big.hits].sort((a, b) => b[1] - a[1] || a[0] - b[0])[0] || [], aa = pos && qSeq ? qSeq[pos - 1] : '';   // the most contacted residue of the largest site
      if (pos) { const ex = aa ? `${pos} or ${aa}${pos}${aa === 'A' ? 'G' : 'A'}` : `${pos}`, box = $('#res-q'); box.placeholder = ex;
        box.closest('label').title = `a residue number or a variant (${ex.replace(' or ', ', ')}): which predictions, partners and sites contact it; the link keeps it`;
        const eg = $('#res-eg'); if (eg && aa) eg.textContent = ex.split(' or ')[1]; } }
    const nMajor = sites.filter((st) => major(st, np)).length, collapse = !sitesAll && nMajor > 0 && nMajor < K;
    const shown = collapse ? lanesOrder.filter((st) => major(st, np) || (!allOn() && ACTIVE.has(st.c))) : lanesOrder;   // a site chosen elsewhere stays in view
    drawSitesMap(map, shown, doms); renderResLook(sites); attachExport('sites-canvas', `atlas_${P.gene}_sites`, () => drawSitesMap(map, shown, doms));
    { const b = bgOf(), note = $('#sites-bg'); if (note) note.innerHTML = !b ? '' : b.ok
      ? `<b>Background</b>, in gray: how often ${esc(P.gene)}'s ${fmtInt(b.n)} predictions below the cutoff contact each residue, scaled to its maximum; surfaces the model places any partner on. A site on a high background may be such a surface rather than a specific interface.`
      : `No background: too few predictions below the cutoff (${fmtInt(b.n)}; at least ${BG_MIN} are needed).`; }
    const more = $('#sites-more'); more.hidden = !(nMajor > 0 && nMajor < K);
    if (!more.hidden) { more.innerHTML = sitesAll ? `All ${K} sites. <button class="more" type="button">Show the ${nMajor} major sites</button>`
      : `The ${nMajor} major sites, each with at least 3 partners and 2% of all partners; ${K - nMajor} minor site${K - nMajor === 1 ? '' : 's'} hidden. <button class="more" type="button">Show all ${K} sites</button>`;
      $('button', more).onclick = () => { sitesAll = !sitesAll; renderSites(); }; }
    $('#sites-domains').innerHTML = doms.length ? doms.map((d, i) => `<span><b>D${i + 1}</b> ${esc(d.name)} (${d.s}–${d.e})</span>`).join('') : '';
  }
  function renderResLook(sites) {   // the residue in the link: its contacts past the cutoff, by site and partner
    const box = $('#res-look'); if (!box) return; const r = RESQ, L = M.plen;
    box.hidden = !r; if (!r) return;
    if (r > L) { box.innerHTML = `Residue ${r} is past the end of the sequence shown here (${fmtInt(L)} residues).`; return; }
    const hit = M.preds.map((p, i) => i).filter((i) => M.fingerprints[i].includes(r)), best = new Map();
    for (const i of hit) { const k = who(M.preds[i].partner).key; if (k !== P.key) best.set(k, Math.max(best.get(k) || 0, M.preds[i].iLIS || 0)); }
    const inS = sites.map((st) => ({ st, h: st.hits.get(r) || 0 })).filter((x) => x.h).sort((a, b) => b.h - a.h);
    const top = [...best].sort((a, b) => b[1] - a[1]), shown = top.slice(0, 12), aa = qSeq && qSeq[r - 1] ? qSeq[r - 1] : '';
    const off = RESV.wt && aa && RESV.wt !== aa ? `<span class="res-warn">This sequence has ${esc(aa)} at ${r}, not ${esc(RESV.wt)}; check the numbering (isoform or construct).</span> ` : '';
    box.innerHTML = `<b>${esc(RESV.wt && !aa ? RESV.wt : aa)}${r}${RESV.mut ? esc(RESV.mut) : ''}</b>: ` + off + (!hit.length ? `no prediction past the ${cut}% FPR cutoff contacts it.`
      : `contacted in ${fmtInt(hit.length)} of ${fmtInt(M.preds.length)} predictions past the ${cut}% FPR cutoff, by ${fmtInt(best.size)} partner${best.size === 1 ? '' : 's'}`
        + `, in ${inS.map((x) => `<span style="color:${clusterColor(x.st.c, M.k)}">●</span> ${clusterLabel(x.st.c)} (${x.h} of ${x.st.n}${x.st.foot.includes(r) ? ', in its footprint' : ''})`).join(', ')}.`
        + ` Partners: ${shown.map(([k]) => `<a href="#/${sp.id}/${P.key}/${k}${scopeQ}">${esc(gname(k))}</a>`).join(', ')}${top.length > shown.length ? ` and ${top.length - shown.length} more` : ''}.`)
      + ` <button class="more" type="button" id="res-copy">Copy link</button>`
      + (RESV.mut ? `<span class="res-note">A contact in a prediction says where partners bind, not what the substitution does to them; the page does not predict its effect.</span>` : '');
    $('#res-copy').onclick = (e) => { navigator.clipboard && navigator.clipboard.writeText(location.href).then(() => { e.target.textContent = 'Copied'; }, () => {}); };
    if (!resScrolled) { resScrolled = true; $('#c-sites').scrollIntoView({ block: 'start' }); }
  }
  function setRes(v) {   // keep the residue in the link without reloading the page
    RESV = parseRes(v); RESQ = RESV.n; resScrolled = true;
    const [path, qs] = location.hash.split('?'), u = new URLSearchParams(qs || ''); RESQ ? u.set('res', resLabel(RESV)) : u.delete('res');
    history.replaceState(null, '', `${path}${u.toString() ? '?' + u : ''}`);
    if (clustered()) renderSites();
  }
  function drawSitesMap(host, sites, doms) {   // one lane per site on the query's residues: solid = footprint, shade = contact frequency within the site
    if (!$('#sites-canvas', host)) host.innerHTML = '<canvas id="sites-canvas"></canvas>';
    const cv = $('#sites-canvas', host), L = M.plen, W = host.clientWidth, bw = (W - AXL - AXR) / L, xOf = (r) => AXL + (r - 1) * bw;
    doms.forEach((d, i) => { d.idx = i + 1; });
    const nL = doms.length ? lanes(doms, (d) => xOf(Math.max(1, d.s)), (d) => Math.max(xOf(Math.max(1, d.s)) + 2, xOf(Math.min(L, d.e) + 1))) : 0;
    const bgd = bgOf(), BGH = bgd && bgd.ok ? 24 : 0;   // the gray background lane under the sites
    const DRH = 15, LH = sites.length > 14 ? 14 : 20, top = 4 + nL * DRH + (nL ? 6 : 0), H = top + sites.length * LH + BGH + 26, g = canvasCtx(cv, W, H), one = !allOn();
    for (const d of doms) { const x0 = xOf(Math.max(1, d.s)), x1 = Math.max(x0 + 2, xOf(Math.min(L, d.e) + 1)), y = 4 + d.lane * DRH;
      g.fillStyle = '#E3E9F1'; g.fillRect(x0, y, x1 - x0, 12); g.fillStyle = '#4A596F'; g.font = '10px "IBM Plex Sans", system-ui, sans-serif'; g.textBaseline = 'middle';
      if (x1 - x0 > 18) g.fillText('D' + d.idx, x0 + 3, y + 6.5); }
    sites.forEach((st, j) => { const y = top + j * LH, col = clusterColor(st.c, M.k);
      g.fillStyle = '#F3F6F9'; g.fillRect(AXL, y + 3, W - AXL - AXR, LH - 6);
      if (one && ACTIVE.has(st.c)) { g.strokeStyle = '#E67E22'; g.lineWidth = 1.5; g.strokeRect(AXL - 0.5, y + 2.5, W - AXL - AXR + 1, LH - 5); }
      g.globalAlpha = one && !ACTIVE.has(st.c) ? 0.35 : 1;
      g.fillStyle = '#4A596F'; g.font = '11px "IBM Plex Sans", system-ui, sans-serif'; g.textBaseline = 'middle'; g.textAlign = 'right'; g.fillText(clusterLabel(st.c), AXL - 8, y + LH / 2);
      g.textAlign = 'left';
      const dim = one && !ACTIVE.has(st.c) ? 0.35 : 1;
      for (const [r, h] of st.hits) { if (r < 1 || r > L) continue; const f = h / st.n; g.fillStyle = col; g.globalAlpha = dim * (f >= SITE_FRAC ? 1 : 0.12 + 0.5 * f / SITE_FRAC); g.fillRect(xOf(r), y + 3, Math.max(1, bw), LH - 6); }
      g.globalAlpha = 1; });
    if (BGH) { const y0 = top + sites.length * LH + 3, hh = BGH - 6;
      g.fillStyle = '#F3F6F9'; g.fillRect(AXL, y0, W - AXL - AXR, hh); g.fillStyle = '#8A96A6';
      for (let r = 1; r <= L; r++) { const v = bgd.f[r] / (bgd.max || 1); if (v > 0) g.fillRect(xOf(r), y0 + hh * (1 - v), Math.max(1, bw), hh * v); }
      g.fillStyle = '#4A596F'; g.font = '9.5px "IBM Plex Sans", system-ui, sans-serif'; g.textBaseline = 'middle'; g.textAlign = 'right'; g.fillText('Background', AXL - 6, y0 + hh / 2); g.textAlign = 'left'; }   // smaller: the margin is shared with the other plots
    drawTicks(g, resTicks(L, W - AXL - AXR, xtWant()), top + sites.length * LH + BGH + 2, (r) => xOf(r) + bw / 2, W);
    if (RESQ >= 1 && RESQ <= L) { const x = Math.round(xOf(RESQ) + bw / 2) + 0.5;   // the looked-up residue
      g.strokeStyle = '#17263A'; g.lineWidth = 1.5; g.beginPath(); g.moveTo(x, top); g.lineTo(x, top + sites.length * LH); g.stroke();   // dark: no cluster has this color
      g.fillStyle = '#17263A'; g.beginPath(); g.moveTo(x - 4, top - 5); g.lineTo(x + 4, top - 5); g.lineTo(x, top); g.closePath(); g.fill(); }
    cv.onmousemove = (e) => { const b = cv.getBoundingClientRect(), x = e.clientX - b.left, y = e.clientY - b.top, j = Math.floor((y - top) / LH), r = Math.floor((x - AXL) / bw) + 1;
      const hd = y < top && doms.find((d) => r >= d.s && r <= d.e && y >= 4 + d.lane * DRH && y <= 16 + d.lane * DRH);
      if (hd) { cv.style.cursor = ''; return showTip(domTip(hd), e.clientX, e.clientY); }   // a domain box: its name and span
      if (BGH && j === sites.length && r >= 1 && r <= L) { cv.style.cursor = ''; return showTip(`<b>Background</b> · ${qSeq && qSeq[r - 1] ? qSeq[r - 1] : ''}${r}: contacted in ${(100 * bgd.f[r]).toFixed(1)}% of ${fmtInt(bgd.n)} predictions below the cutoff`, e.clientX, e.clientY); }
      const st = sites[j]; if (!st || r < 1 || r > L) { cv.style.cursor = ''; return hideTip(); } const h = st.hits.get(r) || 0;
      cv.style.cursor = 'pointer';
      showTip(`<b>${clusterLabel(st.c)}</b> · ${fmtInt(st.partners.length)} partner${st.partners.length === 1 ? '' : 's'}, ${fmtInt(st.n)} prediction${st.n === 1 ? '' : 's'}<br>${qSeq && qSeq[r - 1] ? qSeq[r - 1] : ''}${r}: ${h} of ${st.n} contact it (${Math.round(100 * h / st.n)}%)${(() => { const q = bgRatio(st), rt = bgRates(st); return q == null || !rt ? '' : `<br>footprint: background ${q.toFixed(1)}× the protein mean; its partners' models contact it in ${Math.round(100 * rt.own)}% of cases, predictions below the cutoff in ${Math.round(100 * rt.bg)}%`; })()} · click to show only this site`, e.clientX, e.clientY); };
    cv.onmouseleave = () => { cv.style.cursor = ''; hideTip(); };
    cv.onclick = (e) => { const b = cv.getBoundingClientRect(), j = Math.floor((e.clientY - b.top - top) / LH); if (sites[j]) { hideTip(); toggleCluster(sites[j].c); } };
  }
  // What "reported" means for this protein: the share of all its predicted partners, and of those below the cutoff, that
  // BioGRID reports (well-studied proteins have many reports whatever the prediction)
  function kbBase() {
    const all = B.partners.filter((x) => x.id !== P.key && !x.rep), low = all.filter((x) => x.best < CUT[cut]);
    const pct = (xs) => (xs.length ? `${(100 * xs.filter((x) => kbOf(x.id)).length / xs.length).toFixed(0)}%` : '—');
    return `<span>for comparison, ${pct(all)} of all ${fmtInt(all.length)} predicted partners of ${esc(P.gene)} are reported, and ${pct(low)} of the ${fmtInt(low.length)} below the cutoff; a report is evidence of an interaction, not of this interface</span>`;
  }
  function renderClusterInfo() {
    const box = $('#cluster-info'); if (!box) return;
    if (!clustered()) { box.innerHTML = `<p class="muted" style="margin:0">${M ? 'Nothing to cluster at this cutoff.' : 'Clustering…'}</p>`; return; }
    const mem = {}, cnt = {};
    M.preds.forEach((p, i) => { const c = M.labels[i]; (mem[c] ||= new Set()).add(who(p.partner).key); cnt[c] = (cnt[c] || 0) + 1; });
    const cap = 60, siteBy = new Map(sitesOf().map((x) => [x.c, x]));
    box.innerHTML = range(M.k).map((c) => {
      const rows = [...(mem[c] || [])].map((x) => sp.byKey.get(x) || { key: x, gene: x }).sort((a, b) => a.gene.toLowerCase().localeCompare(b.gene.toLowerCase()));
      const show = infoOpen.has(c) ? rows : rows.slice(0, cap);
      return `<div class="cl-row"><i style="background:${clusterColor(c, M.k)}"></i><div><b${(() => { const st = siteBy.get(c), q = st && bgRatio(st); return q != null ? ` title="background over this site's footprint: ${q.toFixed(1)}× the protein mean"` : ''; })()}>${clusterLabel(c)}</b> <span class="muted">(${rows.length} / ${cnt[c]})${KB ? ` · <span title="partners with a physical or genetic interaction reported in BioGRID ${KB.release}">${rows.filter((r) => kbOf(r.key)).length} of ${rows.length} reported</span>` : ''}</span>: ${show.map((r) => { const kb = kbOf(r.key);
        return `<a href="#/${sp.id}/${P.key}/${r.key}${scopeQ}"${kb ? ` class="${kb.ph ? 'kb-p' : ''}${kb.ph && kb.ge ? ' ' : ''}${kb.ge ? 'kb-g' : ''}" title="${kbTip(kb)}"` : ''}>${esc(r.gene)}</a>`; }).join(', ')}${rows.length > cap ? ` <button class="more" data-c="${c}">${infoOpen.has(c) ? 'show fewer' : `+${rows.length - cap} more`}</button>` : ''}</div></div>`; }).join('');
    box.querySelectorAll('.more').forEach((b) => b.onclick = () => { const c = +b.dataset.c; infoOpen.has(c) ? infoOpen.delete(c) : infoOpen.add(c); renderClusterInfo(); });
    const leg = $('#info-kb'), all = [...new Set(Object.values(mem).flatMap((x) => [...x]))].map(kbOf).filter(Boolean);
    const nP = all.filter((x) => x.ph).length, nG = all.filter((x) => x.ge).length;
    leg.hidden = !KB; if (KB) leg.innerHTML = `<span>reported in BioGRID ${KB.release}: <span class="kb-mark kb-p">physical</span> (${fmtInt(nP)} partner${nP === 1 ? '' : 's'}), <span class="kb-mark kb-g">genetic</span> (${fmtInt(nG)}), <span class="kb-mark kb-p kb-g">both</span></span>` + (KB ? kbBase() : '');
  }

  /* Interaction Residues: any partner, any screen, any model */
  const partnersByScore = [...B.partners].sort((a, b) => b.best - a.best);
  const PLACE = (key) => sp.byKey.get(key) || { key, gene: key, acc: '', clen: 0, len: 0, occ: [] };
  function fillPartners() {
    const sel = $('#res-partner'), keep = sel.value || (partnersByScore[0] && partnersByScore[0].id);
    sel.innerHTML = partnersByScore.map((p) => { const c = partnerCluster.get(p.id);
      return `<option value="${esc(p.id)}">${esc(gname(p.id))} — iLIS ${p.best.toFixed(3)}${c ? ` (${clusterLabel(c, true)})` : ''}</option>`; }).join('');
    if (keep) sel.value = keep;
    if (!$('#res-body').childElementCount) pickPartner();
  }
  function pickPartner() {
    const id = $('#res-partner').value, part = B.partners.find((p) => p.id === id); if (!part) return;
    const best = [...part.preds].sort((a, b) => b.iLIS - a.iLIS)[0];
    $('#res-rank').innerHTML = part.preds.map((p, i) => `<option value="${i}">${esc(runLabel(sp, B, p.run, P))}${p.rep ? ' (repeat)' : ''} · rank ${p.rank} · iLIS ${fmtNum(p.iLIS, 3)}</option>`).join('');
    $('#res-rank').value = String(part.preds.indexOf(best));
    drawResidues();
  }
  function drawResidues() {
    const id = $('#res-partner').value, part = B.partners.find((p) => p.id === id); if (!part) return;
    const pred = part.preds[+$('#res-rank').value] || part.preds[0];
    ifaceView($('#res-body'), { sp, P, O: PLACE(id), pred, B, canvasId: 'res-canvas' });
    if (sp.viruses) { const pid = $('#res-partner').value, R2 = sp.byKey.get(pid);   // a virus pair: its structure, where the Atlas can reach it
      virusStruct(sp, P.i, R2 ? R2.i : null).then((x) => { if (!gone() && $('#res-partner') && $('#res-partner').value === pid) structLink($('#res-struct'), x, 'btn'); }); }
  }
  $('#res-partner').onchange = pickPartner;
  $('#res-rank').onchange = drawResidues;
  $('#res-find').oninput = (e) => { const q = e.target.value.trim().toLowerCase(); if (!q) return;
    const hit = partnersByScore.find((p) => gname(p.id).toLowerCase() === q) || partnersByScore.find((p) => gname(p.id).toLowerCase().startsWith(q)) || partnersByScore.find((p) => p.id.toLowerCase().includes(q));
    if (hit && hit.id !== $('#res-partner').value) { $('#res-partner').value = hit.id; pickPartner(); } };

  /* 3D structure: the AlphaFold DB model in LIVIA's Mol* page, colored like clip.html */
  function colorComponents() {
    const N = S.len || P.clen || 1, GRAY = '#e6e6e6', at = new Array(N + 1).fill(GRAY);
    if (V.mode === 'plddt' && S.plddt) { for (let sr = 1; sr <= N; sr++) { const v = S.plddt.get('A:' + sr); if (v != null) at[sr] = plddtCol(v); } }
    else if (V.mode === 'am' && S.am) { for (let sr = 1; sr <= N; sr++) { const v = S.am[sr]; if (v != null && v >= V.amCut) at[sr] = amCol(v); } }
    else if (clustered() && S.mapOK) {   // cluster consensus: a residue takes a cluster's color when ≥ thr of its predictions contact it
      const byRes = {}, size = {};
      for (let i = 0; i < M.fingerprints.length; i++) { const lab = M.labels[i]; size[lab] = (size[lab] || 0) + 1;
        for (const r of M.fingerprints[i]) { const sr = toStruct(r); if (!sr || sr < 1) continue; (byRes[sr] ||= {})[lab] = (byRes[sr][lab] || 0) + 1; } }
      for (const res in byRes) { const r = +res; if (r < 1 || r > N) continue; let dom = null, dh = -1, dsz = -1; const cc = byRes[res];
        for (const c in cc) { if (!ACTIVE.has(+c)) continue; if (ACTIVE.size !== 1 && (size[c] || 0) < V.minC) continue; const hits = cc[c]; if (hits / (size[c] || 1) < V.thr) continue; const sz = size[c] || 0;
          if (hits > dh || (hits === dh && sz > dsz) || (hits === dh && sz === dsz && dom != null && +c < dom)) { dh = hits; dsz = sz; dom = +c; } }
        if (dom != null) at[r] = clusterColor(dom, M.k); }
    }
    const byColor = {}; let r = 1;
    while (r <= N) { const c = at[r]; let e = r; while (e < N && at[e + 1] === c) e++; (byColor[c] ||= []).push({ start: r, end: e }); r = e + 1; }
    return Object.entries(byColor).map(([color, ranges]) => ({ chain: 'A', ranges, color }));
  }
  function legend3D() {
    const box = $('#legend-3d'); if (!box) return;
    if (V.mode === 'plddt') box.innerHTML = [['#0053D6', 'very high (&gt; 90)'], ['#65CBF3', 'confident (70–90)'], ['#FFDB13', 'low (50–70)'], ['#FF7D45', 'very low (&lt; 50)']].map(([c, l]) => `<span><i style="background:${c}"></i>${l}</span>`).join('');
    else if (V.mode === 'am') box.innerHTML = `<span>benign</span><span class="amgrad"></span><span>pathogenic</span><span class="muted">mean AlphaMissense pathogenicity over all substitutions at the residue</span>`;
    else if (!clustered()) box.innerHTML = '<span class="muted">No clusters at this cutoff.</span>';
    else if (!S.mapOK) box.innerHTML = `<span class="muted">${esc(S.mapNote || '')}</span>`;
    else { const size = {}; for (const l of M.labels) size[l] = (size[l] || 0) + 1;
      const shown = range(M.k).filter((c) => ACTIVE.has(c) && (ACTIVE.size === 1 || size[c] >= V.minC));
      box.innerHTML = shown.map((c) => `<span><i style="background:${clusterColor(c, M.k)}"></i>${clusterLabel(c)}</span>`).join('') + '<span><i style="background:#e6e6e6"></i>not a consensus contact</span>'; }
  }
  function show3D() {
    V.want = true; if (S.state !== 'ready' || V.shown) return;
    const frame = $('#viewer3d-frame'); if (!frame) return;
    // A color update sent while the viewer page is still booting is lost (clustering can finish in that window), so
    // once the viewer says it is ready, the colors are sent again if they changed since it was built.
    const comps = colorComponents(), built = JSON.stringify(comps);
    const onReady = (ev) => { if (!ev.data || ev.data.type !== 'molstarReady') return; const f = $('#viewer3d-frame');
      if (gone() || !f || ev.source === f.contentWindow) window.removeEventListener('message', onReady);
      if (!gone() && f && ev.source === f.contentWindow && JSON.stringify(colorComponents()) !== built) recolor3D(); };
    window.addEventListener('message', onReady);
    frame.src = URL.createObjectURL(new Blob([buildMolstarPage(S.text, 'mmcif', comps, LIVIA)], { type: 'text/html' }));
    V.shown = true; $('#v3d-msg').hidden = true; legend3D();
  }
  function recolor3D() { legend3D(); if (V.shown) applyColorsToMolstarFrame('viewer3d-frame', colorComponents(), 'mmcif'); }
  function mapStruct() {   // the clustered construct onto the AlphaFold DB model (the UniProt sequence)
    const es = S.entrySeq; if (!es) return;
    const cLen = CQ.qLen || P.clen;
    if (qSeq && qSeq.length === cLen) {
      if (qSeq === es) { S.map = null; S.mapOK = true; S.mapNote = 'sequence 1:1'; }
      else { const mi = CLIPResolver.alignMap(qSeq, es); S.map = mi.map; S.mapOK = mi.covered > 0; S.mapNote = `remapped (${mi.method}, ${Math.round(100 * mi.covered / qSeq.length)}% matched)`; }
    } else if (es.length === cLen) { S.map = null; S.mapOK = true; S.mapNote = 'same length'; }
    else { S.map = null; S.mapOK = false; S.mapNote = `The clustered construct (${fmtInt(cLen)} aa) differs from the model (${fmtInt(es.length)} aa), so clusters are not placed on it.`; }
    $('#struct-badge').innerHTML = `AlphaFold DB <a href="https://alphafold.ebi.ac.uk/entry/${esc(P.acc)}" target="_blank" rel="noopener">${esc(P.acc)}</a> · ${esc(S.mapOK ? S.mapNote : 'not mapped')}`;
  }
  async function loadStructure() {
    const msg = (t) => { const m = $('#v3d-msg'); if (m) { m.hidden = false; m.innerHTML = t; } };
    if (!P.acc) { S.state = 'none'; msg('No UniProt accession for this protein, so there is no AlphaFold DB model to show.'); return; }
    try { await liviaReady(); } catch (e) { if (!gone()) { S.state = 'none'; msg(esc(e.message)); } return; }
    if (gone()) return;
    domainsOf(P.acc).then((d) => { if (gone()) return; S.domains = d; drawFreq(); if (clustered()) renderSites(); });
    const entry = await afdbEntry(P.acc);
    if (gone()) return;
    if (entry && entry.failed) { S.state = 'none'; msg(`The AlphaFold DB model of ${esc(P.acc)} could not be loaded (no answer from alphafold.ebi.ac.uk). Reload the page to try again.`); $('#struct-badge').textContent = ''; return; }
    if (!entry) { S.state = 'none'; msg(`No AlphaFold DB model for ${esc(P.acc)}${(P.clen || P.len) > 2700 ? ' (the database has none for proteins longer than 2,700 residues)' : ''}.`); $('#struct-badge').textContent = ''; return; }
    S.uniSeq = S.entrySeq = entry.seq || ''; mapStruct();
    try { S.text = await (await fetch(entry.cifUrl)).text(); } catch (e) { if (!gone()) { S.state = 'none'; msg('The AlphaFold DB model could not be downloaded.'); } return; }
    if (gone()) return;
    S.len = entry.seq.length || P.clen; S.plddt = parseBfactorsPerResidue(S.text, 'cif'); S.state = 'ready';
    if (!S.mapOK) { V.mode = 'plddt'; app.querySelectorAll('#cmode button').forEach((b) => b.classList.toggle('on', b.dataset.m === 'plddt')); }
    drawFreq(); msg('The 3D viewer loads when this card scrolls into view.');
    if (V.want) show3D();
    if (entry.amUrl) alphaMissense(entry.amUrl).then((a) => { if (!a || gone()) return; S.am = a; $('#cm-am').hidden = false; drawFreq(); });
  }
  $('#cmode').onclick = (e) => { const b = e.target.closest('button'); if (!b) return; V.mode = b.dataset.m; V.auto = false;
    app.querySelectorAll('#cmode button').forEach((x) => x.classList.toggle('on', x === b)); $('#am-cut-wrap').hidden = V.mode !== 'am'; recolor3D(); };
  $('#commonality').onchange = (e) => {   // a preset (10–90%), or "custom…" for any value from 1 to 100%
    if (e.target.value === 'custom') { $('#commonality-custom').hidden = false; $('#commonality-pct').value = +(V.thr * 100).toFixed(1); $('#commonality-pct').focus(); $('#commonality-pct').select(); return; }   // selected: typing replaces the value
    $('#commonality-custom').hidden = true; V.thr = +e.target.value; recolor3D(); };
  $('#commonality-pct').onchange = (e) => { const pct = Math.min(100, Math.max(1, parseFloat(e.target.value) || 50)); e.target.value = pct; V.thr = pct / 100; recolor3D(); };
  $('#min-cluster').onchange = (e) => { V.minC = Math.max(1, +e.target.value || 1); recolor3D(); };
  $('#am-cutoff').onchange = (e) => { V.amCut = +e.target.value; recolor3D(); };
  { const io3 = new IntersectionObserver((ents) => { if (ents.some((x) => x.isIntersecting)) { io3.disconnect(); show3D(); } }, { rootMargin: '300px' }); io3.observe($('#c-3d')); }

  /* Partners by score: one point per partner (its best model; dot size = iLIS average) or per prediction, any two scores,
     cluster colors — on a canvas, so it exports */
  const present = Object.keys(METRICS).filter((k) => k === '_rank' || B.counted.some((p) => Number.isFinite(p[k]) && p[k] !== 0));
  $('#sc-y').innerHTML = present.map((k) => `<option value="${k}">${METRICS[k]}</option>`).join('');
  $('#sc-x').innerHTML = present.map((k) => `<option value="${k}">${METRICS[k]}</option>`).join('');
  $('#sc-y').value = 'iLIS'; $('#sc-x').value = present.includes('iLISA') ? 'iLISA' : 'ipTM';
  const suggestItems = () => {   // each partner once, by its best iLIS, in name order; the dot is its cluster's color
    const by = new Map(); for (const p of B.partners) { const n = gname(p.id), cur = by.get(n); if (!cur || p.best > cur.best) by.set(n, p); }
    return [...by].sort((a, b) => a[0].localeCompare(b[0], undefined, { numeric: true })).map(([n, p]) => { const c = partnerCluster.get(p.id), k = M ? M.k : 1;   // UL2 before UL10
      return { name: n, best: p.best, dot: c ? clusterColor(c, k) : '#DDE3EA', sty: bandSty(FPR.iLIS, p.best), tip: `iLIS ${p.best.toFixed(3)}: ${bandLabel[bandIn(FPR.iLIS, p.best)]}${c ? ` · ${clusterLabel(c)}` : ''}` }; }); };
  ['#sc-find', '#fp-find', '#res-find'].forEach((q) => partnerSuggest($(q), suggestItems));   // all three Find partner boxes suggest partners with their iLIS
  ['#sc-x', '#sc-y', '#sc-pts'].forEach((s) => { $(s).onchange = () => drawScatter(); });
  $('#fp-find').oninput = () => { if (clustered()) drawHeatmap(); };
  $('#res-q').onchange = (e) => setRes(e.target.value);
  $('#res-q').onkeydown = (e) => { if (e.key === 'Enter') { e.preventDefault(); setRes(e.target.value); } };
  $('#sc-find').onchange = () => drawScatter();
  $('#sc-find').onkeydown = (e) => { if (e.key === 'Enter') { e.preventDefault(); drawScatter(); } };
  $('#sc-find').oninput = (e) => { if (!e.target.value) drawScatter(); };
  function drawScatter() {
    const cv = $('#scatter-canvas'); if (!cv) return;
    const xK = $('#sc-x').value, yK = $('#sc-y').value, mode = $('#sc-pts').value, k = M ? M.k : 1;
    const rankBy = yK !== '_rank' ? yK : xK !== '_rank' ? xK : 'iLIS';
    let grank = null;
    if (xK === '_rank' || yK === '_rank') { grank = new Map(); B.counted.map((p, i) => [p, i]).sort((a, b) => ((b[0][rankBy] || 0) - (a[0][rankBy] || 0)) || a[1] - b[1]).forEach(([p], i) => grank.set(p, i + 1)); }   // partners: their best model's rank
    const val = (p, key) => (key === '_rank' ? grank.get(p) : p[key]);
    const per = mode === 'partner', pts = [];
    const src = per ? B.partners.filter((t) => t.id !== P.key && !t.rep).map((t) => ({ t, p: t.counted.reduce((a, b) => ((b.iLIS || 0) > (a.iLIS || 0) ? b : a), t.counted[0]) }))
      : B.counted.map((p) => ({ t: null, p }));
    for (const { t, p } of src) { const x = val(p, xK), y = val(p, yK); if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
      const c = per ? partnerCluster.get(t.id) : predCluster.get(p.label + '|' + p.rank); pts.push({ p, t, x, y, c: c != null && (per || mode === 'all' || p.rank === 1) ? c : 0 }); }
    pts.sort((a, b) => (a.c ? 1 : 0) - (b.c ? 1 : 0));
    const W = cv.parentElement.clientWidth, H = 450, m = { l: 62, r: 18, t: 24, b: 46 };
    const UNIT = new Set(['iLIS', 'ipTM', 'pTM', 'LIS', 'cLIS', 'ipSAE', 'actifpTM']);   // scores bounded by 1 keep their whole 0–1 range
    const dom = (key, vs) => (key === '_rank' ? [0, Math.max(1, vs.length)] : UNIT.has(key) ? [0, 1] : key === 'qPl' || key === 'pPl' ? [0, 100]
      : [Math.min(0, d3.min(vs) ?? 0), Math.max(1e-6, (d3.max(vs) ?? 1) * 1.04)]);
    const xs = d3.scaleLinear(dom(xK, pts.map((q) => q.x)), [m.l, W - m.r]).nice(), ys = d3.scaleLinear(dom(yK, pts.map((q) => q.y)), [H - m.b, m.t]).nice();
    const g = canvasCtx(cv, W, H);
    const title = (key) => (key === '_rank' ? `global rank (by ${METRICS[rankBy]})` : METRICS[key] + (per && !ONE ? ' (best model)' : ''));
    canvasAxes(g, xs, ys, m, W, H, title(xK), title(yK));
    g.save(); g.lineWidth = 1; g.font = '10.5px "IBM Plex Mono", ui-monospace, monospace';
    const dash = (x0, y0, x1, y1) => { const L = Math.hypot(x1 - x0, y1 - y0); g.beginPath(); for (let t = 0; t < L; t += 8) { const a = t / L, b = Math.min(L, t + 4) / L; g.moveTo(x0 + (x1 - x0) * a, y0 + (y1 - y0) * a); g.lineTo(x0 + (x1 - x0) * b, y0 + (y1 - y0) * b); } g.stroke(); };   // dashes as segments: canvas2svg has no setLineDash
    [10, 5, 1].forEach((f, j) => {   // benchmarked cutoffs of the metric on each axis (single models)
      g.strokeStyle = BAND[f]; g.fillStyle = BAND[f];
      const vy = FPR[yK] && FPR[yK][j], vx = FPR[xK] && FPR[xK][j];
      if (vy != null && vy <= ys.domain()[1]) { const y = Math.round(ys(vy)) + 0.5; dash(m.l, y, W - m.r, y); g.textAlign = 'right'; g.textBaseline = 'bottom'; g.fillText(`${f}% FPR (${vy})`, W - m.r - 2, y - 3); }
      if (vx != null && vx <= xs.domain()[1]) { const x = Math.round(xs(vx)) + 0.5; dash(x, m.t, x, H - m.b); g.textAlign = 'center'; g.textBaseline = 'bottom'; g.fillText(`${f}%`, x, m.t - 5); }
    });
    g.restore();
    const rad = (q) => (per ? 2.2 + 8 * Math.min(1, (q.t.avg || 0) / 0.8) : q.c ? 4.3 : 2.5);
    for (const q of pts) { g.beginPath(); g.arc(xs(q.x), ys(q.y), rad(q), 0, 2 * Math.PI); g.fillStyle = q.c ? clusterColor(q.c, k) : '#CDD3DB';
      g.globalAlpha = per ? (q.c ? 0.8 : 0.4) : 1; g.fill(); g.globalAlpha = 1; if (q.c || per) { g.lineWidth = 0.7; g.strokeStyle = '#fff'; g.stroke(); } }
    const reportedPts = KB ? pts.filter((q) => (kbOf(q.t ? q.t.id : q.p.partner) || {}).ph) : [], passing = (q) => (q.p.iLIS || 0) >= CUT[10];
    const ringed = scKnown === 'all' ? reportedPts : scKnown === 'pass' ? reportedPts.filter(passing) : [];   // reported in BioGRID (physical): a translucent dark ring
    for (const q of ringed) { g.beginPath(); g.arc(xs(q.x), ys(q.y), rad(q) + 1.3, 0, 2 * Math.PI);
      g.lineWidth = passing(q) ? 1.4 : 1; g.strokeStyle = passing(q) ? 'rgba(23,38,58,0.6)' : 'rgba(23,38,58,0.3)'; g.stroke(); }
    if (per) {   // label the top partners where a label fits, never on a cutoff label or another label
      const placed = [], shown = [];
      [10, 5, 1].forEach((f, j) => { const vy = FPR[yK] && FPR[yK][j], vx = FPR[xK] && FPR[xK][j];
        if (vy != null) placed.push([W - m.r - 118, ys(vy) - 15, W - m.r, ys(vy)]); if (vx != null) placed.push([xs(vx) - 16, 0, xs(vx) + 16, m.t]); });
      g.font = '600 11.5px "IBM Plex Sans", system-ui, sans-serif';
      for (const q of [...pts].sort((a, b) => (b.t.best || 0) - (a.t.best || 0)).slice(0, 10)) {
        const cx = xs(q.x), cy = ys(q.y), rr = rad(q), name = gname(q.t.id), w = g.measureText(name).width + 4, h = 13;
        for (const [tx, ty, anchor] of [[cx + rr + 4, cy + 4, 'left'], [cx - rr - 4, cy + 4, 'right'], [cx, cy - rr - 5, 'center'], [cx, cy + rr + 13, 'center']]) {
          const bx = anchor === 'left' ? tx : anchor === 'right' ? tx - w : tx - w / 2, box = [bx - 3, ty - h, bx + w + 3, ty + 4];   // a little room around each label
          if (box[0] < m.l || box[2] > W - m.r || box[1] < m.t - 2 || box[3] > H - m.b) continue;
          if (placed.some((o) => box[0] < o[2] && box[2] > o[0] && box[1] < o[3] && box[3] > o[1])) continue;
          placed.push(box); shown.push([name, tx, ty, anchor]); break; } }
      for (const [name, tx, ty, anchor] of shown) { g.textAlign = anchor; g.textBaseline = 'alphabetic'; g.lineWidth = 3; g.strokeStyle = 'rgba(255,255,255,0.92)'; g.strokeText(name, tx, ty); g.fillStyle = '#17263A'; g.fillText(name, tx, ty); }
    }
    const find = $('#sc-find').value.trim().toLowerCase();
    if (find) {
      const gl = (p) => gname(p.id).toLowerCase(), hit = B.partners.find((p) => gl(p) === find) || B.partners.find((p) => gl(p).startsWith(find)) || B.partners.find((p) => gl(p).includes(find));
      $('#sc-find').classList.toggle('nf', !hit);
      if (hit) for (const q of pts.filter((q) => (q.t ? q.t.id : q.p.partner) === hit.id)) { const x = xs(q.x), y = ys(q.y);
        g.beginPath(); g.arc(x, y, 8.5, 0, 2 * Math.PI); g.lineWidth = 2.5; g.strokeStyle = '#17263A'; g.stroke();
        g.font = '600 12px "IBM Plex Sans", system-ui, sans-serif'; g.textAlign = 'center'; g.textBaseline = 'alphabetic'; g.lineWidth = 3; g.strokeStyle = 'rgba(255,255,255,0.92)'; const t = per ? gname(hit.id) : `${gname(hit.id)} (R${q.p.rank})`;
        g.strokeText(t, x, y - 13); g.fillStyle = '#17263A'; g.fillText(t, x, y - 13); }
    } else $('#sc-find').classList.remove('nf');
    const qt = d3.quadtree(pts, (q) => xs(q.x), (q) => ys(q.y));
    const near = (e) => { const b = cv.getBoundingClientRect(); return qt.find(e.clientX - b.left, e.clientY - b.top, 10); };
    const fmtM = (key, v) => (Number.isFinite(v) ? (['iLIA', 'iLISA', 'LIA', 'cLIA', 'qPl', 'pPl'].includes(key) ? v.toFixed(1) : v.toFixed(3)) : '–');
    cv.onmousemove = (e) => { const q = near(e); if (!q) { cv.style.cursor = ''; return hideTip(); } cv.style.cursor = 'pointer'; const p = q.p;
      const keys = [...new Set(['iLIS', 'iLISA', 'iLIA', 'ipTM', yK, xK])].filter((x) => x !== '_rank' && present.includes(x));
      showTip(`<b>${esc(gname(p.partner))}</b> · ${per && !ONE ? `best of ${fmtInt(q.t.counted.length)} model${q.t.counted.length === 1 ? '' : 's'} (${esc(runLabel(sp, B, p.run, P))}, rank ${p.rank}) · iLIS average ${fmtNum(q.t.avg, 3)}` : `${esc(runLabel(sp, B, p.run, P))} · rank ${p.rank}`}${q.c ? ` · <span style="color:${clusterColor(q.c, k)}">●</span> ${clusterLabel(q.c)}` : ''}<br>${keys.map((x) => `${METRICS[x]} ${fmtM(x, p[x])}`).join(' · ')}`, e.clientX, e.clientY); };
    cv.onmouseleave = hideTip;
    cv.onclick = (e) => { const q = near(e); if (!q) return; hideTip(); location.hash = `#/${sp.id}/${P.key}/${q.p.partner}${scopeQ}`; };
    const r1 = pearson(pts.map((q) => q.x), pts.map((q) => q.y)), rho = spearman(pts.map((q) => q.x), pts.map((q) => q.y));
    { const d = $('#sc-def'), used = [xK, yK].some((x) => x === 'iLISA' || x === 'iLIA'); d.hidden = !used;
      if (used) d.textContent = 'iLISA = iLIS × iLIA; iLIA = √(LIA × cLIA), where LIA and cLIA count the residue pairs that enter LIS and cLIS'; }
    $('#sc-rho').textContent = `${fmtInt(pts.length)} ${per ? (ONE ? 'partners, one model each' : 'partners, the best model of each pair') : mode === 'rank1' ? 'predictions; rank-1 models colored by cluster' : 'predictions'} · ${Number.isFinite(r1) ? `Pearson r = ${r1.toFixed(3)} · ` : ''}${Number.isFinite(rho) ? `Spearman ρ = ${rho.toFixed(3)}` : ''}`;
    const cnt = {}; let other = 0; for (const q of pts) q.c ? (cnt[q.c] = (cnt[q.c] || 0) + 1) : other++;
    const unit = per ? 'partners' : 'predictions';
    $('#sc-legend').innerHTML = (per ? `<span><i style="background:#A7B2BF;border-radius:50%"></i>dot size: ${ONE ? 'iLIS' : 'iLIS average over the models'}</span>` : '')
      + Object.keys(cnt).map(Number).sort((a, b) => a - b).map((c) => `<span><i style="background:${clusterColor(c, k)};border-radius:50%"></i>${clusterLabel(c)} (${cnt[c]})</span>`).join('')
      + `<span><i style="background:#CDD3DB;border-radius:50%"></i>not clustered (${fmtInt(other)} ${unit})</span>`
      + (KB ? `<label class="kb-toggle"><i class="kb-ring"></i>reported in BioGRID ${KB.release}, physical<select id="sc-known" aria-label="Which reported partners to ring">${[
        ['pass', `past 10% FPR (${fmtInt(reportedPts.filter(passing).length)})`], ['all', `all (${fmtInt(reportedPts.length)})`], ['off', 'none']].map(([v, l]) => `<option value="${v}"${v === scKnown ? ' selected' : ''}>${l}</option>`).join('')}</select></label>` : '');
    if (KB) $('#sc-known').onchange = (e) => { scKnown = e.target.value; drawScatter(); drawTopList(); };
    attachExport('scatter-canvas', `atlas_${P.gene}_partners`, drawScatter);
  }

  /* top partners beside the scatter: the best model of each pair */
  function drawTopList() {
    const list = B.partners.filter((p) => p.id !== P.key && !p.rep).map((p) => ({ ...p, gene: gname(p.id), c: partnerCluster.get(p.id) || 0 })), k = M ? M.k : 1;
    const tip = (cuts, v, what) => `${what}: ${bandLabel[bandIn(cuts, v)]} (cutoffs ${cuts.join(' / ')})`;
    // the steadiest calls first: the share of a pair's models past the 10% cutoff (3 of 3 = 5 of 5, whatever the screen), then average iLIS
    const past = (p) => p.counted.filter((x) => x.iLIS >= CUT[10]).length / Math.max(1, p.counted.length);
    $('#toplist').innerHTML = [...list].sort((a, b) => past(b) - past(a) || b.avg - a.avg || b.best - a.best).slice(0, 12).map((p) => { const xs = partnerIsos(p); return `<li><span class="tl-name"><a href="#/${sp.id}/${P.key}/${p.id}${scopeQ}" title="${esc(p.gene)}">${esc(p.gene)}</a>${KB && scKnown !== 'off' && (kbOf(p.id) || {}).ph ? `<i class="kb-ring" title="${kbTip(kbOf(p.id))}"></i>` : ''}${xs ? `<span class="iso-tag" title="${esc(isoTip(p, xs))}">×${xs.length}</span>` : ''}</span>
      <span class="tl-c" title="${p.c ? clusterLabel(p.c) : 'not clustered at this cutoff'}"><span class="mdot" style="background:${p.c ? clusterColor(p.c, k) : '#DDE3EA'}"></span>${p.c ? clusterLabel(p.c, true) : '—'}</span>
      <span class="num" style="${bandSty(FPR.iLIS, p.best)}" title="${tip(FPR.iLIS, p.best, ONE ? 'iLIS' : 'iLIS best')}">${p.best.toFixed(3)}</span>
      ${ONE ? '' : `<span class="num" style="${bandSty(FPR_AVG.iLIS, p.avg)}" title="${tip(FPR_AVG.iLIS, p.avg, 'iLIS average')}">${p.avg.toFixed(3)}</span>`}
      <span class="num" style="${bandSty(FPR.ipTM, p.iptmBest)}" title="${tip(FPR.ipTM, p.iptmBest, ONE ? 'ipTM' : 'ipTM best')}">${p.iptmBest.toFixed(2)}</span>
      ${ONE ? '' : `<span class="num" style="${bandSty(FPR_AVG.ipTM, p.iptmAvg)}" title="${tip(FPR_AVG.ipTM, p.iptmAvg, 'ipTM average')}">${p.iptmAvg.toFixed(2)}</span>`}
      ${ONE ? '' : `<span class="num tl-past" title="models past the 10% FPR cutoff (iLIS ≥ ${CUT[10]}), of ${p.counted.length}">${p.counted.filter((x) => x.iLIS >= CUT[10]).length}/${p.counted.length}</span>`}</li>`; }).join('');
  }
  reportedOf(sp, P.i).then((k) => { if (gone() || !k) return; KB = k; drawScatter(); drawTopList(); drawTable(); if (clustered()) renderClusterInfo(); }).catch(() => {});   // this protein's shard
  $('#cut-seg').onclick = (e) => { const f = e.target.dataset.f; if (!f) return; cut = +f; [...$('#cut-seg').children].forEach((b) => b.classList.toggle('on', b.dataset.f === f)); cluster(); };

  /* partner table */
  const T = { sort: 'best', asc: false, page: 0, band: P.pos10 ? 10 : 0, filter: '', src: 0 };
  // a virus protein: the virus page's pairs table, one row per partner (residues, the model in LIVIA, every score), from its virus's pairs file
  const VX = sp.id === 'virus' && sp.viruses ? new Map() : null;
  if (VX) { const v = sp.viruses.find((x) => (x.members || []).includes(P.i));
    if (v) virusRows(sp, v.taxid).then((rows) => { if (gone()) return; for (const x of rows) { if (x.a === P.i) VX.set(x.b, x); else if (x.b === P.i) VX.set(x.a, x); } drawTable(); }); }
  if (VX) T.sort = 'm:iLIS';
  const cols = VX ? [['gene', 'Partner'], ['c', 'Cluster'], ['pair', 'Pair'], ['d3', '3D'], ...VMCOL.map((k) => ['m:' + k, k]), ['contacts', 'Contacts']]
    : ONE ? [['gene', 'Partner'], ['c', 'Cluster'], ['src', 'Source'], ['name', 'Protein'], ['best', 'iLIS'], ['iptmBest', 'ipTM'], ['contacts', 'Contacts']]
    : [['gene', 'Partner'], ['c', 'Cluster'], ['src', 'Source'], ['name', 'Protein'], ['best', 'iLIS best'], ['avg', 'iLIS avg'], ['iptmBest', 'ipTM best'], ['iptmAvg', 'ipTM avg'], ['contacts', 'Contacts'], ['pass', 'Models past']];
  function drawTable() {
    let list = B.partners.map((p) => { const r = sp.byKey.get(p.id); return { ...p, gene: r ? r.gene : p.id, name: r ? r.name : '', c: partnerCluster.get(p.id) || 0, pass: p.counted.filter((x) => x.iLIS >= CUT[10]).length, ...(VX ? (() => { const x = r && VX.get(r.i), o = { vx: x };
      for (const k of VMCOL) o['m:' + k] = x && Number.isFinite(x.m[k]) ? x.m[k] : NaN;
      o['m:iLIS'] = p.best; o['m:ipTM'] = p.iptmBest; return o; })() : {}) }; });
    if (T.src) list = list.filter((p) => (SETS ? p.sets.includes(T.src) : p.src & T.src));
    if (T.band) list = list.filter((p) => p.best >= CUT[T.band]);
    if (T.filter) { const f = T.filter.toLowerCase(); list = list.filter((p) => p.gene.toLowerCase().includes(f) || (p.name || '').toLowerCase().includes(f) || p.id.toLowerCase().includes(f)); }
    const key = T.sort; const nv = (x) => (Number.isFinite(x) ? x : -Infinity); list.sort((a, b) => (typeof a[key] === 'string' ? a[key].localeCompare(b[key]) : nv(a[key]) - nv(b[key])) * (T.asc ? 1 : -1));
    const per = 40, pages = Math.max(1, Math.ceil(list.length / per)); T.page = Math.min(T.page, pages - 1);
    const view = list.slice(T.page * per, T.page * per + per), k = M ? M.k : 1;
    $('#pt-note').innerHTML = `${fmtInt(list.length)} shown · ${fmtInt(B.partners.filter((x) => !x.rep).length)} predicted${KB ? ` · reported in BioGRID ${esc(KB.release)}: <span class="kb-mark kb-p">physical</span> <span class="kb-mark kb-g">genetic</span> <span class="kb-mark kb-p kb-g">both</span>` : ''}`;
    $('#pt').innerHTML = `<thead><tr>${cols.map(([c, l]) => `<th data-c="${c}" class="${T.sort === c ? 'sorted' + (T.asc ? ' asc' : '') : ''}${(['best', 'avg', 'iptmBest', 'iptmAvg', 'contacts', 'pass'].includes(c) || c.startsWith('m:')) ? ' n' : ''}">${l}</th>`).join('')}</tr></thead><tbody>${view.map((p) => {
      const b = bandOf(p.best), xs = partnerIsos(p), open = xs && isoOpen.has(p.id);
      const tag = xs ? ` <button type="button" class="iso-tag" data-iso="${esc(p.id)}" aria-expanded="${!!open}" title="${esc(isoTip(p, xs))}">${xs.length} ${isoWord(xs)} ${open ? '▾' : '▸'}</button>` : '';
      const subs = open ? xs.map((x) => { const c = isoCluster(x), bb = bandOf(x.best);
        return `<tr class="iso-sub"><td class="g">${esc(x.label)}</td><td>${c ? `<span class="mdot" style="background:${clusterColor(c, k)}"></span>${clusterLabel(c, true)}` : '<span class="muted">—</span>'}</td>
          <td></td><td class="nm">${fmtInt(x.n)} models</td><td class="n v" style="color:${BAND_TXT[bb]};font-weight:${BAND_W[bb]}" title="${bandLabel[bb]}">${x.best.toFixed(3)}</td>
          ${ONE ? '' : `<td class="n v" style="${bandSty(FPR_AVG.iLIS, x.avg)}">${x.avg.toFixed(3)}</td>`}<td class="n v" style="${bandSty(FPR.ipTM, x.iptmBest)}">${x.iptmBest.toFixed(2)}</td>
          ${ONE ? '' : `<td class="n v" style="${bandSty(FPR_AVG.ipTM, x.iptmAvg)}">${x.iptmAvg.toFixed(2)}</td>`}<td class="n">${fmtInt(x.contacts)}</td>${ONE ? '' : `<td class="n">${x.pass} / ${x.n}</td>`}</tr>`; }).join('') : '';
      if (VX) { const x = p.vx, href = `#/${sp.id}/${P.key}/${p.id}${scopeQ}`, kb = kbOf(p.id);
        const same = p.rep ? ` <span class="muted" title="The same two sequences as ${esc(gname(p.repOf || ''))} (UniProt holds this sequence under several accessions): counted once, under that partner">same as ${esc(gname(p.repOf || ''))}</span>` : '';
        return `<tr${p.rep ? ' class="rep"' : ''}><td class="g"><a href="#/${sp.id}/${p.id}"${kb ? ` class="${kb.ph ? 'kb-p' : ''}${kb.ph && kb.ge ? ' ' : ''}${kb.ge ? 'kb-g' : ''}" title="${kbTip(kb)}"` : ` title="${esc(p.name)}"`}>${esc(p.gene)}</a>${same}</td>
          <td>${p.c ? `<span class="mdot" style="background:${clusterColor(p.c, k)}"></span>${clusterLabel(p.c, true)}` : '<span class="muted">—</span>'}</td>
          <td><a href="${href}" title="interaction residues">residues</a></td>
          <td>${x && x.shown && x.model ? `<a href="${LIVIA}dimer.html?id=${encodeURIComponent(x.model)}" target="_blank" rel="noopener" title="${esc(x.model)} in LIVIA, from the AlphaFold Database">LIVIA ↗</a>`
            : x && x.addr && x.model ? `<a href="#" class="arch" data-m="${esc(x.model)}" title="${esc(x.model)} in LIVIA, read from the release archive at EBI (${((x.addr.cif_len + x.addr.pae_len) / 1048576).toFixed(1)} MB; not displayed by the AlphaFold Database)">LIVIA ↗</a>`
            : '<span class="muted" title="not displayed by the AlphaFold Database (below its filter); the structure opens from the release archive once its part is indexed">—</span>'}</td>
          ${VMCOL.map((kk) => { const val = p['m:' + kk]; return kk === 'iLIS' ? `<td class="n v"><a href="${href}" style="color:${BAND_TXT[b]};font-weight:${BAND_W[b]}" title="${bandLabel[b]}">${vmfmt(kk, val)}</a></td>`
            : kk === 'ipTM' ? `<td class="n v" style="${bandSty(FPR.ipTM, val)}">${vmfmt(kk, val)}</td>` : `<td class="n">${vmfmt(kk, val)}</td>`; }).join('')}
          <td class="n">${fmtInt(p.contacts)}</td></tr>`; }
      const same = p.rep ? ` <span class="muted" title="The same two sequences as ${esc(gname(p.repOf || ''))}: counted once, under that partner">same as ${esc(gname(p.repOf || ''))}</span>` : '';
      return `<tr${p.rep ? ' class="rep"' : ''}><td class="g"><a href="#/${sp.id}/${P.key}/${p.id}${scopeQ}"${(() => { const kb = kbOf(p.id); return kb ? ` class="${kb.ph ? 'kb-p' : ''}${kb.ph && kb.ge ? ' ' : ''}${kb.ge ? 'kb-g' : ''}" title="${kbTip(kb)}"` : ''; })()}>${esc(p.gene)}</a>${tag}${same}</td>
        <td>${p.c ? `<span class="mdot" style="background:${clusterColor(p.c, k)}"></span>${clusterLabel(p.c, true)}` : '<span class="muted">—</span>'}</td>
        <td class="srcc"${(() => { const t = overlapNote(sp, B, p.preds, P); return t ? ` title="${esc(t)}"` : ''; })()}>${SETS ? setBadges(p.sets) : srcBadges(sp, p.src)}</td><td class="nm" title="${esc(p.name)}">${esc(short(p.name))}</td>
        <td class="n v" style="color:${BAND_TXT[b]};font-weight:${BAND_W[b]}" title="${bandLabel[b]}">${p.best.toFixed(3)}</td>
        ${ONE ? '' : `<td class="n v" style="${bandSty(FPR_AVG.iLIS, p.avg)}" title="${bandLabel[bandIn(FPR_AVG.iLIS, p.avg)]} (average-model cutoffs)">${p.avg.toFixed(3)}</td>`}
        <td class="n v" style="${bandSty(FPR.ipTM, p.iptmBest)}" title="${bandLabel[bandIn(FPR.ipTM, p.iptmBest)]}">${p.iptmBest.toFixed(2)}</td>
        ${ONE ? '' : `<td class="n v" style="${bandSty(FPR_AVG.ipTM, p.iptmAvg)}" title="${bandLabel[bandIn(FPR_AVG.ipTM, p.iptmAvg)]} (average-model cutoffs)">${p.iptmAvg.toFixed(2)}</td>`}
        <td class="n">${fmtInt(p.contacts)}</td>${ONE ? '' : `<td class="n" title="models past the 10% FPR cutoff${p.preds.length > p.counted.length ? ` (${p.preds.length - p.counted.length} more in repeat runs, not counted)` : ''}">${p.pass} / ${p.counted.length}</td>`}</tr>${subs}`; }).join('')}</tbody>`;
    if (VX) $('#pt').querySelectorAll('a.arch').forEach((a) => a.onclick = (e) => { e.preventDefault(); const x = [...VX.values()].find((y) => y.model === a.dataset.m); if (x && x.addr) openFromArchive(x.model, x.addr, a); });
    $('#pt').querySelectorAll('[data-iso]').forEach((btn) => btn.onclick = () => { const id = btn.dataset.iso; if (isoOpen.has(id)) isoOpen.delete(id); else isoOpen.add(id); drawTable(); });
    $('#pt').querySelectorAll('th').forEach((th) => th.onclick = () => { const c = th.dataset.c; T.asc = T.sort === c ? !T.asc : (c === 'gene' || c === 'name' || c === 'c'); T.sort = c; drawTable(); });
    $('#pager').innerHTML = pages > 1 ? `<button class="btn" id="pp" ${T.page ? '' : 'disabled'}>Previous</button><span>Page ${T.page + 1} of ${pages}</span><button class="btn" id="pn" ${T.page < pages - 1 ? '' : 'disabled'}>Next</button>` : '';
    if (pages > 1) { $('#pp').onclick = () => { T.page--; drawTable(); }; $('#pn').onclick = () => { T.page++; drawTable(); }; }
  }
  $('#pt-band').value = String(T.band);
  $('#pt-band').onchange = (e) => { T.band = +e.target.value; T.page = 0; drawTable(); };
  if (SETS) $('#pt-src').innerHTML = `<option value="">all categories</option>${SETS.map((x) => `<option>${esc(x)}</option>`).join('')}`;   // one screen: filter by its categories
  $('#pt-src').onchange = (e) => { T.src = SETS ? e.target.value : +e.target.value; T.page = 0; drawTable(); };
  $('#pt-filter').oninput = (e) => { T.filter = e.target.value; T.page = 0; drawTable(); };

  /* network: edge color = what BioGRID reports for the pair, edge width = average iLIS */
  const netBox = $('#net');
  ['#net-shade', '#net-kb', '#net-ev'].forEach((q) => { $(q).onchange = () => { $('#net-ev').disabled = !$('#net-kb').checked; if (NET && NET.restyle) NET.restyle(); }; });
  { const w = d3.select('#net-w'); [[0.1, 10], [0.4, 95], [0.7, 180]].forEach(([a, x0]) => { w.append('line').attr('x1', x0).attr('x2', x0 + 44).attr('y1', 10).attr('y2', 10).attr('stroke', '#50637A').attr('stroke-width', EWID(a)).attr('stroke-linecap', 'round');
      w.append('text').attr('x', x0 + 22).attr('y', 27).attr('text-anchor', 'middle').attr('font-size', 10.5).attr('font-family', 'IBM Plex Mono').attr('fill', '#5A697C').text(a.toFixed(1)); }); }
  const io = new IntersectionObserver(async (ents) => { if (!ents.some((e) => e.isIntersecting)) return; io.disconnect(); await drawNet(); }, { rootMargin: '200px' });
  io.observe(netBox);
  $('#net-n').onchange = drawNet; $('#net-cut').onchange = drawNet;
  async function drawNet() {
    let E, K; try { [E, K] = await Promise.all([edges(sp, B.setId), reported(sp)]); } catch (e) { if (!gone()) netBox.innerHTML = `<div class="empty">${esc(e.message)}</div>`; return; }
    if (gone()) return;
    const n = +$('#net-n').value, c = CUT[+$('#net-cut').value];
    const qAdj = ISO ? new Map(B.partners.filter((x) => x.id !== P.key && !x.rep && sp.byKey.has(x.id)).map((x) => [sp.byKey.get(x.id).i, { best: x.best, avg: x.avg, src: x.src }]))   // this isoform's partners
      : E.adj.get(P.i) || new Map();
    const nb = [...qAdj].filter(([j, e]) => j !== P.i && e.best >= c).sort((a, b) => b[1].best - a[1].best).slice(0, n);   // the pair with itself is drawn as a ring, not as a partner
    const homoOf = (i) => { const e = (E.adj.get(i) || new Map()).get(i); return !!(e && e.best >= c); };
    netBox.innerHTML = '<svg></svg>';
    if (!nb.length) { netBox.innerHTML = '<div class="empty">No partners past this cutoff.</div>'; return; }
    const nodes = [{ id: P.i, row: P, q: true }, ...nb.map(([j, e]) => ({ id: j, row: sp.rows[j], e }))];
    const links = nb.map(([j, e]) => ({ source: P.i, target: j, best: e.best, avg: e.avg, q: true }));
    for (let a = 1; a < nodes.length; a++) { const mm = E.adj.get(nodes[a].id); if (!mm) continue;
      for (let b = a + 1; b < nodes.length; b++) { const e = mm.get(nodes[b].id); if (e && e.best >= c) links.push({ source: nodes[a].id, target: nodes[b].id, best: e.best, avg: e.avg }); } }
    const deg = new Map(); links.forEach((l) => { deg.set(l.source, (deg.get(l.source) || 0) + 1); deg.set(l.target, (deg.get(l.target) || 0) + 1); l.pubs = K ? K.pubs(l.source, l.target) : 0; l.gen = K ? K.gen(l.source, l.target) : 0; });
    if (!K) { $('#net-kb').checked = false; $('#net-kb').disabled = true; $('#net-ev').disabled = true; }
    const svg = d3.select(netBox).select('svg'), W = netBox.clientWidth, H = netBox.clientHeight;
    svg.attr('width', W).attr('height', H);
    const g = svg.append('g');
    const zoom = d3.zoom().scaleExtent([0.2, 4]).on('zoom', (ev) => g.attr('transform', ev.transform));
    svg.call(zoom);
    const link = g.append('g').selectAll('line').data(links).join('line').attr('stroke-width', (d) => EWID(d.avg)).attr('stroke-linecap', 'round');
    link.filter((d) => !d.q).raise();
    const restyle = () => { const st = edgeStyle('net', K);   // recolor in place: no new layout
      link.attr('stroke', st.color).attr('stroke-opacity', (d) => (st.hit(d) ? 0.95 : d.q ? 0.55 : 0.8)); link.filter(st.hit).raise();   // reported pairs on top
      $('#net-kbkey').innerHTML = edgeKey(st, K, links); };
    restyle();
    link.on('mousemove', (ev, d) => showTip(`<b>${esc(sp.rows[typeof d.source === 'object' ? d.source.id : d.source].gene)}</b> × <b>${esc(sp.rows[typeof d.target === 'object' ? d.target.id : d.target].gene)}</b> · iLIS ${ONE ? d.best.toFixed(3) : `best ${d.best.toFixed(3)} · average ${d.avg.toFixed(3)}`}${d.pubs || d.gen ? ` · reported in BioGRID (${kbPubs(d.pubs, d.gen)})` : ''}`, ev.clientX, ev.clientY)).on('mouseleave', hideTip);
    const r = (d) => (d.q ? 18 : 6 + Math.min(10, Math.sqrt(deg.get(d.id) || 1) * 1.6));
    const node = g.append('g').selectAll('g').data(nodes).join('g').style('cursor', 'pointer')
      .call(d3.drag().on('start', (ev, d) => { if (!ev.active) sim.alphaTarget(0.25).restart(); d.fx = d.x; d.fy = d.y; })
        .on('drag', (ev, d) => { d.fx = ev.x; d.fy = ev.y; }).on('end', (ev, d) => { if (!ev.active) sim.alphaTarget(0); if (!d.q) { d.fx = null; d.fy = null; } }));
    const homo = new Set(nodes.filter((d) => homoOf(d.id)).map((d) => d.id));
    const circle = node.append('circle').attr('r', r).attr('stroke', (d) => (homo.has(d.id) ? HOMO_RING : '#fff')).attr('stroke-width', (d) => (homo.has(d.id) ? 3.2 : d.q ? 2.5 : 1.8));
    node.append('text').text((d) => d.row.gene).attr('text-anchor', 'middle').attr('dy', (d) => -r(d) - 6)   // every label, the center one too: dark text on a white edge
      .attr('font-family', 'IBM Plex Sans, sans-serif').attr('font-size', (d) => (d.q ? 14 : 11)).attr('font-weight', (d) => (d.q ? 700 : 600)).attr('fill', '#17263A')
      .attr('paint-order', 'stroke').attr('stroke', 'rgba(255,255,255,0.92)').attr('stroke-width', (d) => (d.q ? 4 : 3)).attr('stroke-linejoin', 'round');
    node.filter((d) => d.q).raise();
    const placeLabels = liftLabels(g, node);
    NET = { link, restyle, recolor() {
      const k = M ? M.k : 1;
      circle.attr('fill', (d) => { if (d.q) return '#1A5276'; const cl = partnerCluster.get(d.row.key); return cl ? clusterColor(cl, k) : '#C3CCD6'; });
      const used = [...new Set(nodes.filter((d) => !d.q).map((d) => partnerCluster.get(d.row.key)).filter(Boolean))].sort((a, b) => a - b);
      $('#net-legend').innerHTML = (used.length ? used.map((cl) => `<span><i style="background:${clusterColor(cl, k)};border-radius:50%"></i>${clusterLabel(cl)}</span>`).join('') + '<span><i style="background:#C3CCD6;border-radius:50%"></i>not clustered at this cutoff</span>' : '') + (homo.size ? homoKey(homo.size) : '');
    } };
    NET.recolor();
    node.on('mousemove', (ev, d) => showTip(d.q ? `<b>${esc(d.row.gene)}</b> · ${fmtInt(d.row.pos10)} partners past 10% FPR` : `<b>${esc(d.row.gene)}</b> · iLIS ${ONE ? d.e.best.toFixed(3) : `best ${d.e.best.toFixed(3)} · average ${d.e.avg.toFixed(3)}`} with ${esc(P.gene)}<br>${srcBadges(sp, d.e.src)}`, ev.clientX, ev.clientY))
      .on('mouseleave', hideTip)
      .on('click', (ev, d) => { if (!d.q && !ev.defaultPrevented) location.hash = `#/${sp.id}/${d.row.key}${scopeQ}`; });
    const q = nodes[0]; q.fx = W / 2; q.fy = H / 2;
    const R = Math.min(W, H) * 0.42;
    const sim = d3.forceSimulation(nodes).force('link', d3.forceLink(links).id((d) => d.id).distance((l) => (l.q ? R * (1.15 - 0.55 * Math.min(1, l.best)) : 60)).strength((l) => (l.q ? 0.5 : 0.35)))
      .force('charge', d3.forceManyBody().strength(-340)).force('collide', d3.forceCollide().radius((d) => r(d) + 14)).force('x', d3.forceX(W / 2).strength(0.04)).force('y', d3.forceY(H / 2).strength(0.05))
      .on('tick', () => { placeLabels(); for (const sel of [link]) sel.attr('x1', (d) => d.source.x).attr('y1', (d) => d.source.y).attr('x2', (d) => d.target.x).attr('y2', (d) => d.target.y); node.attr('transform', (d) => `translate(${d.x},${d.y})`); });
    let fitted = false;   // once the layout settles, zoom so every node and label fits the box (the zoom stays free afterwards)
    sim.on('end', () => { if (fitted) return; fitted = true;
      const xs = nodes.map((d) => d.x), ys = nodes.map((d) => d.y);
      const x0 = Math.min(...xs) - 48, x1 = Math.max(...xs) + 48, y0 = Math.min(...ys) - 34, y1 = Math.max(...ys) + 24;
      const s = Math.min(1, 0.96 * Math.min(W / (x1 - x0), H / (y1 - y0)));
      svg.transition().duration(450).call(zoom.transform, d3.zoomIdentity.translate(W / 2 - s * (x0 + x1) / 2, H / 2 - s * (y0 + y1) / 2).scale(s)); });
    svgExport($('#net-x'), `atlas_${P.gene}_network`, () => $('svg', netBox));
  }

  drawTopList(); drawTable(); fillPartners(); drawScatter(); drawFreq(); drawHeatmap(); renderClusterInfo(); legend3D();
  loadStructure();
  cluster();
  let rsz; window.onresize = () => { clearTimeout(rsz); rsz = setTimeout(() => { if (clustered()) renderSites(); drawFreq(); drawHeatmap(); drawScatter(); const rb = $('#res-body'); if (rb && rb._redraw) rb._redraw(); }, 150); };
}

/* ── pair page ───────────────────────────────────────────────────────────────────────────────────────── */
// A pair reported in BioGRID: a badge in the pair's header that jumps to a card listing each publication (PubMed title,
// authors and year from NCBI where the record is a PubMed id), from the same BioGRID release as the counts elsewhere.
async function pairRefs(sp, P, O, gone) {
  const R = sp.byKey.get(O.key); if (P.i == null || !R || R.i == null) return;
  let K; try { K = await reportedOf(sp, P.i); } catch (e) { return; }
  if (gone() || !K) return;
  const ph = K.pubs(P.i, R.i), ge = K.gen(P.i, R.i), refs = K.refs && K.refs(P.i, R.i); if (!ph && !ge) return;
  const head = app.querySelector('.phead .srcs');
  if (head) head.insertAdjacentHTML('afterend', `<a class="kb-badge" href="#c-refs" title="${esc(kbTip0(K, ph, ge))}">Reported in BioGRID ${esc(K.release)} · ${[ph ? `physical, ${ph} publication${ph === 1 ? '' : 's'}` : '', ge ? `genetic, ${ge} publication${ge === 1 ? '' : 's'}` : ''].filter(Boolean).join(' · ')} ↓</a>`);
  const badge = app.querySelector('.kb-badge'); if (badge) badge.onclick = (e) => { e.preventDefault(); $('#c-refs').scrollIntoView({ behavior: 'auto', block: 'start' }); };
  const rows = refs ? [...refs.p.map((x) => ['physical', x]), ...refs.g.map((x) => ['genetic', x])] : [];
  const link = (x) => { const [k, v] = [x.slice(0, x.indexOf(':')), x.slice(x.indexOf(':') + 1)];
    return k === 'PUBMED' ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${esc(v)}/" target="_blank" rel="noopener">PMID ${esc(v)} ↗</a>` : k === 'DOI' ? `<a href="https://doi.org/${esc(v)}" target="_blank" rel="noopener">doi:${esc(v)} ↗</a>` : esc(x); };
  app.insertAdjacentHTML('beforeend', `<div class="card" id="c-refs"><div class="card-head"><h2>Reported interactions <span class="muted">BioGRID ${esc(K.release)}</span></h2>
    <a class="btn" href="https://thebiogrid.org/" target="_blank" rel="noopener">BioGRID ↗</a></div>
    <p class="muted" style="margin:2px 0 10px">Publications that BioGRID lists for ${esc(P.gene)} and ${esc(O.gene)}: a reported interaction, not a validation of this prediction.</p>
    ${rows.length ? `<div class="tbl-wrap"><table class="pt ref-tbl"><thead><tr><th>Evidence</th><th>Publication</th><th>Title</th><th>Authors</th><th class="n">Year</th></tr></thead><tbody>${rows.map(([ty, x]) =>
      `<tr data-ref="${esc(x)}"${rows.length > 25 && rows.indexOf(rows.find((z) => z[1] === x)) >= 25 ? ' hidden' : ''}><td>${ty}</td><td>${link(x)}</td><td class="t muted">…</td><td class="a"></td><td class="n y"></td></tr>`).join('')}</tbody></table></div>
      ${rows.length > 25 ? `<div class="pager"><button class="more" type="button" id="refs-all">show all ${fmtInt(rows.length)}</button></div>` : ''}`
      : `<p class="muted">${fmtInt(ph + ge)} publication${ph + ge === 1 ? '' : 's'}; this page's file does not list them. Search the pair on BioGRID.</p>`}</div>`);
  const all = $('#refs-all'); if (all) all.onclick = () => { app.querySelectorAll('#c-refs tr[hidden]').forEach((t) => { t.hidden = false; }); all.remove(); };
  const ids = [...new Set(rows.map(([, x]) => x).filter((x) => x.startsWith('PUBMED:')).map((x) => x.slice(7)))];
  for (let i = 0; i < ids.length; i += 150) {   // NCBI E-utilities esummary: title, first author, year
    let d; try { d = await (await fetch(`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=${ids.slice(i, i + 150).join(',')}`)).json(); } catch (e) { d = null; }
    if (gone()) return;
    for (const tr of app.querySelectorAll('#c-refs tr[data-ref^="PUBMED:"]')) { const r = d && d.result && d.result[tr.dataset.ref.slice(7)]; if (!r && d) continue;
      tr.querySelector('.t').textContent = r ? r.title || '' : ''; tr.querySelector('.t').classList.remove('muted');
      if (r) { const au = r.authors || []; tr.querySelector('.a').textContent = au.length ? au[0].name + (au.length > 1 ? ' et al.' : '') : ''; tr.querySelector('.y').textContent = (r.pubdate || '').slice(0, 4); } }
  }
  for (const td of app.querySelectorAll('#c-refs td.t.muted')) td.textContent = '';
}
const kbTip0 = (K, ph, ge) => `reported in BioGRID ${K.release}: ${[ph ? `physical, ${ph} publication${ph === 1 ? '' : 's'}` : '', ge ? `genetic, ${ge} publication${ge === 1 ? '' : 's'}` : ''].filter(Boolean).join('; ')}`;
async function viewPair(spId, q1, q2, setId = '') {   // setId: the scope the pair was opened from (a screen or a thematic set)
  const gen = ROUTE, sp = await species(spId), P = resolveRow(sp, q1), O0 = resolveRow(sp, q2);
  if (stale(gen)) return;
  if (!P) { app.innerHTML = `<div class="empty">No protein “${esc(q1)}”.</div>`; return; }
  let scope = '', label = '';
  if (setId && sp.dsIds.includes(setId)) { scope = setId; label = (await regDataset(setId)).short; }
  else if (setId) for (const id of sp.dsIds) { const TS = await setsOf(await dataset(id)).catch(() => null), S = TS && TS.byId.get(setId); if (S) { scope = setId; label = S.short; break; } }
  if (stale(gen)) return;
  const scopeQ = scope ? `?set=${encodeURIComponent(scope)}` : '';
  if (P.key !== q1 || (O0 && O0.key !== q2)) { location.replace(`#/${sp.id}/${P.key}/${O0 ? O0.key : q2}${scopeQ}`); return; }
  const O = O0 || { key: q2, gene: q2, name: '', acc: '', clen: 0, len: 0, occ: [], id: q2 };
  document.title = `${P.gene} · ${O.gene} · LIVIA Atlas`;
  const crumbs = `<div class="crumbs"><a href="#/">Atlas</a> / <a href="#/${sp.id}">${esc(sp.reg.label)}</a> / <a href="#/${sp.id}/${P.key}${scopeQ}">${esc(P.gene)}</a> / ${O0 ? `<a href="#/${sp.id}/${O.key}">${esc(O.gene)}</a>` : esc(O.gene)}</div>`;
  app.innerHTML = crumbs + '<div class="loading">Loading…</div>';
  let B;
  try { B = await merged(sp, P, scope, true); } catch (e) { B = { partners: [] }; }   // every model of the pair: all of a split gene's files
  if (stale(gen)) return;
  const part = B.partners.find((p) => p.id === O.key);
  if (!part) { app.innerHTML = crumbs + `<div class="empty">${esc(P.gene)} and ${esc(O.gene)} were not predicted together${scope ? ` in ${esc(label)}. <a href="#/${sp.id}/${P.key}/${O.key}">Every screen</a>` : ' in these screens'}.</div>`; return; }
  const best = [...part.preds].sort((a, b) => b.iLIS - a.iLIS)[0], b = bandOf(part.best), oLen = O.clen || part.preds[0].pLen, one = part.counted.length <= 1;   // one model: best = average
  const who = (R, len, col) => `<div class="who"><b style="color:${col}">${esc(R.gene)}</b> <span>${esc(short(R.name) || '')}</span>
      <div class="ids">${R.acc ? uniprotLink(R.acc) : '<span>no UniProt entry</span>'}<span>${esc(R.id)}</span>${len ? `<span>${fmtInt(len)} aa</span>` : ''}</div></div>`;
  const num = (v, d) => `<td class="n">${fmtNum(v, d)}</td>`;
  const band = (cuts, v, d) => `<td class="n" style="${bandSty(cuts, v)}">${fmtNum(v, d)}</td>`;   // a score in the color of the FPR band it passes
  const lab = (p) => runLabel(sp, B, p.run, P) + (p.rep ? ` · repeat of ${runLabel(sp, B, B.runs.get(p.run).repeatOf, P)}` : '');
  app.innerHTML = `${crumbs}
    <div class="phead"><div><h1><span style="color:var(--query-t)">${esc(P.gene)}</span> <span style="color:var(--ink-3);font-weight:600">×</span> <span style="color:var(--partner-t)">${esc(O.gene)}</span></h1>
        <div class="pairwho">${who(P, P.clen, 'var(--query-t)')}${who(O, oLen, 'var(--partner-t)')}</div>
        <div class="srcs">Predicted in ${part.sets.length ? setBadges(part.sets) : srcBadges(sp, part.src)}${(() => { const t = srcPapers(sp, part.src); return t ? ` <span class="muted">· ${t}</span>` : ''; })()}${scope ? ` <span class="muted">· only ${esc(label)} models shown · <a href="#/${sp.id}/${P.key}/${O.key}">every model</a></span>` : ''}</div>${(() => { const t = overlapNote(sp, B, part.preds, P); return t ? `<div class="srcs muted ovl">${esc(t)}</div>` : ''; })()}
        <div class="actions"><a class="btn" href="#/${sp.id}/${P.key}">${esc(P.gene)} page</a>${O0 ? `<a class="btn" href="#/${sp.id}/${O.key}">${esc(O.gene)} page</a>` : ''}${citeBtn(`${P.gene} × ${O.gene} (${sp.reg.label})`, sp.dsIds.filter((_, i) => part.src & (1 << i)))}<span id="pair-struct"></span></div></div>
      <div class="kpis"><div class="kpi"><b style="color:${BAND_TXT[b]}">${part.best.toFixed(3)}</b><span>iLIS${one ? '' : ' best'} · ${bandLabel[b]}</span></div>
        ${one ? '' : `<div class="kpi"><b style="color:${bandCol(FPR_AVG.iLIS, part.avg)}">${part.avg.toFixed(3)}</b><span>iLIS average · ${bandLabel[bandIn(FPR_AVG.iLIS, part.avg)]}</span></div>`}
        <div class="kpi"><b style="color:${bandCol(FPR.ipTM, part.iptmBest)}">${part.iptmBest.toFixed(2)}</b><span>ipTM${one ? '' : ' best'} · ${bandLabel[bandIn(FPR.ipTM, part.iptmBest)]}</span></div>
        ${one ? '' : `<div class="kpi"><b style="color:${bandCol(FPR_AVG.ipTM, part.iptmAvg)}">${part.iptmAvg.toFixed(2)}</b><span>ipTM average · ${bandLabel[bandIn(FPR_AVG.ipTM, part.iptmAvg)]}</span></div>`}</div></div>
    <div class="card"><div class="card-head"><h2>Ranked models</h2><span class="muted">every model of every screen · rank = the prediction's own model order, ipTM-based, not the iLIS order · click one to show its interaction residues${part.preds.some((p) => p.rep) ? ' · a repeat folded the same two sequences again: listed, not counted' : ''}</span></div>
      <p class="legend-text">Scores are colored by the false-positive-rate band they pass, each metric by its own benchmarked cutoffs (listed on the About page).</p><div class="legend" style="margin:0 0 10px">
        ${[1, 5, 10, 0].map((f) => `<span><i style="background:${BAND[f]}"></i>${f ? f + '% FPR' : 'below 10% FPR'}</span>`).join('')}</div>
      <div class="tbl-wrap"><table class="pt models"><thead>
        <tr><th rowspan="2">Source</th><th rowspan="2" title="The model's rank within its prediction: the predictor's own ipTM-based order, not the iLIS order">Rank</th><th rowspan="2" class="n">iLIS</th><th rowspan="2" class="n">iLISA</th><th rowspan="2" class="n">ipTM</th><th rowspan="2" class="n">LIS</th><th rowspan="2" class="n">cLIS</th>
          <th colspan="2" class="grp" title="PAE ≤ 12 Å">Local interaction residues (LIR)</th><th colspan="2" class="grp" title="the interaction residues that also have Cβ ≤ 8 Å">Contact residues (cLIR)</th></tr>
        <tr><th class="n sub q">${esc(P.gene)}</th><th class="n sub p">${esc(O.gene)}</th><th class="n sub q">${esc(P.gene)}</th><th class="n sub p">${esc(O.gene)}</th></tr></thead>
        <tbody>${part.preds.map((p, i) => `<tr data-i="${i}"><td><span class="src" style="--c:${runColor(sp, p)}">${esc(lab(p))}</span></td><td>${p.rank}</td>
          <td class="n">${fmtNum(p.iLIS, 3)} <span class="band b${bandOf(p.iLIS)}">${bandLabel[bandOf(p.iLIS)]}</span></td>
          ${band(FPR.iLISA, p.iLISA, 1)}${band(FPR.ipTM, p.ipTM, 2)}${band(FPR.LIS, p.LIS, 3)}${band(FPR.cLIS, p.cLIS, 3)}${num(p.qLIR, 0)}${num(p.pLIR, 0)}${num(p.qcLIR, 0)}${num(p.pcLIR, 0)}</tr>`).join('')}</tbody></table></div></div>
    <div class="card"><div class="card-head"><h2>Interaction Residues</h2><select id="model-pick" aria-label="Model" style="font:13px var(--sans);padding:5px 8px;border:1px solid var(--line);border-radius:7px">${part.preds.map((p, i) =>
        `<option value="${i}">${esc(lab(p))} · rank ${p.rank} · iLIS ${fmtNum(p.iLIS, 3)}</option>`).join('')}</select></div>
      <div class="legend" style="margin:2px 0 12px"><span><i style="background:#E0E0E0"></i>not an interaction residue</span><span><i style="background:#80CBC4"></i><i style="background:#FFAB91;margin-left:-2px"></i>interaction residue (LIR: PAE ≤ 12 Å)</span><span><i style="background:#00897B"></i><i style="background:#E64A19;margin-left:-2px"></i>contact (cLIR: also Cβ ≤ 8 Å)</span></div>
      <div id="iface"></div></div>`;
  const pick = (i) => {
    $('#model-pick').value = String(i);
    app.querySelectorAll('.models tbody tr').forEach((x) => x.classList.toggle('on', +x.dataset.i === i));
    ifaceView($('#iface'), { sp, P, O, pred: part.preds[i], B, canvasId: 'iface-canvas' });
  };
  $('#model-pick').onchange = (e) => pick(+e.target.value);
  app.querySelectorAll('.models tbody tr').forEach((tr) => tr.onclick = () => pick(+tr.dataset.i));
  pick(part.preds.indexOf(best));
  pairRefs(sp, P, O, () => stale(gen));
  if (sp.viruses) { const R2 = sp.byKey.get(O.key); virusStruct(sp, P.i, R2 ? R2.i : null).then((x) => { if (!stale(gen)) structLink($('#pair-struct'), x, 'btn'); }); }   // a virus pair: its model in LIVIA
  let rsz; window.onresize = () => { clearTimeout(rsz); rsz = setTimeout(() => { const f = $('#iface'); if (f && f._redraw) f._redraw(); }, 150); };
}

/* species page: the merged index — its screens, counts and most connected proteins */
async function viewSpecies(spId) {
  const gen = ROUTE, sp = await species(spId), c = sp.manifest.counts, hubs = [...sp.rows].sort((a, b) => b.pos10 - a.pos10).slice(0, 24);
  const TSs = await Promise.all(sp.dsIds.map(async (id) => { try { return await setsOf(await dataset(id)); } catch (e) { return null; } }));
  if (stale(gen)) return;
  document.title = `${sp.reg.label} · LIVIA Atlas`;
  app.innerHTML = `<div class="crumbs"><a href="#/">Atlas</a> / <a href="#/${sp.id}">${esc(sp.reg.label)}</a></div>
    <div class="dshead"><h1>${esc(sp.reg.heading || sp.reg.label + ' protein interactions')}</h1><div class="pname"><i>${esc(sp.reg.name)}</i> · ${sp.dsIds.length === 1 ? 'one screen' : sp.dsIds.length + ' screens'}, one page per ${sp.manifest.keyedBy ? 'gene' : 'protein'}</div></div>
    ${kpiRow(c)}
    <div class="card"><h2>Search</h2><div id="sp-search" style="margin-top:10px"></div></div>
    ${sp.viruses ? `<div class="card" id="vir-card"><div class="card-head"><div><h2>Viruses</h2><div class="muted">${fmtInt(sp.viruses.length)} viruses; in each, every pair of its proteins was folded. Open one for its network.</div></div>
      <input type="search" id="vir-filter" placeholder="Filter by name or family" aria-label="Filter viruses by name, family, genus or species" style="width:220px"></div>
      <div class="tbl-wrap"><table class="pt" id="vir-t"></table></div><p class="muted" id="vir-note" style="margin:8px 0 0"></p></div>` : ''}
    <div class="card"><div class="card-head"><div><h2>Network of your ${sp.manifest.keyedBy ? 'genes' : 'proteins'} of interest</h2><div class="muted">name a few ${sp.manifest.keyedBy ? 'genes' : 'proteins'}; see the predicted pairs among them and, if you like, the partners they share</div></div>
      <a class="btn" href="#/${sp.id}/network">Build a network →</a></div></div>
    <div class="card"><h2>Screens</h2><div class="screens">${sp.manifest.datasets.map((d, di) => `<div class="screen"><span class="src" style="--c:${sp.dsColor[di]}">${esc(d.short)}</span>
      <div><a href="#/datasets/${d.id}"><b>${esc(d.title)}</b></a><div class="muted">${fmtInt(d.counts.proteins)} proteins · ${fmtInt(d.counts.pairs)} pairs · ${d.counts.runs ? `${fmtInt(d.counts.runs)} predictions · ` : ''}${fmtInt(d.counts.predictions)} models</div>
      <div class="cite">${d.url ? `<a href="${esc(d.url)}" target="_blank" rel="noopener">${esc(d.citation)} ↗</a>` : esc(d.citation)}</div></div></div>`).join('')}</div>
      ${sp.dsIds.length > 1 ? `<p class="muted" style="margin:10px 0 0">${fmtInt(c.pairsInSeveral)} pairs were predicted in more than one screen; their pages keep every model with its source.</p>` : ''}</div>
    <div class="card"><h2>Most connected proteins <span class="muted">partners past the 10% FPR cutoff, any screen</span></h2>
      <div class="chips">${hubs.map((r) => `<a class="chip" href="#/${sp.id}/${r.key}">${esc(r.gene)} <span class="num" style="color:var(--ink-3)">${fmtInt(r.pos10)}</span></a>`).join('')}</div></div>
    ${TSs.map(setsCard).join('')}`;
  mountSearch($('#sp-search'), { spId });
  if (sp.viruses) {   // the virus list: by family, then name (or any column clicked); filtered by name, family, genus or species
    const VS = { key: 'family', asc: true }, name = (v) => v.name.toLowerCase();
    const COLS = [['name', 'Virus', name], ['family', 'Family', (v) => (v.family || '\uffff').toLowerCase()], ['host', 'Host', (v) => (v.host || '\uffff').toLowerCase()], ['n', 'Proteins', (v) => v.n],
      ['pairs', 'Pairs folded', (v) => v.het + v.hom], ['p10', `Past 10% FPR${cutNote(10)}`, (v) => v.hpos + v.mpos], ['p5', `Past 5% FPR${cutNote(5)}`, (v) => v.pos5 || 0],
      ['p1', `Past 1% FPR${cutNote(1)}`, (v) => v.pos1 || 0], ['pwp', 'Proteins with a partner', (v) => v.pwp]];
    const draw = () => { const f = $('#vir-filter').value.trim().toLowerCase(), al = VALIAS[f] || [], [, , val] = COLS.find(([k]) => k === VS.key);
      const all = sp.viruses.filter((v) => !f || name(v).includes(f) || al.includes(name(v)) || [v.family, v.genus, v.species].some((x) => x && x.toLowerCase().includes(f)))
        .sort((a, b) => { const x = val(a), y = val(b), d = typeof x === 'string' ? x.localeCompare(y) : x - y; return (VS.asc ? d : -d) || name(a).localeCompare(name(b)); }), show = all.slice(0, 60);
      const n = (x) => (x == null ? '<span class="muted">—</span>' : fmtInt(x));
      $('#vir-t').innerHTML = `<thead><tr>${COLS.map(([k, l]) => `<th data-k="${k}" class="${['name', 'family', 'host'].includes(k) ? '' : 'n'}${VS.key === k ? ' sorted' + (VS.asc ? ' asc' : '') : ''}">${l}</th>`).join('')}</tr></thead><tbody>${show.map((v) =>
        `<tr><td class="g"><a href="#/${sp.id}/taxon/${v.taxid}">${esc(v.name)}</a></td><td title="${esc(v.genus ? 'genus ' + v.genus : '')}">${esc(v.family || '—')}</td><td class="nm">${esc(v.host || '—')}</td><td class="n">${fmtInt(v.n)}</td><td class="n">${fmtInt(v.het + v.hom)}</td><td class="n">${fmtInt(v.hpos + v.mpos)}</td><td class="n">${n(v.pos5)}</td><td class="n">${n(v.pos1)}</td><td class="n">${fmtInt(v.pwp)}</td></tr>`).join('')}</tbody>`;
      $('#vir-t').querySelectorAll('th').forEach((th) => th.onclick = () => { const k = th.dataset.k; VS.asc = VS.key === k ? !VS.asc : ['name', 'family', 'host'].includes(k); VS.key = k; draw(); });
      $('#vir-note').textContent = (all.length > show.length ? `${fmtInt(show.length)} of ${fmtInt(all.length)} shown; type to narrow.` : `${fmtInt(all.length)} shown.`)
        + ' Family and host from the ICTV Virus Metadata Resource (MSL40); click a column to sort.'; };
    $('#vir-filter').oninput = draw; draw();
  }
}

/* ── virus page: one virus of the viral dataset, whose proteins were folded all against all. Its network (every protein
   a node, pairs past the cutoff as edges colored by FPR band, a ring for a predicted homodimer), its pairs and its proteins.
   Reads only the species index and edge list on the site. ───────────────────────────────────────────────────────── */
async function viewVirus(spId, taxid) {
  const gen = ROUTE, sp = await species(spId), v = (sp.viruses || []).find((x) => String(x.taxid) === String(taxid));
  if (stale(gen)) return;
  if (!v) {   // a species taxid (NCBI): its strain here, or the list of them
    const kin = (sp.viruses || []).filter((x) => x.lineage.some((t) => String(t) === String(taxid)));   // any NCBI taxon between the family and the virus
    if (kin.length === 1) { location.replace(`#/${sp.id}/taxon/${kin[0].taxid}`); return; }
    app.innerHTML = kin.length ? `<div class="empty">Taxon ${esc(taxid)} holds ${fmtInt(kin.length)} viruses here: ${kin.slice(0, 60).map((x) => `<a href="#/${sp.id}/taxon/${x.taxid}">${esc(x.name)}</a>`).join(', ')}${kin.length > 60 ? ` and ${fmtInt(kin.length - 60)} more` : ''}.</div>`
      : `<div class="empty">No virus with taxon ${esc(taxid)} here. <a href="#/${sp.id}">All viruses</a></div>`; return; }
  document.title = `${v.name} · LIVIA Atlas`;
  // this virus's pairs at iLIS >= 0.223 (homodimers included) with every score lis.py wrote: data/species/virus/pairs/<taxid>.tsv
  let text; try { text = await getText(sp.base + `pairs/${v.taxid}.tsv`); } catch (e) { if (!stale(gen)) app.innerHTML = `<div class="empty">${esc(e.message)}</div>`; return; }
  if (stale(gen)) return;
  const q = hashPath().q; let cut = [10, 5, 1].includes(+q.get('cut')) ? +q.get('cut') : 5, topk = q.has('top') ? Math.max(0, +q.get('top') || 0) : 0;   // defaults: 5% FPR, every edge
  const R = (i) => sp.rows[i], lines = text.trim().split('\n'), head = lines[0].split('\t');
  const all = lines.slice(1).map((l) => { const t = l.split('\t'), m = {}, ad = {}; let model = '', shown = false; head.forEach((h, k) => { if (h === 'model') model = t[k] || ''; else if (h === 'afdb') shown = t[k] === '1';
    else if (ARCH_COLS.includes(h)) { if (t[k] !== '' && t[k] != null) ad[h] = +t[k]; } else if (k > 1) m[h] = t[k] === '' ? NaN : +t[k]; });
    const addr = ARCH_COLS.every((h) => Number.isFinite(ad[h])) ? ad : null;   // its place in the release archive, once its chunk is indexed
    return { a: +t[0], b: +t[1], best: m.iLIS, avg: m.iLIS, iptm: m.ipTM, m, model, shown, addr }; }).sort((x, y) => y.best - x.best);   // one model per pair: best = average
  const MCOL = head.slice(2).filter((h) => h !== 'model' && h !== 'afdb' && !ARCH_COLS.includes(h)),   // scores only: the archive address columns drive the LIVIA link
    mfmt = (k, x) => (!Number.isFinite(x) ? '–' : k === 'iLISA' ? x.toFixed(1) : k === 'ipTM' ? x.toFixed(2) : x.toFixed(3));
  let sortKey = 'iLIS', sortAsc = false;
  const lab = (r) => {   // a short node label: the gene, or the code that names a polyprotein product ("Serine protease NS3" → NS3)
    const g = r.gene; if (g.length <= 12) return g;
    const ns = g.match(/^Non[ -]?structural protein\s+(\w+)$/i); if (ns) return 'NS' + ns[1];
    const last = g.split(/\s+/).pop(); if (/^(NS\d+[A-Z]?|nsP\d|VP\d+[a-z]?|[A-Za-z]{1,3}\d*[A-Z]?|2[Kk]|\d[A-C])$/.test(last) && last.length <= 5) return last;
    return g.slice(0, 11) + '…'; };
  app.innerHTML = `<div class="crumbs"><a href="#/">Atlas</a> / <a href="#/${sp.id}">${esc(sp.reg.label)}</a> / <a href="#/${sp.id}/taxon/${v.taxid}">${esc(v.name)}</a></div>
    <div class="phead vh"><div><h1>${esc(v.name)}</h1>
      <div class="pname">every pair of its ${fmtInt(v.n)} proteins folded with AlphaFold-Multimer (one model each) and scored with lis.py</div>
      ${v.family ? `<div class="pname">${esc(v.family)}${v.genus ? ` · <i>${esc(v.genus)}</i>` : ''}${v.species ? ` · species <i>${esc(v.species)}</i>` : ''}${v.host ? ` · host: ${esc(v.host)}` : ''} <span class="muted">(ICTV VMR MSL40)</span></div>` : ''}
      ${(() => { const d = ((REG && REG.datasets) || []).find((x) => x.id === 'viral-dimers-afdb'); return d ? `<div class="pname vsrc">Predictions from the AlphaFold Database release of viral protein complexes (EMBL-EBI, Google DeepMind, NVIDIA and collaborators; models CC BY 4.0): <a href="${esc(d.paper)}" target="_blank" rel="noopener">Han, Narain et al. 2026 ↗</a> · data: <a href="https://doi.org/10.5281/zenodo.${recOf(d) || ''}" target="_blank" rel="noopener">LIVIA Atlas on Zenodo ↗</a></div>` : ''; })()}
      <div class="ids"><a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=${v.taxid}" target="_blank" rel="noopener">NCBI taxon ${v.taxid}</a><span>${fmtInt(v.het)} heterodimers · ${fmtInt(v.hom)} homodimers folded</span></div>
      <div class="actions">${citeBtn(`${v.name} (virus)`, sp.dsIds)}</div></div>
      <div class="kpis"><div class="kpi"><b>${fmtInt(v.n)}</b><span>proteins</span></div><div class="kpi f10"><b>${fmtInt(v.hpos)}</b><span>heterodimers past 10% FPR${cutNote(10)}</span></div>
        <div class="kpi"><b>${fmtInt(v.mpos)}</b><span>homodimers past 10% FPR${cutNote(10)}</span></div><div class="kpi"><b>${fmtInt(v.pwp)}</b><span>proteins with a partner</span></div></div></div>
    <div class="card" id="vn-card"><div class="card-head"><h2>Network</h2><div class="controls" style="margin:0"><label class="ctl" title="a large virus is easier to read with each protein's strongest pairs only; the pairs table lists every pair">Edges<select id="vn-top"><option value="0">all</option><option value="5">top 5 per protein</option><option value="3">top 3 per protein</option><option value="1">top 1 per protein</option></select></label>
        <label class="ctl" title="find groups of proteins predicted to bind each other more than the rest (Louvain communities, iLIS-weighted) and lay each group out apart"><input type="checkbox" id="vn-comm" checked> Group by community</label>
        <div class="ctl"><span>Cutoff</span><div class="seg" id="vn-cut">${[10, 5, 1].map((f) => `<button data-f="${f}" class="${f === cut ? 'on' : ''}">${f}% FPR · ${CUT[f].toFixed(3)}</button>`).join('')}</div></div></div></div>
      <p class="muted" style="margin:2px 0 12px">Every protein of the virus is a node; an edge joins two proteins whose pair passed the cutoff, shaded in gray by its iLIS (darker is higher), its width the iLIS. A black ring marks a protein predicted to form a homodimer; gray nodes have no partner at this cutoff. Grouped by community, proteins predicted to bind each other more than the rest sit together, one color per group. Click a protein for its page, an edge for the pair.</p>
      <div class="net" id="vn-net"></div><div class="legend" id="vn-legend"></div><div id="vn-x"></div></div>
    <div class="card"><div class="card-head"><h2>Pairs <span class="muted" id="vp-note"></span></h2><div class="controls" style="margin:0"><label class="ctl" title="list each protein's pair with itself as well"><input type="checkbox" id="vp-homo"> Homodimers</label><input type="search" id="vp-find" placeholder="Filter pairs" aria-label="Filter pairs by protein" style="width:170px"><button class="btn" id="vp-csv" type="button">↓ CSV</button></div></div>
      <p class="legend-text">iLIS and ipTM are colored by the false-positive-rate band they pass, each by its own benchmarked cutoffs (every cutoff is on the About page).</p><div class="legend" style="margin:0 0 10px">
        ${[1, 5, 10, 0].map((f) => `<span><i style="background:${BAND[f]}"></i>${f ? f + '% FPR' : 'below 10% FPR'}</span>`).join('')}<span class="muted">Click a column to sort.</span></div>
      <div class="tbl-wrap"><table class="pt compact" id="vp"></table></div><div class="pager" id="vp-more"></div></div>
    <div class="card"><div class="card-head"><h2>Proteins <span class="muted">${fmtInt(v.n)}</span></h2></div><div class="tbl-wrap"><table class="pt" id="vprot"></table></div></div>`;
  let shownPairs = 60;
  const pairHref = (x) => `#/${sp.id}/${R(x.a).key}/${R(x.b).key}`;
  function tables() {
    const c = CUT[cut], vq = (($('#vp-find') || {}).value || '').trim().toLowerCase(), hit = (x) => !vq || [R(x.a), R(x.b)].some((r) => [r.gene, r.name, r.key].some((t) => String(t || '').toLowerCase().includes(vq)));
    const P0 = all.filter((x) => x.best >= c), P = vq ? P0.filter(hit) : P0, deg = new Map(), homo = new Set();   // counts from every pair at the cutoff; the filter narrows the list only
    for (const x of P0) if (x.a === x.b) homo.add(x.a); else { deg.set(x.a, (deg.get(x.a) || 0) + 1); deg.set(x.b, (deg.get(x.b) || 0) + 1); }
    const withHomo = !!($('#vp-homo') || {}).checked, nh = P.filter((x) => x.a === x.b).length;   // homodimers listed only when asked for
    $('#vp-note').textContent = `${fmtInt(P.length)} past iLIS ${c} (${cut}% FPR), ${fmtInt(nh)} of them homodimers${withHomo || !nh ? '' : ', not listed'}`;
    const PS = [...(withHomo ? P : P.filter((x) => x.a !== x.b))].sort((x, y) => { const d = (x.m[sortKey] ?? -Infinity) - (y.m[sortKey] ?? -Infinity); return (sortAsc ? d : -d) || y.best - x.best; });
    $('#vp').innerHTML = `<thead><tr><th>Protein</th><th>Partner</th><th title="the interaction residues of the pair">Pair</th><th title="the model in LIVIA, for the pairs the AlphaFold Database displays (its filter: ipSAE ≥ 0.60 and pDockQ2 ≥ 0.23)">3D</th>${MCOL.map((k) => `<th class="n${k === sortKey ? ' sorted' + (sortAsc ? ' asc' : '') : ''}" data-k="${k}" style="cursor:pointer">${k}</th>`).join('')}</tr></thead><tbody>${PS.slice(0, shownPairs).map((x) => { const b = bandOf(x.best);
      return `<tr><td class="g"><a href="#/${sp.id}/${R(x.a).key}">${esc(R(x.a).gene)}</a></td><td class="g">${x.a === x.b ? '<span class="muted" title="homodimer: the protein with itself">self</span>' : `<a href="#/${sp.id}/${R(x.b).key}">${esc(R(x.b).gene)}</a>`}</td>
        <td><a href="${pairHref(x)}" title="interaction residues">residues</a></td>
        <td>${x.shown && x.model ? `<a href="${LIVIA}dimer.html?id=${encodeURIComponent(x.model)}" target="_blank" rel="noopener" title="${esc(x.model)} in LIVIA, from the AlphaFold Database">LIVIA ↗</a>`
          : x.addr && x.model ? `<a href="#" class="arch" data-m="${esc(x.model)}" title="${esc(x.model)} in LIVIA, read from the release archive at EBI (${((x.addr.cif_len + x.addr.pae_len) / 1048576).toFixed(1)} MB; not displayed by the AlphaFold Database)">LIVIA ↗</a>`
          : '<span class="muted" title="not displayed by the AlphaFold Database (below its filter); the structure opens from the release archive once its part is indexed">—</span>'}</td>
        ${MCOL.map((k) => (k === 'iLIS' ? `<td class="n v"><a href="${pairHref(x)}" style="color:${BAND_TXT[b]};font-weight:${BAND_W[b]}">${mfmt(k, x.m[k])}</a></td>` : k === 'ipTM' ? `<td class="n v" style="${bandSty(FPR.ipTM, x.m[k])}">${mfmt(k, x.m[k])}</td>` : `<td class="n">${mfmt(k, x.m[k])}</td>`)).join('')}</tr>`; }).join('')}</tbody>`;
    $('#vp').querySelectorAll('th[data-k]').forEach((th) => th.onclick = () => { const k = th.dataset.k; if (k === sortKey) sortAsc = !sortAsc; else { sortKey = k; sortAsc = false; } tables(); });
    $('#vp').querySelectorAll('a.arch').forEach((a) => a.onclick = (e) => { e.preventDefault(); const x = all.find((y) => y.model === a.dataset.m); if (x && x.addr) openFromArchive(x.model, x.addr, a); });
    $('#vp-more').innerHTML = PS.length > shownPairs ? `<button class="more" type="button">show ${fmtInt(Math.min(200, PS.length - shownPairs))} more</button>` : '';
    const mb = $('#vp-more button'); if (mb) mb.onclick = () => { shownPairs += 200; tables(); };
    const prot = [...v.members].sort((a, b) => (deg.get(b) || 0) - (deg.get(a) || 0) || R(a).gene.localeCompare(R(b).gene));
    $('#vprot').innerHTML = `<thead><tr><th>Protein</th><th>Name</th><th class="n">Length</th><th class="n">Partners in the virus</th><th>Homodimer</th></tr></thead><tbody>${prot.map((i) =>
      `<tr><td class="g"><a href="#/${sp.id}/${R(i).key}">${esc(R(i).gene)}</a></td><td class="nm" title="${esc(R(i).name)}">${esc(short(R(i).name))}</td><td class="n">${fmtInt(R(i).clen || R(i).len)}</td>
        <td class="n">${fmtInt(deg.get(i) || 0)}</td><td>${homo.has(i) ? 'predicted' : '<span class="muted">—</span>'}</td></tr>`).join('')}</tbody>`;
    return { P, deg, homo };
  }
  function draw() {
    const { P, deg, homo } = tables(), box = $('#vn-net'); box.innerHTML = '<svg></svg>';
    const W = box.clientWidth, H = box.clientHeight, svg = d3.select(box).select('svg').attr('width', W).attr('height', H), g = svg.append('g');
    const zoom = d3.zoom().scaleExtent([0.1, 8]).on('zoom', (ev) => { g.attr('transform', ev.transform); relabel(ev.transform.k); }); svg.call(zoom);
    const het = P.filter((x) => x.a !== x.b), k = topk; $('#vn-top').value = String(k);
    let shown = het;
    if (k) { const keep = new Set(), per = new Map(); for (const x of het) { for (const i of [x.a, x.b]) { if (!per.has(i)) per.set(i, []); per.get(i).push(x); } }
      for (const xs of per.values()) xs.sort((a, b) => b.best - a.best).slice(0, k).forEach((x) => keep.add(x)); shown = het.filter((x) => keep.has(x)); }   // an edge stays if it is among either end's k strongest
    const nodes = v.members.map((i) => ({ id: i, row: R(i) })), links = shown.map((x) => ({ source: x.a, target: x.b, best: x.best, avg: x.avg, x }));
    const r = (d) => 5 + Math.min(9, Math.sqrt(deg.get(d.id) || 0) * 1.7);
    // communities on every pair past the cutoff (not only the edges drawn); groups of three or more get a color and a place of their own
    const idx = new Map(nodes.map((d, n) => [d.id, n])), cm = louvain(nodes.length, het.map((x) => [idx.get(x.a), idx.get(x.b), x.best]));
    const csize = new Map(); nodes.forEach((d, n) => { d.c = cm[n]; csize.set(d.c, (csize.get(d.c) || 0) + 1); });
    const big = [...csize].filter(([c, sz]) => sz >= 3).map(([c]) => c).sort((a, b) => a - b), grouped = $('#vn-comm').checked;
    const ccol = (d) => (!deg.get(d.id) ? '#C3CCD6' : grouped && big.includes(d.c) && big.indexOf(d.c) < TAB10.length ? TAB10[big.indexOf(d.c)] : '#1A5276');
    const center = new Map();   // each group's place: groups on a ring around the largest, spaced by their size
    if (grouped) { const R0 = Math.min(W, H) * 0.42; big.forEach((c, n) => { if (n === 0) center.set(c, [W / 2, H / 2]); else { const t = (n - 1) / Math.max(1, big.length - 1) * 2 * Math.PI; center.set(c, [W / 2 + R0 * Math.cos(t), H / 2 + R0 * 0.8 * Math.sin(t)]); } }); }
    const home = (d) => center.get(d.c) || [W / 2, H / 2];
    if (grouped) nodes.forEach((d, n) => { const [cx, cy] = home(d); d.x = cx + 30 * Math.cos(n); d.y = cy + 30 * Math.sin(n); });   // each node starts at its group's place
    const rank = new Map([...nodes].sort((a, b) => (deg.get(b.id) || 0) - (deg.get(a.id) || 0) || (homo.has(b.id) - homo.has(a.id))).map((d, n) => [d.id, n]));   // label priority: most partners first
    const link = g.append('g').selectAll('line').data(links).join('line').attr('stroke', (d) => EGRAY(d.best)).attr('stroke-width', (d) => EWID(d.avg)).attr('stroke-opacity', 0.85).attr('stroke-linecap', 'round').style('cursor', 'pointer');
    link.on('mousemove', (ev, d) => showTip(`<b>${esc(d.x && R(d.x.a).gene)}</b> × <b>${esc(d.x && R(d.x.b).gene)}</b> · iLIS ${d.best.toFixed(3)} · ${bandLabel[bandOf(d.best)]}${Number.isFinite(d.x.iptm) ? ` · ipTM ${d.x.iptm.toFixed(2)}` : ''}<br>click for the pair`, ev.clientX, ev.clientY))
      .on('mouseleave', hideTip).on('click', (ev, d) => { hideTip(); location.hash = pairHref(d.x); });
    const node = g.append('g').selectAll('g').data(nodes).join('g').style('cursor', 'pointer')
      .call(d3.drag().on('start', (ev, d) => { if (!ev.active) sim.alphaTarget(0.25).restart(); d.fx = d.x; d.fy = d.y; })
        .on('drag', (ev, d) => { d.fx = ev.x; d.fy = ev.y; }).on('end', (ev, d) => { if (!ev.active) sim.alphaTarget(0); d.fx = null; d.fy = null; }));
    node.append('circle').attr('r', r).attr('fill', ccol).attr('stroke', (d) => (homo.has(d.id) ? HOMO_RING : '#fff')).attr('stroke-width', (d) => (homo.has(d.id) ? 3 : 1.5));
    const label = node.append('text').text((d) => lab(d.row)).attr('text-anchor', 'middle')
      .attr('font-family', 'IBM Plex Sans, sans-serif').attr('font-weight', 600).attr('fill', '#17263A')
      .attr('paint-order', 'stroke').attr('stroke', 'rgba(255,255,255,0.92)').attr('stroke-linejoin', 'round');
    const placeLabels = liftLabels(g, node);
    function relabel(z) {   // zoomed in, more labels; each at the same size on screen. About 40 labels at the fitted view of a large virus, every label for a small one
      const budget = nodes.length <= 60 ? Infinity : 40 * z * z / Math.max(0.2, fitK * fitK);
      label.attr('display', (d) => (rank.get(d.id) < budget ? null : 'none')).attr('font-size', 10.5 / z).attr('stroke-width', 3 / z).attr('dy', (d) => -r(d) - 5 / z);
    }
    let fitK = 1; relabel(1);
    node.on('mousemove', (ev, d) => showTip(`<b>${esc(d.row.gene)}</b>${d.row.name && d.row.name !== d.row.gene ? ` · ${esc(short(d.row.name))}` : ''}<br>${fmtInt(deg.get(d.id) || 0)} partner${(deg.get(d.id) || 0) === 1 ? '' : 's'} in the virus at this cutoff${homo.has(d.id) ? ' · predicted homodimer' : ''}`, ev.clientX, ev.clientY))
      .on('mouseleave', hideTip).on('click', (ev, d) => { if (!ev.defaultPrevented) location.hash = `#/${sp.id}/${d.row.key}`; });
    const same = (l) => (typeof l.source === 'object' ? l.source.c === l.target.c : true);
    const sim = d3.forceSimulation(nodes).force('link', d3.forceLink(links).id((d) => d.id).distance((l) => 60 + 60 * (1 - Math.min(1, l.best))).strength((l) => (grouped && !same(l) ? 0.02 : 0.5)))
      .force('charge', d3.forceManyBody().strength(nodes.length > 150 ? -90 : -200)).force('collide', d3.forceCollide().radius((d) => r(d) + 8))
      .force('x', d3.forceX((d) => home(d)[0]).strength((d) => (grouped && center.has(d.c) ? 0.4 : 0.07))).force('y', d3.forceY((d) => home(d)[1]).strength((d) => (grouped && center.has(d.c) ? 0.4 : 0.09)))
      .on('tick', () => { placeLabels(); link.attr('x1', (d) => d.source.x).attr('y1', (d) => d.source.y).attr('x2', (d) => d.target.x).attr('y2', (d) => d.target.y); node.attr('transform', (d) => `translate(${d.x},${d.y})`); });
    let fitted = false;
    sim.on('end', () => { if (fitted) return; fitted = true; const xs = nodes.map((d) => d.x), ys = nodes.map((d) => d.y), pad = 40;
      const x0 = Math.min(...xs) - pad, x1 = Math.max(...xs) + pad, y0 = Math.min(...ys) - pad - 12, y1 = Math.max(...ys) + pad, k = Math.min(2, W / (x1 - x0), H / (y1 - y0));
      fitK = k; svg.transition().duration(400).call(zoom.transform, d3.zoomIdentity.translate(W / 2 - k * (x0 + x1) / 2, H / 2 - k * (y0 + y1) / 2).scale(k)); });
    $('#vn-legend').innerHTML = (k ? `<span class="muted">${fmtInt(links.length)} of ${fmtInt(het.length)} edges drawn: each protein's ${k} strongest</span>` : '')
      + '<span><i style="background:linear-gradient(90deg, #C5CCD4, #1E2A38);height:4px;width:60px;border-radius:2px"></i>iLIS, 0.223 to 0.85+</span>'
      + (grouped && big.length ? big.slice(0, TAB10.length).map((c, n) => `<span><i style="background:${TAB10[n]};border-radius:50%"></i>community ${n + 1} (${fmtInt(csize.get(c))})</span>`).join('') + (big.length > TAB10.length ? `<span class="muted">+${big.length - TAB10.length} smaller communities in dark blue</span>` : '') : '')
      + `${homoKey(homo.size)}<span><i style="background:#C3CCD6;border-radius:50%"></i>no partner at this cutoff (${fmtInt(v.members.filter((i) => !deg.get(i)).length)})</span>`;
    svgExport($('#vn-x'), `atlas_virus_${v.taxid}_network`, () => $('svg', box));
  }
  $('#vn-cut').onclick = (e) => { const f = e.target.dataset.f; if (!f) return; cut = +f; [...$('#vn-cut').children].forEach((b) => b.classList.toggle('on', b.dataset.f === f));
    const [path, qs] = location.hash.split('?'), u = new URLSearchParams(qs || ''); cut === 5 ? u.delete('cut') : u.set('cut', cut); history.replaceState(null, '', `${path}${u.toString() ? '?' + u : ''}`); shownPairs = 60; draw(); };
  $('#vn-comm').onchange = () => draw();
  $('#vn-top').onchange = (e) => { topk = +e.target.value; const [path, qs] = location.hash.split('?'), u = new URLSearchParams(qs || ''); topk ? u.set('top', topk) : u.delete('top'); history.replaceState(null, '', `${path}${u.toString() ? '?' + u : ''}`); draw(); };
  $('#vp-homo').onchange = () => { shownPairs = 60; tables(); };
  $('#vp-find').oninput = () => { shownPairs = 60; tables(); };   // filter by either protein's gene, name or accession
  $('#vp-csv').onclick = () => { const csvq = (x) => (/[",\n]/.test(x) ? `"${String(x).replace(/"/g, '""')}"` : x), c = CUT[cut];
    const text = `Protein_1,Protein_2,accession_1,accession_2,${MCOL.join(',')},band,AFDB_model,displayed_by_AFDB\n` + all.filter((x) => x.best >= c).map((x) => [csvq(R(x.a).gene), csvq(R(x.b).gene), R(x.a).acc || R(x.a).key, R(x.b).acc || R(x.b).key, ...MCOL.map((k) => (Number.isFinite(x.m[k]) ? x.m[k] : '')), bandLabel[bandOf(x.best)], x.model, x.shown ? 'yes' : 'no'].join(',')).join('\n') + '\n';
    const u = URL.createObjectURL(new Blob([text], { type: 'text/csv' })), a = document.createElement('a'); a.href = u; a.download = `atlas_virus_${v.taxid}_pairs.csv`; a.click(); setTimeout(() => URL.revokeObjectURL(u), 3000); };
  draw();
}

/* ── network builder: the proteins a reader names and the predicted pairs among them. It reads only the species index
   and edge list already on the site (no bundle); a protein or an edge opens its page. ─────────────────────────────── */
/* ── nested network: baits, then the candidates that join round by round. Round 1 takes a candidate with a predicted pair
   to any bait past the cutoff; from round 2 a candidate needs at least k connections into the network accepted so far.
   Every candidate carries a p-value: the chance of that many connections by chance, P(X ≥ k | m, p_c), where m is the
   accepted network's size (interim: the coverage index that gives the partners each candidate was tested against is not
   built yet) and p_c the candidate's own share of partners past the cutoff (proteins.json, the same best-over-models iLIS
   the edges use), shrunk toward the benchmark rate so a protein with few partners is not over-read. ───────────────── */
function binomTail(m, p, k) {   // P(X ≥ k), X ~ Binomial(m, p)
  if (k <= 0) return 1; if (m < k || p <= 0) return 0; if (p >= 1) return 1;
  let s = 0, t = Math.pow(1 - p, m); for (let i = 0; i < k; i++) { s += t; t *= ((m - i) / (i + 1)) * (p / (1 - p)); } return Math.min(1, Math.max(0, 1 - s));
}
function bhQ(ps) {   // Benjamini–Hochberg q-values, in the input order
  const o = ps.map((p, i) => [p, i]).sort((a, b) => a[0] - b[0]), q = new Array(ps.length); let run = 1;
  for (let r = o.length - 1; r >= 0; r--) { run = Math.min(run, (o[r][0] * o.length) / (r + 1)); q[o[r][1]] = run; } return q;
}
function nestedRun(E, baits, cands, o) {   // o: { c: cutoff, f: its FPR in %, k, rounds, strict, pc(i), m1(i) } → rounds, groups and each candidate's support
  // mt (optional): how many of the accepted network a candidate was folded with (the tested-pair index), every round;
  // a candidate folded with fewer than the rule needs cannot be judged. m1 (optional, without mt): the baits a candidate was folded with (from the baits' own predictions). Round 1 then tests only those,
  // with m = that count; a candidate folded with no bait cannot be judged in round 1. Without it, m is the network's size.
  const acc = new Map(baits.map((i) => [i, 0])), info = new Map(), rounds = [], nb = (i) => E.adj.get(i) || new Map();
  let cumFP = 0;
  for (let r = 1; r <= o.rounds; r++) {
    const kReq = r === 1 ? 1 : o.k, m = acc.size, rest = cands.filter((i) => !acc.has(i)), rows = []; let untested = 0;
    for (const i of rest) { const mi = o.mt ? o.mt(i, acc) : r === 1 && o.m1 ? o.m1(i) : m; if (mi < kReq) { untested++; continue; }   // too few tested partners in the network to reach the rule: cannot be judged
      const sup = [...nb(i)].filter(([j, e]) => acc.has(j) && e.best >= o.c).map(([j]) => j), p0 = o.pc(i);
      const nsup = o.fam ? new Set(sup.map(o.fam)).size : sup.length;   // paralogs counted once: one per family
      rows.push({ i, sup, n: nsup, m: mi, p: binomTail(mi, p0, Math.max(1, nsup)), e: binomTail(mi, p0, kReq) }); }
    const q = bhQ(rows.map((x) => x.p)); rows.forEach((x, n) => { x.q = q[n]; });
    const pass = rows.filter((x) => x.n >= kReq && (!o.strict || x.q <= 0.05)), eFP = rows.reduce((s, x) => s + x.e, 0); cumFP += eFP;
    const hidden = r === 1 ? [] : rows.filter((x) => x.n >= 1 && x.n < kReq);
    for (const x of [...pass, ...hidden]) info.set(x.i, { ...x, r, hidden: x.n < kReq || (o.strict && x.q > 0.05) });
    for (const x of pass) acc.set(x.i, r);
    rounds.push({ r, m, kReq, tested: rows.length, untested: o.mt || (r === 1 && o.m1) ? untested : null, accepted: pass.length, hidden: hidden.length, eFP, cumFP });
    if (!pass.length) break;
  }
  return { acc, info, rounds };
}

async function viewNested(spId, q) {
  const gen = ROUTE, sp = await species(spId); if (stale(gen)) return;
  document.title = `Nested network · ${sp.reg.label} · LIVIA Atlas`;
  const S = { baits: q.get('baits') || '', ids: q.get('ids') || '', cut: [10, 5, 1].includes(+q.get('cut')) ? +q.get('cut') : 1, k: [1, 2, 3].includes(+q.get('k')) ? +q.get('k') : 2,
    rounds: [1, 2, 3].includes(+q.get('rounds')) ? +q.get('rounds') : 3, strict: q.get('strict') === '1', para: q.get('para') === '1', hidden: false, col: null };
  const G = sp.manifest.keyedBy ? 'genes' : 'proteins';
  app.innerHTML = `<div class="crumbs"><a href="#/">Atlas</a> / <a href="#/${sp.id}">${esc(sp.reg.label)}</a> / <a href="#/${sp.id}/network">Network</a> / <a href="${esc(location.hash)}">Nested</a></div>
    <div class="dshead"><h1>Nested network</h1><div class="pname">${esc(sp.reg.label)} · baits, then the candidates predicted to join them, round by round</div></div>
    <div class="card"><div class="card-head"><h2>Baits and candidates</h2><span class="muted">names as in the network builder: symbols, accessions, older names, or a table</span></div>
      <label class="nlab">Baits<textarea class="ids" id="ns-baits" rows="2" spellcheck="false" placeholder="one or a few ${G}">${esc(S.baits.split(',').join(', '))}</textarea></label>
      <label class="nlab">Candidates <span class="muted">(IP-MS preys, screen hits, GWAS or proteomics hits)</span><textarea class="ids" id="ns-ids" rows="4" spellcheck="false">${esc(S.ids.split(',').join(', '))}</textarea></label>
      <div class="controls" style="margin-top:8px"><button class="btn" id="ns-filebtn" type="button">Load candidates from a file</button><input type="file" id="ns-file" accept=".txt,.csv,.tsv,.tab,text/plain,text/csv,text/tab-separated-values" hidden>
        <span class="ex-row" id="ns-ex"></span></div>
      <div class="controls" style="margin-top:10px"><div class="ctl"><span>Cutoff</span><div class="seg" id="ns-cut">${[10, 5, 1].map((f) => `<button data-f="${f}" class="${f === S.cut ? 'on' : ''}">${f}% FPR · ${CUT[f].toFixed(3)}</button>`).join('')}</div></div>
        <label title="from round 2, a candidate joins with at least this many connections into the network accepted so far">At least <select id="ns-k">${[1, 2, 3].map((k) => `<option value="${k}"${k === S.k ? ' selected' : ''}>${k}</option>`).join('')}</select> connections from round 2</label>
        <label>Rounds <select id="ns-rounds">${[1, 2, 3].map((k) => `<option value="${k}"${k === S.rounds ? ' selected' : ''}>up to ${k}</option>`).join('')}</select></label>
        <label title="accept only candidates whose chance p-value passes Benjamini–Hochberg q ≤ 0.05 within their round"><input type="checkbox" id="ns-strict"${S.strict ? ' checked' : ''}> Strict (BH q ≤ 0.05)</label>
        <label id="ns-para-wrap" title="two paralogs in the network count as one connection (MMseqs2 families, 30% identity over half of both sequences)"><input type="checkbox" id="ns-para"${S.para ? ' checked' : ''}> Count paralogs once</label>
        <button class="btn" id="ns-go" type="button">Build the nested network</button></div>
      <p class="muted" id="ns-status" style="margin:10px 0 0"></p></div>
    <div class="card" id="ns-card" hidden><div class="card-head"><h2>Rounds</h2></div><div class="tbl-wrap"><table class="pt compact" id="ns-rounds-t"></table></div><p class="muted" id="ns-base" style="margin:8px 0 0"></p>
      <p class="legend-text">p: the chance of at least that many connections by chance, from the candidate's own share of partners past the cutoff and the accepted network's size (m). Expected false positives: the sum of those chances over every candidate tested in the round; since each candidate's share of partners past the cutoff includes its real partners, this errs high (compare the random networks below). The cutoffs are benchmarked on the top-ranked model's iLIS while an edge takes the best iLIS over every model, so the nominal false positive rate is a lower bound. m is how many of the accepted network the candidate was folded with (from the index of tested pairs; for a species without it, the baits in round 1 and the network's size after).</p></div>
    <div class="card" id="ns-ncard" hidden><div class="card-head"><h2>Network</h2><div class="controls" style="margin:0"><label><input type="checkbox" id="ns-hidden"> show hidden candidates</label>
      <button class="btn" id="ns-csv" type="button">↓ CSV</button><button class="btn" id="ns-graphml" type="button">↓ GraphML</button></div></div>
      <p class="muted" style="margin:2px 0 10px">Columns by round: baits at the left, then each round's accepted ${G}. An edge is a predicted pair past the cutoff, colored by the round of the later ${G.slice(0, -1)}. Click a ${G.slice(0, -1)} for its page, an edge for the pair.</p>
      <div class="net" id="ns-net"></div><div class="legend" id="ns-legend"></div><div id="ns-x"></div></div>
    <div class="card" id="ns-tcard" hidden><div class="card-head"><h2>Accepted ${G}</h2></div><div class="tbl-wrap"><table class="pt compact" id="ns-t"></table></div></div>`;
  const status = (t) => { $('#ns-status').innerHTML = t; };
  const RC = ['#1A5276', '#E67E22', '#27AE60', '#8E44AD'], gname = (i) => sp.rows[i].gene;
  let res = null, last = null;
  $('#ns-cut').onclick = (e) => { const f = e.target.dataset.f; if (!f) return; S.cut = +f; [...$('#ns-cut').children].forEach((b) => b.classList.toggle('on', b.dataset.f === f)); if (res) run(); };
  if (!(sp.manifest.files || {}).paralogs) $('#ns-para-wrap').style.display = 'none';
  ['#ns-k', '#ns-rounds', '#ns-strict', '#ns-para'].forEach((s) => { $(s).onchange = () => { if (res) run(); }; });
  $('#ns-hidden').onchange = () => { if (last) drawNet(last); };
  $('#ns-go').onclick = () => run();
  $('#ns-filebtn').onclick = () => $('#ns-file').click();
  $('#ns-file').onchange = async (e) => { const f = e.target.files[0]; if (!f) return; $('#ns-ids').value = await f.text(); e.target.value = ''; run(); };
  examplePicker($('#ns-ex'), sp.id, (t) => {
    const R = readIdInput(sp, t), bait = R.T ? R.T.rows.find((r) => r.includes('bait')) : null;   // an example with a bait names it in the role column
    if (!bait) { status('This example has no bait; name one in the Baits box, then Build.'); $('#ns-ids').value = t; return; }
    const bn = [bait[R.col], ...bait].find((x) => x && resolveHow(sp, String(x).trim(), true)) || ''; $('#ns-baits').value = bn; $('#ns-ids').value = t.split('\n').filter((l) => !l.split('\t').includes('bait')).join('\n'); run(); });
  async function run() {
    const Rb = readIdInput(sp, $('#ns-baits').value), Rc = readIdInput(sp, $('#ns-ids').value, S.col); S.col = Rc.col;
    const res1 = (toks) => { const ok = [], miss = []; for (const t of toks) { const h = resolveHow(sp, t, true); if (h) { if (!ok.includes(h.row.i)) ok.push(h.row.i); } else miss.push(t); } return { ok, miss }; };
    const B = res1(Rb.toks), C = res1(Rc.toks); S.k = +$('#ns-k').value; S.rounds = +$('#ns-rounds').value; S.strict = $('#ns-strict').checked; S.para = $('#ns-para').checked;
    const cands = C.ok.filter((i) => !B.ok.includes(i));
    const qs = new URLSearchParams({ baits: B.ok.map((i) => sp.rows[i].key).join(','), cut: S.cut, k: S.k, rounds: S.rounds }); if (S.strict) qs.set('strict', '1'); if (S.para) qs.set('para', '1'); if (Rc.toks.length <= 400) qs.set('ids', Rc.toks.join(','));
    history.replaceState(null, '', `#/${sp.id}/nested?${qs}`);
    if (!B.ok.length) { status(`Name at least one bait the ${esc(sp.reg.label)} screens hold${B.miss.length ? ` (not found: ${esc(B.miss.join(', '))})` : ''}.`); return; }
    if (!cands.length) { status('Give the candidates to test against the baits.'); return; }
    status('Reading the edge list…'); let E, K; try { [E, K] = await Promise.all([edges(sp), reported(sp)]); } catch (e) { status(esc(e.message)); return; }
    if (stale(gen)) return;
    const c = CUT[S.cut], f = S.cut / 100, posKey = { 10: 'pos10', 5: 'pos5', 1: 'pos1' }[S.cut], A = 20;
    const pc = (i) => { const r = sp.rows[i]; return r.partners ? ((r[posKey] || 0) + A * f) / (r.partners + A) : f; };   // shrunk toward the benchmark rate
    const [TP, PA] = await Promise.all([tested(sp), S.para ? paralogs(sp) : null]); if (stale(gen)) return;
    const fam = PA ? PA.of : null;
    const mt = TP ? (i, acc) => { if (!fam) { let n = 0; for (const j of acc.keys()) if (TP.has(i, j)) n++; return n; }
      const s = new Set(); for (const j of acc.keys()) if (TP.has(i, j)) s.add(fam(j)); return s.size; } : null;   // families, when paralogs count once
    if (!TP) status('Reading the baits’ predictions (which candidates each was folded with)…');
    const tb = TP ? [] : await Promise.all(B.ok.map(async (b) => { try { const all = await merged(sp, sp.rows[b], '', true); return new Set(all.preds.map((x) => sp.byKey.get(x.partner)).filter(Boolean).map((x) => x.i)); } catch (e) { return null; } }));
    if (stale(gen)) return;
    const m1 = !TP && tb.every(Boolean) ? (i) => tb.filter((s) => s.has(i)).length : null;   // a bait whose predictions did not load: fall back to the network's size
    const o = { c, f, k: S.k, rounds: S.rounds, strict: S.strict, pc, m1, mt, fam };
    res = nestedRun(E, B.ok, cands, o);
    // random networks of this size: the same candidates and rules, baits drawn at random from the proteins with predictions
    const pool = sp.rows.filter((r) => r.partners > 0 && !cands.includes(r.i)).map((r) => r.i), base = [];
    const ob = { ...o, m1: null };   // random baits: the tested-pair index when there is one; else their predictions are not read and round 1 uses the network's size
    for (let t = 0; t < 200; t++) { const rb = new Set(); while (rb.size < B.ok.length && rb.size < pool.length) rb.add(pool[Math.floor(Math.random() * pool.length)]);
      base.push(nestedRun(E, [...rb], cands, ob).acc.size - rb.size); }
    base.sort((a, b) => a - b);
    const nAcc = res.acc.size - B.ok.length, never = cands.filter((i) => (sp.rows[i].partners || 0) < S.k).length;
    const kbShare = (r) => { const ids = [...res.acc].filter(([, rr]) => rr === r).map(([i]) => i); if (!K || !ids.length) return '';
      const before = new Set([...res.acc].filter(([, rr]) => rr < r).map(([i]) => i)); const n = ids.filter((i) => [...before].some((j) => K.pubs(i, j) > 0)).length; return `${fmtInt(n)} of ${fmtInt(ids.length)}`; };
    $('#ns-card').hidden = false;
    $('#ns-rounds-t').innerHTML = `<thead><tr><th>Round</th><th class="n" title="the network accepted before this round">m</th><th>rule</th><th class="n">tested</th><th class="n" title="candidates folded with fewer of the network than the rule needs: not tested, not negatives">cannot be judged</th><th class="n">accepted</th><th class="n" title="at least one connection, but fewer than the rule needs">hidden</th><th class="n">expected false positives</th><th class="n">cumulative</th><th title="accepted ${G} with a pair reported in BioGRID (physical) into the network before the round; the file has no low- or high-throughput flag">BioGRID support</th></tr></thead><tbody>`
      + res.rounds.map((x) => `<tr><td>${x.r}</td><td class="n">${fmtInt(x.m)}</td><td>≥ ${x.kReq} connection${x.kReq === 1 ? '' : 's'}${S.strict ? ', q ≤ 0.05' : ''}</td><td class="n">${fmtInt(x.tested)}</td><td class="n">${x.untested == null ? '<span class="muted">–</span>' : fmtInt(x.untested)}</td><td class="n"><b>${fmtInt(x.accepted)}</b></td><td class="n">${fmtInt(x.hidden)}</td><td class="n">${x.eFP.toFixed(1)}</td><td class="n">${x.cumFP.toFixed(1)}</td><td>${kbShare(x.r) || '<span class="muted">–</span>'}</td></tr>`).join('') + '</tbody>';
    $('#ns-base').innerHTML = `Random networks of this size (${fmtInt(B.ok.length)} random bait${B.ok.length === 1 ? '' : 's'} drawn from the ${G} with predictions, the same candidates and rules, 200 runs) accept ${d3.mean(base).toFixed(1)} candidates on average (95th percentile ${fmtInt(Math.ceil(d3.quantile(base, 0.95)))}); this network accepts <b>${fmtInt(nAcc)}</b>.`
      + (never ? ` ${fmtInt(never)} candidate${never === 1 ? '' : 's'} can never reach ${S.k} connections: fewer than ${S.k} partners folded in all.` : '')
      + (TP ? ` Every round counts only what was folded: m is how many of the network each candidate was folded with (the index of tested pairs), and a candidate folded with fewer than the rule needs cannot be judged.`
        : ` Round 1 counts only the candidates folded with a bait (from the baits’ own predictions); later rounds need the index of tested pairs, which this species does not have yet, so m there is the network’s size.`);
    status(`${fmtInt(B.ok.length)} bait${B.ok.length === 1 ? '' : 's'} · ${fmtInt(cands.length)} candidates${fam ? ' · paralogs counted once' : ''} · ${fmtInt(nAcc)} accepted over ${res.rounds.filter((x) => x.accepted).length} round${res.rounds.filter((x) => x.accepted).length === 1 ? '' : 's'} at iLIS ${c} (${S.cut}% FPR)`
      + `${B.miss.length || C.miss.length ? ` · not found: ${esc([...B.miss, ...C.miss].slice(0, 30).join(', '))}` : ''}`);
    last = { E, K, B: B.ok, o }; drawNet(last); table();
  }
  function nodesOf(showHidden) {   // baits, the accepted, and (when asked) the hidden; at most 400, kept by p then connections
    let ids = [...res.acc.keys()]; if (showHidden) ids = ids.concat([...res.info].filter(([, x]) => x.hidden).map(([i]) => i));
    if (ids.length > 400) { const B = new Set(last.B); ids = [...B, ...ids.filter((i) => !B.has(i)).sort((a, b) => (res.info.get(a).p - res.info.get(b).p) || (res.info.get(b).n - res.info.get(a).n)).slice(0, 400 - B.size)]; }
    return ids;
  }
  function drawNet({ E, o }) {
    $('#ns-ncard').hidden = false; const showH = $('#ns-hidden').checked, ids = nodesOf(showH), set = new Set(ids);
    const col = (i) => (res.acc.has(i) ? res.acc.get(i) : res.info.get(i).r), box = $('#ns-net'); box.innerHTML = '<svg></svg>';
    const W = box.clientWidth, H = box.clientHeight, nc = Math.max(...ids.map(col)) + 1, svg = d3.select(box).select('svg').attr('width', W).attr('height', H), g = svg.append('g');
    svg.call(d3.zoom().scaleExtent([0.2, 6]).on('zoom', (ev) => g.attr('transform', ev.transform)));
    const links = []; for (const a of ids) for (const [b, e] of E.adj.get(a) || []) if (b > a && set.has(b) && e.best >= o.c) links.push({ a, b, e, r: Math.max(col(a), col(b)) });
    const byCol = d3.range(nc).map((k) => ids.filter((i) => col(i) === k)), y = new Map();
    byCol.forEach((xs) => xs.forEach((i, n) => y.set(i, (n + 1) / (xs.length + 1))));
    const nbs = new Map(ids.map((i) => [i, []])); for (const l of links) { nbs.get(l.a).push(l.b); nbs.get(l.b).push(l.a); }
    for (let sweep = 0; sweep < 4; sweep++) for (let k = 1; k < nc; k++) {   // order each column by the mean height of its links to earlier columns
      const xs = byCol[k].map((i) => { const up = nbs.get(i).filter((j) => col(j) < k); return [i, up.length ? d3.mean(up, (j) => y.get(j)) : y.get(i)]; }).sort((a, b) => a[1] - b[1]);
      xs.forEach(([i], n) => y.set(i, (n + 1) / (xs.length + 1))); byCol[k] = xs.map(([i]) => i); }
    const X = (i) => 70 + (col(i) / Math.max(1, nc - 1)) * (W - 180), Y = (i) => 20 + y.get(i) * (H - 40);
    g.append('g').selectAll('line').data(links).join('line').attr('x1', (d) => X(d.a)).attr('y1', (d) => Y(d.a)).attr('x2', (d) => X(d.b)).attr('y2', (d) => Y(d.b))
      .attr('stroke', (d) => RC[Math.min(d.r, RC.length - 1)]).attr('stroke-opacity', (d) => (col(d.a) === col(d.b) ? 0.25 : 0.55)).attr('stroke-width', (d) => EWID(d.e.avg)).style('cursor', 'pointer')
      .on('mousemove', (ev, d) => showTip(`<b>${esc(gname(d.a))}</b> × <b>${esc(gname(d.b))}</b> · iLIS ${d.e.best.toFixed(3)}`, ev.clientX, ev.clientY)).on('mouseleave', hideTip)
      .on('click', (ev, d) => { hideTip(); location.hash = `#/${sp.id}/${sp.rows[d.a].key}/${sp.rows[d.b].key}`; });
    const many = ids.length > 90, node = g.append('g').selectAll('g').data(ids).join('g').attr('transform', (i) => `translate(${X(i)},${Y(i)})`).style('cursor', 'pointer');
    node.append('circle').attr('r', (i) => (col(i) === 0 ? 8 : 5.5)).attr('fill', (i) => (res.info.get(i) && res.info.get(i).hidden ? '#fff' : RC[Math.min(col(i), RC.length - 1)]))
      .attr('stroke', (i) => RC[Math.min(col(i), RC.length - 1)]).attr('stroke-width', 1.6).attr('stroke-dasharray', (i) => (res.info.get(i) && res.info.get(i).hidden ? '2 2' : null));
    node.filter((i) => !many || col(i) === 0 || (res.info.get(i) && res.info.get(i).p < 0.01)).append('text').text(gname).attr('x', 9).attr('dy', '0.32em').attr('font-size', 10.5)
      .attr('font-family', 'IBM Plex Sans, sans-serif').attr('font-weight', (i) => (col(i) === 0 ? 700 : 600)).attr('fill', '#17263A').attr('paint-order', 'stroke').attr('stroke', 'rgba(255,255,255,0.92)').attr('stroke-width', 3);
    node.on('mousemove', (ev, i) => { const x = res.info.get(i); showTip(`<b>${esc(gname(i))}</b> · ${col(i) === 0 ? 'bait' : `${x.hidden ? 'hidden, ' : ''}round ${x.r}`}${x ? `<br>${fmtInt(x.n)} connection${x.n === 1 ? '' : 's'} into the network · p ${x.p.toExponential(1)} · q ${x.q.toFixed(3)}` : ''}`, ev.clientX, ev.clientY); })
      .on('mouseleave', hideTip).on('click', (ev, i) => { location.hash = `#/${sp.id}/${sp.rows[i].key}`; });
    $('#ns-legend').innerHTML = RC.slice(0, nc).map((cc, k) => `<span><i style="background:${cc};border-radius:50%"></i>${k === 0 ? 'baits' : `round ${k}`} (${fmtInt(byCol[k].filter((i) => !(res.info.get(i) || {}).hidden).length)})</span>`).join('')
      + (showH ? '<span><i style="background:#fff;border:1.5px dashed #5B6573;border-radius:50%"></i>hidden: below the rule</span>' : '') + (ids.length >= 400 ? '<span class="muted">capped at 400, kept by p, then connections</span>' : '');
    svgExport($('#ns-x'), `atlas_${sp.id}_nested`, () => $('svg', box));
  }
  function rowsOut() { return [...res.acc.keys(), ...[...res.info].filter(([, x]) => x.hidden).map(([i]) => i)].map((i) => { const x = res.info.get(i), r = res.acc.has(i) ? res.acc.get(i) : x.r;
    return { i, grp: r === 0 ? 'bait' : x.hidden ? 'hidden' : `round ${r}`, r, x }; }); }
  function table() {
    $('#ns-tcard').hidden = false;
    $('#ns-t').innerHTML = `<thead><tr><th>${G.slice(0, -1)}</th><th>group</th><th class="n">connections</th><th>supported by</th><th class="n" title="the chance of at least that many connections by chance">p</th><th class="n">q</th><th class="n" title="the chance of meeting this round's rule by chance; summed over the round's candidates it gives the expected false positives">chance of the rule</th></tr></thead><tbody>`
      + rowsOut().filter((d) => d.grp !== 'hidden').sort((a, b) => (a.r - b.r) || ((a.x ? a.x.p : 0) - (b.x ? b.x.p : 0))).map((d) => `<tr><td class="g"><a href="#/${sp.id}/${sp.rows[d.i].key}">${esc(gname(d.i))}</a></td><td>${d.grp}</td><td class="n">${d.x ? fmtInt(d.x.n) : '–'}</td>
        <td>${d.x ? esc(d.x.sup.map(gname).slice(0, 12).join(', ')) + (d.x.sup.length > 12 ? ` +${d.x.sup.length - 12}` : '') : '<span class="muted">–</span>'}</td><td class="n">${d.x ? d.x.p.toExponential(1) : '–'}</td><td class="n">${d.x ? d.x.q.toFixed(3) : '–'}</td><td class="n">${d.x ? d.x.e.toExponential(1) : '–'}</td></tr>`).join('') + '</tbody>';
  }
  const dl = (text, name, type) => { const u = URL.createObjectURL(new Blob([text], { type })), a = document.createElement('a'); a.href = u; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(u), 3000); };
  $('#ns-csv').onclick = () => { if (!res) return; const q2 = (v) => (/[",\n]/.test(v) ? `"${String(v).replace(/"/g, '""')}"` : v);
    dl('protein,key,group,round,connections,supported_by,p,q,chance_of_rule\n' + rowsOut().map((d) => [q2(gname(d.i)), sp.rows[d.i].key, d.grp, d.r, d.x ? d.x.n : '', q2(d.x ? d.x.sup.map(gname).join(' ') : ''), d.x ? d.x.p : '', d.x ? d.x.q : '', d.x ? d.x.e : ''].join(',')).join('\n') + '\n', `atlas_${sp.id}_nested.csv`, 'text/csv'); };
  $('#ns-graphml').onclick = () => { if (!res || !last) return; const x = (v) => String(v).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
    const R = rowsOut(), set = new Set(R.map((d) => d.i)), col = new Map(R.map((d) => [d.i, d.r])), E = last.E;
    const keys = [['gene', 'node', 'gene', 'string'], ['key', 'node', 'key', 'string'], ['group', 'node', 'group', 'string'], ['round', 'node', 'round', 'int'], ['conn', 'node', 'connections', 'int'],
      ['p', 'node', 'p', 'double'], ['q', 'node', 'q', 'double'], ['ilis', 'edge', 'iLIS_best', 'double'], ['eround', 'edge', 'round', 'int']];
    const nodes = R.map((d) => `<node id="n${d.i}"><data key="gene">${x(gname(d.i))}</data><data key="key">${x(sp.rows[d.i].key)}</data><data key="group">${d.grp}</data><data key="round">${d.r}</data>${d.x ? `<data key="conn">${d.x.n}</data><data key="p">${d.x.p}</data><data key="q">${d.x.q}</data>` : ''}</node>`);
    const eds = []; for (const a of set) for (const [b, e] of E.adj.get(a) || []) if (b > a && set.has(b) && e.best >= last.o.c) eds.push(`<edge source="n${a}" target="n${b}"><data key="ilis">${e.best}</data><data key="eround">${Math.max(col.get(a), col.get(b))}</data></edge>`);
    dl(`<?xml version="1.0" encoding="UTF-8"?>\n<graphml xmlns="http://graphml.graphdrawing.org/xmlns">\n${keys.map(([id, f, nm, t]) => `<key id="${id}" for="${f}" attr.name="${nm}" attr.type="${t}"/>`).join('\n')}\n<graph id="atlas_${sp.id}_nested" edgedefault="undirected">\n${nodes.join('\n')}\n${eds.join('\n')}\n</graph>\n</graphml>\n`, `atlas_${sp.id}_nested.graphml`, 'application/xml'); };
  if (S.baits && S.ids) run(); else $('#ns-baits').focus();
}

// The network builder's example tables: published IP-MS hit lists (the paper's own cutoff), under CC BY 4.0
// Example tables for the network views (data/examples/index.json): published hit lists and gene sets under CC BY or CC0,
// each at its paper's own cutoff; a picker lists a species' examples, and the row under it cites the one chosen.
let EXIDX = null;
const examplesFor = async (spId) => { if (!EXIDX) EXIDX = getJSON('data/examples/index.json').catch(() => []); return (await EXIDX).filter((x) => x.species === spId); };
function examplePicker(host, spId, onLoad) {   // host: an element; onLoad(text, entry) when the reader loads one
  examplesFor(spId).then((xs) => { if (!xs.length || !host.isConnected) return;
    host.innerHTML = `<label>Example <select class="ex-pick" aria-label="Example table">${xs.map((x, n) => `<option value="${n}">${esc(x.what)} · ${esc(x.type)}</option>`).join('')}</select></label>
      <button class="btn ex-load" type="button">Load</button><span class="muted ex-cite"></span>`;
    const sel = host.querySelector('.ex-pick'), cite = () => { const x = xs[+sel.value];
      host.querySelector('.ex-cite').innerHTML = `<a href="https://doi.org/${esc(x.doi)}" target="_blank" rel="noopener">${esc(x.cite)}</a>, ${esc(x.license)} · <a href="data/examples/${esc(x.file)}" download>download</a>`; };
    sel.onchange = cite; cite();
    host.querySelector('.ex-load').onclick = async () => { const x = xs[+sel.value]; let t; try { t = await getText(`data/examples/${x.file}`); } catch (e) { return; } onLoad(t, x); }; });
}
async function viewNetwork(spId, q) {
  const gen = ROUTE, sp = await species(spId);
  const TSs = await Promise.all(sp.dsIds.map(async (id) => { try { return await setsOf(await dataset(id)); } catch (e) { return null; } }));
  if (stale(gen)) return;
  document.title = `Network · ${sp.reg.label} · LIVIA Atlas`;
  const scopes = [['', sp.dsIds.length > 1 ? 'every screen' : sp.dsShort[0]], ...(sp.dsIds.length > 1 ? sp.manifest.datasets.map((d) => [d.id, d.short]) : []),
    ...TSs.filter(Boolean).flatMap((T) => T.list.map((x) => [x.id, x.short]))];
  const S = { ids: q.get('ids') || '', add: ['none', 'shared', 'top'].includes(q.get('add')) ? q.get('add') : 'none', k: Math.max(1, Math.min(50, +q.get('k') || 5)),
    cut: q.get('cut') === 'c' ? 'c' : [10, 5, 1].includes(+q.get('cut')) ? +q.get('cut') : 10, cutv: Math.min(1, Math.max(CUT[10], +q.get('cutv') || 0.4)), iptm: Math.min(1, Math.max(0, +q.get('iptm') || 0)),
    set: scopes.some(([id]) => id === (q.get('set') || '')) ? q.get('set') || '' : '',
    mind: [0, 1, 2, 3].includes(+q.get('mind')) ? +q.get('mind') : 0, grp: q.get('grp') || '', exp: (q.get('exp') || '').split(',').filter(Boolean), click: q.get('click') === 'open' ? 'open' : 'add', col: null, data: [], ncol: '' };   // exp: proteins expanded by a click (keys), in the order clicked
  const eg = [...sp.rows].sort((a, b) => b.pos10 - a.pos10).slice(0, 5).map((r) => r.gene).join(', ');
  app.innerHTML = `<div class="crumbs"><a href="#/">Atlas</a> / <a href="#/${sp.id}">${esc(sp.reg.label)}</a> / <a href="${esc(location.hash)}">Network</a></div>
    <div class="dshead"><h1>Network of your ${sp.manifest.keyedBy ? 'genes' : 'proteins'} of interest</h1><div class="pname"><a href="#/${sp.id}/nested">Nested network ↗</a> · baits and candidates, accepted round by round</div><div class="pname">${esc(sp.reg.label)} · the predicted pairs among the ${sp.manifest.keyedBy ? 'genes' : 'proteins'} you name</div></div>
    <div class="card"><div class="card-head"><h2>Proteins</h2><span class="muted">gene symbols, UniProt accessions (isoforms too)${sp.manifest.keyedBy ? ', FlyBase IDs, CG numbers' : ''} or older names · commas, spaces or new lines, or a table</span></div>
      <textarea class="ids" id="nw-ids" rows="3" spellcheck="false" placeholder="for example: ${esc(eg)}">${esc(S.ids.split(',').join(', '))}</textarea>
      <div class="controls" style="margin-top:8px"><button class="btn" id="nw-filebtn" type="button" title="a list or a table of names: txt, csv or tsv; a table's name column is found for you">Load a file</button><input type="file" id="nw-file" accept=".txt,.csv,.tsv,.tab,text/plain,text/csv,text/tab-separated-values" hidden>
        <span class="ex-row" id="nw-ex"></span>
        <label id="nw-col-wrap" hidden>Names in <select id="nw-col"></select></label><span class="muted" id="nw-table"></span></div>
      <div class="controls" style="margin-top:10px"><label>Show <select id="nw-add"><option value="none">only these proteins</option><option value="shared">+ partners that two or more share</option><option value="top">+ each one's top partners</option></select></label>
        <label title="hide proteins with fewer predicted pairs than this in the drawing (repeated until every protein left has at least this many)">Min. pairs <select id="nw-mind">${[0, 1, 2, 3].map((n) => `<option value="${n}"${n === S.mind ? ' selected' : ''}>${n ? `${n}+` : 'any'}</option>`).join('')}</select></label>
        <span id="nw-k-wrap"><input type="number" id="nw-k" min="1" max="50" value="${S.k}" style="width:56px" aria-label="partners per protein"> per protein</span>
        <div class="ctl"><span>Cutoff</span><div class="seg" id="nw-cut">${[10, 5, 1].map((f) => `<button data-f="${f}" class="${f === S.cut ? 'on' : ''}">${f}% FPR · ${CUT[f].toFixed(3)}</button>`).join('')}<button data-f="c" class="${S.cut === 'c' ? 'on' : ''}" title="an iLIS cutoff of your own">custom</button></div>
          <input type="number" id="nw-cutv" min="${CUT[10]}" max="1" step="0.01" value="${S.cutv}" style="width:70px${S.cut === 'c' ? '' : ';display:none'}" aria-label="custom iLIS cutoff" title="iLIS from ${CUT[10]} (the edge list holds the pairs past 10% FPR) to 1"></div>
        <label title="also require this ipTM (best over models); blank or 0 for none">ipTM ≥ <input type="number" id="nw-iptm" min="0" max="1" step="0.05" value="${S.iptm || ''}" placeholder="any" style="width:64px"></label>
        <label>In <select id="nw-set">${scopes.map(([id, l]) => `<option value="${esc(id)}"${id === S.set ? ' selected' : ''}>${esc(l)}</option>`).join('')}</select></label>
        <button class="btn" id="nw-go" type="button">Draw the network</button></div>
      <p class="muted" id="nw-status" style="margin:10px 0 0"></p></div>
    <div class="card" id="nw-card" hidden><div class="card-head"><h2>Network</h2><div class="controls" style="margin:0">${edgeCtl('nw')}
        <label class="ctl" id="nw-ncol-wrap" style="display:none" title="color the proteins by a column of the table you gave">Color proteins by<select id="nw-ncol"></select></label>
        <label class="ctl" title="lay the network out in groups, each in its own area with an outline">Group by<select id="nw-grp"></select></label>
        <label title="also draw pairs reported in BioGRID that were not predicted past the cutoff"><input type="checkbox" id="nw-unpred"> + reported, not predicted</label>
        <button class="btn" id="nw-check" type="button" style="display:none" title="for each dashed pair: was it folded and scored below the cutoff, or never folded in these screens? Reads each protein's predictions">Check which were folded</button>
        <button class="btn" id="nw-link" type="button" title="copy a link that opens this network">Copy link</button><button class="btn" id="nw-copyids" type="button" title="copy every protein in this network (yours and the added partners), comma separated">Copy proteins</button><button class="btn" id="nw-useids" type="button" title="put every protein in this network into the input box and draw it again as the proteins of interest">Use as input</button><button class="btn" id="nw-csv" type="button">↓ CSV</button><button class="btn" id="nw-graphml" type="button" title="the network for Cytoscape, Gephi or yEd: node group, edge iLIS, ipTM, screens and BioGRID publications">↓ GraphML</button>
        <button class="btn" id="nw-livia" type="button" title="the same network in LIVIA's network page: Leiden communities, layouts, Cytoscape export">Open in LIVIA Network ↗</button></div></div>
      <div class="controls" style="margin:2px 0 6px"><label class="ctl" title="find a protein in this network: it is centered and marked">Find <input type="search" id="nw-find" placeholder="a protein in the network" aria-label="Find a protein in the network" style="width:190px"></label><div class="ctl"><span>Click a protein to</span><div class="seg" id="nw-click"><button data-m="add" class="${S.click === 'add' ? 'on' : ''}" title="add its top partners past the cutoff, in place (the number per protein above); ⌘ or Ctrl-click opens its page in a new tab">add its partners</button><button data-m="open" class="${S.click === 'open' ? 'on' : ''}">open its page</button></div></div>
        <button class="btn" id="nw-unexp" type="button" style="display:none" title="remove the partners added by clicks">Undo added partners</button></div>
      <p class="muted" style="margin:2px 0 12px">${sp.manifest.keyedBy ? 'Genes' : 'Proteins'} of interest are large and dark, added partners small and light; a ring of dashes marks a protein you expanded. Click an edge for the interaction residues of the pair. Drag to move, scroll to zoom.</p>
      <div class="net" id="nw-net"></div>

      <div class="netkey"><div class="kbkey" id="nw-key"></div><div class="kbkey" id="nw-nkey"></div>
        <div><span>Edge width · ${sp.one ? 'iLIS' : 'average iLIS'}</span><svg id="nw-w" width="260" height="30" aria-hidden="true"></svg></div></div><div id="nw-x"></div>
<div id="nw-heat-wrap" hidden><div class="card-head" style="margin-top:14px"><h3 style="margin:0">Pairs as a matrix</h3><div class="controls" style="margin:0"><span class="muted" id="nw-heat-note"></span>
        <label class="ctl" title="every folded pair with its best iLIS (below the cutoff in grays), or only the pairs past the cutoff">Show <select id="nw-heat-show"><option value="all">every folded pair</option><option value="pos">pairs past the cutoff only</option></select></label>
        <label class="ctl">Colors <select id="nw-heat-cs">${['Blues', 'Viridis', 'YlGnBu', 'Reds', 'Grays', 'Cividis'].map((k) => `<option>${k}</option>`).join('')}</select></label></div></div>
        <div id="nw-heat" style="width:100%"></div><div class="legend" id="nw-heat-key"></div></div></div>`;
  $('#nw-add').value = S.add;
  const showK = () => { $('#nw-k-wrap').hidden = $('#nw-add').value !== 'top' && S.click !== 'add'; }; showK();   // the number also sets how many a click adds
  const redraw = () => { if (!$('#nw-card').hidden) draw(); };   // once a network is drawn, every option redraws it at once
  $('#nw-cut').onclick = (e) => { const f = e.target.dataset.f; if (!f) return; S.cut = f === 'c' ? 'c' : +f; [...$('#nw-cut').children].forEach((b) => b.classList.toggle('on', b.dataset.f === f));
    $('#nw-cutv').style.display = S.cut === 'c' ? '' : 'none'; redraw(); };
  { let t = 0; $('#nw-cutv').oninput = (e) => { clearTimeout(t); t = setTimeout(() => { const v = parseFloat(e.target.value); if (!Number.isFinite(v)) return;
      S.cutv = Math.min(1, Math.max(CUT[10], v)); if (S.cutv !== v) e.target.value = S.cutv; redraw(); }, 400); };
    let u = 0; $('#nw-iptm').oninput = (e) => { clearTimeout(u); u = setTimeout(() => { const v = parseFloat(e.target.value); S.iptm = Number.isFinite(v) ? Math.min(1, Math.max(0, v)) : 0; redraw(); }, 400); }; }
  $('#nw-add').onchange = () => { showK(); redraw(); };
  $('#nw-set').onchange = redraw;
  $('#nw-mind').onchange = () => { S.mind = +$('#nw-mind').value; redraw(); };
  { let t = 0; $('#nw-k').oninput = () => { clearTimeout(t); t = setTimeout(redraw, 350); }; }
  const cutV = () => (S.cut === 'c' ? S.cutv : CUT[S.cut]), cutLab = () => (S.cut === 'c' ? 'custom' : `${S.cut}% FPR`);   // the iLIS cutoff in force and its name
  const passE = (e) => e.best >= cutV() && (!S.iptm || (Number.isFinite(e.iptm) && e.iptm >= S.iptm));   // a pair past the cutoff and, when asked, the ipTM
  let net = null, KBN = null, EB = null, TPN = null;   // TPN: the species' tested pairs (every screen), when the file is there   // KBN: this species' BioGRID pairs; EB: its edges (homodimer rings), for the drawing
  // A reported pair that was not predicted past the cutoff is one of two things: folded and scored below the cutoff, or never
  // folded in these screens. edges.tsv holds only pairs past 10% FPR, so the answer comes from one end's own predictions
  // (its bundle lists every partner it was folded with). FOLD: scope|a,b → { st: 'low', best, iptm } | { st: 'none' } | { st: 'err' }
  const FOLD = new Map(), fkey = (a, b) => `${S.set}|${Math.min(a, b)},${Math.max(a, b)}`;
  async function foldedWith(i) {   // → Map(partner row → best iLIS over every model, best ipTM), or null when the protein is in no screen of this scope
    let all; try { all = await merged(sp, sp.rows[i], S.set, true); } catch (e) { if (/no predictions|No interaction data/i.test(e.message)) return null; throw e; }   // true: every isoform file too
    const m = new Map(); for (const p of all.preds) { const r = sp.byKey.get(p.partner); if (!r) continue; const x = m.get(r.i);
      if (!x || (p.iLIS || 0) > x.best) m.set(r.i, { best: p.iLIS || 0, iptm: Number.isFinite(p.ipTM) ? Math.max(p.ipTM, x ? x.iptm || 0 : 0) : x ? x.iptm : NaN }); }
    return m; }
  async function checkFolded(pairs, tick) {   // reads as few proteins as it can: each time the protein that ends the most unchecked pairs
    let left = pairs.filter(([a, b]) => !FOLD.has(fkey(a, b))), n = 0;
    while (left.length) { const cnt = new Map(); for (const [a, b] of left) for (const v of [a, b]) cnt.set(v, (cnt.get(v) || 0) + 1);
      const v = [...cnt].sort((x, y) => y[1] - x[1])[0][0]; let m, bad = false; try { m = await foldedWith(v); } catch (e) { bad = true; }
      for (const [a, b] of left.filter(([a, b]) => a === v || b === v)) { const o = a === v ? b : a;
        FOLD.set(fkey(a, b), bad ? { st: 'err' } : m && m.has(o) ? { st: 'low', ...m.get(o) } : { st: 'none' }); }
      left = left.filter(([a, b]) => a !== v && b !== v); if (tick) tick(++n, n + new Set(left.flat()).size); } }
  // Proteins colored by a column of the reader's table: numbers on a scale (centered on 0 when the column has both signs;
  // a p-value column as -log10), or up to 12 categories. Proteins the table does not give a value for stay light gray.
  const PCOL = /(^|[^a-z])(p|q)[-_. ]?(val|value)|padj|p\.adj|adj[-_. ]?p|fdr|bfdr/i;
  // The column a table's proteins are colored by until the reader picks one: a fold change or score, else a p-value, else the
  // first number column, else a column of groups; none for a plain list
  const FCOL = /fold|log2|lfc|\bfc\b|ratio|enrich|score|effect|\bbeta\b|\bz\b|abundance|intensity/i;
  const autoColor = (D) => { const num = D.filter((x) => x.kind === 'num'), pick = num.find((x) => FCOL.test(x.name) && !PCOL.test(x.name)) || num.find((x) => PCOL.test(x.name)) || num[0] || D.find((x) => x.kind === 'cat');
    return pick ? pick.name : ''; };
  function nodeColors() {
    const D = (S.data || []).find((x) => x.name === S.ncol), base = (d) => (d.q ? '#1A5276' : '#AEBBCA');
    if (!D) return { fill: base, key: '' };
    const none = '#DDE3EA';
    if (D.kind === 'cat') { const cats = [...new Set(D.vals.values())].sort(), col = (v) => TAB10[cats.indexOf(v) % TAB10.length];
      return { fill: (d) => (D.vals.has(d.id) ? col(D.vals.get(d.id)) : none),
        key: `<div class="kbrow"><span class="muted">proteins by ${esc(D.name)}:</span>${cats.map((v) => `<span><i style="background:${col(v)};border-radius:50%"></i>${esc(v)}</span>`).join('')}<span><i style="background:${none};border-radius:50%"></i>no value</span></div>` }; }
    const raw = [...D.vals.values()].filter(Number.isFinite), isP = PCOL.test(D.name) || (raw.length && raw.every((v) => v > 0 && v <= 1) && d3.median(raw) < 0.2);
    const pos = raw.filter((v) => v > 0), top = isP && pos.length ? Math.max(...pos.map((v) => -Math.log10(v))) + 1 : 1, zeros = isP && raw.some((v) => v <= 0);
    const isRatio = !isP && /fold|ratio|\bfc\b/i.test(D.name) && !/log/i.test(D.name) && pos.length === raw.length && raw.length && Math.max(...raw) / Math.min(...raw) > 10;   // a raw fold change spanning a decade or more: shown as log2
    const tf = (v) => (isP ? (v > 0 ? -Math.log10(v) : top) : isRatio ? Math.log2(v) : v), vs = raw.map(tf).sort((a, b) => a - b), lo = d3.quantile(vs, 0.02), hi = d3.quantile(vs, 0.98);   // a p-value of 0 sits one step past the column's largest -log10
    let sc, lab = isP ? `−log10 ${D.name}` : isRatio ? `log2 ${D.name}` : D.name, ends;
    if (lo < 0 && hi > 0) { const m = Math.max(-lo, hi); sc = d3.scaleDiverging((t) => d3.interpolateRdBu(1 - t)).domain([-m, 0, m]).clamp(true); ends = [-m, m]; }
    else { sc = d3.scaleSequential(d3.interpolateYlOrRd).domain([lo, hi === lo ? lo + 1 : hi]).clamp(true); ends = [lo, hi]; }
    const stops = d3.range(0, 1.01, 0.1).map((t) => sc(ends[0] + t * (ends[1] - ends[0]))).join(', '), f2 = (v) => (Math.abs(v) >= 100 ? v.toFixed(0) : v.toFixed(Math.abs(v) >= 10 ? 1 : 2));
    return { fill: (d) => { const v = D.vals.get(d.id); return Number.isFinite(v) ? sc(tf(v)) : none; },
      key: `<div class="kbrow"><span class="muted">proteins by ${esc(lab)}:</span><span><i class="kb-grad" style="background:linear-gradient(90deg, ${stops})"></i>${f2(ends[0])} to ${zeros ? '≥ ' : ''}${f2(ends[1])}</span>${zeros ? '<span class="muted">values of 0 shown at the top, one step past the smallest nonzero value</span>' : ''}<span><i style="background:${none};border-radius:50%"></i>no value</span></div>` };
  }
  const foldText = (f) => (!f ? 'not checked yet: hover again in a moment' : f.st === 'low' && !Number.isFinite(f.best) ? `folded, scored below the cutoff ${cutV()}` : f.st === 'low' ? `folded: best iLIS ${f.best < CUT[10] ? '≈' + f.best.toFixed(2) : f.best.toFixed(3)}${Number.isFinite(f.iptm) ? `, ipTM ${f.iptm.toFixed(2)}` : ''}, below the cutoff ${CUT[S.cut]}`
    : f.st === 'none' ? 'never folded in these screens: not tested, not a negative' : 'could not be checked (the prediction files did not load)');
  { const w = d3.select('#nw-w'); [[0.1, 10], [0.4, 95], [0.7, 180]].forEach(([a, x0]) => { w.append('line').attr('x1', x0).attr('x2', x0 + 44).attr('y1', 10).attr('y2', 10).attr('stroke', '#50637A').attr('stroke-width', EWID(a)).attr('stroke-linecap', 'round');
    w.append('text').attr('x', x0 + 22).attr('y', 27).attr('text-anchor', 'middle').attr('font-size', 10.5).attr('font-family', 'IBM Plex Mono').attr('fill', '#5A697C').text(a.toFixed(1)); }); }
  const status = (t) => { $('#nw-status').innerHTML = t; };
  const gname = (i) => sp.rows[i].gene, screens = (m) => sp.dsShort.filter((_, di) => m & (1 << di)).join(' + ');
  let tableNote = '';
  function readInput() {   // the names to draw, through readIdInput; the column picker and the note follow it
    const R = readIdInput(sp, $('#nw-ids').value, S.col), wrap = $('#nw-col-wrap'); S.col = R.col; S.data = R.data; tableNote = R.note;
    if (R.T) { $('#nw-col').innerHTML = (R.strongN >= 2 ? `<option value="all"${S.col === 'all' ? ' selected' : ''}>all name columns</option>` : '') + R.T.cols.map((x) => `<option value="${x.c}"${x.c === S.col ? ' selected' : ''}${x.num ? ' disabled' : ''}>${esc(x.name)}${x.num ? ' (numbers)' : ` (${fmtInt(x.hits)} found)`}</option>`).join(''); wrap.hidden = R.named < 2; }
    else wrap.hidden = true;
    $('#nw-table').textContent = tableNote; return R.toks;
  }
  async function draw(seed = null) {   // seed: the positions and zoom to keep when a click adds partners
    const toks = readInput();
    const ns = $('#nw-ncol'); if (S.ncol && !(S.data || []).some((x) => x.name === S.ncol)) S.ncol = '';
    if (!S.ncolUser) S.ncol = autoColor(S.data || []);   // chosen for the reader until they pick one
    ns.innerHTML = `<option value="">of interest · added partners</option>` + (S.data || []).map((x) => `<option value="${esc(x.name)}"${x.name === S.ncol ? ' selected' : ''}>${esc(x.name)}${x.kind === 'cat' ? ' (groups)' : ''}</option>`).join('');
    $('#nw-ncol-wrap').style.display = (S.data || []).length ? '' : 'none';
    const cats = (S.data || []).filter((x) => x.kind === 'cat'); if (S.grp.startsWith('col:') && !cats.some((x) => 'col:' + x.name === S.grp)) S.grp = '';
    $('#nw-grp').innerHTML = [['', 'none'], ['comm', 'communities (predicted pairs)'], ['q', 'of interest · added partners'], ...cats.map((x) => ['col:' + x.name, x.name])].map(([v, l]) => `<option value="${esc(v)}"${v === S.grp ? ' selected' : ''}>${esc(l)}</option>`).join('');
    S.ids = toks.join(','); S.add = $('#nw-add').value; S.k = Math.max(1, Math.min(50, +$('#nw-k').value || 5)); S.set = $('#nw-set').value;
    const writeURL = () => { const qs = new URLSearchParams({ ids: S.ids, add: S.add }); if (S.add === 'top' || S.exp.length) qs.set('k', S.k); qs.set('cut', S.cut); if (S.cut === 'c') qs.set('cutv', S.cutv); if (S.iptm) qs.set('iptm', S.iptm); if (S.mind) qs.set('mind', S.mind); if (S.grp) qs.set('grp', S.grp); if (S.set) qs.set('set', S.set);
      if (S.exp.length) qs.set('exp', S.exp.join(',')); if (S.click === 'open') qs.set('click', 'open');
      history.replaceState(null, '', `#/${sp.id}/network?${qs}`); };
    writeURL();
    if (!toks.length) { status('Name at least one protein.'); $('#nw-card').hidden = true; return; }
    const found = [], missing = [], how = new Map();
    for (const t of toks) { const h = resolveHow(sp, t, true); if (!h) { missing.push(t); continue; } how.set(h.how, (how.get(h.how) || 0) + 1); if (!found.includes(h.row)) found.push(h.row); }
    if (!found.length) { status(`None of these names is in the ${esc(sp.reg.label)} screens: ${esc(missing.join(', '))}.`); $('#nw-card').hidden = true; return; }
    const nFound = [...how.values()].reduce((s, n) => s + n, 0), twice = nFound - found.length;   // before the cap trims the list
    const CAP = 400, over = found.length > CAP ? found.length : 0; if (over) found.length = CAP;
    if (toks.length > CAP) { S.ids = found.map((r) => r.key).join(','); writeURL(); }   // a long file: the link carries the proteins drawn, not every line
    const howText = `${fmtInt(nFound)} name${nFound === 1 ? '' : 's'} found, by ${[...how].sort((a, b) => b[1] - a[1]).map(([h, n]) => `${h} ${fmtInt(n)}`).join(', ')}${twice ? ` (${fmtInt(twice)} named the same protein twice)` : ''}`;
    status('Reading the edge list…');
    let E, K, TP; try { [E, K, TP] = await Promise.all([edges(sp, S.set), reported(sp), S.set ? null : tested(sp)]); } catch (e) { status(esc(e.message)); return; }
    KBN = K; EB = E; TPN = TP;
    if (stale(gen)) return;
    const c = cutV(), Q = new Set(found.map((r) => r.i)), keep = new Set(Q), grew = new Set(), expd = new Set();
    const nb = (i) => [...(E.adj.get(i) || new Map())].filter(([, e]) => passE(e));
    if (S.add === 'shared') { const n = new Map(); for (const i of Q) for (const [j] of nb(i)) if (!Q.has(j)) n.set(j, (n.get(j) || 0) + 1);
      [...n].filter(([, k]) => k >= 2).sort((a, b) => b[1] - a[1] || gname(a[0]).localeCompare(gname(b[0]))).slice(0, CAP - keep.size).forEach(([j]) => keep.add(j)); }
    else if (S.add === 'top') for (const i of Q) { for (const [j] of nb(i).sort((a, b) => b[1].best - a[1].best).slice(0, S.k)) { if (keep.size >= CAP) break; keep.add(j); } }
    let expNote = '', lastFresh = 0, lastAdded = 0;   // what the latest click could add, and did
    for (const key of S.exp) {   // clicks: each expanded protein's top partners, in the order clicked, while under the cap
      const r = sp.byKey.get(key); if (!r || !keep.has(r.i)) continue; expd.add(r.i);
      const fresh = nb(r.i).filter(([j]) => !keep.has(j)).sort((a, b) => b[1].best - a[1].best).slice(0, S.k);
      if (!fresh.length && seed && seed.at === r.i) { expNote = ` · ${esc(r.gene)} has no other partner past this cutoff`; }
      let added = 0; for (const [j] of fresh) { if (keep.size >= CAP) break; keep.add(j); grew.add(j); added++; }
      if (seed && seed.at === r.i) { lastFresh = fresh.length; lastAdded = added; } }
    if (seed && !expNote && lastFresh && !lastAdded) expNote = ` · the network is at ${CAP} proteins; a click adds no more`;
    if (expNote && S.exp.length) { S.exp.pop(); writeURL(); }   // a click that added nothing is not kept in the link
    $('#nw-unexp').style.display = S.exp.length ? '' : 'none';
    let links = [];
    for (const a of keep) for (const [b, e] of E.adj.get(a) || []) if (b > a && keep.has(b) && passE(e)) links.push({ source: a, target: b, ...e, pubs: K ? K.pubs(a, b) : 0, gen: K ? K.gen(a, b) : 0 });
    let pruned = 0;   // Min. pairs: drop proteins with fewer pairs in the drawing, again until none is left below it (a k-core)
    if (S.mind) for (let again = true; again;) { again = false; const dg = new Map(); for (const l of links) { dg.set(l.source, (dg.get(l.source) || 0) + 1); dg.set(l.target, (dg.get(l.target) || 0) + 1); }
      for (const i of [...keep]) if ((dg.get(i) || 0) < S.mind) { keep.delete(i); pruned++; again = true; }
      if (again) links = links.filter((l) => keep.has(l.source) && keep.has(l.target)); }
    const nRep = links.filter((l) => l.pubs || l.gen).length, extra = [];
    if (K && $('#nw-unpred').checked) { const ks = [...keep], has = new Set(links.map((l) => l.source * sp.rows.length + l.target));
      for (let x = 0; x < ks.length; x++) for (let y = x + 1; y < ks.length; y++) { const a = Math.min(ks[x], ks[y]), b = Math.max(ks[x], ks[y]), n = K.pubs(a, b), gg = K.gen(a, b);
        if ((n || gg) && !has.has(a * sp.rows.length + b)) extra.push({ source: a, target: b, pubs: n, gen: gg, unpred: true }); } }
    const Qin = [...Q].filter((i) => keep.has(i)), added = keep.size - Qin.length, alone = Qin.filter((i) => !links.some((l) => l.source === i || l.target === i)).map(gname);
    status(`${fmtInt(Qin.length)} ${sp.manifest.keyedBy ? 'genes' : 'proteins'} of interest${added ? `, ${fmtInt(added)} added partner${added === 1 ? '' : 's'}` : ''} · ${fmtInt(links.length)} predicted pair${links.length === 1 ? '' : 's'} past iLIS ${c} (${cutLab()})${S.iptm ? ` with ipTM ≥ ${S.iptm}` : ''}`
      + `${K ? ` · ${fmtInt(nRep)} of them reported in BioGRID ${K.release}${extra.length ? `, plus ${fmtInt(extra.length)} reported pair${extra.length === 1 ? '' : 's'} not predicted (dashed)` : ''}` : ' · no BioGRID records for this species'}`
      + `${missing.length ? ` · not found: ${esc(missing.slice(0, 30).join(', '))}${missing.length > 30 ? ` and ${fmtInt(missing.length - 30)} more` : ''}${missing.some((t) => /^ENS[A-Z]*[GTP]\d{6,}/i.test(t) || /^\d+$/.test(t)) ? ' (Ensembl and Entrez IDs are not in the Atlas index; use gene symbols or UniProt accessions)' : ''}` : ''}${alone.length && links.length ? ` · no pair here for ${esc(alone.join(', '))}` : ''}${keep.size >= CAP ? ` · capped at ${CAP} proteins; open it in LIVIA Network for more` : ''}${expNote}`
      + `<br><span class="muted">${esc(howText)}${over ? ` · the first ${CAP} of ${fmtInt(over)} proteins drawn` : ''}${tableNote ? ` · ${esc(tableNote)}` : ''}${pruned ? ` · ${fmtInt(pruned)} hidden with fewer than ${S.mind} pair${S.mind === 1 ? '' : 's'} in the drawing` : ''}</span>`);
    $('#nw-card').hidden = false;
    if (!links.length) { $('#nw-net').innerHTML = '<div class="empty">No predicted pair among these proteins at this cutoff. Try + partners, or a lower cutoff.</div>'; net = null; $('#nw-heat-wrap').hidden = true; return; }
    if (!K) { $('#nw-kb').checked = false; $('#nw-kb').disabled = true; $('#nw-ev').disabled = true; }
    if (TPN) for (const d of extra) { const k = fkey(d.source, d.target); if (!FOLD.has(k)) FOLD.set(k, TPN.has(d.source, d.target) ? { st: 'low', best: TPN.score ? TPN.score(d.source, d.target) : NaN } : { st: 'none' }); }   // the index says folded or not; the score comes on hover
    graph([...keep].map((i) => ({ id: i, row: sp.rows[i], q: Q.has(i), grew: grew.has(i), expd: expd.has(i) })), links, extra, seed);
  }
  // The drawn proteins as a matrix: a predicted pair past the cutoff colored by its best iLIS; folded but below the cutoff light
  // gray; never folded white (from the index of tested pairs; without it, every other cell is gray). Rows follow the grouping.
  let heatArgs = null;
  // Matrix color scales: light at the cutoff to dark at 0.85 and above, starting well clear of the gray of pairs below the cutoff
  const HEATCS = (k) => { const f = { Blues: (t) => d3.interpolateBlues(0.35 + 0.65 * t), Viridis: (t) => d3.interpolateViridis(0.95 - 0.95 * t), YlGnBu: (t) => d3.interpolateYlGnBu(0.3 + 0.7 * t),
    Reds: (t) => d3.interpolateReds(0.35 + 0.65 * t), Grays: (t) => d3.interpolateGreys(0.45 + 0.55 * t), Cividis: (t) => d3.interpolateCividis(0.95 - 0.95 * t) }[k] || ((t) => d3.interpolateBlues(0.35 + 0.65 * t));
    return d3.range(0, 1.0001, 0.1).map((t) => [+t.toFixed(2), f(t)]); };
  // Plotly, loaded the first time a matrix is drawn (zoom, pan, every label on hover and when zoomed)
  const plotly = () => (window.Plotly ? Promise.resolve(window.Plotly) : (window.__plotlyP = window.__plotlyP || new Promise((ok, no) => {
    const s = document.createElement('script'); s.src = 'https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js'; s.onload = () => ok(window.Plotly); s.onerror = () => no(new Error('Plotly did not load')); document.head.appendChild(s); })));
  // The drawn proteins as a matrix: a predicted pair past the cutoff colored by its best iLIS; folded but below the cutoff gray;
  // never folded white (from the index of tested pairs; without it, every other cell is gray). Rows follow the grouping.
  async function heatmap(nodes, links, gk, groups, gname2, deg) {
    heatArgs = [nodes, links, gk, groups, gname2, deg];
    const wrap = $('#nw-heat-wrap'), box = $('#nw-heat'), n = nodes.length; wrap.hidden = n < 2; if (n < 2) return;
    const gi = (d) => { const x = groups.indexOf(gk.get(d.id)); return x < 0 ? groups.length : x; };
    const ord = [...nodes].sort((a, b) => (gi(a) - gi(b)) || (b.q - a.q) || ((deg.get(b.id) || 0) - (deg.get(a.id) || 0)) || a.row.gene.localeCompare(b.row.gene));
    const M = new Map(), key = (a, b) => (a < b ? a + ',' + b : b + ',' + a);
    for (const l of links) { const a = typeof l.source === 'object' ? l.source.id : l.source, b = typeof l.target === 'object' ? l.target.id : l.target; M.set(key(a, b), l); }
    const TP = S.set ? null : TPN, cut = cutV(), names = ord.map((d) => d.row.gene);
    const showAll = $('#nw-heat-show').value !== 'pos', E0 = EB && !S.set ? EB : null;   // below the cutoff: the edge list's score (0.223 and up), else the index's
    const zi = [], zg = [], tx = []; let folded = 0, never = 0;
    for (let i = 0; i < n; i++) { const ri = [], rg = [], rt = [];
      for (let j = 0; j < n; j++) { const a = ord[i].id, b = ord[j].id, l = i === j ? null : M.get(key(a, b)), f = i !== j && !l && (!TP || TP.has(a, b));
        let v = NaN; if (f) { const e = E0 && (E0.adj.get(a) || new Map()).get(b); v = e ? e.best : TP && TP.score ? TP.score(a, b) : NaN; }
        ri.push(l ? l.best : null); rg.push(f && showAll ? (Number.isFinite(v) ? v : 0) : null);
        const vs = Number.isFinite(v) ? `best iLIS ${v < CUT[10] ? '≈' + v.toFixed(2) : v.toFixed(3)}, below the cutoff` : TP ? 'folded, below the cutoff' : 'not past the cutoff';
        rt.push(i === j ? names[i] : l ? `${names[i]} × ${names[j]}<br>iLIS ${l.best.toFixed(3)}${Number.isFinite(l.iptm) ? ` · ipTM ${l.iptm.toFixed(2)}` : ''}` : f ? `${names[i]} × ${names[j]}<br>${vs}` : `${names[i]} × ${names[j]}<br>never folded in these screens`);
        if (j > i && !l) { if (f) folded++; else never++; } }
      zi.push(ri); zg.push(rg); tx.push(rt); }
    $('#nw-heat-note').textContent = `${fmtInt(n)} × ${fmtInt(n)} · ${fmtInt(links.length)} of ${fmtInt(n * (n - 1) / 2)} pairs past iLIS ${cut}` + (TP ? ` · ${fmtInt(folded)} folded below it · ${fmtInt(never)} never folded` : '');
    $('#nw-heat-key').innerHTML = `${showAll ? `<span><i class="kb-grad" style="background:linear-gradient(90deg, #F1F3F5, #8E98A3)"></i>${TP ? `folded, below the cutoff: best iLIS 0 to ${cut}` : 'below the cutoff or not folded'}</span>` : '<span class="muted">pairs below the cutoff not shown</span>'}${TP ? '<span><i style="background:#fff;border:1px solid #C3CCD6"></i>never folded</span>' : ''}<span class="muted">scroll or drag a box to zoom · double-click to reset · click a colored cell for the pair</span>`;
    let P; try { P = await plotly(); } catch (e) { box.innerHTML = `<div class="empty">${esc(e.message)}</div>`; return; }
    const shapes = []; ord.forEach((d, k) => { if (k && gi(d) !== gi(ord[k - 1])) for (const s of [{ x0: k - 0.5, x1: k - 0.5, y0: -0.5, y1: n - 0.5 }, { y0: k - 0.5, y1: k - 0.5, x0: -0.5, x1: n - 0.5 }]) shapes.push({ type: 'line', ...s, line: { color: '#5B6573', width: 1 } }); });
    const side = Math.max(360, Math.min(900, 18 * n + 160)), tick = Math.max(6, Math.min(11, 520 / n));
    const cs = $('#nw-heat-cs').value || 'Blues';
    P.react(box, [
      { type: 'heatmap', z: zg, x: names, y: names, text: tx, hovertemplate: '%{text}<extra></extra>', colorscale: [[0, '#F1F3F5'], [1, '#8E98A3']], zmin: 0, zmax: cut, showscale: false, xgap: 1, ygap: 1, hoverongaps: false },   // below the cutoff: grays, darker toward it
      { type: 'heatmap', z: zi, x: names, y: names, text: tx, hovertemplate: '%{text}<extra></extra>', colorscale: HEATCS(cs), zmin: CUT[10], zmax: 0.85,
        colorbar: { title: { text: 'best iLIS', side: 'right' }, thickness: 12, len: 0.6 }, xgap: 1, ygap: 1, hoverongaps: false }],
      { width: Math.min(box.clientWidth || 900, side + 120), height: side, margin: { l: 110, r: 20, t: 110, b: 20 }, plot_bgcolor: '#FFFFFF', paper_bgcolor: 'rgba(0,0,0,0)', shapes,
        xaxis: { side: 'top', tickangle: -60, tickfont: { size: tick, family: 'IBM Plex Sans, sans-serif' }, automargin: true, showgrid: false, constrain: 'domain' },
        yaxis: { autorange: 'reversed', tickfont: { size: tick, family: 'IBM Plex Sans, sans-serif' }, automargin: true, showgrid: false, scaleanchor: 'x' }, dragmode: 'zoom' },
      { displaylogo: false, responsive: true, scrollZoom: true, toImageButtonOptions: { filename: `atlas_${sp.id}_matrix`, format: 'svg' }, modeBarButtonsToRemove: ['select2d', 'lasso2d'] });
    box.removeAllListeners && box.removeAllListeners('plotly_click'); box.removeAllListeners && box.removeAllListeners('plotly_relayout');
    box.on('plotly_relayout', (ev) => {   // zoomed in, the labels grow with the room each one has (6 to 14 px)
      const xr = ev['xaxis.range[0]'] != null ? [ev['xaxis.range[0]'], ev['xaxis.range[1]']] : Array.isArray(ev['xaxis.range']) ? ev['xaxis.range'] : ev['xaxis.autorange'] ? [-0.5, n - 0.5] : null;
      if (!xr) return; const shown = Math.max(1, Math.abs(xr[1] - xr[0])), f = Math.max(6, Math.min(14, 520 / shown));
      if (Math.abs(f - (box.__tick || tick)) < 0.5) return; box.__tick = f; P.relayout(box, { 'xaxis.tickfont.size': f, 'yaxis.tickfont.size': f }); });
    box.on('plotly_click', (ev) => { const p = ev.points && ev.points.find((q) => q.curveNumber === 1); if (!p || p.z == null) return; const a = ord[p.pointIndex[0]], b = ord[p.pointIndex[1]], u = `#/${sp.id}/${a.row.key}/${b.row.key}`, e2 = ev.event || {}; if (e2.metaKey || e2.ctrlKey) window.open(u, '_blank'); else location.hash = u; });
  }
  function graph(nodes, links, extra = [], seed = null) {
    const box = $('#nw-net'); box.innerHTML = '<svg></svg>';
    const W = box.clientWidth, H = box.clientHeight, svg = d3.select(box).select('svg').attr('width', W).attr('height', H), g = svg.append('g');
    let label = null, rank = null, fitK = seed ? seed.t.k : 1;   // labels: set once the nodes are drawn; relabel() shows more as the view zooms in
    const zoom = d3.zoom().scaleExtent([0.1, 8]).on('zoom', (ev) => { g.attr('transform', ev.transform); relabel(ev.transform.k); }); svg.call(zoom);
    if (seed) { svg.call(zoom.transform, seed.t);   // a click added partners: the old nodes stay put, the new ones start at the protein clicked
      const [ax, ay] = seed.pos.get(seed.at) || [W / 2, H / 2];
      nodes.forEach((d, n) => { const p = seed.pos.get(d.id); if (p) { [d.x, d.y] = p; } else { d.x = ax + 25 * Math.cos(n); d.y = ay + 25 * Math.sin(n); } }); }
    const deg = new Map(); links.forEach((l) => { deg.set(l.source, (deg.get(l.source) || 0) + 1); deg.set(l.target, (deg.get(l.target) || 0) + 1); });
    const r = (d) => (d.q ? 10 : 4) + Math.min(8, Math.sqrt(deg.get(d.id) || 0) * 1.4);
    // Group by: a key per protein; groups of two or more get a place of their own on a ring and an outline
    let gkey = null, gname2 = new Map();
    if (S.grp === 'comm') { const ix = new Map(nodes.map((d, n) => [d.id, n])), cm = louvain(nodes.length, links.map((l) => [ix.get(typeof l.source === 'object' ? l.source.id : l.source), ix.get(typeof l.target === 'object' ? l.target.id : l.target), l.best]));
      const size = new Map(); cm.forEach((k) => size.set(k, (size.get(k) || 0) + 1)); const order = [...size].filter(([, s]) => s >= 3).sort((a, b) => b[1] - a[1]).map(([k]) => k);
      gkey = (d) => { const k = cm[ix.get(d.id)]; return order.includes(k) ? 'c' + k : null; }; order.forEach((k, n) => gname2.set('c' + k, `community ${n + 1}`)); }
    else if (S.grp === 'q') { gkey = (d) => (d.q ? 'q' : 'a'); gname2.set('q', `${sp.manifest.keyedBy ? 'genes' : 'proteins'} of interest`); gname2.set('a', 'added partners'); }
    else if (S.grp.startsWith('col:')) { const D = (S.data || []).find((x) => 'col:' + x.name === S.grp); if (D) { gkey = (d) => (D.vals.has(d.id) ? 'v' + D.vals.get(d.id) : null); for (const v of new Set(D.vals.values())) gname2.set('v' + v, String(v)); } }
    const gk = new Map(), gsize = new Map(); if (gkey) for (const d of nodes) { const k = gkey(d); if (k != null) { gk.set(d.id, k); gsize.set(k, (gsize.get(k) || 0) + 1); } }
    const groups = [...gsize].filter(([, s]) => s >= 2).sort((a, b) => b[1] - a[1]).map(([k]) => k), gcenter = new Map();
    groups.forEach((k, n) => { if (groups.length === 1) { gcenter.set(k, [W / 2, H / 2]); return; } const t = (n / groups.length) * 2 * Math.PI - Math.PI / 2, R = Math.min(W, H) * 0.36;
      gcenter.set(k, [W / 2 + R * Math.cos(t) * 1.25, H / 2 + R * Math.sin(t)]); });
    const home = (d) => gcenter.get(gk.get(d.id));
    if (groups.length && !seed) nodes.forEach((d, n) => { const h = home(d); if (h) { d.x = h[0] + 20 * Math.cos(n); d.y = h[1] + 20 * Math.sin(n); } });
    const catOrder = S.grp.startsWith('col:') ? [...gname2.keys()].sort((a, b) => a.localeCompare(b)) : null;   // the same order nodeColors gives a category column
    const hullG = g.append('g').attr('class', 'nw-hulls').style('pointer-events', 'none'), GCOL = (n) => TAB10[(catOrder ? catOrder.indexOf(groups[n]) : n) % TAB10.length];
    const hulls = hullG.selectAll('g').data(groups).join('g');
    hulls.append('path').attr('fill', (k, n) => GCOL(n)).attr('fill-opacity', 0.07).attr('stroke', (k, n) => GCOL(n)).attr('stroke-opacity', 0.45).attr('stroke-width', 1.5).attr('stroke-linejoin', 'round');
    hulls.append('text').text((k) => `${gname2.get(k) || k} (${fmtInt(gsize.get(k))})`).attr('font-family', 'IBM Plex Sans, sans-serif').attr('font-size', 12).attr('font-weight', 700).attr('fill', (k, n) => GCOL(n)).attr('text-anchor', 'middle');
    const drawHulls = () => { if (!groups.length) return; hulls.each(function (k) { const pts = []; for (const d of nodes) if (gk.get(d.id) === k && d.x != null) { const rr = r(d) + 14; for (let a = 0; a < 6; a++) pts.push([d.x + rr * Math.cos(a * Math.PI / 3), d.y + rr * Math.sin(a * Math.PI / 3)]); }
      const h = pts.length >= 3 ? d3.polygonHull(pts) : null; const el = d3.select(this);
      el.select('path').attr('d', h ? d3.line().curve(d3.curveCatmullRomClosed.alpha(0.6))(h) : null);
      if (h) { const top = h.reduce((a, b) => (b[1] < a[1] ? b : a)); const cx = d3.mean(h, (p) => p[0]); el.select('text').attr('x', cx).attr('y', top[1] - 6); } }); };
    const dash = g.append('g').selectAll('line').data(extra).join('line').attr('stroke', kbCol).attr('stroke-opacity', 0.75).attr('stroke-width', 1.4).attr('stroke-dasharray', '5 4').style('cursor', 'pointer');
    const link = g.append('g').selectAll('line').data(links).join('line').attr('stroke-width', (d) => EWID(d.avg)).attr('stroke-linecap', 'round').style('cursor', 'pointer');
    const homo = new Set(nodes.filter((d) => { const e = (EB && EB.adj.get(d.id) || new Map()).get(d.id); return e && passE(e); }).map((d) => d.id));
    const ends = (d) => [typeof d.source === 'object' ? d.source.id : d.source, typeof d.target === 'object' ? d.target.id : d.target];
    const K = KBN, restyle = () => { const st = edgeStyle('nw', K);   // recolor in place: no new layout
      link.attr('stroke', st.color).attr('stroke-opacity', (d) => (st.hit(d) ? 0.95 : 0.8)); link.filter(st.hit).raise();   // reported pairs on top
      const shownExtra = st.kb ? extra.filter((d) => kbHit(d, st.ev)) : []; dash.attr('display', (d) => (shownExtra.includes(d) ? null : 'none'));
      const fo = (d) => FOLD.get(fkey(...ends(d)));
      dash.attr('stroke-dasharray', (d) => ((fo(d) || {}).st === 'none' ? '1.2 4.5' : '5 4')).attr('stroke-opacity', (d) => ((fo(d) || {}).st === 'none' ? 0.55 : 0.75));
      const nf = (s) => shownExtra.filter((d) => (fo(d) || {}).st === s).length, nck = shownExtra.filter((d) => fo(d)).length;
      const foldKey = shownExtra.length && nck ? `<div class="kbrow"><span class="muted">reported, not predicted (${fmtInt(shownExtra.length)}):</span><span><i class="kb-dash"></i>folded, scored below the cutoff (${fmtInt(nf('low'))})</span><span><i class="kb-dot"></i>never folded in these screens (${fmtInt(nf('none'))})</span>${nck < shownExtra.length ? `<span class="muted">not checked (${fmtInt(shownExtra.length - nck)})</span>` : ''}${nf('err') ? `<span class="muted">could not be checked (${fmtInt(nf('err'))})</span>` : ''}</div>` : '';
      $('#nw-key').innerHTML = edgeKey(st, K, links, shownExtra) + foldKey + (homo.size ? `<div class="kbrow">${homoKey(homo.size)}</div>` : '');
      const cb = $('#nw-check'); cb.style.display = shownExtra.length ? '' : 'none'; if (!cb.dataset.busy) cb.textContent = nck && nck === shownExtra.length ? 'Checked: which were folded' : 'Check which were folded'; };
    restyle();
    const pubsText = (d) => (d.pubs || d.gen ? ` · reported in BioGRID (${kbPubs(d.pubs, d.gen)})` : '');
    link.on('mousemove', (ev, d) => { const [a, b] = ends(d); showTip(`<b>${esc(gname(a))}</b> × <b>${esc(gname(b))}</b> · iLIS ${sp.one ? d.best.toFixed(3) : `best ${d.best.toFixed(3)} · average ${d.avg.toFixed(3)}`}${Number.isFinite(d.iptm) ? ` · ipTM ${d.iptm.toFixed(2)}` : ''}${pubsText(d)}<br>${srcBadges(sp, d.src)} · click for the interaction residues`, ev.clientX, ev.clientY); })
      .on('mouseleave', hideTip).on('click', (ev, d) => { hideTip(); const [a, b] = ends(d), u = `#/${sp.id}/${sp.rows[a].key}/${sp.rows[b].key}${S.set ? `?set=${encodeURIComponent(S.set)}` : ''}`; if (ev.metaKey || ev.ctrlKey) window.open(u, '_blank'); else location.hash = u; });
    dash.on('mousemove', (ev, d) => { const [a, b] = ends(d), f = FOLD.get(fkey(a, b));
        showTip(`<b>${esc(gname(a))}</b> × <b>${esc(gname(b))}</b>${pubsText(d)}<br>not predicted past iLIS ${cutV()}${S.iptm ? ` with ipTM ≥ ${S.iptm}` : ''}: ${esc(foldText(f))}<br>click for the pair's predictions`, ev.clientX, ev.clientY);
        if ((!f || (f.st === 'low' && !Number.isFinite(f.best))) && !d.checking) { d.checking = true; if (f) FOLD.delete(fkey(a, b)); checkFolded([[a, b]]).then(() => { d.checking = false; restyle(); }); } })
      .on('mouseleave', hideTip).on('click', (ev, d) => { hideTip(); const [a, b] = ends(d), u = `#/${sp.id}/${sp.rows[a].key}/${sp.rows[b].key}`; if (ev.metaKey || ev.ctrlKey) window.open(u, '_blank'); else location.hash = u; });
    const node = g.append('g').selectAll('g').data(nodes).join('g').style('cursor', 'pointer')
      .call(d3.drag().on('start', (ev, d) => { if (!ev.active) sim.alphaTarget(0.25).restart(); d.fx = d.x; d.fy = d.y; })
        .on('drag', (ev, d) => { d.fx = ev.x; d.fy = ev.y; }).on('end', (ev, d) => { if (!ev.active) sim.alphaTarget(0); d.fx = null; d.fy = null; }));
    node.append('circle').attr('class', 'nfill').attr('r', r).attr('fill', (d) => (d.q ? '#1A5276' : '#AEBBCA')).attr('stroke', (d) => (homo.has(d.id) ? HOMO_RING : '#fff')).attr('stroke-width', (d) => (homo.has(d.id) ? 3 : d.q ? 2.2 : 1.5));
    node.filter((d) => d.expd).append('circle').attr('r', (d) => r(d) + 4.5).attr('fill', 'none').attr('stroke', '#1A5276').attr('stroke-width', 1.4).attr('stroke-dasharray', '3 2.5');   // expanded by a click
    // every protein gets a label; the fitted view shows the proteins of interest and the best-connected partners, and zooming in shows
    // more (the budget grows with the square of the zoom), each label at the same size on screen
    rank = new Map([...nodes].sort((a, b) => (b.q - a.q) || ((deg.get(b.id) || 0) - (deg.get(a.id) || 0))).map((d, n) => [d.id, n]));
    label = node.append('text').text((d) => d.row.gene).attr('text-anchor', 'middle')
      .attr('font-family', 'IBM Plex Sans, sans-serif').attr('font-weight', (d) => (d.q ? 700 : 600)).attr('fill', '#17263A')
      .attr('paint-order', 'stroke').attr('stroke', 'rgba(255,255,255,0.92)').attr('stroke-linejoin', 'round');
    relabel(d3.zoomTransform(svg.node()).k);
    // Labels in priority order (proteins of interest, then the best connected), each kept only if it does not overlap one already
    // placed, on screen. Below the fitted zoom they shrink with the view; zoomed in they keep their size, so more fit.
    function relabel(k) { if (!label) return; const z = Math.max(k, fitK), t = d3.zoomTransform(svg.node()), boxes = [];
      const fs = (d) => (d.q ? 13 : 10.5) * (k / z);   // font size on screen
      const shown = new Set(), order = [...nodes].sort((a, b) => rank.get(a.id) - rank.get(b.id));
      const partners = nodes.length <= 80 || k >= fitK * 1.3;   // added partners get labels once the view is zoomed in past the fit
      for (const d of order) { if (d.x == null || (!d.q && !partners)) continue; const f = fs(d); if (f < 6.5) continue;
        const w = String(d.row.gene).length * f * 0.62 + 8, h = f + 5, cx = t.applyX(d.x), cy = t.applyY(d.y) - (r(d) * k + 5 * (k / z)) - f / 2;
        const b = [cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2];
        if (boxes.some((o) => b[0] < o[2] && b[2] > o[0] && b[1] < o[3] && b[3] > o[1])) continue;
        boxes.push(b); shown.add(d.id); }
      label.attr('display', (d) => (shown.has(d.id) ? null : 'none')).attr('font-size', (d) => (d.q ? 13 : 10.5) / z).attr('stroke-width', 3 / z).attr('dy', (d) => -r(d) - 5 / z); }
    node.filter((d) => d.q).raise();
    const recolor = () => { const nc = nodeColors(); node.select('circle.nfill').attr('fill', nc.fill); $('#nw-nkey').innerHTML = nc.key; }; recolor();
    const placeLabels = liftLabels(g, node);
    node.on('mousemove', (ev, d) => showTip(`<b>${esc(d.row.gene)}</b>${d.row.name ? ` · ${esc(short(d.row.name))}` : ''}<br>${fmtInt(deg.get(d.id) || 0)} pair${(deg.get(d.id) || 0) === 1 ? '' : 's'} in this network · ${fmtInt(d.row.pos10)} partners past 10% FPR in the Atlas${(S.data || []).filter((x) => x.vals.has(d.id)).map((x) => `<br>${esc(x.name)}: ${esc(String(x.vals.get(d.id)))}`).join('')}<br>${S.click === 'add' ? `click to add its top ${S.k} partners · ⌘ or Ctrl-click opens its page in a new tab` : 'click for its page'}`, ev.clientX, ev.clientY))
      .on('mouseleave', hideTip).on('click', (ev, d) => { if (ev.defaultPrevented) return; hideTip();
        if (ev.metaKey || ev.ctrlKey) { window.open(`#/${sp.id}/${d.row.key}`, '_blank'); return; }   // a new tab, as for a link
        if (S.click === 'open') { location.hash = `#/${sp.id}/${d.row.key}`; return; }
        S.exp.push(d.row.key);   // a second click on the same protein adds its next partners
        draw({ at: d.id, t: d3.zoomTransform(svg.node()), pos: new Map(nodes.map((x) => [x.id, [x.x, x.y]])) }); });
    const sim = d3.forceSimulation(nodes).force('link', d3.forceLink([...links, ...extra]).id((d) => d.id).distance((l) => 70 + 60 * (1 - Math.min(1, l.best || 0))).strength((l) => (l.unpred ? 0 : 0.4)))
      .force('charge', d3.forceManyBody().strength(groups.length ? -140 : -260)).force('collide', d3.forceCollide().radius((d) => r(d) + 10))
      .force('x', d3.forceX((d) => (home(d) || [W / 2])[0]).strength((d) => (home(d) ? 0.35 : 0.05))).force('y', d3.forceY((d) => (home(d) || [0, H / 2])[1]).strength((d) => (home(d) ? 0.35 : 0.06)))
      .on('tick', () => { placeLabels(); drawHulls(); for (const sel of [dash, link]) sel.attr('x1', (d) => d.source.x).attr('y1', (d) => d.source.y).attr('x2', (d) => d.target.x).attr('y2', (d) => d.target.y); node.attr('transform', (d) => `translate(${d.x},${d.y})`); });
    if (seed) sim.alpha(0.35);   // settle the new nodes without reshuffling the rest
    let fitted = !!seed;   // once the layout settles, zoom so every node and label fits (the zoom stays free afterwards); kept as it was after a click
    sim.on('end', () => { if (fitted) return; fitted = true; const xs = nodes.map((d) => d.x), ys = nodes.map((d) => d.y);
      const x0 = Math.min(...xs) - 48, x1 = Math.max(...xs) + 48, y0 = Math.min(...ys) - 34, y1 = Math.max(...ys) + 24, sc = Math.min(1.4, 0.96 * Math.min(W / (x1 - x0), H / (y1 - y0)));
      fitK = sc; svg.transition().duration(450).call(zoom.transform, d3.zoomIdentity.translate(W / 2 - sc * (x0 + x1) / 2, H / 2 - sc * (y0 + y1) / 2).scale(sc)); });
    sim.on('end.labels', () => relabel(d3.zoomTransform(svg.node()).k));
    svgExport($('#nw-x'), `atlas_${sp.id}_network`, () => $('svg', box));
    const focus = (id) => { const d = nodes.find((n) => n.id === id); if (!d || d.x == null) return;   // center the protein, zoomed in, and mark it
      const k = Math.max(d3.zoomTransform(svg.node()).k, fitK * 2.2, 1.2);
      svg.transition().duration(550).call(zoom.transform, d3.zoomIdentity.translate(W / 2 - k * d.x, H / 2 - k * d.y).scale(k));
      g.selectAll('circle.nw-found').remove();
      const ring = g.append('circle').attr('class', 'nw-found').attr('cx', d.x).attr('cy', d.y).attr('r', r(d) + 6).attr('fill', 'none').attr('stroke', '#E4572E').attr('stroke-width', 3).style('pointer-events', 'none');
      ring.transition().delay(550).duration(900).attr('r', r(d) + 16).attr('stroke-opacity', 0.2).transition().duration(600).attr('r', r(d) + 6).attr('stroke-opacity', 1); };
    net = { link, nodes, links, extra, restyle, recolor, focus };
    heatmap(nodes, links, gk, groups, gname2, deg);
  }
  $('#nw-go').onclick = () => draw();
  $('#nw-filebtn').onclick = () => $('#nw-file').click();
  examplePicker($('#nw-ex'), sp.id, (t) => { $('#nw-ids').value = t;
    S.col = null; S.exp = []; S.ncolUser = false; $('#nw-add').value = 'top'; $('#nw-k').value = 3; showK(); draw(); });   // an example opens with each hit's top partners, so the list grows outward
  $('#nw-file').onchange = async (e) => { const f = e.target.files[0]; if (!f) return; $('#nw-ids').value = await f.text(); S.col = null; S.exp = []; S.ncolUser = false; e.target.value = ''; draw(); };
  $('#nw-col').onchange = (e) => { S.col = e.target.value === 'all' ? 'all' : +e.target.value; S.exp = []; S.ncolUser = false; draw(); };
  $('#nw-ids').oninput = () => { S.col = null; S.ncolUser = false; };
  $('#nw-click').onclick = (e) => { const m = e.target.dataset.m; if (!m) return; S.click = m; [...$('#nw-click').children].forEach((b) => b.classList.toggle('on', b.dataset.m === m)); showK();
    const [path, qs] = location.hash.split('?'), u = new URLSearchParams(qs || ''); m === 'open' ? u.set('click', 'open') : u.delete('click'); history.replaceState(null, '', `${path}?${u}`); };
  $('#nw-unexp').onclick = () => { S.exp = []; draw(); };
  $('#nw-grp').onchange = (e) => { S.grp = e.target.value; draw(); };
  $('#nw-heat-cs').onchange = () => { if (heatArgs) heatmap(...heatArgs); };
  $('#nw-heat-show').onchange = () => { if (heatArgs) heatmap(...heatArgs); };
  $('#nw-ncol').onchange = (e) => { S.ncol = e.target.value; S.ncolUser = true; if (net && net.recolor) net.recolor(); };
  $('#nw-check').onclick = async () => { if (!net || !net.extra.length) return; const b = $('#nw-check'); if (b.dataset.busy) return; b.dataset.busy = '1';
    const e2 = (l) => [typeof l.source === 'object' ? l.source.id : l.source, typeof l.target === 'object' ? l.target.id : l.target];
    b.textContent = 'Checking…'; await checkFolded(net.extra.map(e2), (n, of) => { b.textContent = `Checking ${n} of ${of} proteins…`; });
    delete b.dataset.busy; if (net && net.restyle) net.restyle(); };
  ['#nw-shade', '#nw-kb', '#nw-ev'].forEach((q) => { $(q).onchange = () => { $('#nw-ev').disabled = !$('#nw-kb').checked; if (net && net.restyle) net.restyle(); }; }); $('#nw-unpred').onchange = () => { if ($('#nw-unpred').checked && !$('#nw-kb').disabled) { $('#nw-kb').checked = true; $('#nw-ev').disabled = false; } draw(); };   // those lines only mean something with the crimson layer
  $('#nw-ids').onkeydown = (e) => { if (e.key === 'Enter' && (e.metaKey || e.ctrlKey)) { e.preventDefault(); draw(); } };
  const drawnIds = () => (net ? [...net.nodes].sort((a, b) => (b.q - a.q) || a.row.gene.localeCompare(b.row.gene)).map((d) => d.row.gene) : []);   // yours first, then the added, by name
  // Find: the proteins of the drawn network, with the iLIS of each one's strongest pair here; picking one centers it
  partnerSuggest($('#nw-find'), () => { if (!net) return []; const best = new Map();
    for (const l of net.links) { const a = typeof l.source === 'object' ? l.source.id : l.source, b = typeof l.target === 'object' ? l.target.id : l.target; for (const v of [a, b]) best.set(v, Math.max(best.get(v) || 0, l.best)); }
    return [...net.nodes].map((d) => ({ name: d.row.gene, best: best.get(d.id) || 0, dot: d.q ? '#1A5276' : '#AEBBCA', sty: '', tip: d.q ? 'of interest' : 'added partner' }))
      .sort((a, b) => a.name.localeCompare(b.name)); });
  $('#nw-find').onchange = () => { if (!net) return; const q = $('#nw-find').value.trim().toLowerCase(); const d = net.nodes.find((n) => n.row.gene.toLowerCase() === q) || net.nodes.find((n) => n.row.gene.toLowerCase().startsWith(q));
    if (d) net.focus(d.id); };
  $('#nw-copyids').onclick = async () => { const b = $('#nw-copyids'), t = drawnIds().join(', '); if (!t) return;
    try { await navigator.clipboard.writeText(t); b.textContent = `Copied ${fmtInt(drawnIds().length)}`; } catch (e) { window.prompt('Copy these proteins:', t); } setTimeout(() => { b.textContent = 'Copy proteins'; }, 1500); };
  $('#nw-useids').onclick = () => { const ids = drawnIds(); if (!ids.length) return; $('#nw-ids').value = ids.join(', '); S.col = null; S.exp = []; S.ncolUser = false;
    $('#nw-add').value = 'none'; showK(); window.scrollTo({ top: $('#nw-ids').getBoundingClientRect().top + scrollY - 90, behavior: 'smooth' }); draw(); };
  $('#nw-link').onclick = async () => { const b = $('#nw-link'); try { await navigator.clipboard.writeText(location.href); b.textContent = 'Copied'; } catch (e) { window.prompt('Copy this link:', location.href); } setTimeout(() => { b.textContent = 'Copy link'; }, 1500); };
  const rowsOut = () => (net ? net.links.map((l) => { const a = typeof l.source === 'object' ? l.source.id : l.source, b = typeof l.target === 'object' ? l.target.id : l.target; return { a, b, l }; }) : []);
  $('#nw-csv').onclick = () => { const csvq = (v) => (/[",\n]/.test(v) ? `"${String(v).replace(/"/g, '""')}"` : v);
    const text = 'Protein_1,Protein_2,iLIS_best,iLIS_avg,ipTM_best,screens,BioGRID_physical,BioGRID_genetic\n' + rowsOut().map(({ a, b, l }) => [csvq(gname(a)), csvq(gname(b)), l.best.toFixed(3), l.avg.toFixed(3), Number.isFinite(l.iptm) ? l.iptm.toFixed(2) : '', csvq(screens(l.src)), l.pubs || 0, l.gen || 0].join(',')).join('\n') + '\n';
    const u = URL.createObjectURL(new Blob([text], { type: 'text/csv' })), d = document.createElement('a'); d.href = u; d.download = `atlas_${sp.id}_network.csv`; d.click(); setTimeout(() => URL.revokeObjectURL(u), 3000); };
  $('#nw-graphml').onclick = () => { if (!net) return;
    const x = (v) => String(v).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
    const deg = new Map(); rowsOut().forEach(({ a, b }) => { deg.set(a, (deg.get(a) || 0) + 1); deg.set(b, (deg.get(b) || 0) + 1); });
    const keys = [['gene', 'node', 'string'], ['key', 'node', 'string'], ['accession', 'node', 'string'], ['group', 'node', 'string'], ['degree', 'node', 'int'],
      ...(S.data || []).map((x, n) => [`data${n}`, 'node', x.kind === 'num' ? 'double' : 'string', x.name]),
      ['type', 'edge', 'string'], ['iLIS_best', 'edge', 'double'], ['iLIS_avg', 'edge', 'double'], ['ipTM_best', 'edge', 'double'], ['screens', 'edge', 'string'], ['BioGRID_physical', 'edge', 'int'], ['BioGRID_genetic', 'edge', 'int']];
    const d = (k, v) => (v === '' || v == null || (typeof v === 'number' && !Number.isFinite(v)) ? '' : `<data key="${k}">${x(v)}</data>`);
    const nodesX = net.nodes.map((n) => `<node id="n${n.id}">${d('gene', n.row.gene)}${d('key', n.row.key)}${d('accession', n.row.acc || '')}${d('group', n.q ? 'query' : n.grew ? 'added by a click' : 'added partner')}${d('degree', deg.get(n.id) || 0)}${(S.data || []).map((c2, k) => (c2.vals.has(n.id) ? d(`data${k}`, c2.vals.get(n.id)) : '')).join('')}</node>`);
    const ends = (l) => [typeof l.source === 'object' ? l.source.id : l.source, typeof l.target === 'object' ? l.target.id : l.target];
    const edgesX = [...rowsOut().map(({ a, b, l }) => `<edge source="n${a}" target="n${b}">${d('type', 'predicted')}${d('iLIS_best', +l.best.toFixed(3))}${d('iLIS_avg', +l.avg.toFixed(3))}${d('ipTM_best', Number.isFinite(l.iptm) ? +l.iptm.toFixed(2) : '')}${d('screens', screens(l.src))}${d('BioGRID_physical', l.pubs || 0)}${d('BioGRID_genetic', l.gen || 0)}</edge>`),
      ...net.extra.map((l) => { const [a, b] = ends(l); const f = FOLD.get(fkey(a, b)); return `<edge source="n${a}" target="n${b}">${d('type', `reported, not predicted${!f ? '' : f.st === 'low' ? '; folded, below the cutoff' : f.st === 'none' ? '; never folded in these screens' : ''}`)}${f && f.st === 'low' ? d('iLIS_best', +f.best.toFixed(3)) : ''}${d('BioGRID_physical', l.pubs || 0)}${d('BioGRID_genetic', l.gen || 0)}</edge>`; })];
    const text = `<?xml version="1.0" encoding="UTF-8"?>\n<graphml xmlns="http://graphml.graphdrawing.org/xmlns">\n${keys.map(([id, f, t, nm]) => `<key id="${id}" for="${f}" attr.name="${x(nm || id)}" attr.type="${t}"/>`).join('\n')}\n`
      + `<graph id="atlas_${sp.id}_network" edgedefault="undirected">\n${nodesX.join('\n')}\n${edgesX.join('\n')}\n</graph>\n</graphml>\n`;
    const u = URL.createObjectURL(new Blob([text], { type: 'application/xml' })), a = document.createElement('a'); a.href = u; a.download = `atlas_${sp.id}_network.graphml`; a.click(); setTimeout(() => URL.revokeObjectURL(u), 3000); };
  $('#nw-livia').onclick = () => { if (!net) return; const w = window.open(`${LIVIA}network.html?post=1`, '_blank'); if (!w) return;
    const text = 'name,Symbol_1,Symbol_2,iLIS,ipTM\n' + rowsOut().map(({ a, b, l }) => `${gname(a)}___${gname(b)},${gname(a)},${gname(b)},${l.best.toFixed(3)},${Number.isFinite(l.iptm) ? l.iptm.toFixed(2) : ''}`).join('\n') + '\n';   // the pair table LIVIA's network reads (names already split)
    const keys = {}; for (const n of net.nodes) { keys[n.row.gene] = n.row.key; keys[String(n.row.gene).toUpperCase()] = keys[String(n.row.gene).toUpperCase()] || n.row.key; }   // LIVIA's cLIP opens each protein's Atlas page
    handTo(w, { type: 'livia-load', name: `atlas_${sp.id}_network.csv`, data: new Blob([text], { type: 'text/csv' }), ilis: cutV(), atlas: { url: location.href.split('#')[0], sp: sp.id, keys } }); };
  if (S.ids) draw(); else $('#nw-ids').focus();
}

/* ── router ─────────────────────────────────────────────────────────────────────────────────────────── */
// GoatCounter counts page loads only; the atlas routes by #/…, so each view is counted by hand, on the atlas's own
// counter (livia-atlas.goatcounter.com). A view is counted by its species only ('/human', '/virus', …; every other page as
// '/'): the counter shows how much each species is used, never which protein, pair or virus a reader opened.
const COUNT_SP = new Set(['human', 'fly', 'worm', 'zebrafish', 'yeast', 'virus']);
function trackView() {
  const sp = location.hash.replace(/^#\/?/, '').split(/[/?]/)[0], path = COUNT_SP.has(sp) ? `/${sp}` : '/';
  const count = () => window.goatcounter && window.goatcounter.count && window.goatcounter.count({ path, title: path === '/' ? 'LIVIA Atlas' : `${sp} · LIVIA Atlas` });
  if (window.goatcounter && window.goatcounter.count) count(); else window.addEventListener('load', () => setTimeout(count, 0), { once: true });
}
const hashPath = () => { const [path, q] = location.hash.replace(/^#\/?/, '').split('?'); return { parts: path.split('/').filter(Boolean).map(decodeURIComponent), q: new URLSearchParams(q || '') }; };
let LAST_PATH = null, holdTimer = null;
// Every page with three or more cards gets the section bar the protein page has: sticky, the card in view marked, a
// click jumps to it. Built after the page renders; the protein page keeps its own.
function mountSections() {
  if (app.querySelector('.subnav')) return;
  const cards = [...app.querySelectorAll(':scope > .card, :scope > .reading > .card')].filter((c) => !c.hidden && c.querySelector('h2'));
  if (cards.length < 3) return;
  const items = cards.map((c, i) => { if (!c.id) c.id = `sec-${i}`; const h = c.querySelector('h2'); return [c.id, ((h.childNodes[0] && h.childNodes[0].textContent) || h.textContent).trim()]; });
  const bar = el(`<nav class="subnav" aria-label="Sections">${items.map(([t, l]) => `<button data-t="${t}">${esc(l)}</button>`).join('')}</nav>`);
  cards[0].before(bar);
  bar.querySelectorAll('button').forEach((b) => b.onclick = () => { const t = document.getElementById(b.dataset.t); if (t) window.scrollTo({ top: t.getBoundingClientRect().top + window.scrollY - 112, behavior: 'auto' }); });
  let cur = null, raf = 0;
  const spy = () => { raf = 0; if (!bar.isConnected) { window.removeEventListener('scroll', onScroll); return; }
    const lim = bar.getBoundingClientRect().bottom + 24; let on = null;
    for (const b of bar.querySelectorAll('button')) { const t = document.getElementById(b.dataset.t); if (t && !t.hidden && t.getBoundingClientRect().top <= lim) on = b; }
    if (on === cur) return; cur = on; bar.querySelectorAll('button').forEach((b) => b.classList.toggle('on', b === on));
    if (on && bar.scrollWidth > bar.clientWidth + 1) bar.scrollLeft = Math.max(0, on.offsetLeft - (bar.clientWidth - on.offsetWidth) / 2); };
  const onScroll = () => { if (!raf) raf = requestAnimationFrame(spy); };
  window.addEventListener('scroll', onScroll, { passive: true }); requestAnimationFrame(spy);
}
// Any table without its own controls: a click on a column sorts by it (numbers as numbers), and a box above filters the
// rows by text. Tables that sort and filter themselves (header cells with data-k or data-c) are left alone.
function enhanceTables(root = app) {
  for (const t of root.querySelectorAll('table.pt, table.sets')) {
    if (t.dataset.enh || t.querySelector('th[data-k], th[data-c]') || !t.tBodies[0] || t.tBodies[0].rows.length < 2) continue;
    t.dataset.enh = '1';
    const wrap = t.closest('.tbl-wrap') || t, box = el('<input type="search" class="tbl-find" placeholder="Filter rows" aria-label="Filter the table">');
    wrap.before(box);
    box.oninput = () => { const q = box.value.trim().toLowerCase(); for (const r of t.tBodies[0].rows) r.hidden = !!q && !r.textContent.toLowerCase().includes(q); };
    if (!t.tHead || t.tHead.rows.length !== 1) continue;   // several header rows: filter only
    const num = (s) => { const v = parseFloat(String(s).replace(/[,%×]/g, '').replace(/^[^\d.-]+/, '')); return Number.isFinite(v) ? v : null; };
    [...t.tHead.rows[0].cells].forEach((th, k) => { th.classList.add('sortable'); th.title = th.title || 'click to sort';
      th.onclick = () => { const asc = th.dataset.dir !== 'asc'; [...t.tHead.rows[0].cells].forEach((x) => { delete x.dataset.dir; x.classList.remove('sorted', 'asc'); });
        th.dataset.dir = asc ? 'asc' : 'desc'; th.classList.add('sorted'); if (asc) th.classList.add('asc');
        const rows = [...t.tBodies[0].rows], val = (r) => (r.cells[k] ? r.cells[k].textContent.trim() : ''), allNum = rows.every((r) => val(r) === '' || val(r) === '—' || num(val(r)) != null);
        rows.sort((a, b) => { const x = val(a), y = val(b); const d = allNum ? (num(x) ?? -Infinity) - (num(y) ?? -Infinity) : x.localeCompare(y); return asc ? d : -d; });
        for (const r of rows) t.tBodies[0].appendChild(r); }; });
  }
}
async function route() {
  const gen = ++ROUTE, { parts, q } = hashPath(), setId = q.get('set') || '', here = location.hash;
  // The same page with another isoform or scope (?iso=, ?set=) keeps the reader where they were: the page holds its height
  // while it redraws, then returns to the same place. Any other page starts at the top.
  const path = parts.join('/'), stay = path === LAST_PATH, keepY = window.scrollY; LAST_PATH = path;
  clearTimeout(holdTimer);
  if (stay) app.style.minHeight = `${app.offsetHeight}px`; else { app.style.minHeight = ''; window.scrollTo(0, 0); }
  hideTip(); window.onresize = null;
  document.title = 'LIVIA Atlas';
  stopClip();
  try {
    await registry();
    if (stale(gen)) return;
    if (!parts.length) await viewHome();
    else if (parts[0] === 'datasets') await (parts[2] ? viewSet(parts[1], parts[2]) : parts[1] ? viewDataset(parts[1]) : viewDatasets());
    else if (parts[0] === 'about') viewAbout();
    else if (parts[0] === 'themes' && parts[1]) await viewTheme(parts[1]);
    else if (await regSpecies(parts[0])) { if (parts.length === 1) await viewSpecies(parts[0]); else if (parts[1] === 'network' && parts.length === 2) await viewNetwork(parts[0], q); else if (parts[1] === 'nested' && parts.length === 2) await viewNested(parts[0], q); else if (parts[1] === 'taxon' && parts.length === 3) await viewVirus(parts[0], parts[2]); else if (parts.length === 2) await viewProtein(parts[0], parts[1], setId, q.get('iso')); else await viewPair(parts[0], parts[1], parts[2], setId); }
    else if (await regDataset(parts[0])) {   // links from before the species pages: #/<screen>/<name>[/<name>]
      const d = await regDataset(parts[0]);
      if (parts.length === 1 || !d.species) { location.replace(`#/datasets/${d.id}`); return; }
      const sp = await species(d.species), k = (n) => { const r = sp.byName.get(n); return r ? r.key : n; };
      location.replace(`#/${sp.id}/${parts.slice(1).map(k).join('/')}`); return;
    }
    else app.innerHTML = `<div class="empty">Nothing at “${esc(parts.join('/'))}”. <a href="#/">Go to the atlas home</a></div>`;
  } catch (e) { if (!stale(gen)) app.innerHTML = `<div class="empty">${esc(e.message || e)}</div>`; console.error(e); }
  if (stale(gen)) return;
  if (stay) { window.scrollTo({ top: keepY, behavior: 'instant' }); holdTimer = setTimeout(() => { app.style.minHeight = ''; }, 4000); }
  if (location.hash === here) trackView();   // not for a view that redirected
  if (!stale(gen)) { mountSections(); enhanceTables(); setTimeout(() => { if (!stale(gen)) { mountSections(); enhanceTables(); } }, 2500); }   // cards and tables that fill in later
  if (!$('#top-search').firstChild) mountSearch($('#top-search'));
}
window.addEventListener('hashchange', route);
route();
// Column headers stay in view while a long table scrolls under the page header (and a protein page's section bar). The
// wrapper scrolls sideways, which stops CSS sticky, so the header cells are moved down instead, never past the last row.
{
  let raf = 0;
  const pin = () => {
    raf = 0;
    const line = Math.max(0, ...['.top', '.subnav'].map((q) => { const e = document.querySelector(q); return e ? e.getBoundingClientRect().bottom : 0; }));
    for (const t of document.querySelectorAll('table.pt')) {
      const r = t.getBoundingClientRect(), h = t.tHead ? t.tHead.offsetHeight : 0;
      const y = h && r.top < line ? Math.max(0, Math.min(line - r.top, r.height - 2 * h)) : 0;
      t.style.setProperty('--hy', `${y}px`); t.classList.toggle('pinned', y > 0);
    }
  };
  const later = () => { if (!raf) raf = requestAnimationFrame(pin); };
  addEventListener('scroll', later, { passive: true }); addEventListener('resize', later);
}
