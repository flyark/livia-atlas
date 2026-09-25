'use strict';
/* LIVIA cLIP Atlas — search AlphaFold-predicted protein interactions by protein.
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
const BAND = { 1: '#6D4FD1', 5: '#16956A', 10: '#C78B00', 0: '#A7B2BF' };
const bandOf = (v) => (v >= CUT[1] ? 1 : v >= CUT[5] ? 5 : v >= CUT[10] ? 10 : 0);
const bandLabel = { 1: '1% FPR', 5: '5% FPR', 10: '10% FPR', 0: 'below' };
// Benchmarked cutoffs at 10 / 5 / 1% FPR for single models and for the average over a pair's models: AFM-LIS
// thresholds_data_yfh_lipdockq.xlsx ("total group"; Y2H reference sets in yeast, fly and human — Kim et al. 2026, FlyPredictome).
const FPR = { iLIS: [0.223, 0.339, 0.551], ipTM: [0.48, 0.59, 0.72], iLIA: [620.3, 1247.4, 3078.8], iLISA: [143.9, 360.6, 1241.0],
  LIS: [0.168, 0.257, 0.439], cLIS: [0.298, 0.449, 0.716], ipSAE: [0.165, 0.363, 0.615], actifpTM: [0.745, 0.880, 0.963] };
const FPR_AVG = { iLIS: [0.072, 0.120, 0.268], ipTM: [0.292, 0.336, 0.442] };
const bandIn = (cuts, v) => (v >= cuts[2] ? 1 : v >= cuts[1] ? 5 : v >= cuts[0] ? 10 : 0);
const bandCol = (cuts, v) => BAND[bandIn(cuts, v)];   // a value's color = its FPR band under its own metric's cutoff
const ARCHIVE = { doi: '10.5281/zenodo.22964479', url: 'https://doi.org/10.5281/zenodo.22964479' };   // the atlas's data record: the concept DOI, always the latest version
const REF = {
  livia: ['Kim & Perrimon (2026) LIVIA, bioRxiv', '10.64898/2026.05.01.721633'],
  flypredictome: ['Kim et al. (2026) FlyPredictome, bioRxiv', '10.64898/2026.04.14.718529'],
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
async function screenFile(ds, rel) {
  const z = ds.reg.zip;
  if (!z) { const res = await fetch(new URL(rel.split('/').map(encodeURIComponent).join('/'), ds.manifest.bundleBase || ds.base).href); return res.ok ? res : null; }
  if (!OFFS.has(ds.id)) OFFS.set(ds.id, fetch(new URL(z.offsets, location.href).href).then((r) => (r.ok ? r.json() : {})).catch(() => ({})));
  const at = (await OFFS.get(ds.id))[rel]; if (!at) return null;
  const res = await fetch(DEV && z.dev ? new URL(z.dev, location.href).href : z.url, { headers: { Range: `bytes=${at[0]}-${at[0] + at[1] - 1}` } });
  if (res.status !== 206) { try { if (res.body) res.body.cancel(); } catch (e) { /* nothing to cancel */ } return null; }
  return res;
}
async function datasetRows(ds) {
  if (ds.rows) return ds.rows;
  const prot = await getJSON(ds.base + 'proteins.json'), c = Object.fromEntries(prot.columns.map((k, i) => [k, i]));
  ds.rows = prot.rows.map((r) => ({ id: r[c.id], gene: r[c.gene], acc: r[c.acc], partners: r[c.partners], pos10: r[c.pos10], pos5: r[c.pos5], pos1: r[c.pos1] }));
  return ds.rows;
}
async function species(id) {   // the species index: one row per protein over every screen, keyed by UniProt accession (fly: FlyBase gene)
  if (SPC[id]) return SPC[id];
  const reg = await regSpecies(id);
  if (!reg) throw new Error(`Unknown species “${id}”.`);
  const base = new URL(reg.base, location.href).href;
  const [manifest, prot] = await Promise.all([getJSON(base + 'manifest.json'), getJSON(base + 'proteins.json')]);
  const c = Object.fromEntries(prot.columns.map((k, i) => [k, i]));
  const rows = prot.rows.map((r, i) => {
    const occ = String(r[c.occ] || '').split(' ').filter(Boolean).map((t) => { const k = t.indexOf(':'); return { di: +t.slice(0, k), name: t.slice(k + 1) }; });
    return { i, key: r[c.key], gene: r[c.gene], acc: r[c.acc], name: r[c.name], syn: r[c.syn], len: r[c.len], clen: r[c.clen], status: r[c.status], occ,
      id: occ.length ? occ[0].name : r[c.key], partners: r[c.partners], pos10: r[c.pos10], pos5: r[c.pos5], pos1: r[c.pos1], best: r[c.bestIlis],
      src: occ.reduce((m, o) => m | (1 << o.di), 0) };
  });
  const byKey = new Map(rows.map((r) => [r.key, r])), byName = new Map(), byGene = new Map();
  for (const r of rows) { for (const o of r.occ) byName.set(o.name, r); if (!byGene.has(r.gene)) byGene.set(r.gene, r); }
  const keys = rows.map((r) => ({ gene: r.gene.toLowerCase(), key: r.key.toLowerCase(), acc: (r.acc || '').toLowerCase(), ids: r.occ.map((o) => o.name.toLowerCase()),
    syn: (r.syn || '').toLowerCase().split(/\s+/).filter(Boolean), name: (r.name || '').toLowerCase() }));
  const dsIds = manifest.datasets.map((d) => d.id), regs = await Promise.all(dsIds.map(regDataset));
  SPC[id] = { id, reg, base, manifest, rows, byKey, byName, byGene, keys, dsIds, dsShort: manifest.datasets.map((d) => d.short),
    dsColor: regs.map((r) => (r && r.color) || '#5B6B7F'), edges: null, cache: new Map() };
  return SPC[id];
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
    const lines = (await getText(url)).split('\n');
    const head = lines[0].split('\t'), cb = head.indexOf('iLIS_best'), ca = head.indexOf('iLIS_avg'), cs = head.indexOf('src');
    const adj = new Map();
    for (let k = 1; k < lines.length; k++) {
      if (!lines[k]) continue;
      const t = lines[k].split('\t'), i = map ? map[+t[0]] : +t[0], j = map ? map[+t[1]] : +t[1], e = { best: +t[cb], avg: +t[ca], src: cs >= 0 ? +t[cs] : 1 };
      if (i < 0 || j < 0 || i == null || j == null) continue;
      if (!adj.has(i)) adj.set(i, new Map()); if (!adj.has(j)) adj.set(j, new Map());
      adj.get(i).set(j, e); adj.get(j).set(i, e);
    }
    return { adj };
  })().catch((e) => { sp.edgesBy.delete(setId); throw e; }));   // a failed read is not kept: the next draw tries again
  return sp.edgesBy.get(setId);
}
function parseCSV(text) {
  const out = []; let row = [], f = '', q = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (q) { if (ch === '"') { if (text[i + 1] === '"') { f += '"'; i++; } else q = false; } else f += ch; }
    else if (ch === '"') q = true;
    else if (ch === ',') { row.push(f); f = ''; }
    else if (ch === '\n') { row.push(f); out.push(row); row = []; f = ''; }
    else if (ch !== '\r') f += ch;
  }
  if (f || row.length) { row.push(f); out.push(row); }
  return out;
}
const shiftRanges = (s, k) => (!k || !s || s === '[]' ? s : '[' + String(s).replace(/[\[\]\s]/g, '').split(',').filter(Boolean).map((t) => t.split('-').map((x) => +x + k).join('-')).join(',') + ']');
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
      const t = parseCSV(await zip.file(csvName).async('string'));
      return { rel, bytes, header: t[0], H: Object.fromEntries(t[0].map((k, i) => [k, i])), rows: t.slice(1).filter((r) => r.length > 10),
        seqs: parseFasta(faName ? await zip.file(faName).async('string') : ''), cons: conName ? parseCons(await zip.file(conName).async('string')) : null };
    })();
    ds.raw.set(rel, job); job.catch(() => ds.raw.delete(rel));
  }
  return ds.raw.get(rel);
}
// Every prediction of one protein across the screens of its species. A pair predicted in two screens, or both ways
// round in one, keeps every model; each keeps its screen and chain order (its "run") and rank. A gene-keyed bundle
// (FlyPredictome) holds every construct of the gene: its construct table says which construct belongs to which gene,
// and each prediction also keeps its screen category (set) and prediction run (batch).
function merged(sp, P, scope = '') {   // scope: one screen of the species (its dataset id) or one thematic set
  const ck = P.key + (scope ? '?' + scope : ''), onlyDi = scope ? sp.dsIds.indexOf(scope) : -1, setId = onlyDi >= 0 ? '' : scope;
  if (!sp.cache.has(ck)) {
    const job = (async () => {
      const parts = (await Promise.all(P.occ.filter((o) => onlyDi < 0 || o.di === onlyDi).map(async (o) => {   // a screen that cannot be reached is left out
        try { const ds = await dataset(sp.dsIds[o.di]); return { di: o.di, name: o.name, ds, raw: await bundleRaw(ds, o.name), TS: await setsOf(ds) }; } catch (e) { return null; }
      }))).filter(Boolean);
      if (!parts.length) throw new Error(`No interaction data for ${P.gene}.`);
      const preds = [], runs = new Map(), seqs = new Map(), cons = new Map(), sets = [], TS = (parts.find((x) => x.TS) || {}).TS || null;
      for (const part of parts) {
        const C = part.raw.cons;
        if (C) for (const [nm, c] of C) cons.set(nm, c);
        const keyOf = (nm) => { const c = C && C.get(nm); if (c) return c.key; const r = sp.byName.get(nm); return r ? r.key : nm; };
        for (const [nm, s] of part.raw.seqs) seqs.set(C ? nm : keyOf(nm), s);   // gene-keyed: the reference under the gene key, constructs by name
        const H = part.raw.H, num = (r, c) => (H[c] == null || r[H[c]] === '' ? NaN : +r[H[c]]);
        for (const r of part.raw.rows) {
          const nm = r[H.name], at = nm.indexOf('___'), a = nm.slice(0, at), b = nm.slice(at + 3), qi = C ? keyOf(a) === P.key : a === part.name;
          if (!qi && (C ? keyOf(b) !== P.key : b !== part.name)) continue;
          const tags = part.TS && H.batch != null ? part.TS.runSets.get(+r[H.batch]) || [] : null;   // the run's thematic sets
          if (setId && !(tags && tags.some((t) => t.id === setId))) continue;
          const ps = primarySet(tags), qc = qi ? a : b, pc = qi ? b : a, key = keyOf(pc), set = ps ? ps.short : H.set != null ? r[H.set] : '';
          const rid = part.di + '|' + nm + (H.batch != null ? '|' + r[H.batch] : '');
          if (!runs.has(rid)) runs.set(rid, { id: rid, di: part.di, qi, key, qc, pc, set, tags: tags ? tags.map((t) => t.id) : [] });
          if (set && !sets.includes(set)) sets.push(set);
          const s = (x, y) => (qi ? x : y);
          const p = { partner: key, run: rid, di: part.di, qi, qc, pc, set, tags: runs.get(rid).tags, rank: num(r, 'rank'), iLIS: num(r, 'iLIS'), iLIA: num(r, 'iLIA'), iLISA: num(r, 'iLISA'), ipTM: num(r, 'ipTM'),
            pTM: num(r, 'pTM'), LIS: num(r, 'LIS'), cLIS: num(r, 'cLIS'), LIA: num(r, 'LIA'), cLIA: num(r, 'cLIA'), ipSAE: num(r, 'ipSAE'), actifpTM: num(r, 'actifpTM'),
            qPl: num(r, s('pLDDT_i', 'pLDDT_j')), pPl: num(r, s('pLDDT_j', 'pLDDT_i')), qLIR: num(r, s('LIR_i', 'LIR_j')), pLIR: num(r, s('LIR_j', 'LIR_i')),
            qcLIR: num(r, s('cLIR_i', 'cLIR_j')), pcLIR: num(r, s('cLIR_j', 'cLIR_i')), qLen: num(r, s('len_i', 'len_j')), pLen: num(r, s('len_j', 'len_i')),
            qL: r[H[s('LIR_indices_i', 'LIR_indices_j')]], pL: r[H[s('LIR_indices_j', 'LIR_indices_i')]],
            qC: r[H[s('cLIR_indices_i', 'cLIR_indices_j')]], pC: r[H[s('cLIR_indices_j', 'cLIR_indices_i')]], row: r, hdr: part.raw.header };
          if (!Number.isFinite(p.iLISA)) p.iLISA = (p.iLIS || 0) * (p.iLIA || 0);
          preds.push(p);
        }
      }
      const byP = new Map();
      for (const p of preds) { if (!byP.has(p.partner)) byP.set(p.partner, []); byP.get(p.partner).push(p); }
      const partners = [...byP].map(([key, ps]) => {
        const ids = [...new Set(ps.map((p) => p.run))].sort((x, y) => runs.get(x).di - runs.get(y).di || (runs.get(y).qi ? 1 : 0) - (runs.get(x).qi ? 1 : 0));
        ps.sort((x, y) => ids.indexOf(x.run) - ids.indexOf(y.run) || x.rank - y.rank);
        const il = ps.map((p) => p.iLIS || 0), ip = ps.map((p) => p.ipTM || 0);
        return { id: key, row: sp.byKey.get(key) || null, preds: ps, runs: ids, src: ps.reduce((m, p) => m | (1 << p.di), 0), best: Math.max(...il), avg: mean(il),
          ilisaBest: Math.max(...ps.map((p) => p.iLISA || 0)), iptmBest: Math.max(...ip), iptmAvg: mean(ip), contacts: Math.max(...ps.map((p) => p.qcLIR || 0)),
          sets: [...new Set(ps.map((p) => p.set).filter(Boolean))] };
      });
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
        choices.sort((x, y) => y.n - x.n || (x.id ? 1 : -1));
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
        const rows = [], aside = new Map(), nameN = new Map();
        for (const p of preds) {
          if (!inClip(p)) { const c = cons.get(p.qc), k = c ? p.qc : 'screen ' + p.di; if (!aside.has(k)) aside.set(k, { di: p.di, len: p.qLen, con: c || null, n: 0 }); aside.get(k).n++; continue; }
          nameN.set(p.qc, (nameN.get(p.qc) || 0) + 1);
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
        return { choice, rows, qLen, qName, aside };
      };
      for (const p of preds) if (!(p.iLIS >= CUT[10])) { delete p.row; delete p.hdr; }   // only rows past the lowest cutoff are ever clustered
      const C0 = clipFor(choices.length ? choices[0].id : '');
      if (scope && !preds.length) throw new Error(`${P.gene} has no predictions in this ${onlyDi >= 0 ? 'screen' : 'set'}.`);
      return { parts, preds, partners, runs, seqs, cons, sets, TS, setId: scope, choices, clipFor, C0, clipRows: C0.rows, qLabel: P.key, labels, qLen: C0.qLen, qName: C0.qName, aside: C0.aside };
    })();
    sp.cache.set(ck, job); job.catch(() => sp.cache.delete(ck));
  }
  return sp.cache.get(ck);
}
const runLabel = (sp, B, rid, P) => { const ru = B.runs.get(rid); if (!ru) return ''; if (ru.label) return ru.label;
  const O = sp.byKey.get(ru.key), og = O ? O.gene : ru.key;
  return sp.dsShort[ru.di] + (ru.twin ? ` · ${ru.qi ? P.gene + '–' + og : og + '–' + P.gene}` : ''); };
// FlyPredictome's source categories (the paper's), as badges
const SET_COL = { 'Ligand–receptor': '#2B7A78', Kinase: '#B04A73', 'Kinase–TF': '#9A4A8F', 'Literature-derived': '#8C6D31', 'Large-scale proteomics': '#2B5F8E',
  'Subcellular organelle': '#5E7D2B', 'Inferred from orthologs': '#B5543C', Other: '#5B6B7F' };
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
  if (!AFDB.has(acc)) AFDB.set(acc, fetch(`https://alphafold.ebi.ac.uk/api/prediction/${encodeURIComponent(acc)}`).then((r) => (r.ok ? r.json() : null)).then((d) => {
    const list = Array.isArray(d) ? d : d ? [d] : [];
    const e = list.find((x) => (x.uniprotAccession === acc) && (x.uniprotStart || x.sequenceStart || 1) === 1) || list[0];
    return e && e.cifUrl ? { cifUrl: e.cifUrl, seq: e.uniprotSequence || e.sequence || '', amUrl: e.amAnnotationsUrl || null } : null;
  }).catch(() => null));
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
    <span>W</span><input list="lm-exp-dim" placeholder="auto"><span>H</span><input list="lm-exp-dim" placeholder="auto"><span>font ×</span><input list="lm-exp-font" value="1"></div>`);
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
    if (s) scored.push([s + Math.min(10, Math.log10(1 + sp.rows[i].pos10) * 3), i]);
  }
  scored.sort((a, b) => b[0] - a[0]);
  return scored.slice(0, limit).map(([score, i]) => ({ row: sp.rows[i], score }));
}
const findProteins = (sp, q, limit = 10) => scoreProteins(sp, q, limit).map((h) => h.row);
function mountSearch(host, { big = false, spId = null, autofocus = false, only = null, set = '' } = {}) {   // spId null: every species; only: a set's keys
  host.innerHTML = `<div class="search ${big ? 'big' : ''}">
      <svg class="glass" width="${big ? 20 : 17}" height="${big ? 20 : 17}" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round"><circle cx="11" cy="11" r="7"/><path d="m20 20-3.6-3.6"/></svg>
      <input type="search" placeholder="${big ? (spId ? 'Gene, accession or protein name' : 'Gene, UniProt accession, FlyBase ID or protein name') : 'Search a protein'}" aria-label="Search proteins" autocomplete="off" spellcheck="false">
      <div class="suggest" hidden></div></div>`;
  const input = $('input', host), box = $('.suggest', host);
  let items = [], on = -1, timer = null;
  const go = (it) => { box.hidden = true; input.value = ''; location.hash = `#/${it.sp.id}/${it.row.key}${set ? '?set=' + encodeURIComponent(set) : ''}`; };
  const paint = () => { [...box.children].forEach((c, k) => c.classList.toggle('on', k === on)); };
  async function update() {
    const q = input.value, sps = spId ? [await species(spId)] : await Promise.all(((await registry()).species || []).map((x) => species(x.id).catch(() => null)));
    const many = sps.filter(Boolean).length > 1;
    items = sps.filter(Boolean).flatMap((sp) => scoreProteins(sp, q, only ? 400 : 10).filter((h) => !only || only.has(h.row.key)).map((h) => ({ ...h, sp })))
      .sort((x, y) => y.score - x.score).slice(0, 10);
    on = items.length ? 0 : -1;
    box.innerHTML = items.map(({ row: r, sp }) => `<div class="sg"><b>${esc(r.gene)}</b><span class="nm">${esc(short(r.name))}</span><span class="ct">${fmtInt(r.pos10)} / ${fmtInt(r.partners)}</span>
        <span class="sub">${many ? `<span class="sp-tag">${esc(sp.reg.label)}</span>` : ''}${esc(r.acc || '—')} · ${esc(r.id)}</span></div>`).join('')
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
  document.addEventListener('click', (e) => { if (!host.contains(e.target)) box.hidden = true; });
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
  g.save(); g.fillStyle = '#6B7A8D'; g.strokeStyle = '#B9C4CF'; g.lineWidth = 1; g.font = '10.5px "IBM Plex Mono", ui-monospace, monospace'; g.textBaseline = 'top'; g.textAlign = 'left';
  for (const r of ticks) { const x = Math.round(xc(r)) + 0.5, s = String(r), w = g.measureText(s).width;
    g.beginPath(); g.moveTo(x, y); g.lineTo(x, y + 4); g.stroke(); g.fillText(s, Math.max(1, Math.min(W - w - 1, x - w / 2)), y + 6); }
  g.restore();
}
function lanes(items, x0Of, x1Of) {   // stack overlapping boxes into lanes; returns lane count
  const end = []; for (const d of items) { const x0 = x0Of(d); let l = 0; while (l < end.length && end[l] > x0 + 0.5) l++; d.lane = l; end[l] = x1Of(d); } return end.length;
}
// Interface of one prediction on both proteins: gray = not in the interface, light = LIR, dark = cLIR; UniProt domains above.
function drawIfaceTracks(cv, tracks) {
  const W = cv.parentElement.clientWidth, pad = 8, gap = 18;
  for (const t of tracks) { t.x = (r) => pad + (r - 1) / t.len * (W - 2 * pad); t.nl = t.doms.length ? lanes(t.doms, (d) => t.x(d.start), (d) => t.x(d.end + 1)) : 0; }
  const bh = (t) => 22 + t.nl * 15 + (t.nl ? 4 : 0) + 22 + 20;
  const H = tracks.reduce((s, t) => s + bh(t), 0) + gap * (tracks.length - 1);
  const g = canvasCtx(cv, W, H);
  let y = 0; cv._blocks = [];
  for (const t of tracks) {
    const L = t.len, bw = Math.max(1.3, (W - 2 * pad) / L), name = t.label || t.gene;
    g.font = '600 13.5px "IBM Plex Sans", system-ui, sans-serif'; g.fillStyle = t.col.clir; g.textBaseline = 'alphabetic'; g.textAlign = 'left'; g.fillText(name, pad, y + 15);
    const gw = g.measureText(name).width; g.font = '12px "IBM Plex Mono", ui-monospace, monospace'; g.fillStyle = '#6B7A8D';
    const extent = t.span ? `residues ${fmtInt(t.span[0])}–${fmtInt(t.span[1])} of ${fmtInt(L)}` : t.own ? `${fmtInt(L)} aa construct, its own numbering` : `${fmtInt(L)} aa`;
    g.fillText(`${extent} · ${t.lir.length} interface · ${t.clir.length} contact`, pad + gw + 10, y + 15);
    y += 22;
    for (const d of t.doms) { const x0 = t.x(Math.max(1, d.start)), x1 = Math.max(x0 + 2, t.x(Math.min(L, d.end) + 1)), yy = y + d.lane * 15;
      g.fillStyle = '#E3E9F1'; g.fillRect(x0, yy, x1 - x0, 13); g.strokeStyle = '#9FB0C4'; g.lineWidth = 0.6; g.strokeRect(x0 + 0.3, yy + 0.3, x1 - x0 - 0.6, 12.4);
      g.font = '10.5px "IBM Plex Sans", system-ui, sans-serif'; g.fillStyle = '#34445A'; let s = d.name; while (s.length > 2 && g.measureText(s).width > x1 - x0 - 6) s = s.slice(0, -2) + '…';
      if (s.length > 2 && g.measureText(s).width <= x1 - x0 - 6) g.fillText(s, x0 + 3, yy + 10); }
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
    const blk = cv._blocks.find((k) => yy >= k.y0 - 4 && yy <= k.y1 + 4); if (!blk) return hideTip();
    const t = blk.t, r = Math.floor((x - pad) / (W - 2 * pad) * t.len) + 1; if (r < 1 || r > t.len) return hideTip();
    const st = t.span && (r < t.span[0] || r > t.span[1]) ? 'not in the folded construct' : t.clirSet.has(r) ? 'contact (cLIR)' : t.lirSet.has(r) ? 'interface (LIR)' : 'not in the interface';
    const dom = t.doms.filter((d) => r >= d.start && r <= d.end).map((d) => d.name).join(', ');
    showTip(`<b>${esc(t.gene)}</b> ${t.seq && t.seq[r - 1] ? t.seq[r - 1] : ''}${r} · ${st}${dom ? `<br>${esc(dom)}` : ''}`, e.clientX, e.clientY); };
  cv.onmouseleave = hideTip;
}
// One sequence as a continuous, searchable flow in 10-residue groups (position numbers are CSS, not text), as clip.html.
function seqPanel(label, seq, lir, clir, col, len, span) {   // span: the folded part of a longer gene; the rest is dimmed
  const ext = span ? `residues ${fmtInt(span[0])}–${fmtInt(span[1])} of ${fmtInt(len)}` : `${fmtInt(len || seq.length)} residues`;
  const head = `<div class="seqh"><b style="color:${col.clir}">${esc(label)}</b> <span class="muted">${ext} · ${lir.length} interface · ${clir.length} contact</span></div>`;
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
  host.innerHTML = `${none ? `<p class="note">This model has no residue pair with PAE ≤ 12 Å between the two proteins, so it has no interface.</p>` : ''}
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
  // UniProt domain coordinates fit only the UniProt sequence: a full-length human construct; for fly, a FlyBase
  // reference that is UniProt's sequence (checked against the AlphaFold DB entry)
  const domOK = async (R, s, seq) => { if (!R || !R.acc || s.own) return false;
    if (!sp.manifest.keyedBy) return !!(R.len && R.len === s.len);
    const e = await afdbEntry(R.acc); return !!(e && seq && e.seq === seq); };
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
  g.fillStyle = '#6B7A8D'; g.font = '11px "IBM Plex Mono", ui-monospace, monospace';
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
const TRY = { human: ['TP53', 'MDM2', 'CTNNB1', 'MAPK3', 'SMAD4', 'KRAS'], fly: ['arm', 'dsh', 'Ras85D', 'N', 'yki', 'Akt'] };   // example proteins on the home page
async function viewHome() {
  const gen = ROUTE, reg = await registry(), SPX = reg.species[0], sp = await species(SPX.id);
  const all = await Promise.all((reg.species || []).map((x) => species(x.id)));   // totals over every species and screen: they grow as screens are added
  if (stale(gen)) return;
  const tot = (k) => all.reduce((s, x) => s + (x.manifest.counts[k] || 0), 0), nScreens = all.reduce((s, x) => s + x.dsIds.length, 0);
  app.innerHTML = `
    <section class="hero hero-center">
      <h1>Every predicted <em>partner</em>, down to the <em>residue</em>.</h1>
      <p class="lede">Search AlphaFold-Multimer interaction screens by protein. Every prediction is scored with iLIS and its interface is resolved
        to residues: which partners a protein is predicted to bind, how confidently, and where.</p>
      <div id="home-search" class="hero-search"></div>
      <div class="totals"><span><b>${fmtInt(tot('runs'))}</b> predictions</span><span><b>${fmtInt(tot('predictions'))}</b> models</span><span><b>${fmtInt(tot('pairs'))}</b> protein pairs</span>
        <span><b>${fmtInt(tot('proteins'))}</b> proteins</span><span><b>${fmtInt(nScreens)}</b> screen${nScreens === 1 ? '' : 's'}</span></div>
      <div class="chips"><span class="lbl">Try</span>${(reg.species || []).map((x, i) => (reg.species.length > 1 ? `<span class="lbl${i ? ' sp' : ''}">${esc(x.label)}</span>` : '')
        + (TRY[x.id] || []).map((g) => `<a class="chip" data-sp="${x.id}" data-g="${g}">${g}</a>`).join('')).join('')}</div>
      <a class="showcase" id="showcase" href="#/${sp.id}/P04637" style="color:inherit;text-decoration:none">
        <div class="sc-top"><b>TP53</b><span id="sc-note">clustering its partners…</span></div>
        <canvas id="sc-cv" height="150"></canvas>
        <div class="sc-foot"><span>Contact residue frequency, colored by cluster</span><span style="color:var(--blue)">Open TP53 →</span></div></a>
    </section>
    ${(REG.themes || []).length ? '<h2 class="section-h">Themes</h2><div class="datasets live-row" id="themes"></div>' : ''}
    <h2 class="section-h">Datasets</h2>
    <div class="datasets live-row">${reg.datasets.filter((d) => d.status !== 'planned').map(dsCard).join('')}</div>
    <h2 class="section-h">How it works</h2>
    <div class="steps">
      <div class="step"><h4>1 · Predict</h4><p>A screen of protein pairs folded with AlphaFold-Multimer (or any predictor that reports PAE).</p></div>
      <div class="step"><h4>2 · Score</h4><p>lis.py scores each prediction over its confident residue pairs only: iLIS, with cutoffs benchmarked at 10%, 5% and 1% false-positive rate (${cite('flypredictome')}).</p></div>
      <div class="step"><h4>3 · Resolve</h4><p>The contact residues of every interface are kept, so partners can be compared by where they bind and clustered by their interaction fingerprints with LIVIA cLIP.</p></div>
    </div>`;
  mountSearch($('#home-search'), { big: true, autofocus: true });
  showcase(sp, sp.byKey.get('P04637'));
  app.querySelectorAll('.chip[data-g]').forEach((a) => a.onclick = async () => { const s2 = await species(a.dataset.sp), r = resolveRow(s2, a.dataset.g); if (r) location.hash = `#/${s2.id}/${r.key}`; });
  fillDsStats(); fillThemes();
}
async function fillThemes() {   // home: each theme's species and totals, from its members' counts
  const reg = await registry(), box = $('#themes'); if (!box) return;
  const cards = await Promise.all((reg.themes || []).map(async (T) => { const rows = await themeMembers(T), sum = (k) => rows.reduce((a, r) => a + (r.counts[k] || 0), 0);
    return `<div class="ds live"><span class="badge on">Theme</span><h3><a href="#/themes/${T.id}">${esc(T.title)}</a></h3>
      <div class="sp">${rows.map((r) => esc(r.S.label)).join(' · ')}</div>
      <div class="stats"><div><b>${fmtInt(sum('proteins'))}</b><span>proteins</span></div><div><b>${fmtInt(sum('pairs'))}</b><span>pairs</span></div><div><b>${fmtInt(sum('pairsFpr10'))}</b><span>past 10% FPR</span></div></div>
      <div class="src-line">${esc(T.about || '')}</div></div>`; }));
  box.innerHTML = cards.join('');
}
function dsCard(d) {
  const live = d.status === 'live', spx = ((REG && REG.species) || []).find((s) => s.id === d.species);
  const head = `<span class="badge ${live ? 'on' : 'soon'}">${live ? 'Searchable' : d.status === 'external' ? 'Separate site' : 'Planned'}</span>
    <h3>${live ? `<a href="#/datasets/${d.id}">${esc(d.title)}</a>` : d.url ? `<a href="${esc(d.url)}" target="_blank" rel="noopener">${esc(d.title)}</a>` : esc(d.title)}</h3>
    <div class="sp">${esc(spx ? spx.name : d.species)}${live && d.short ? ` · <span class="src" style="--c:${d.color}">${esc(d.short)}</span>` : ''}</div>`;
  const src = d.paper ? `<a href="${esc(d.paper)}" target="_blank" rel="noopener">${esc(d.source)} ↗</a>` : esc(d.source);   // every paper reference links to the paper
  return `<div class="ds ${live ? 'live' : ''}" data-ds="${d.id}">${head}${live ? '<div class="stats"></div>' : ''}
    <div class="src-line">${src}${d.note ? ` — ${esc(d.note)}` : ''}</div></div>`;
}
async function fillDsStats() {
  for (const d of (await registry()).datasets.filter((x) => x.status === 'live')) {
    let s; try { s = await dataset(d.id); } catch (e) { continue; }
    const k = s.manifest.counts, box = app.querySelector(`[data-ds="${d.id}"] .stats`);
    if (box) box.innerHTML = `<div><b>${fmtInt(k.proteins)}</b><span>proteins</span></div><div><b>${fmtInt(k.pairs)}</b><span>pairs</span></div><div><b>${fmtInt(k.pairsFpr10)}</b><span>past 10% FPR</span></div>`;
    const TS = await setsOf(s).catch(() => null), scr = TS ? TS.list.filter((x) => x.type === 'screen') : [];   // its named screens, as thematic sets
    if (box && scr.length && !box.parentElement.querySelector('.setchips')) box.insertAdjacentHTML('afterend', `<div class="setchips"><span>Sets</span>${scr.map((x) =>
      `<a class="src" style="--c:${x.color}" href="#/datasets/${d.id}/${x.id}">${esc(x.short)}</a>`).join('')}<a href="#/datasets/${d.id}">all ${TS.list.length} ›</a></div>`);
  }
}
async function showcase(sp, P) {
  if (!P) return;
  const gen = ROUTE;
  try {
    const B = await merged(sp, P), M = await runClip(B.clipRows, B.qLabel, CUT[10]);
    const cv = $('#sc-cv'); if (!cv || stale(gen)) return;
    const L = M.plen || P.clen, n = M.fingerprints.length, tot = new Uint16Array(L + 2), byC = new Map();
    for (let i = 0; i < n; i++) for (const r of M.fingerprints[i]) { if (r < 1 || r > L) continue; tot[r]++; const m = byC.get(r) || {}; m[M.labels[i]] = (m[M.labels[i]] || 0) + 1; byC.set(r, m); }
    let max = 1; for (let r = 1; r <= L; r++) max = Math.max(max, tot[r]);
    const W = cv.clientWidth, H = 150, g = canvasCtx(cv, W, H);
    const x = (r) => (r - 1) / L * W, bw = Math.max(1, W / L), base = H - 18;
    g.fillStyle = '#F4F7FA'; g.fillRect(0, 0, W, base);
    for (let r = 1; r <= L; r++) { if (!tot[r]) continue; let dom = 1, b = -1; const cc = byC.get(r); for (const c in cc) if (cc[c] > b) { b = cc[c]; dom = +c; }
      const h = tot[r] / max * (base - 6); g.fillStyle = clusterColor(dom, M.k); g.fillRect(x(r), base - h, bw, h); }
    g.fillStyle = '#8593A5'; g.font = '10.5px "IBM Plex Mono", monospace'; g.textAlign = 'left'; g.fillText('1', 0, H - 4); g.textAlign = 'right'; g.fillText(String(L), W, H - 4);
    const parts = new Set(M.preds.map((p) => (B.labels.get(p.partner) || {}).key)).size;
    $('#sc-note').textContent = `${M.k} clusters · ${fmtInt(parts)} partners past 10% FPR`;
  } catch (e) { const n = !stale(gen) && $('#sc-note'); if (n) n.textContent = ''; }
}

async function viewDatasets() {
  const gen = ROUTE, reg = await registry();
  if (stale(gen)) return;
  app.innerHTML = `<div class="crumbs"><a href="#/">Atlas</a> / Datasets</div>
    <h2 class="section-h" style="margin-top:4px">Datasets</h2>
    <div class="datasets live-row">${reg.datasets.filter((d) => d.status !== 'planned').map(dsCard).join('')}</div>`;
  fillDsStats();
}

function viewAbout() {
  app.innerHTML = `<div class="reading"><div class="crumbs"><a href="#/">Atlas</a> / About</div>
    <div class="card" style="margin-top:6px"><h2>What this is</h2>
      <p style="max-width:78ch">LIVIA cLIP Atlas makes large AlphaFold-Multimer interaction screens searchable at the level of residues. Every prediction is scored with
      <b>iLIS</b>, the integrated local interaction score computed by lis.py over residue pairs with predicted aligned error of at most 12 Å (LIS), and over those that
      are also in contact, Cβ–Cβ distance of at most 8 Å (cLIS): iLIS = √(LIS × cLIS). Benchmarked cutoffs mark predictions at a 10%, 5% and 1% false-positive rate
      (iLIS ≥ 0.223, 0.339, 0.551; ${cite('flypredictome')}).</p>
      <table class="cuts"><caption>Benchmarked cutoffs, from the AFM-LIS benchmark on Y2H reference sets in yeast, fly and human
        (${cite('flypredictome')}; <a href="https://github.com/flyark/AFM-LIS" target="_blank" rel="noopener">AFM-LIS</a>)</caption>
        <thead><tr><th></th><th>10% FPR</th><th>5% FPR</th><th>1% FPR</th></tr></thead>
        <tbody>${[['iLIS, best model', FPR.iLIS, 3], ['iLIS, average over models', FPR_AVG.iLIS, 3], ['ipTM, best model', FPR.ipTM, 2], ['ipTM, average over models', FPR_AVG.ipTM, 3]]
          .map(([l, c, d]) => `<tr><th>${l}</th>${c.map((v, j) => `<td style="color:${BAND[[10, 5, 1][j]]}">≥ ${v.toFixed(d)}</td>`).join('')}</tr>`).join('')}</tbody></table>
      <p style="max-width:78ch">Each protein has one page per species that gathers its predictions from every screen. A pair predicted in two screens, or both ways round,
      keeps every model with its source. The interface residues on both proteins are kept for every prediction, and each protein page runs
      <a href="${LIVIA}clip.html" target="_blank" rel="noopener">LIVIA cLIP</a> in the browser: partners are clustered by their interaction fingerprints, and the clusters
      are mapped onto the AlphaFold DB structure with pLDDT and, for human proteins, AlphaMissense.</p>
      <p style="max-width:78ch">Fly pages are organized by FlyBase gene. FlyPredictome folded many constructs of a gene: isoforms, fragments, phosphosite windows,
      point mutants. Every name is resolved to its gene, and a construct is placed on the gene's reference sequence when its sequence is known, so its
      contacts are drawn in the gene's residue numbering. Clustering uses the full-length construct; the others are listed with their partners.
      The table of construct names and their genes is on the FlyPredictome page.</p></div>
    <div class="card"><h2>Cite</h2>
      <p style="max-width:78ch">Kim, A.-R. &amp; Perrimon, N. (2026). LIVIA: a browser-based tool for assessing and visualizing predicted protein interactions. bioRxiv.
      <a href="https://doi.org/${REF.livia[1]}" target="_blank" rel="noopener">doi.org/${REF.livia[1]}</a></p>
      <p style="max-width:78ch">iLIS cutoffs: Kim, A.-R. et al. (2026). FlyPredictome: a structural atlas of predicted protein-protein interactions in <i>Drosophila</i>. bioRxiv.
      <a href="https://doi.org/${REF.flypredictome[1]}" target="_blank" rel="noopener">doi.org/${REF.flypredictome[1]}</a></p>
      <p class="muted" style="font-size:14px">Each dataset page lists the screen it comes from; please cite that source too.</p></div></div>`;
}

// The counts of a screen, a species or a set, in one row: a prediction is a pair folded once (its ranked models are
// counted as models); protein pairs are unique pairs, a protein with itself left out
const kpiRow = (k) => `<div class="kpirow"><div class="kpi"><b>${fmtInt(k.proteins)}</b><span>proteins</span></div><div class="kpi"><b>${fmtInt(k.pairs)}</b><span>protein pairs</span></div>
  ${k.runs ? `<div class="kpi"><b>${fmtInt(k.runs)}</b><span>predictions</span></div>` : ''}<div class="kpi"><b>${fmtInt(k.predictions)}</b><span>models</span></div>
  <div class="kpi f10"><b>${fmtInt(k.pairsFpr10)}</b><span>pairs past 10% FPR</span></div><div class="kpi f5"><b>${fmtInt(k.pairsFpr5)}</b><span>past 5% FPR</span></div>
  <div class="kpi f1"><b>${fmtInt(k.pairsFpr1)}</b><span>past 1% FPR</span></div></div>`;
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
  document.title = `${T.title} · LIVIA cLIP Atlas`;
  const rows = await themeMembers(T), sum = (k) => rows.reduce((a, r) => a + (r.counts[k] || 0), 0);
  if (stale(gen)) return;
  const cite = (src) => (src && src.url ? `<a href="${esc(src.url)}" target="_blank" rel="noopener">${esc(shortCite(src))} ↗</a>` : '');
  app.innerHTML = `<div class="crumbs"><a href="#/">Atlas</a> / Themes / ${esc(T.title)}</div>
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
  document.title = `${S.short} · ${ds.reg.short} · LIVIA cLIP Atlas`;
  const prot = await getJSON(new URL(S.files.proteins, ds.base).href), c = Object.fromEntries(prot.columns.map((k, i) => [k, i]));
  if (stale(gen)) return;
  const rows = prot.rows.map((r) => ({ key: r[c.key], pos10: r[c.pos10] })), keys = new Set(rows.map((r) => r.key));
  const hubs = [...rows].sort((a, b) => b.pos10 - a.pos10).slice(0, 24), k = S.counts;
  const link = (u, t) => `<a href="${esc(u)}" target="_blank" rel="noopener">${esc(t)} ↗</a>`;
  app.innerHTML = `<div class="crumbs"><a href="#/">Atlas</a> / <a href="#/datasets">Datasets</a> / <a href="#/datasets/${ds.id}">${esc(ds.reg.title)}</a> / ${esc(S.short)}</div>
    <div class="dshead"><h1>${esc(S.title)}</h1>
      <div class="pname"><span class="src" style="--c:${S.color}">${esc(S.short)}</span> ${S.type === 'screen' ? 'A screen' : 'A source category'} within ${esc(ds.reg.title)} · <i>${esc(m.species.name)}</i></div>
      ${S.source ? `<div class="cite">${link(S.source.url, `${S.source.citation} doi:${S.source.doi}`)}</div>` : ''}
      <div class="cite">Part of ${m.source.url ? link(m.source.url, m.source.citation) : esc(m.source.citation)}</div></div>
    ${kpiRow(k)}
    <div class="card"><h2>Search this set</h2><p class="muted" style="margin:2px 0 10px">Protein pages opened from here show only this set's predictions, with a switch to all of ${esc(ds.reg.short)}.</p><div id="set-search"></div></div>
    <div class="card"><h2>Most connected proteins in this set <span class="muted">partners past the 10% FPR cutoff</span></h2>
      <div class="chips">${hubs.map((r) => { const R = sp.byKey.get(r.key); return `<a class="chip" href="#/${sp.id}/${r.key}?set=${S.id}">${esc(R ? R.gene : r.key)} <span class="num" style="color:var(--ink-3)">${fmtInt(r.pos10)}</span></a>`; }).join('')}</div></div>
    <div class="card"><h2>Data</h2><p class="muted" style="font-size:14px;margin:4px 0 0">This set is part of the ${esc(ds.reg.title)} download in the
      <a href="${ARCHIVE.url}" target="_blank" rel="noopener">LIVIA cLIP Atlas record on Zenodo ↗</a>;
      <a href="https://github.com/flyark/livia-atlas/blob/main/tools/extract_set.py" target="_blank" rel="noopener">extract_set.py ↗</a> pulls out just this set as a table.</p></div>`;
  mountSearch($('#set-search'), { spId: sp.id, only: keys, set: S.id });
}
async function viewDataset(dsId) {   // one screen: what it is, its counts and files; proteins link to their species pages
  const gen = ROUTE, ds = await dataset(dsId), m = ds.manifest, k = m.counts, sp = await species(ds.reg.species), TS = await setsOf(ds);
  const rows = await datasetRows(ds), hubs = [...rows].sort((a, b) => b.pos10 - a.pos10).slice(0, 24);
  if (stale(gen)) return;
  const scopeQ = sp.dsIds.length > 1 ? `?set=${encodeURIComponent(ds.id)}` : '';   // one of several screens: its protein pages open in its scope
  app.innerHTML = `<div class="crumbs"><a href="#/">Atlas</a> / <a href="#/datasets">Datasets</a> / ${esc(ds.reg.title)}</div>
    <div class="dshead"><h1>${esc(ds.reg.title)}</h1><div class="pname"><i>${esc(m.species.name)}</i> · ${esc(m.source.method)} · ${esc(m.analysis.tool)}, PAE ≤ ${m.analysis.paeCutoff} Å, Cβ ≤ ${m.analysis.cbCutoff} Å</div>
      <div class="cite">${m.source.url ? `<a href="${esc(m.source.url)}" target="_blank" rel="noopener">${esc(m.source.citation)}${m.source.doi ? ` doi:${esc(m.source.doi)}` : ''} ↗</a>` : esc(m.source.citation)}</div></div>
    ${kpiRow(k)}
    <div class="card"><h2>Search</h2><p class="muted" style="margin:2px 0 10px">${sp.dsIds.length > 1 ? `Protein pages opened from here show only this screen, <span class="src" style="--c:${ds.reg.color}">${esc(ds.reg.short)}</span>, with a switch to every ${esc(sp.reg.label.toLowerCase())} screen.`
      : `Protein pages show every prediction of this screen${TS ? ', with a switch to each of its thematic sets' : ''}.`}</p><div id="ds-search"></div></div>
    <div class="card"><h2>Most connected proteins in this screen <span class="muted">partners past the 10% FPR cutoff</span></h2>
      <div class="chips">${hubs.map((r) => { const R = sp.byName.get(r.id); return `<a class="chip" href="#/${sp.id}/${R ? R.key : r.id}${scopeQ}">${esc(r.gene)} <span class="num" style="color:var(--ink-3)">${fmtInt(r.pos10)}</span></a>`; }).join('')}</div></div>
    ${setsCard(TS)}
    <div class="card"><h2>Data</h2><p class="muted" style="font-size:14px;margin:4px 0 0">Every file of this screen is in the
      <a href="${ARCHIVE.url}" target="_blank" rel="noopener">LIVIA cLIP Atlas record on Zenodo (doi:${ARCHIVE.doi}) ↗</a>. A protein page's
      <b>Data</b> menu downloads that protein's predictions.${m.files.identity ? ` <a href="${ds.base}${m.files.identity}" download>Construct names and their FlyBase genes</a> (table).` : ''}</p></div>`;
  mountSearch($('#ds-search'), sp.dsIds.length > 1 ? { spId: sp.id, set: ds.id, only: new Set(rows.map((r) => { const R = sp.byName.get(r.id); return R ? R.key : r.id; })) } : { spId: sp.id });
}

/* ── protein page: LIVIA cLIP, natively, over every screen, with a partner overview, a network and a partner table ── */
let CLIPW = null, clipSeq = 0; const clipWait = new Map();
function runClip(rows, gene, cut) {
  if (!CLIPW) { CLIPW = new Worker('clipworker.js?v=20260925h'); CLIPW.onmessage = (e) => { const w = clipWait.get(e.data.id); if (w) { clipWait.delete(e.data.id); e.data.ok ? w.resolve(e.data) : w.reject(new Error(e.data.message)); } }; }
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
const ECOL = () => d3.scaleLinear().domain([CUT[10], CUT[5], CUT[1], 0.85]).range(['#E0AE2E', '#16956A', '#6D4FD1', '#2A1B7A']).interpolate(d3.interpolateLab).clamp(true);
const EWID = (a) => 0.5 + 4.3 * Math.max(0, Math.min(1, a / 0.8));
function resolveRow(sp, q) {   // a key, any screen's name, an accession, a gene symbol (exact case first), a CG number or an older name
  if (sp.byKey.has(q)) return sp.byKey.get(q);
  if (sp.byName.has(q)) return sp.byName.get(q);
  if (sp.byGene.has(q)) return sp.byGene.get(q);
  let i = sp.rows.findIndex((r) => (' ' + (r.syn || '') + ' ').includes(' ' + q + ' '));   // an older name, exact case
  if (i >= 0) return sp.rows[i];
  const l = q.toLowerCase(); i = sp.keys.findIndex((k) => k.gene === l || k.acc === l || k.ids.includes(l));
  if (i < 0) i = sp.keys.findIndex((k) => k.syn.includes(l));
  return i >= 0 ? sp.rows[i] : null;
}

async function viewProtein(spId, q, setId = '') {   // setId: show only that thematic set's predictions
  const gen = ROUTE, gone = () => stale(gen);   // after every await: stop if the reader has moved to another page
  const sp = await species(spId), P = resolveRow(sp, q), qs = setId ? `?set=${encodeURIComponent(setId)}` : '';
  if (gone()) return;
  if (!P) { app.innerHTML = `<div class="empty">No protein “${esc(q)}” in the ${esc(sp.reg.label.toLowerCase())} screens. <a href="#/${sp.id}">Search ${esc(sp.reg.label.toLowerCase())} proteins</a></div>`; return; }
  if (P.key !== q) { location.replace(`#/${sp.id}/${P.key}${qs}`); return; }
  document.title = `${P.gene} · LIVIA cLIP Atlas`;
  const flags = [];
  if (P.status === 'renamed') flags.push(`<span class="flag">named ${esc(P.occ.map((o) => o.name).filter((n, i, a) => a.indexOf(n) === i).join(' / '))} in the screens; UniProt renamed it</span>`);
  if (P.status === 'unreviewed') flags.push('<span class="flag">unreviewed UniProt entry</span>');
  if (P.status === 'other species') flags.push('<span class="flag">not a fly protein: folded as a partner of fly proteins</span>');
  if (P.status === 'construct') flags.push('<span class="flag">an engineered construct or a retired gene, kept under its screen name</span>');
  if (P.status === 'obsolete') flags.push('<span class="flag">UniProt has since retired this entry; the sequence is the one the screen folded</span>');
  const fbLink = /^FBgn\d{7}$/.test(P.key) ? `<a href="https://flybase.org/reports/${P.key}" target="_blank" rel="noopener">FlyBase ${P.key}</a>` : '';
  const nav = [['c-overview', 'Overview'], ['c-freq', 'Frequency'], ['c-fp', 'Fingerprint'], ['c-info', 'Clusters'], ['c-res', 'Residues'], ['c-3d', '3D structure'], ['c-scatter', 'Scatter'], ['c-net', 'Network'], ['c-pt', 'Partners']];
  const chips = '<div class="chips cl-chips" data-chips></div>';
  const xticks = '<label class="xt">x-ticks <input type="number" class="xticks" min="2" max="40" placeholder="auto"></label>';
  const occ = (await Promise.all(P.occ.map(async (o) => { try { return { ...o, ds: await dataset(sp.dsIds[o.di]) }; } catch (e) { return null; } }))).filter(Boolean);
  if (gone()) return;
  const dataMenu = occ.map((o, i) => { const u = bundleUrl(o.ds, o.name), label = `${sp.dsShort[o.di]}${occ.filter((x) => x.di === o.di).length > 1 ? ' · ' + o.name : ''}`;
    if (o.ds.reg.zip) return `<div class="dm-row"><span class="src" style="--c:${sp.dsColor[o.di]}">${esc(label)}</span><a href="#" data-dl="${i}">Download .zip</a></div>`;   // inside the archive: saved from the read
    return `<div class="dm-row"><span class="src" style="--c:${sp.dsColor[o.di]}">${esc(label)}</span><a href="${LIVIA}clip.html?data=${encodeURIComponent(u)}&gene=${encodeURIComponent(P.gene)}" target="_blank" rel="noopener">Open in LIVIA cLIP ↗</a><a href="${u}" download>Download .zip</a></div>`; }).join('');
  app.innerHTML = `<div class="crumbs"><a href="#/">Atlas</a> / <a href="#/${sp.id}">${esc(sp.reg.label)}</a> / ${esc(P.gene)}</div>
    <div class="phead"><div><h1>${esc(P.gene)}</h1><div class="pname" title="${esc(P.name)}">${esc(short(P.name) || P.id)}</div>
      <div class="ids">${fbLink}${uniprotLink(P.acc)}${fbLink ? '' : `<span>${esc(P.id)}</span>`}${P.clen ? `<span>${fmtInt(P.clen)} aa</span>` : ''}</div>
      <div class="srcs scope" id="scope">In ${srcBadges(sp, P.src)}</div>
      <div class="flags" id="flags">${flags.join('')}</div>
      <div class="actions"><details class="dmenu"><summary class="btn">Data &amp; LIVIA cLIP ▾</summary><div class="dm-pop">${dataMenu}</div></details></div></div>
      <div class="kpis"><div class="kpi"><b id="kp-all">${fmtInt(P.partners)}</b><span>partners predicted</span></div><div class="kpi f10"><b id="kp-10">${fmtInt(P.pos10)}</b><span>past 10% FPR</span></div>
        <div class="kpi f5"><b id="kp-5">${fmtInt(P.pos5)}</b><span>past 5% FPR</span></div><div class="kpi f1"><b id="kp-1">${fmtInt(P.pos1)}</b><span>past 1% FPR</span></div></div></div>
    <div class="setbar" id="setbar" hidden></div>
    <nav class="subnav" aria-label="Sections">${nav.map(([t, l]) => `<button data-t="${t}">${l}</button>`).join('')}</nav>
    <div class="card" id="c-overview"><div class="card-head"><h2>Partners by score</h2><span class="muted" id="sc-sub"></span></div>
      <div class="overview"><div><div class="scat" id="scat"></div>
          <div class="legend"><span><i style="background:#A7B2BF;border-radius:50%"></i>dot size: iLIS average over the models</span><span>color: cluster (gray: not clustered)</span></div></div>
        <div><h3>Top partners <span class="muted">by iLIS of the best model</span></h3>
          <div class="tl-head"><span></span><span>Partner</span><span>Cluster</span><span>iLIS<br>best</span><span>iLIS<br>avg</span><span>ipTM<br>best</span><span>ipTM<br>avg</span></div>
          <ol class="toplist" id="toplist"></ol>
          <div class="legend tl-key"><span>FPR band, each value by its own benchmarked cutoff</span><span><i style="background:#6D4FD1"></i>1%</span><span><i style="background:#16956A"></i>5%</span><span><i style="background:#C78B00"></i>10%</span><span><i style="background:#A7B2BF"></i>below</span></div></div></div></div>
    <div class="card" id="c-clip"><div class="card-head"><div><h2 id="clip-title">${esc(P.gene)} — interactome</h2><div class="muted" id="clip-sub">Loading the predictions…</div></div>
        <div class="clip-ctl"><span class="muted">iLIS cutoff</span><div class="seg" id="cut-seg">${[10, 5, 1].map((f) => `<button data-f="${f}" class="${f === 10 ? 'on' : ''}">${f}% FPR</button>`).join('')}</div></div></div>
      <div class="stat4"><div><b id="s-partners">–</b><span>partners</span></div><div><b id="s-preds">–</b><span>models</span></div><div><b id="s-k">–</b><span>clusters</span></div><div><b id="s-len">–</b><span>query length</span></div></div>
      <p class="explain">LIVIA cLIP, run in your browser on the predictions shown: the residues of ${esc(P.gene)} that a prediction past the cutoff contacts
        (cLIR: PAE ≤ 12 Å and Cβ ≤ 8 Å) form its interaction fingerprint. Fingerprints are compared by cosine distance and joined by average linkage, and the number
        of clusters is chosen by silhouette. Clusters are numbered by size: Cluster 1 is the largest.</p><p class="note" id="clip-aside" hidden></p></div>
    <div class="card" id="c-freq"><div class="card-head"><h2>Contact residue frequency</h2><div class="card-tools"><span class="muted">predictions contacting each residue, colored by their most frequent cluster</span>${xticks}</div></div>
      ${chips}<div class="plot" id="freq-wrap"></div><div class="domlegend" id="freq-domains"></div><div class="hot" id="hot"></div></div>
    <div class="card" id="c-fp"><div class="card-head"><h2>Clustered interaction fingerprint</h2><div class="card-tools"><span class="muted">one row per prediction, in dendrogram order · hover for the partner, click to open the pair</span>${xticks}</div></div>
      ${chips}<div class="plot" id="fp-wrap"></div>
      <div class="legend"><span><i style="background:#08306B"></i>contact residue (cLIR)</span><span><i style="background:#F7FBFF;box-shadow:inset 0 0 0 1px #C9D6E3"></i>no contact</span><span>left: dendrogram and cluster of each prediction</span></div></div>
    <div class="card" id="c-info"><div class="card-head"><h2>Cluster info</h2><span class="muted">Cluster n (proteins / predictions) · largest first</span></div><div class="clinfo" id="cluster-info"></div></div>
    <div class="card" id="c-res"><div class="card-head"><h2>Interaction Residues</h2>
        <div class="controls" style="margin:0"><select id="res-partner" aria-label="Partner" style="max-width:300px"></select><input type="search" id="res-find" placeholder="Find partner" style="width:130px"><select id="res-rank" aria-label="Model" style="max-width:280px"></select></div></div>
      <div class="legend" style="margin:2px 0 12px"><span><i style="background:#E0E0E0"></i>not in the interface</span><span><i style="background:#80CBC4"></i><i style="background:#FFAB91;margin-left:-2px"></i>interface (LIR: PAE ≤ 12 Å)</span><span><i style="background:#00897B"></i><i style="background:#E64A19;margin-left:-2px"></i>contact (cLIR: also Cβ ≤ 8 Å)</span></div>
      <div id="res-body"></div></div>
    <div class="card" id="c-3d"><div class="card-head"><h2>3D structure</h2><span class="muted" id="struct-badge"></span></div>
      <p class="muted" style="margin:2px 0 6px">${esc(P.gene)} as predicted alone in the AlphaFold Database, residues colored by the cluster that consensus-contacts them · click clusters to isolate.
        <b>C<i>n</i> (N)</b>: N = predictions (AlphaFold ranks) in that cluster.</p>
      <div class="controls"><span>Highlight residues contacted by ≥</span><select id="commonality">${[0.5, 0.6, 0.7, 0.8, 0.9].map((v) => `<option value="${v}">${Math.round(v * 100)}%</option>`).join('')}</select>
        <span>of a cluster's members · color only clusters with ≥</span><input type="number" id="min-cluster" min="1" value="5" style="width:60px"><span>predictions</span></div>
      ${chips}
      <div class="controls"><span>Color by</span><div class="seg" id="cmode"><button data-m="cluster" class="on">cluster</button><button data-m="plddt">pLDDT</button><button data-m="am" id="cm-am" hidden>AlphaMissense</button></div>
        <span id="am-cut-wrap" hidden>AM pathogenicity average ≥ <select id="am-cutoff">${[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9].map((v) => `<option value="${v}">${Math.round(v * 100)}%</option>`).join('')}</select></span></div>
      <div class="viewer3d"><iframe id="viewer3d-frame" title="3D structure viewer"></iframe><div class="v3d-msg" id="v3d-msg">Loading the AlphaFold DB model…</div></div>
      <div class="legend" id="legend-3d"></div></div>
    <div class="card" id="c-scatter"><div class="card-head"><h2>Interaction scatter plot</h2>
        <div class="controls" style="margin:0"><label>Y <select id="sc-y"></select></label><label>X <select id="sc-x"></select></label>
          <label>show <select id="sc-pts"><option value="all">every prediction</option><option value="rank1">rank-1 per pair</option></select></label>
          <input type="search" id="sc-find" list="sc-list" placeholder="Find partner" style="width:140px"><datalist id="sc-list"></datalist></div></div>
      <div class="muted" id="sc-rho" style="margin:2px 0 8px"></div><div class="plot" id="scat2"><canvas id="scatter-canvas"></canvas></div><div class="legend" id="sc-legend"></div></div>
    <div class="card" id="c-net"><div class="card-head"><h2>Network</h2>
      <div class="controls" style="margin:0"><span>Partners</span><select id="net-n"><option>30</option><option>60</option><option selected>100</option><option>200</option></select>
        <span>Cutoff</span><select id="net-cut"><option value="10">10% FPR</option><option value="5">5% FPR</option><option value="1">1% FPR</option></select></div></div>
      <p class="muted" style="margin:2px 0 12px">${esc(P.gene)} at the center; partners sit closer the higher their iLIS and are filled with their cluster color. Edges between partners join partners
        predicted to bind each other, in any screen. Drag to move, scroll to zoom, click to open.</p>
      <div class="net" id="net"><div class="loading" style="padding:20px">Loading the network…</div></div>
      <div class="netkey"><div><span>Edge color · best iLIS over the models</span><div class="grad" id="net-grad"></div><div class="gl" id="net-gl"></div></div>
        <div><span>Edge width · average iLIS</span><svg id="net-w" width="260" height="30" aria-hidden="true"></svg></div></div>
      <div class="legend" id="net-legend"></div><div id="net-x"></div></div>
    <div class="card" id="c-pt"><div class="card-head"><h2>Partners <span class="muted" id="pt-note"></span></h2>
      <div class="controls" style="margin:0"><select id="pt-src"><option value="0">all screens</option>${sp.dsIds.map((_, di) => (P.src & (1 << di) ? `<option value="${1 << di}">${esc(sp.dsShort[di])}</option>` : '')).join('')}</select>
        <select id="pt-band"><option value="10">past 10% FPR</option><option value="5">past 5% FPR</option><option value="1">past 1% FPR</option><option value="0">all predicted</option></select>
        <input type="search" id="pt-filter" placeholder="Filter partners" style="width:180px"></div></div>
      <div class="tbl-wrap"><table class="pt" id="pt"></table></div><div class="pager" id="pager"></div></div>`;
  app.querySelectorAll('.subnav button').forEach((b) => b.onclick = () => { const t = document.getElementById(b.dataset.t); if (t) window.scrollTo({ top: t.getBoundingClientRect().top + window.scrollY - 112, behavior: 'smooth' }); });
  { const dm = $('.dmenu'); document.addEventListener('click', (e) => { if (dm && dm.open && !dm.contains(e.target)) dm.open = false; }); }
  app.querySelectorAll('[data-dl]').forEach((a) => a.onclick = async (e) => {   // a bundle inside a screen archive: read it, save it
    e.preventDefault(); const o = occ[+a.dataset.dl];
    try { const raw = await bundleRaw(o.ds, o.name), u = URL.createObjectURL(new Blob([raw.bytes], { type: 'application/zip' }));
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
  {
    const n = new Map(); for (const p of B0.preds) { const d = sp.dsIds[p.di]; n.set(d, (n.get(d) || 0) + 1); for (const t of p.tags || []) n.set(t, (n.get(t) || 0) + 1); }
    const screens = occ.length > 1 ? occ.map((o) => ({ id: o.ds.id, short: o.ds.reg.short, color: o.ds.reg.color, title: o.ds.reg.title })) : [];
    const sets = TS0 ? TS0.list.filter((x) => n.get(x.id)).sort((x, y) => (x.type === 'screen' ? 0 : 1) - (y.type === 'screen' ? 0 : 1)) : [];
    const chip = (x) => `<a class="src scope-chip${SET && SET.id === x.id ? ' on' : ''}" style="--c:${x.color}" href="#/${sp.id}/${P.key}?set=${encodeURIComponent(x.id)}"
      title="${esc(x.title)}: ${fmtInt(n.get(x.id) || 0)} models">${esc(x.short)} <span class="n">${fmtInt(n.get(x.id) || 0)}</span></a>`;
    if (screens.length || sets.length) $('#scope').innerHTML = `<span class="lbl">In</span><a class="src scope-chip all${SET ? '' : ' on'}" href="#/${sp.id}/${P.key}"
      title="every screen and set">All <span class="n">${fmtInt(B0.preds.length)}</span></a>${screens.map(chip).join('')}${screens.length && sets.length ? '<span class="sep"></span>' : ''}${sets.map(chip).join('')}`;
  }
  if (SET) {
    const who = SET.source && SET.source.citation ? `${SET.source.citation.split(' ')[0]} et al. ${(SET.source.citation.match(/\((\d{4})\)/) || [])[1] || ''}` : '';
    const what = SET.type === 'dataset' ? `the ${SET.title.charAt(0).toLowerCase() + SET.title.slice(1)}` : SET.type === 'screen'
      ? `the ${SET.title.charAt(0).toLowerCase() + SET.title.slice(1)} of ${TS0.ds.reg.short}` : `the ${SET.title} category of ${TS0.ds.reg.short}`;
    const about = SET.type === 'dataset' ? `#/datasets/${SET.id}` : `#/datasets/${TS0.ds.id}/${SET.id}`;
    const bar = $('#setbar'); bar.hidden = false;
    bar.innerHTML = `<span class="src" style="--c:${SET.color}">${esc(SET.short)}</span><span>Only ${esc(what)}${who ? ` (<a href="${esc(SET.source.url)}" target="_blank" rel="noopener">${esc(who)}</a>)` : ''}:
      ${fmtInt(B.preds.length)} of ${fmtInt(B0.preds.length)} models.</span><span><a href="#/${sp.id}/${P.key}">Show every prediction</a> · <a href="${about}">about this ${SET.type === 'dataset' ? 'screen' : 'set'}</a></span>`;
    const others = B.partners.filter((x) => x.id !== P.key), cnt = (c) => others.filter((x) => x.best >= c).length;
    $('#kp-all').textContent = fmtInt(others.length); $('#kp-10').textContent = fmtInt(cnt(CUT[10])); $('#kp-5').textContent = fmtInt(cnt(CUT[5])); $('#kp-1').textContent = fmtInt(cnt(CUT[1]));
    document.title = `${P.gene} · ${SET.short} · LIVIA cLIP Atlas`;
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
    $('#s-len').textContent = fmtInt(CQ.qLen || P.clen || 0);
    const a = $('#clip-aside');
    if (B.cons.size) {   // a gene-keyed screen: its other constructs by kind
      const KIND = { phosphosite: ['phosphosite window', 'phosphosite windows'], fragment: ['fragment', 'fragments'], mutant: ['point mutant', 'point mutants'],
        variant: ['variant', 'variants'], peptide: ['peptide', 'peptides'], isoform: ['isoform', 'isoforms'], gene: ['construct', 'constructs'] };
      const pickable = new Set((B.choices || []).map((c) => c.id)), say = (list) => { const byKind = new Map();
        for (const x of list) { const k = x.con ? x.con.kind : 'construct'; if (!byKind.has(k)) byKind.set(k, []); byKind.get(k).push(x); }
        return [...byKind].map(([k, xs]) => { const [one, many] = KIND[k] || ['construct', 'constructs'], ex = xs.slice(0, 3).map((x) => (x.con ? x.con.label : '')).filter(Boolean);
          return `${fmtInt(xs.length)} ${xs.length === 1 ? one : many}${ex.length ? ` (${ex.join(', ')}${xs.length > 3 ? ', …' : ''})` : ''}`; }).join(', '); };
      const own = [...CQ.aside.entries()].filter(([k]) => pickable.has(k)).map(([, x]) => x), rest = [...CQ.aside.entries()].filter(([k]) => !pickable.has(k)).map(([, x]) => x);
      const uses = CQ.qName ? `${(B.cons.get(CQ.qName) || {}).label || CQ.qName}, in its own numbering` : `the ${fmtInt(CQ.qLen)} aa reference, with every construct placed or mapped on it`;
      a.textContent = CQ.qName ? `Clustering uses ${uses}. The ${P.gene} reference and its other constructs are in the menu above.` : `Clustering uses ${uses}.`
        + (own.length ? ` ${say(own)} ${own.length === 1 ? 'shares' : 'share'} too little sequence with ${CQ.qName ? 'it' : 'the reference'} to be drawn on it; choose one in the menu above to cluster it on its own.` : '')
        + (rest.length ? ` ${say(rest)} ${rest.length === 1 ? 'is' : 'are'} listed with ${rest.length === 1 ? 'its' : 'their'} partners and on the pair pages, not clustered.` : '');
      a.hidden = !(own.length || rest.length || CQ.qName);
    } else if (CQ.aside.size) { a.hidden = false;
      a.textContent = [...CQ.aside.values()].map((x) => `${sp.dsShort[x.di]} predicted ${P.gene} as a ${fmtInt(x.len)} aa construct`).join('; ') + `, so those models are listed but left out of the clustering, which uses the ${fmtInt(CQ.qLen)} aa construct.`;
    } else a.hidden = true;
  };
  constructNote();
  if (B.choices && B.choices.length > 1) {   // which construct to cluster: the reference, or another construct with models of its own
    $('.clip-ctl').insertAdjacentHTML('afterbegin', `<select id="clip-construct" class="setpick" aria-label="Construct to cluster">${B.choices.map((c) =>
      `<option value="${esc(c.id)}">${esc(c.label)} · ${fmtInt(c.n)} models</option>`).join('')}</select>`);
    $('#clip-construct').value = CQ.choice;
    $('#clip-construct').onchange = (e) => { CQ = B.clipFor(e.target.value); setQSeq(); constructNote(); mapStruct(); cluster(); };
  }
  app.querySelectorAll('.xticks').forEach((i) => { i.oninput = () => { app.querySelectorAll('.xticks').forEach((o) => { if (o !== i) o.value = i.value; }); drawFreq(); drawHeatmap(); }; });

  /* cLIP ─ clustering + everything drawn from it */
  async function cluster() {
    $('#clip-sub').textContent = 'Clustering…';
    let m = null, err = null;
    try { m = await runClip(CQ.rows, B.qLabel, CUT[cut]); } catch (e) { err = e; }
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
    $('#s-partners').textContent = fmtInt(partnerCluster.size); $('#s-preds').textContent = fmtInt(n); $('#s-k').textContent = k || '–';
    $('#clip-title').innerHTML = `${esc(P.gene)} — interactome <span class="muted">(${fmtInt(partnerCluster.size)} partners, ${fmtInt(n)} predictions, ${k} clusters)</span>`;
    if (M) $('#clip-sub').textContent = clustered() ? `query length ${fmtInt(M.plen)} aa · cosine + average linkage · silhouette-optimal k · iLIS ≥ ${CUT[cut]} (${cut}% FPR)`
      : `${n ? 'Only one prediction' : 'No predictions'} past the ${cut}% FPR cutoff, so there is nothing to cluster.`;
    const want = !clustered() && V.mode === 'cluster' ? 'plddt' : clustered() && V.auto && S.mapOK ? 'cluster' : null;   // no clusters: show pLDDT until there are
    if (want) { V.auto = want === 'plddt'; V.mode = want; app.querySelectorAll('#cmode button').forEach((b) => b.classList.toggle('on', b.dataset.m === want)); }
    renderChips(); drawFreq(); drawHeatmap(); renderClusterInfo(); recolor3D(); drawScatter(); drawOverview(); drawTable(); fillPartners(); if (NET) NET.recolor();
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
    paintChips(); drawFreq(); drawHeatmap(); recolor3D();
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
    const yTicks = (vals, fmt, yOf) => { g.fillStyle = '#6B7A8D'; g.font = '10.5px "IBM Plex Mono", ui-monospace, monospace'; g.textAlign = 'right'; g.textBaseline = 'middle'; for (const v of vals) g.fillText(fmt(v), AXL - 6, yOf(v)); };
    const frame = (y0, y1) => { g.fillStyle = '#F6F8FB'; g.fillRect(AXL, y0, W - AXL - AXR, y1 - y0); g.strokeStyle = '#D5DDE6'; g.lineWidth = 1; g.beginPath(); g.moveTo(AXL + 0.5, y0); g.lineTo(AXL + 0.5, y1 + 0.5); g.lineTo(W - AXR, y1 + 0.5); g.stroke(); };
    for (const d of doms) { const x0 = xOf(Math.max(1, d.s)), x1 = Math.max(x0 + 2, xOf(Math.min(L, d.e) + 1)), y = 4 + d.lane * DRH;   // domains D1, D2 … (names in the legend), as LIVIA cLIP
      g.fillStyle = '#E3E9F1'; g.fillRect(x0, y, x1 - x0, 12); g.strokeStyle = '#9FB0C4'; g.lineWidth = 0.6; g.strokeRect(x0 + 0.3, y + 0.3, x1 - x0 - 0.6, 11.4);
      g.fillStyle = '#34445A'; g.font = '10px "IBM Plex Sans", system-ui, sans-serif'; g.textAlign = 'center'; g.textBaseline = 'middle'; const lb = 'D' + d.idx; if (g.measureText(lb).width < x1 - x0 - 2) g.fillText(lb, (x0 + x1) / 2, y + 6.5); }
    frame(padT, fBase);   // # partners (predictions contacting the residue); bar color = the residue's most frequent cluster
    g.strokeStyle = '#E3E8EE'; g.lineWidth = 1; for (const f of [0.5, 1]) { const y = Math.round(fBase - f * fH) + 0.5; g.beginPath(); g.moveTo(AXL, y); g.lineTo(W - AXR, y); g.stroke(); }
    for (let r = 1; r <= L; r++) { if (!tot[r]) continue; let dom = 1, b = -1; const cc = byC[r]; for (const c in cc) if (cc[c] > b) { b = cc[c]; dom = +c; }
      const h = tot[r] / max * fH; g.fillStyle = clusterColor(dom, k); g.fillRect(xOf(r), fBase - h, Math.max(1, bw), h); }
    yTicks([...new Set([0, Math.round(max / 2), max])], String, (v) => fBase - v / max * fH); yLabel('# partners', padT, fBase, '#34445A');
    const line = (vals, yOf, col) => { let prev = null; g.lineWidth = 1.5;
      for (let r = 1; r <= L; r++) { const sr = toStruct(r), v = sr ? vals(sr) : null; if (v == null) { prev = null; continue; } const x = xc(r), y = yOf(v);
        if (prev) { const mx = (prev.x + x) / 2, my = (prev.y + y) / 2; g.strokeStyle = col(prev.v); g.beginPath(); g.moveTo(prev.x, prev.y); g.lineTo(mx, my); g.stroke(); g.strokeStyle = col(v); g.beginPath(); g.moveTo(mx, my); g.lineTo(x, y); g.stroke(); }
        prev = { x, y, v }; } };
    if (pl) { frame(pTop, pBase); const yOf = (v) => pBase - v / 100 * pH; line((sr) => S.plddt.get('A:' + sr), yOf, plddtCol); yTicks([0, 50, 100], String, yOf); yLabel('pLDDT', pTop, pBase, '#2F6FA8'); }
    if (am) { frame(aTop, aBase); const yOf = (v) => aBase - v * aH; line((sr) => S.am[sr], yOf, amCol); yTicks([0, 0.5, 1], (v) => v.toFixed(1), yOf); yLabel('AM path.', aTop, aBase, '#A33'); }
    drawTicks(g, resTicks(L, W - AXL - AXR, xtWant()), aBase + 2, xc, W);
    cv.onmousemove = (e) => { const b = cv.getBoundingClientRect(), r = Math.floor((e.clientX - b.left - AXL) / bw) + 1; if (r < 1 || r > L) return hideTip();
      const cc = byC[r] || {}, parts = Object.keys(cc).map(Number).sort((a, z) => cc[z] - cc[a]).map((c) => `<span style="color:${clusterColor(c, k)}">●</span> ${clusterLabel(c, true)} ${cc[c]}`).join(' · ');
      const sr = toStruct(r), plv = pl && sr ? S.plddt.get('A:' + sr) : null, amv = am && sr ? S.am[sr] : null, dom = doms.filter((d) => r >= d.s && r <= d.e).map((d) => d.name).join(', ');
      showTip(`<b>${qSeq[r - 1] || ''}${r}</b> · ${tot[r]} prediction${tot[r] === 1 ? '' : 's'}${parts ? '<br>' + parts : ''}${dom ? `<br>${esc(dom)}` : ''}${plv != null ? `<br>pLDDT ${plv.toFixed(0)}` : ''}${amv != null ? `${plv != null ? ' · ' : '<br>'}AM ${amv.toFixed(2)}` : ''}`, e.clientX, e.clientY); };
    cv.onmouseleave = hideTip;
    $('#freq-domains').innerHTML = doms.length ? `<b>Domains</b> ${doms.map((d) => `<span><b>D${d.idx}</b> ${esc(d.name)} <span class="muted">(${d.start}–${d.end})</span></span>`).join('')}${P.acc ? uniprotLink(P.acc, `UniProt ${esc(P.acc)} ↗`) : ''}` : '';
    const hot = range(L).filter((r) => tot[r]).sort((a, b) => tot[b] - tot[a]).slice(0, 12);
    $('#hot').innerHTML = hot.length ? `<span class="hot-lbl">Most contacted</span>${hot.map((r) => `<span title="${tot[r]} predictions">${qSeq[r - 1] || ''}${r} · ${tot[r]}</span>`).join('')}` : '';
    attachExport('freq', `atlas_${P.gene}_frequency`, drawFreq);
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

  function renderClusterInfo() {
    const box = $('#cluster-info'); if (!box) return;
    if (!clustered()) { box.innerHTML = `<p class="muted" style="margin:0">${M ? 'Nothing to cluster at this cutoff.' : 'Clustering…'}</p>`; return; }
    const mem = {}, cnt = {};
    M.preds.forEach((p, i) => { const c = M.labels[i]; (mem[c] ||= new Set()).add(who(p.partner).key); cnt[c] = (cnt[c] || 0) + 1; });
    const cap = 60;
    box.innerHTML = range(M.k).map((c) => {
      const rows = [...(mem[c] || [])].map((x) => sp.byKey.get(x) || { key: x, gene: x }).sort((a, b) => a.gene.toLowerCase().localeCompare(b.gene.toLowerCase()));
      const show = infoOpen.has(c) ? rows : rows.slice(0, cap);
      return `<div class="cl-row"><i style="background:${clusterColor(c, M.k)}"></i><div><b>${clusterLabel(c)}</b> <span class="muted">(${rows.length} / ${cnt[c]})</span>: ${show.map((r) =>
        `<a href="#/${sp.id}/${P.key}/${r.key}${scopeQ}">${esc(r.gene)}</a>`).join(', ')}${rows.length > cap ? ` <button class="more" data-c="${c}">${infoOpen.has(c) ? 'show fewer' : `+${rows.length - cap} more`}</button>` : ''}</div></div>`; }).join('');
    box.querySelectorAll('.more').forEach((b) => b.onclick = () => { const c = +b.dataset.c; infoOpen.has(c) ? infoOpen.delete(c) : infoOpen.add(c); renderClusterInfo(); });
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
    $('#res-rank').innerHTML = part.preds.map((p, i) => `<option value="${i}">${esc(runLabel(sp, B, p.run, P))} · rank ${p.rank} · iLIS ${fmtNum(p.iLIS, 3)}</option>`).join('');
    $('#res-rank').value = String(part.preds.indexOf(best));
    drawResidues();
  }
  function drawResidues() {
    const id = $('#res-partner').value, part = B.partners.find((p) => p.id === id); if (!part) return;
    const pred = part.preds[+$('#res-rank').value] || part.preds[0];
    ifaceView($('#res-body'), { sp, P, O: PLACE(id), pred, B, canvasId: 'res-canvas' });
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
    frame.src = URL.createObjectURL(new Blob([buildMolstarPage(S.text, 'mmcif', colorComponents(), LIVIA)], { type: 'text/html' }));
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
    domainsOf(P.acc).then((d) => { if (gone()) return; S.domains = d; drawFreq(); });
    const entry = await afdbEntry(P.acc);
    if (gone()) return;
    if (!entry) { S.state = 'none'; msg(`No AlphaFold DB model for ${esc(P.acc)} (the database has none for proteins longer than 2,700 residues).`); $('#struct-badge').textContent = ''; return; }
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
  $('#commonality').onchange = (e) => { V.thr = +e.target.value; recolor3D(); };
  $('#min-cluster').onchange = (e) => { V.minC = Math.max(1, +e.target.value || 1); recolor3D(); };
  $('#am-cutoff').onchange = (e) => { V.amCut = +e.target.value; recolor3D(); };
  { const io3 = new IntersectionObserver((ents) => { if (ents.some((x) => x.isIntersecting)) { io3.disconnect(); show3D(); } }, { rootMargin: '300px' }); io3.observe($('#c-3d')); }

  /* Interaction scatter plot: every prediction (or rank 1 of each run), cluster colors, any two metrics — on a canvas, so it exports */
  const present = Object.keys(METRICS).filter((k) => k === '_rank' || B.preds.some((p) => Number.isFinite(p[k]) && p[k] !== 0));
  $('#sc-y').innerHTML = present.map((k) => `<option value="${k}">${METRICS[k]}</option>`).join('');
  $('#sc-x').innerHTML = present.map((k) => `<option value="${k}">${METRICS[k]}</option>`).join('');
  $('#sc-y').value = 'iLIS'; $('#sc-x').value = present.includes('iLISA') ? 'iLISA' : 'ipTM';
  $('#sc-list').innerHTML = [...new Set(B.partners.map((p) => gname(p.id)))].sort((a, b) => a.localeCompare(b)).map((g) => `<option value="${esc(g)}"></option>`).join('');
  ['#sc-x', '#sc-y', '#sc-pts'].forEach((s) => { $(s).onchange = () => drawScatter(); });
  $('#sc-find').onchange = () => drawScatter();
  $('#sc-find').onkeydown = (e) => { if (e.key === 'Enter') { e.preventDefault(); drawScatter(); } };
  $('#sc-find').oninput = (e) => { if (!e.target.value) drawScatter(); };
  function drawScatter() {
    const cv = $('#scatter-canvas'); if (!cv) return;
    const xK = $('#sc-x').value, yK = $('#sc-y').value, mode = $('#sc-pts').value, k = M ? M.k : 1;
    const rankBy = yK !== '_rank' ? yK : xK !== '_rank' ? xK : 'iLIS';
    let grank = null;
    if (xK === '_rank' || yK === '_rank') { grank = new Map(); B.preds.map((p, i) => [p, i]).sort((a, b) => ((b[0][rankBy] || 0) - (a[0][rankBy] || 0)) || a[1] - b[1]).forEach(([p], i) => grank.set(p, i + 1)); }
    const val = (p, key) => (key === '_rank' ? grank.get(p) : p[key]);
    const pts = [];
    for (const p of B.preds) { const x = val(p, xK), y = val(p, yK); if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
      const c = predCluster.get(p.label + '|' + p.rank); pts.push({ p, x, y, c: c != null && (mode === 'all' || p.rank === 1) ? c : 0 }); }
    pts.sort((a, b) => (a.c ? 1 : 0) - (b.c ? 1 : 0));
    const W = cv.parentElement.clientWidth, H = 450, m = { l: 62, r: 18, t: 24, b: 46 };
    const UNIT = new Set(['iLIS', 'ipTM', 'pTM', 'LIS', 'cLIS', 'ipSAE', 'actifpTM']);   // scores bounded by 1 keep their whole 0–1 range
    const dom = (key, vs) => (key === '_rank' ? [0, Math.max(1, vs.length)] : UNIT.has(key) ? [0, 1] : key === 'qPl' || key === 'pPl' ? [0, 100]
      : [Math.min(0, d3.min(vs) ?? 0), Math.max(1e-6, (d3.max(vs) ?? 1) * 1.04)]);
    const xs = d3.scaleLinear(dom(xK, pts.map((q) => q.x)), [m.l, W - m.r]).nice(), ys = d3.scaleLinear(dom(yK, pts.map((q) => q.y)), [H - m.b, m.t]).nice();
    const g = canvasCtx(cv, W, H);
    const title = (key) => (key === '_rank' ? `global rank (by ${METRICS[rankBy]})` : METRICS[key]);
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
    for (const q of pts) { g.beginPath(); g.arc(xs(q.x), ys(q.y), q.c ? 4.3 : 2.5, 0, 2 * Math.PI); g.fillStyle = q.c ? clusterColor(q.c, k) : '#CDD3DB'; g.fill(); if (q.c) { g.lineWidth = 0.7; g.strokeStyle = '#fff'; g.stroke(); } }
    const find = $('#sc-find').value.trim().toLowerCase();
    if (find) {
      const gl = (p) => gname(p.id).toLowerCase(), hit = B.partners.find((p) => gl(p) === find) || B.partners.find((p) => gl(p).startsWith(find)) || B.partners.find((p) => gl(p).includes(find));
      $('#sc-find').classList.toggle('nf', !hit);
      if (hit) for (const q of pts.filter((q) => q.p.partner === hit.id)) { const x = xs(q.x), y = ys(q.y);
        g.beginPath(); g.arc(x, y, 8.5, 0, 2 * Math.PI); g.lineWidth = 2.5; g.strokeStyle = '#17263A'; g.stroke();
        g.font = '600 12px "IBM Plex Sans", system-ui, sans-serif'; g.textAlign = 'center'; g.textBaseline = 'alphabetic'; g.lineWidth = 3; g.strokeStyle = 'rgba(255,255,255,0.92)'; const t = `${gname(hit.id)} (R${q.p.rank})`;
        g.strokeText(t, x, y - 13); g.fillStyle = '#17263A'; g.fillText(t, x, y - 13); }
    } else $('#sc-find').classList.remove('nf');
    const qt = d3.quadtree(pts, (q) => xs(q.x), (q) => ys(q.y));
    const near = (e) => { const b = cv.getBoundingClientRect(); return qt.find(e.clientX - b.left, e.clientY - b.top, 10); };
    const fmtM = (key, v) => (Number.isFinite(v) ? (['iLIA', 'iLISA', 'LIA', 'cLIA', 'qPl', 'pPl'].includes(key) ? v.toFixed(1) : v.toFixed(3)) : '–');
    cv.onmousemove = (e) => { const q = near(e); if (!q) { cv.style.cursor = ''; return hideTip(); } cv.style.cursor = 'pointer'; const p = q.p;
      const keys = [...new Set(['iLIS', 'iLISA', 'iLIA', 'ipTM', yK, xK])].filter((x) => x !== '_rank' && present.includes(x));
      showTip(`<b>${esc(gname(p.partner))}</b> · ${esc(runLabel(sp, B, p.run, P))} · rank ${p.rank}${q.c ? ` · <span style="color:${clusterColor(q.c, k)}">●</span> ${clusterLabel(q.c)}` : ''}<br>${keys.map((x) => `${METRICS[x]} ${fmtM(x, p[x])}`).join(' · ')}`, e.clientX, e.clientY); };
    cv.onmouseleave = hideTip;
    cv.onclick = (e) => { const q = near(e); if (!q) return; hideTip(); location.hash = `#/${sp.id}/${P.key}/${q.p.partner}${scopeQ}`; };
    const r1 = pearson(pts.map((q) => q.x), pts.map((q) => q.y)), rho = spearman(pts.map((q) => q.x), pts.map((q) => q.y));
    $('#sc-rho').textContent = `${Number.isFinite(r1) ? `Pearson r = ${r1.toFixed(3)} · ` : ''}${Number.isFinite(rho) ? `Spearman ρ = ${rho.toFixed(3)} ` : ''}(n = ${fmtInt(pts.length)})`;
    const cnt = {}; let other = 0; for (const q of pts) q.c ? (cnt[q.c] = (cnt[q.c] || 0) + 1) : other++;
    $('#sc-legend').innerHTML = Object.keys(cnt).map(Number).sort((a, b) => a - b).map((c) => `<span><i style="background:${clusterColor(c, k)};border-radius:50%"></i>cluster ${c} (${cnt[c]})</span>`).join('')
      + `<span><i style="background:#CDD3DB;border-radius:50%"></i>other predictions (${fmtInt(other)})</span>`;
    attachExport('scatter-canvas', `atlas_${P.gene}_scatter`, drawScatter);
  }

  /* overview: iLIS × ipTM for every partner (best model of any screen), dot size = iLIS average, color = cluster */
  function drawOverview() {
    const host = $('#scat'); if (!host) return;
    const list = B.partners.map((p) => ({ ...p, gene: gname(p.id), c: partnerCluster.get(p.id) || 0 }));
    $('#sc-sub').textContent = `${fmtInt(list.length)} partners · best model of each pair`;
    const W = host.clientWidth, H = 400, m = { l: 48, r: 14, t: 26, b: 42 };
    const x = d3.scaleLinear([0, 1], [m.l, W - m.r]), y = d3.scaleLinear([0, 1], [H - m.b, m.t]), rad = (v) => 2.2 + 8 * Math.min(1, v / 0.8);   // both scores span 0–1: room above the top partners for their labels
    host.innerHTML = '';
    const svg = d3.select(host).append('svg').attr('width', W).attr('height', H);
    svg.append('g').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(5));
    svg.append('g').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5));
    svg.selectAll('.tick text').attr('font-family', 'IBM Plex Mono').attr('fill', '#6B7A8D'); svg.selectAll('.domain, .tick line').attr('stroke', '#D5DDE6');
    svg.append('text').attr('x', (m.l + W - m.r) / 2).attr('y', H - 6).attr('text-anchor', 'middle').attr('font-size', 12).attr('fill', '#34445A').text('ipTM (best model)');
    svg.append('text').attr('transform', `translate(13,${(m.t + H - m.b) / 2}) rotate(-90)`).attr('text-anchor', 'middle').attr('font-size', 12).attr('fill', '#34445A').text('iLIS (best model)');
    [10, 5, 1].forEach((f, j) => {   // benchmarked cutoffs: iLIS across, ipTM down
      const v = FPR.iLIS[j], u = FPR.ipTM[j];
      svg.append('line').attr('x1', m.l).attr('x2', W - m.r).attr('y1', y(v)).attr('y2', y(v)).attr('stroke', BAND[f]).attr('stroke-dasharray', '4 4').attr('opacity', 0.7);
      svg.append('text').attr('x', W - m.r - 2).attr('y', y(v) - 4).attr('text-anchor', 'end').attr('font-size', 10.5).attr('font-family', 'IBM Plex Mono').attr('fill', BAND[f]).text(`${f}% FPR`);
      svg.append('line').attr('x1', x(u)).attr('x2', x(u)).attr('y1', m.t).attr('y2', H - m.b).attr('stroke', BAND[f]).attr('stroke-dasharray', '2 5').attr('opacity', 0.55);
      svg.append('text').attr('x', x(u)).attr('y', m.t - 8).attr('text-anchor', 'middle').attr('font-size', 10).attr('font-family', 'IBM Plex Mono').attr('fill', BAND[f]).text(`${f}%`);
    });
    const k = M ? M.k : 1, color = (p) => (p.c ? clusterColor(p.c, k) : '#B7C2CE');
    const pts = [...list].sort((a, b) => (a.c ? 1 : 0) - (b.c ? 1 : 0) || a.best - b.best);
    svg.append('g').selectAll('circle').data(pts).join('circle').attr('cx', (p) => x(p.iptmBest)).attr('cy', (p) => y(p.best)).attr('r', (p) => rad(p.avg))
      .attr('fill', color).attr('fill-opacity', (p) => (p.c ? 0.8 : 0.35)).attr('stroke', '#fff').attr('stroke-width', 0.8).style('cursor', 'pointer')
      .on('mousemove', (ev, p) => showTip(`<b>${esc(p.gene)}</b> · iLIS ${p.best.toFixed(3)} (avg ${p.avg.toFixed(3)}) · iLISA ${p.ilisaBest.toFixed(1)} · ipTM ${p.iptmBest.toFixed(2)}${p.c ? ` · ${clusterLabel(p.c)}` : ''}<br>${SETS ? setBadges(p.sets) : srcBadges(sp, p.src)}`, ev.clientX, ev.clientY))
      .on('mouseleave', hideTip).on('click', (ev, p) => { hideTip(); location.hash = `#/${sp.id}/${P.key}/${p.id}${scopeQ}`; });
    const placed = [...[CUT[10], CUT[5], CUT[1]].map((v) => [W - m.r - 60, y(v) - 16, W - m.r, y(v)]),   // cutoff labels are taken
      ...FPR.ipTM.map((u) => [x(u) - 16, 0, x(u) + 16, m.t])], labels = [];
    for (const p of [...list].sort((a, b) => b.best - a.best).slice(0, 10)) {   // label the top partners where a label fits, never on another label
      const cx = x(p.iptmBest), cy = y(p.best), rr = rad(p.avg), w = p.gene.length * 7 + 4, h = 13;
      for (const [tx, ty, anchor] of [[cx + rr + 4, cy + 4, 'start'], [cx - rr - 4, cy + 4, 'end'], [cx, cy - rr - 5, 'middle'], [cx, cy + rr + 13, 'middle']]) {
        const bx = anchor === 'start' ? tx : anchor === 'end' ? tx - w : tx - w / 2, box = [bx, ty - h + 2, bx + w, ty + 3];
        if (box[0] < m.l || box[2] > W - m.r || box[1] < m.t - 2) continue;
        if (placed.some((q) => box[0] < q[2] && box[2] > q[0] && box[1] < q[3] && box[3] > q[1])) continue;
        placed.push(box); labels.push({ p, tx, ty, anchor }); break;
      }
    }
    svg.append('g').selectAll('text').data(labels).join('text').attr('x', (d) => d.tx).attr('y', (d) => d.ty).attr('text-anchor', (d) => d.anchor)
      .attr('font-size', 11.5).attr('font-weight', 600).attr('fill', '#17263A').attr('paint-order', 'stroke').attr('stroke', 'rgba(255,255,255,0.92)').attr('stroke-width', 3).text((d) => d.p.gene);
    const tip = (cuts, v, what) => `${what}: ${bandLabel[bandIn(cuts, v)]} (cutoffs ${cuts.join(' / ')})`;
    $('#toplist').innerHTML = [...list].sort((a, b) => b.best - a.best).slice(0, 12).map((p) => `<li><a href="#/${sp.id}/${P.key}/${p.id}${scopeQ}">${esc(p.gene)}</a>
      <span class="tl-c" title="${p.c ? clusterLabel(p.c) : 'not clustered at this cutoff'}"><span class="mdot" style="background:${p.c ? clusterColor(p.c, k) : '#DDE3EA'}"></span>${p.c ? clusterLabel(p.c, true) : '—'}</span>
      <span class="num" style="color:${bandCol(FPR.iLIS, p.best)}" title="${tip(FPR.iLIS, p.best, 'iLIS best')}">${p.best.toFixed(3)}</span>
      <span class="num" style="color:${bandCol(FPR_AVG.iLIS, p.avg)}" title="${tip(FPR_AVG.iLIS, p.avg, 'iLIS average')}">${p.avg.toFixed(3)}</span>
      <span class="num" style="color:${bandCol(FPR.ipTM, p.iptmBest)}" title="${tip(FPR.ipTM, p.iptmBest, 'ipTM best')}">${p.iptmBest.toFixed(2)}</span>
      <span class="num" style="color:${bandCol(FPR_AVG.ipTM, p.iptmAvg)}" title="${tip(FPR_AVG.ipTM, p.iptmAvg, 'ipTM average')}">${p.iptmAvg.toFixed(2)}</span></li>`).join('');
    svgExport(host, `atlas_${P.gene}_partners`, () => $('svg', host));
  }
  $('#cut-seg').onclick = (e) => { const f = e.target.dataset.f; if (!f) return; cut = +f; [...$('#cut-seg').children].forEach((b) => b.classList.toggle('on', b.dataset.f === f)); cluster(); };

  /* partner table */
  const T = { sort: 'best', asc: false, page: 0, band: P.pos10 ? 10 : 0, filter: '', src: 0 };
  const cols = [['gene', 'Partner'], ['c', 'Cluster'], ['src', 'Source'], ['name', 'Protein'], ['best', 'iLIS best'], ['avg', 'iLIS avg'], ['iptmBest', 'ipTM best'], ['iptmAvg', 'ipTM avg'], ['contacts', 'Contacts'], ['pass', 'Models past']];
  function drawTable() {
    let list = B.partners.map((p) => { const r = sp.byKey.get(p.id); return { ...p, gene: r ? r.gene : p.id, name: r ? r.name : '', c: partnerCluster.get(p.id) || 0, pass: p.preds.filter((x) => x.iLIS >= CUT[10]).length }; });
    if (T.src) list = list.filter((p) => (SETS ? p.sets.includes(T.src) : p.src & T.src));
    if (T.band) list = list.filter((p) => p.best >= CUT[T.band]);
    if (T.filter) { const f = T.filter.toLowerCase(); list = list.filter((p) => p.gene.toLowerCase().includes(f) || (p.name || '').toLowerCase().includes(f) || p.id.toLowerCase().includes(f)); }
    const key = T.sort; list.sort((a, b) => (typeof a[key] === 'string' ? a[key].localeCompare(b[key]) : a[key] - b[key]) * (T.asc ? 1 : -1));
    const per = 40, pages = Math.max(1, Math.ceil(list.length / per)); T.page = Math.min(T.page, pages - 1);
    const view = list.slice(T.page * per, T.page * per + per), k = M ? M.k : 1;
    $('#pt-note').textContent = `${fmtInt(list.length)} shown · ${fmtInt(B.partners.length)} predicted`;
    $('#pt').innerHTML = `<thead><tr>${cols.map(([c, l]) => `<th data-c="${c}" class="${T.sort === c ? 'sorted' + (T.asc ? ' asc' : '') : ''}${['best', 'avg', 'iptmBest', 'iptmAvg', 'contacts', 'pass'].includes(c) ? ' n' : ''}">${l}</th>`).join('')}</tr></thead><tbody>${view.map((p) => {
      const b = bandOf(p.best);
      return `<tr><td class="g"><a href="#/${sp.id}/${P.key}/${p.id}${scopeQ}">${esc(p.gene)}</a></td>
        <td>${p.c ? `<span class="mdot" style="background:${clusterColor(p.c, k)}"></span>${clusterLabel(p.c, true)}` : '<span class="muted">—</span>'}</td>
        <td class="srcc">${SETS ? setBadges(p.sets) : srcBadges(sp, p.src)}</td><td class="nm" title="${esc(p.name)}">${esc(short(p.name))}</td>
        <td class="n v" style="color:${BAND[b]}" title="${bandLabel[b]}">${p.best.toFixed(3)}</td>
        <td class="n v" style="color:${bandCol(FPR_AVG.iLIS, p.avg)}" title="${bandLabel[bandIn(FPR_AVG.iLIS, p.avg)]} (average-model cutoffs)">${p.avg.toFixed(3)}</td>
        <td class="n v" style="color:${bandCol(FPR.ipTM, p.iptmBest)}" title="${bandLabel[bandIn(FPR.ipTM, p.iptmBest)]}">${p.iptmBest.toFixed(2)}</td>
        <td class="n v" style="color:${bandCol(FPR_AVG.ipTM, p.iptmAvg)}" title="${bandLabel[bandIn(FPR_AVG.ipTM, p.iptmAvg)]} (average-model cutoffs)">${p.iptmAvg.toFixed(2)}</td>
        <td class="n">${fmtInt(p.contacts)}</td><td class="n" title="models past the 10% FPR cutoff">${p.pass} / ${p.preds.length}</td></tr>`; }).join('')}</tbody>`;
    $('#pt').querySelectorAll('th').forEach((th) => th.onclick = () => { const c = th.dataset.c; T.asc = T.sort === c ? !T.asc : (c === 'gene' || c === 'name' || c === 'c'); T.sort = c; drawTable(); });
    $('#pager').innerHTML = pages > 1 ? `<button class="btn" id="pp" ${T.page ? '' : 'disabled'}>Previous</button><span>Page ${T.page + 1} of ${pages}</span><button class="btn" id="pn" ${T.page < pages - 1 ? '' : 'disabled'}>Next</button>` : '';
    if (pages > 1) { $('#pp').onclick = () => { T.page--; drawTable(); }; $('#pn').onclick = () => { T.page++; drawTable(); }; }
  }
  $('#pt-band').value = String(T.band);
  $('#pt-band').onchange = (e) => { T.band = +e.target.value; T.page = 0; drawTable(); };
  if (SETS) $('#pt-src').innerHTML = `<option value="">all categories</option>${SETS.map((x) => `<option>${esc(x)}</option>`).join('')}`;   // one screen: filter by its categories
  $('#pt-src').onchange = (e) => { T.src = SETS ? e.target.value : +e.target.value; T.page = 0; drawTable(); };
  $('#pt-filter').oninput = (e) => { T.filter = e.target.value; T.page = 0; drawTable(); };

  /* network: edge color = best iLIS over the models, edge width = average iLIS */
  const netBox = $('#net'), ecol = ECOL();
  { const stops = d3.range(0, 1.0001, 0.1).map((t) => { const v = CUT[10] + t * (0.85 - CUT[10]); return `${ecol(v)} ${(t * 100).toFixed(0)}%`; });
    $('#net-grad').style.background = `linear-gradient(90deg, ${stops.join(', ')})`;
    const pos = (v) => ((v - CUT[10]) / (0.85 - CUT[10]) * 100).toFixed(1) + '%';
    $('#net-gl').innerHTML = [[CUT[10], '0.223'], [CUT[5], '0.339'], [CUT[1], '0.551'], [0.85, '0.85+']].map(([v, l]) => `<span style="left:${pos(v)}">${l}</span>`).join('');
    const w = d3.select('#net-w'); [[0.1, 10], [0.4, 95], [0.7, 180]].forEach(([a, x0]) => { w.append('line').attr('x1', x0).attr('x2', x0 + 44).attr('y1', 10).attr('y2', 10).attr('stroke', '#50637A').attr('stroke-width', EWID(a)).attr('stroke-linecap', 'round');
      w.append('text').attr('x', x0 + 22).attr('y', 27).attr('text-anchor', 'middle').attr('font-size', 10.5).attr('font-family', 'IBM Plex Mono').attr('fill', '#6B7A8D').text(a.toFixed(1)); }); }
  const io = new IntersectionObserver(async (ents) => { if (!ents.some((e) => e.isIntersecting)) return; io.disconnect(); await drawNet(); }, { rootMargin: '200px' });
  io.observe(netBox);
  $('#net-n').onchange = drawNet; $('#net-cut').onchange = drawNet;
  async function drawNet() {
    let E; try { E = await edges(sp, B.setId); } catch (e) { if (!gone()) netBox.innerHTML = `<div class="empty">${esc(e.message)}</div>`; return; }
    if (gone()) return;
    const n = +$('#net-n').value, c = CUT[+$('#net-cut').value];
    const nb = [...(E.adj.get(P.i) || new Map())].filter(([, e]) => e.best >= c).sort((a, b) => b[1].best - a[1].best).slice(0, n);
    netBox.innerHTML = '<svg></svg>';
    if (!nb.length) { netBox.innerHTML = '<div class="empty">No partners past this cutoff.</div>'; return; }
    const nodes = [{ id: P.i, row: P, q: true }, ...nb.map(([j, e]) => ({ id: j, row: sp.rows[j], e }))];
    const links = nb.map(([j, e]) => ({ source: P.i, target: j, best: e.best, avg: e.avg, q: true }));
    for (let a = 1; a < nodes.length; a++) { const mm = E.adj.get(nodes[a].id); if (!mm) continue;
      for (let b = a + 1; b < nodes.length; b++) { const e = mm.get(nodes[b].id); if (e && e.best >= c) links.push({ source: nodes[a].id, target: nodes[b].id, best: e.best, avg: e.avg }); } }
    const deg = new Map(); links.forEach((l) => { deg.set(l.source, (deg.get(l.source) || 0) + 1); deg.set(l.target, (deg.get(l.target) || 0) + 1); });
    const svg = d3.select(netBox).select('svg'), W = netBox.clientWidth, H = netBox.clientHeight;
    svg.attr('width', W).attr('height', H);
    const g = svg.append('g');
    const zoom = d3.zoom().scaleExtent([0.2, 4]).on('zoom', (ev) => g.attr('transform', ev.transform));
    svg.call(zoom);
    const link = g.append('g').selectAll('line').data(links).join('line').attr('stroke', (d) => ecol(d.best)).attr('stroke-width', (d) => EWID(d.avg))
      .attr('stroke-opacity', (d) => (d.q ? 0.5 : 0.85)).attr('stroke-linecap', 'round');
    link.filter((d) => !d.q).raise();
    link.on('mousemove', (ev, d) => showTip(`<b>${esc(sp.rows[typeof d.source === 'object' ? d.source.id : d.source].gene)}</b> × <b>${esc(sp.rows[typeof d.target === 'object' ? d.target.id : d.target].gene)}</b> · iLIS best ${d.best.toFixed(3)} · average ${d.avg.toFixed(3)}`, ev.clientX, ev.clientY)).on('mouseleave', hideTip);
    const r = (d) => (d.q ? 18 : 6 + Math.min(10, Math.sqrt(deg.get(d.id) || 1) * 1.6));
    const node = g.append('g').selectAll('g').data(nodes).join('g').style('cursor', 'pointer')
      .call(d3.drag().on('start', (ev, d) => { if (!ev.active) sim.alphaTarget(0.25).restart(); d.fx = d.x; d.fy = d.y; })
        .on('drag', (ev, d) => { d.fx = ev.x; d.fy = ev.y; }).on('end', (ev, d) => { if (!ev.active) sim.alphaTarget(0); if (!d.q) { d.fx = null; d.fy = null; } }));
    const circle = node.append('circle').attr('r', r).attr('stroke', '#fff').attr('stroke-width', (d) => (d.q ? 2.5 : 1.8));
    node.append('text').text((d) => d.row.gene).attr('text-anchor', 'middle').attr('dy', (d) => -r(d) - 6)   // every label, the center one too: dark text on a white edge
      .attr('font-family', 'IBM Plex Sans, sans-serif').attr('font-size', (d) => (d.q ? 14 : 11)).attr('font-weight', (d) => (d.q ? 700 : 600)).attr('fill', '#17263A')
      .attr('paint-order', 'stroke').attr('stroke', 'rgba(255,255,255,0.92)').attr('stroke-width', (d) => (d.q ? 4 : 3)).attr('stroke-linejoin', 'round');
    node.filter((d) => d.q).raise();
    NET = { recolor() {
      const k = M ? M.k : 1;
      circle.attr('fill', (d) => { if (d.q) return '#1A5276'; const cl = partnerCluster.get(d.row.key); return cl ? clusterColor(cl, k) : '#C3CCD6'; });
      const used = [...new Set(nodes.filter((d) => !d.q).map((d) => partnerCluster.get(d.row.key)).filter(Boolean))].sort((a, b) => a - b);
      $('#net-legend').innerHTML = used.length ? used.map((cl) => `<span><i style="background:${clusterColor(cl, k)};border-radius:50%"></i>${clusterLabel(cl)}</span>`).join('') + '<span><i style="background:#C3CCD6;border-radius:50%"></i>not clustered at this cutoff</span>' : '';
    } };
    NET.recolor();
    node.on('mousemove', (ev, d) => showTip(d.q ? `<b>${esc(d.row.gene)}</b> · ${fmtInt(d.row.pos10)} partners past 10% FPR` : `<b>${esc(d.row.gene)}</b> · iLIS best ${d.e.best.toFixed(3)} · average ${d.e.avg.toFixed(3)} with ${esc(P.gene)}<br>${srcBadges(sp, d.e.src)}`, ev.clientX, ev.clientY))
      .on('mouseleave', hideTip)
      .on('click', (ev, d) => { if (!d.q && !ev.defaultPrevented) location.hash = `#/${sp.id}/${d.row.key}${scopeQ}`; });
    const q = nodes[0]; q.fx = W / 2; q.fy = H / 2;
    const R = Math.min(W, H) * 0.42;
    const sim = d3.forceSimulation(nodes).force('link', d3.forceLink(links).id((d) => d.id).distance((l) => (l.q ? R * (1.15 - 0.55 * Math.min(1, l.best)) : 60)).strength((l) => (l.q ? 0.5 : 0.35)))
      .force('charge', d3.forceManyBody().strength(-340)).force('collide', d3.forceCollide().radius((d) => r(d) + 14)).force('x', d3.forceX(W / 2).strength(0.04)).force('y', d3.forceY(H / 2).strength(0.05))
      .on('tick', () => { link.attr('x1', (d) => d.source.x).attr('y1', (d) => d.source.y).attr('x2', (d) => d.target.x).attr('y2', (d) => d.target.y); node.attr('transform', (d) => `translate(${d.x},${d.y})`); });
    let fitted = false;   // once the layout settles, zoom so every node and label fits the box (the zoom stays free afterwards)
    sim.on('end', () => { if (fitted) return; fitted = true;
      const xs = nodes.map((d) => d.x), ys = nodes.map((d) => d.y);
      const x0 = Math.min(...xs) - 48, x1 = Math.max(...xs) + 48, y0 = Math.min(...ys) - 34, y1 = Math.max(...ys) + 24;
      const s = Math.min(1, 0.96 * Math.min(W / (x1 - x0), H / (y1 - y0)));
      svg.transition().duration(450).call(zoom.transform, d3.zoomIdentity.translate(W / 2 - s * (x0 + x1) / 2, H / 2 - s * (y0 + y1) / 2).scale(s)); });
    svgExport($('#net-x'), `atlas_${P.gene}_network`, () => $('svg', netBox));
  }

  drawOverview(); drawTable(); fillPartners(); drawScatter(); drawFreq(); drawHeatmap(); renderClusterInfo(); legend3D();
  loadStructure();
  cluster();
  let rsz; window.onresize = () => { clearTimeout(rsz); rsz = setTimeout(() => { drawFreq(); drawHeatmap(); drawOverview(); drawScatter(); const rb = $('#res-body'); if (rb && rb._redraw) rb._redraw(); }, 150); };
}

/* ── pair page ───────────────────────────────────────────────────────────────────────────────────────── */
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
  document.title = `${P.gene} · ${O.gene} · LIVIA cLIP Atlas`;
  const crumbs = `<div class="crumbs"><a href="#/">Atlas</a> / <a href="#/${sp.id}">${esc(sp.reg.label)}</a> / <a href="#/${sp.id}/${P.key}${scopeQ}">${esc(P.gene)}</a> / ${esc(O.gene)}</div>`;
  app.innerHTML = crumbs + '<div class="loading">Loading…</div>';
  let B;
  try { B = await merged(sp, P, scope); } catch (e) { B = { partners: [] }; }
  if (stale(gen)) return;
  const part = B.partners.find((p) => p.id === O.key);
  if (!part) { app.innerHTML = crumbs + `<div class="empty">${esc(P.gene)} and ${esc(O.gene)} were not predicted together${scope ? ` in ${esc(label)}. <a href="#/${sp.id}/${P.key}/${O.key}">Every screen</a>` : ' in these screens'}.</div>`; return; }
  const best = [...part.preds].sort((a, b) => b.iLIS - a.iLIS)[0], b = bandOf(part.best), oLen = O.clen || part.preds[0].pLen;
  const who = (R, len, col) => `<div class="who"><b style="color:${col}">${esc(R.gene)}</b> <span>${esc(short(R.name) || '')}</span>
      <div class="ids">${R.acc ? uniprotLink(R.acc) : '<span>no UniProt entry</span>'}<span>${esc(R.id)}</span>${len ? `<span>${fmtInt(len)} aa</span>` : ''}</div></div>`;
  const num = (v, d) => `<td class="n">${fmtNum(v, d)}</td>`;
  const lab = (p) => runLabel(sp, B, p.run, P);
  app.innerHTML = `${crumbs}
    <div class="phead"><div><h1><span style="color:var(--query)">${esc(P.gene)}</span> <span style="color:var(--ink-3);font-weight:600">×</span> <span style="color:var(--partner)">${esc(O.gene)}</span></h1>
        <div class="pairwho">${who(P, P.clen, 'var(--query)')}${who(O, oLen, 'var(--partner)')}</div>
        <div class="srcs">Predicted in ${part.sets.length ? setBadges(part.sets) : srcBadges(sp, part.src)}${scope ? ` <span class="muted">· only ${esc(label)} models shown · <a href="#/${sp.id}/${P.key}/${O.key}">every model</a></span>` : ''}</div>
        <div class="actions"><a class="btn" href="#/${sp.id}/${P.key}">${esc(P.gene)} page</a>${O0 ? `<a class="btn" href="#/${sp.id}/${O.key}">${esc(O.gene)} page</a>` : ''}</div></div>
      <div class="kpis"><div class="kpi"><b style="color:${BAND[b]}">${part.best.toFixed(3)}</b><span>iLIS best · ${bandLabel[b]}</span></div><div class="kpi"><b>${part.ilisaBest.toFixed(1)}</b><span>iLISA best</span></div>
        <div class="kpi"><b style="color:${bandCol(FPR_AVG.iLIS, part.avg)}">${part.avg.toFixed(3)}</b><span>iLIS average · ${bandLabel[bandIn(FPR_AVG.iLIS, part.avg)]}</span></div>
        <div class="kpi"><b style="color:${bandCol(FPR.ipTM, part.iptmBest)}">${part.iptmBest.toFixed(2)}</b><span>ipTM best · ${bandLabel[bandIn(FPR.ipTM, part.iptmBest)]}</span></div></div></div>
    <div class="card"><div class="card-head"><h2>Ranked models</h2><span class="muted">every model of every screen · click one to show its interface</span></div>
      <div class="tbl-wrap"><table class="pt models"><thead>
        <tr><th rowspan="2">Source</th><th rowspan="2">Rank</th><th rowspan="2" class="n">iLIS</th><th rowspan="2" class="n">iLISA</th><th rowspan="2" class="n">ipTM</th><th rowspan="2" class="n">LIS</th><th rowspan="2" class="n">cLIS</th>
          <th colspan="2" class="grp">Interface residues (LIR)</th><th colspan="2" class="grp">Contact residues (cLIR)</th></tr>
        <tr><th class="n sub q">${esc(P.gene)}</th><th class="n sub p">${esc(O.gene)}</th><th class="n sub q">${esc(P.gene)}</th><th class="n sub p">${esc(O.gene)}</th></tr></thead>
        <tbody>${part.preds.map((p, i) => `<tr data-i="${i}"><td><span class="src" style="--c:${runColor(sp, p)}">${esc(lab(p))}</span></td><td>${p.rank}</td>
          <td class="n">${fmtNum(p.iLIS, 3)} <span class="band b${bandOf(p.iLIS)}">${bandLabel[bandOf(p.iLIS)]}</span></td>
          ${num(p.iLISA, 1)}<td class="n" style="color:${bandCol(FPR.ipTM, p.ipTM)};font-weight:600">${fmtNum(p.ipTM, 2)}</td>${num(p.LIS, 3)}${num(p.cLIS, 3)}${num(p.qLIR, 0)}${num(p.pLIR, 0)}${num(p.qcLIR, 0)}${num(p.pcLIR, 0)}</tr>`).join('')}</tbody></table></div></div>
    <div class="card"><div class="card-head"><h2>Interaction Residues</h2><select id="model-pick" aria-label="Model" style="font:13px var(--sans);padding:5px 8px;border:1px solid var(--line);border-radius:7px">${part.preds.map((p, i) =>
        `<option value="${i}">${esc(lab(p))} · rank ${p.rank} · iLIS ${fmtNum(p.iLIS, 3)}</option>`).join('')}</select></div>
      <div class="legend" style="margin:2px 0 12px"><span><i style="background:#E0E0E0"></i>not in the interface</span><span><i style="background:#80CBC4"></i><i style="background:#FFAB91;margin-left:-2px"></i>interface (LIR: PAE ≤ 12 Å)</span><span><i style="background:#00897B"></i><i style="background:#E64A19;margin-left:-2px"></i>contact (cLIR: also Cβ ≤ 8 Å)</span></div>
      <div id="iface"></div></div>`;
  const pick = (i) => {
    $('#model-pick').value = String(i);
    app.querySelectorAll('.models tbody tr').forEach((x) => x.classList.toggle('on', +x.dataset.i === i));
    ifaceView($('#iface'), { sp, P, O, pred: part.preds[i], B, canvasId: 'iface-canvas' });
  };
  $('#model-pick').onchange = (e) => pick(+e.target.value);
  app.querySelectorAll('.models tbody tr').forEach((tr) => tr.onclick = () => pick(+tr.dataset.i));
  pick(part.preds.indexOf(best));
  let rsz; window.onresize = () => { clearTimeout(rsz); rsz = setTimeout(() => { const f = $('#iface'); if (f && f._redraw) f._redraw(); }, 150); };
}

/* species page: the merged index — its screens, counts and most connected proteins */
async function viewSpecies(spId) {
  const gen = ROUTE, sp = await species(spId), c = sp.manifest.counts, hubs = [...sp.rows].sort((a, b) => b.pos10 - a.pos10).slice(0, 24);
  const TSs = await Promise.all(sp.dsIds.map(async (id) => { try { return await setsOf(await dataset(id)); } catch (e) { return null; } }));
  if (stale(gen)) return;
  document.title = `${sp.reg.label} · LIVIA cLIP Atlas`;
  app.innerHTML = `<div class="crumbs"><a href="#/">Atlas</a> / ${esc(sp.reg.label)}</div>
    <div class="dshead"><h1>${esc(sp.reg.label)} protein interactions</h1><div class="pname"><i>${esc(sp.reg.name)}</i> · ${sp.dsIds.length === 1 ? 'one screen' : sp.dsIds.length + ' screens'}, one page per ${sp.manifest.keyedBy ? 'gene' : 'protein'}</div></div>
    ${kpiRow(c)}
    <div class="card"><h2>Search</h2><div id="sp-search" style="margin-top:10px"></div></div>
    <div class="card"><h2>Screens</h2><div class="screens">${sp.manifest.datasets.map((d, di) => `<div class="screen"><span class="src" style="--c:${sp.dsColor[di]}">${esc(d.short)}</span>
      <div><a href="#/datasets/${d.id}"><b>${esc(d.title)}</b></a><div class="muted">${fmtInt(d.counts.proteins)} proteins · ${fmtInt(d.counts.pairs)} pairs · ${d.counts.runs ? `${fmtInt(d.counts.runs)} predictions · ` : ''}${fmtInt(d.counts.predictions)} models</div>
      <div class="cite">${d.url ? `<a href="${esc(d.url)}" target="_blank" rel="noopener">${esc(d.citation)} ↗</a>` : esc(d.citation)}</div></div></div>`).join('')}</div>
      ${sp.dsIds.length > 1 ? `<p class="muted" style="margin:10px 0 0">${fmtInt(c.pairsInSeveral)} pairs were predicted in more than one screen; their pages keep every model with its source.</p>` : ''}</div>
    <div class="card"><h2>Most connected proteins <span class="muted">partners past the 10% FPR cutoff, any screen</span></h2>
      <div class="chips">${hubs.map((r) => `<a class="chip" href="#/${sp.id}/${r.key}">${esc(r.gene)} <span class="num" style="color:var(--ink-3)">${fmtInt(r.pos10)}</span></a>`).join('')}</div></div>
    ${TSs.map(setsCard).join('')}`;
  mountSearch($('#sp-search'), { spId });
}

/* ── router ─────────────────────────────────────────────────────────────────────────────────────────── */
// GoatCounter counts page loads only; the atlas routes by #/…, so each view is counted by hand, on the atlas's own
// counter (livia-atlas.goatcounter.com), as /human/P04637, /fly/FBgn0000117 …
function trackView() {
  const count = () => window.goatcounter && window.goatcounter.count && window.goatcounter.count({ path: location.hash.replace(/^#/, '') || '/', title: document.title });
  if (window.goatcounter && window.goatcounter.count) count(); else window.addEventListener('load', () => setTimeout(count, 0), { once: true });
}
const hashPath = () => { const [path, q] = location.hash.replace(/^#\/?/, '').split('?'); return { parts: path.split('/').filter(Boolean).map(decodeURIComponent), q: new URLSearchParams(q || '') }; };
async function route() {
  const gen = ++ROUTE, { parts, q } = hashPath(), setId = q.get('set') || '', here = location.hash;
  window.scrollTo(0, 0); hideTip(); window.onresize = null;
  document.title = 'LIVIA cLIP Atlas';
  stopClip();
  try {
    await registry();
    if (stale(gen)) return;
    if (!parts.length) await viewHome();
    else if (parts[0] === 'datasets') await (parts[2] ? viewSet(parts[1], parts[2]) : parts[1] ? viewDataset(parts[1]) : viewDatasets());
    else if (parts[0] === 'about') viewAbout();
    else if (parts[0] === 'themes' && parts[1]) await viewTheme(parts[1]);
    else if (await regSpecies(parts[0])) { if (parts.length === 1) await viewSpecies(parts[0]); else if (parts.length === 2) await viewProtein(parts[0], parts[1], setId); else await viewPair(parts[0], parts[1], parts[2], setId); }
    else if (await regDataset(parts[0])) {   // links from before the species pages: #/<screen>/<name>[/<name>]
      const d = await regDataset(parts[0]);
      if (parts.length === 1 || !d.species) { location.replace(`#/datasets/${d.id}`); return; }
      const sp = await species(d.species), k = (n) => { const r = sp.byName.get(n); return r ? r.key : n; };
      location.replace(`#/${sp.id}/${parts.slice(1).map(k).join('/')}`); return;
    }
    else app.innerHTML = `<div class="empty">Nothing at “${esc(parts.join('/'))}”. <a href="#/">Go to the atlas home</a></div>`;
  } catch (e) { if (!stale(gen)) app.innerHTML = `<div class="empty">${esc(e.message || e)}</div>`; console.error(e); }
  if (stale(gen)) return;
  if (location.hash === here) trackView();   // not for a view that redirected
  if (!$('#top-search').firstChild) mountSearch($('#top-search'));
}
window.addEventListener('hashchange', route);
route();
