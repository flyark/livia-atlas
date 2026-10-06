'use strict';
// LIVIA cLIP for one protein, computed with LIVIA's own cLIP modules (loaded from the LIVIA site, so the atlas and
// clip.html always cluster the same way): fingerprints of the protein's contact residues per prediction past the
// cutoff → cosine distance → average linkage → silhouette-optimal k → dendrogram order. Clusters are renumbered by
// size as in clip.html: Cluster 1 is the largest (ties: lower original label first).
// The distances and the linkage come from clipfast.js, which builds the same tree as LIVIA's module in a fraction of
// the time (LIVIA's loop rescans every pair at every merge: a 7,000-prediction interactome took it minutes).
// A very large interactome is thinned before clustering, and the page is told: above LIMIT predictions only the
// top-ranked passing model of each pair is clustered (cLIP itself takes every passing model), and above LIMIT pairs
// the pairs with the highest iLIS; the residue frequency and the sites then count pairs, not models.
// Input: lis.py rows as objects (header → value), already merged across screens and named "<query>___<partner>".
// Messages back: { id, stage: 'clustering', n, pairs, thinned } while it runs, then { id, ok: true, … }.
const LIMIT = 4000;
let loaded = false;
self.onmessage = (ev) => {
  const { id, livia, rows, gene, cut } = ev.data;
  try {
    if (!loaded) { importScripts(livia + 'js/clip-pipeline.js?v=20261006a', livia + 'js/clip-clustering.js', 'clipfast.js?v=20261002a'); loaded = true; }
    let use = rows, thinned = null;
    const pairs = (rs) => new Set(rs.map((r) => r.name)).size;
    if (rows.length > LIMIT) {   // one model per pair: its top-ranked model among those that pass (ranks start at 1, or at 0 for AlphaFold 3 output)
      const minRank = new Map(); for (const r of rows) { const k = +r.rank; if (!minRank.has(r.name) || k < minRank.get(r.name)) minRank.set(r.name, k); }
      use = rows.filter((r) => +r.rank === minRank.get(r.name)); thinned = { from: rows.length, how: 'rank1', pairs: pairs(use) };
      if (use.length > LIMIT) { use = [...use].sort((a, b) => +b.iLIS - +a.iLIS).slice(0, LIMIT); thinned.how = 'top'; thinned.pairs = pairs(use); }
      thinned.to = use.length;
    }
    const pipe = CLIPPipeline.runPipelineRows(CLIPPipeline.convertRows(use), gene, { ilisCutoff: cut, topN: Infinity, sortBy: 'iLIS' });
    const n = pipe.fingerprints.length;
    self.postMessage({ id, stage: 'clustering', n, pairs: pairs(pipe.rows), thinned });
    let k = n ? 1 : 0, labels = n ? [1] : [], Z = [], order = n ? [0] : [];
    if (n >= 2) {
      const D = CLIPFast.cosineDistMatrix(pipe.fingerprints);
      Z = CLIPFast.averageLinkage(D, n);
      const optimalK = CLIPCluster.findOptimalK(Z, D, n, 30, 0.05), raw = CLIPCluster.fcluster(Z, n, optimalK);
      const cnt = {}; for (const l of raw) cnt[l] = (cnt[l] || 0) + 1;
      const remap = {}; Object.keys(cnt).map(Number).sort((x, y) => cnt[y] - cnt[x] || x - y).forEach((old, i) => { remap[old] = i + 1; });
      k = optimalK; labels = raw.map((l) => remap[l]); order = CLIPCluster.leafOrder(Z, n);
    }
    self.postMessage({ id, ok: true, k, labels, order, Z, plen: pipe.proteinLen, fingerprints: pipe.fingerprints, thinned,
      preds: pipe.rows.map((r) => ({ partner: String(r.Symbol_2), rank: +r.Rank, iLIS: +r.iLIS })) });
  } catch (e) {
    self.postMessage({ id, ok: false, message: String(e && e.message || e) });
  }
};
