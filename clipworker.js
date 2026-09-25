'use strict';
// LIVIA cLIP for one protein, computed with LIVIA's own cLIP modules (loaded from the LIVIA site, so the atlas and
// clip.html always cluster the same way): fingerprints of the protein's contact residues per prediction past the
// cutoff → cosine distance → average linkage → silhouette-optimal k → dendrogram order. Clusters are renumbered by
// size as in clip.html: Cluster 1 is the largest (ties: lower original label first).
// Input: lis.py rows as objects (header → value), already merged across screens and named "<query>___<partner>".
let loaded = false;
self.onmessage = (ev) => {
  const { id, livia, rows, gene, cut } = ev.data;
  try {
    if (!loaded) { importScripts(livia + 'js/clip-pipeline.js', livia + 'js/clip-clustering.js'); loaded = true; }
    const pipe = CLIPPipeline.runPipelineRows(CLIPPipeline.convertRows(rows), gene, { ilisCutoff: cut, topN: Infinity, sortBy: 'iLIS' });
    const n = pipe.fingerprints.length;
    let k = n ? 1 : 0, labels = n ? [1] : [], Z = [], order = n ? [0] : [];
    if (n >= 2) {
      const cl = CLIPCluster.clusterInteractions(pipe.fingerprints, { maxK: 30, dropThreshold: 0.05 });
      const cnt = {}; for (const l of cl.labels) cnt[l] = (cnt[l] || 0) + 1;
      const remap = {}; Object.keys(cnt).map(Number).sort((x, y) => cnt[y] - cnt[x] || x - y).forEach((old, i) => { remap[old] = i + 1; });
      k = cl.optimalK; labels = cl.labels.map((l) => remap[l]); Z = cl.Z; order = CLIPCluster.leafOrder(cl.Z, n);
    }
    self.postMessage({ id, ok: true, k, labels, order, Z, plen: pipe.proteinLen, fingerprints: pipe.fingerprints,
      preds: pipe.rows.map((r) => ({ partner: String(r.Symbol_2), rank: +r.Rank, iLIS: +r.iLIS })) });
  } catch (e) {
    self.postMessage({ id, ok: false, message: String(e && e.message || e) });
  }
};
