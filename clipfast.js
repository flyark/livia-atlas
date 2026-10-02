'use strict';
// Average-linkage clustering of cLIP fingerprints for large interactomes, giving the same tree as LIVIA's
// CLIPCluster.averageLinkage (clip-clustering.js) in a fraction of the time and memory: the distance matrix is one
// Float32Array row per fingerprint, and the next merge is found from each row's cached nearest neighbor (O(n) per
// merge) instead of rescanning every pair (O(n²) per merge, O(n³) in all). Ties break as LIVIA's loop does: the smallest
// distance, then the lowest id, then the lowest partner id; merged nodes take the ids n, n + 1, … in merge order.
// Exposed as CLIPFast (a worker's global): cosineDistMatrix, averageLinkage.
const CLIPFast = (() => {
  // cosine distance between residue sets, 1 − |A∩B| / (√|A| √|B|), as a square matrix of Float32Array rows
  function cosineDistMatrix(fingerprints) {
    const n = fingerprints.length, D = new Array(n);
    for (let i = 0; i < n; i++) D[i] = new Float32Array(n);
    let L = 0; for (const f of fingerprints) for (const r of f) if (r > L) L = r;
    const W = (L >> 5) + 1, bits = new Uint32Array(n * W), norm = new Float64Array(n);   // each fingerprint as a bitset
    for (let i = 0; i < n; i++) { const o = i * W; for (const r of fingerprints[i]) bits[o + (r >> 5)] |= 1 << (r & 31); norm[i] = Math.sqrt(fingerprints[i].length); }
    const pop = (x) => { x -= (x >>> 1) & 0x55555555; x = (x & 0x33333333) + ((x >>> 2) & 0x33333333); return (((x + (x >>> 4)) & 0x0F0F0F0F) * 0x01010101) >>> 24; };
    for (let i = 0; i < n; i++) {
      const oi = i * W, Di = D[i];
      for (let j = i + 1; j < n; j++) {
        const oj = j * W; let inter = 0;
        for (let w = 0; w < W; w++) { const x = bits[oi + w] & bits[oj + w]; if (x) inter += pop(x); }
        let d = norm[i] === 0 || norm[j] === 0 ? 1 : 1 - inter / (norm[i] * norm[j]);
        if (d < 0) d = 0;
        Di[j] = d; D[j][i] = d;
      }
    }
    return D;
  }
  // average linkage (UPGMA, Lance–Williams) → Z: [[idA, idB, dist, size], …] with merged ids n, n + 1, …
  function averageLinkage(D0, n) {
    if (n < 2) return [];
    const W = D0.map((row) => Float64Array.from(row));   // working distances, Float64 as LIVIA keeps them
    const id = new Int32Array(n), size = new Int32Array(n), alive = new Uint8Array(n).fill(1), nn = new Int32Array(n), nd = new Float64Array(n);
    for (let s = 0; s < n; s++) { id[s] = s; size[s] = 1; }
    const EPS = 1e-12;
    const rescan = (s) => {   // the nearest other slot: the smallest distance, ties to the lowest id
      const row = W[s]; let best = Infinity, bi = -1, bid = Infinity;
      for (let t = 0; t < n; t++) { if (t === s || !alive[t]) continue; const d = row[t];
        if (d < best - EPS || (Math.abs(d - best) <= EPS && id[t] < bid)) { best = d; bi = t; bid = id[t]; } }
      nn[s] = bi; nd[s] = best;
    };
    for (let s = 0; s < n; s++) rescan(s);
    const Z = []; let nextId = n, left = n;
    while (left > 1) {
      let best = Infinity, bs = -1, bid = Infinity;   // the merge: the smallest nd, ties to the lowest id (LIVIA scans ids ascending)
      for (let s = 0; s < n; s++) { if (!alive[s]) continue; const d = nd[s];
        if (d < best - EPS || (Math.abs(d - best) <= EPS && id[s] < bid)) { best = d; bs = s; bid = id[s]; } }
      let a = bs, b = nn[bs];
      if (id[a] > id[b]) { const t = a; a = b; b = t; }   // LIVIA lists the lower id first
      const na = size[a], nb = size[b], ia = id[a], ib = id[b];
      Z.push([ia, ib, best, na + nb]);
      for (let k = 0; k < n; k++) { if (!alive[k] || k === a || k === b) continue; const dn = (na * W[a][k] + nb * W[b][k]) / (na + nb); W[a][k] = dn; W[k][a] = dn; }
      alive[b] = 0; id[a] = nextId++; size[a] = na + nb; left--;
      for (let k = 0; k < n; k++) { if (!alive[k] || k === a) continue;   // rows that pointed at the merged pair look again; the rest check the new node
        if (nn[k] === a || nn[k] === b) rescan(k);
        else { const d = W[k][a]; if (d < nd[k] - EPS || (Math.abs(d - nd[k]) <= EPS && id[a] < id[nn[k]])) { nn[k] = a; nd[k] = d; } } }
      rescan(a);
    }
    return Z;
  }
  return { cosineDistMatrix, averageLinkage };
})();
if (typeof module !== 'undefined') module.exports = CLIPFast;
