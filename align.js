'use strict';
// Pairwise protein alignment for the Atlas: Gotoh's affine-gap global alignment with free end gaps (a shorter or longer
// ortholog costs nothing at the ends), BLOSUM62, open −11, extend −1. Also a windowed identity along the first sequence,
// the confidence strip under an aligned track. Exposed as ALIGN: align(a, b) → { map, score, identity, aligned, cols, ident }, or null when the pair is too long (25M cells).
//   map[i]      (0-based residue of a) → 1-based residue of b, or null; aligned: residues of a with a partner;
//   identity    identical / aligned; cols: the alignment as [[ia|null, ib|null], …] (0-based);
//   ident(i, w) share of identical pairs among the residues of a within ±w of i that are aligned (0 when none).
const ALIGN = (() => {
  const AA = 'ARNDCQEGHILKMFPSTWYV';
  const B62 = [   // BLOSUM62, rows and columns in the order of AA
    [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0], [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],
    [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3], [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],
    [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1], [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],
    [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2], [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],
    [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3], [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],
    [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1], [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],
    [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1], [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2], [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],
    [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0], [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],
    [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1], [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4]];
  const IDX = new Int8Array(128).fill(-1); for (let i = 0; i < AA.length; i++) IDX[AA.charCodeAt(i)] = i;
  const sub = (x, y) => { const i = IDX[x.charCodeAt(0)], j = IDX[y.charCodeAt(0)]; return i < 0 || j < 0 ? (x === y ? 1 : -1) : B62[i][j]; };
  function align(a, b, open = -11, ext = -1) {
    a = a.toUpperCase(); b = b.toUpperCase();
    const n = a.length, m = b.length, W = m + 1, NEG = -1e9;
    if (!n || !m) return { map: new Array(n).fill(null), score: 0, identity: 0, aligned: 0, cols: [], ident: () => 0 };
    if ((n + 1) * (m + 1) > 25e6) return null;   // about 15 bytes a cell: two 5,000-residue proteins are the limit of what a page should hold
    // M: a[i] with b[j]; X: a[i] against a gap (gap in b); Y: b[j] against a gap (gap in a). End gaps are free.
    const M = new Float32Array((n + 1) * W), X = new Float32Array((n + 1) * W), Y = new Float32Array((n + 1) * W), T = new Uint8Array((n + 1) * W * 3);
    M.fill(NEG); X.fill(NEG); Y.fill(NEG); M[0] = 0;
    for (let j = 1; j <= m; j++) Y[j] = 0;   // leading gap in a: free
    for (let i = 1; i <= n; i++) X[i * W] = 0;   // leading gap in b: free
    for (let i = 1; i <= n; i++) {
      const ai = a[i - 1], row = i * W, up = (i - 1) * W;
      for (let j = 1; j <= m; j++) {
        const s = sub(ai, b[j - 1]), k = row + j, d = up + j - 1;
        let best = M[d], t = 0; if (X[d] > best) { best = X[d]; t = 1; } if (Y[d] > best) { best = Y[d]; t = 2; }
        M[k] = best + s; T[k * 3] = t;
        const endB = j === m;   // a trailing gap in b (a's tail unaligned) is free
        const xo = M[up + j] + (endB ? 0 : open), xe = X[up + j] + (endB ? 0 : ext), yo2 = Y[up + j] + (endB ? 0 : open);
        if (xo >= xe && xo >= yo2) { X[k] = xo; T[k * 3 + 1] = 0; } else if (xe >= yo2) { X[k] = xe; T[k * 3 + 1] = 1; } else { X[k] = yo2; T[k * 3 + 1] = 2; }
        const endA = i === n;   // a trailing gap in a (b's tail unaligned) is free
        const yo = M[row + j - 1] + (endA ? 0 : open), ye = Y[row + j - 1] + (endA ? 0 : ext), yx = X[row + j - 1] + (endA ? 0 : open);
        if (yo >= ye && yo >= yx) { Y[k] = yo; T[k * 3 + 2] = 0; } else if (ye >= yx) { Y[k] = ye; T[k * 3 + 2] = 2; } else { Y[k] = yx; T[k * 3 + 2] = 1; }   // codes as everywhere: 0 M, 1 X, 2 Y
      }
    }
    const last = n * W + m; let state = 0, score = M[last]; if (X[last] > score) { score = X[last]; state = 1; } if (Y[last] > score) { score = Y[last]; state = 2; }
    const map = new Array(n).fill(null), cols = []; let i = n, j = m, identical = 0, aligned = 0;
    while (i > 0 || j > 0) {
      if (i === 0) { cols.push([null, j - 1]); j--; continue; } if (j === 0) { cols.push([i - 1, null]); i--; continue; }
      const k = i * W + j;
      if (state === 0) { map[i - 1] = j; aligned++; if (a[i - 1] === b[j - 1]) identical++; cols.push([i - 1, j - 1]); state = T[k * 3]; i--; j--; }
      else if (state === 1) { cols.push([i - 1, null]); state = T[k * 3 + 1]; i--; }
      else { cols.push([null, j - 1]); state = T[k * 3 + 2]; j--; }
    }
    cols.reverse();
    const same = new Uint8Array(n); for (let p = 0; p < n; p++) if (map[p] != null && a[p] === b[map[p] - 1]) same[p] = 1;
    const ident = (p, w = 10) => { let s = 0, t = 0; for (let q = Math.max(0, p - w); q <= Math.min(n - 1, p + w); q++) if (map[q] != null) { t++; s += same[q]; } return t ? s / t : 0; };
    return { map, score, identity: aligned ? identical / aligned : 0, aligned, cols, ident };
  }
  return { align, sub };
})();
if (typeof module !== 'undefined') module.exports = ALIGN;
