'use strict';
// Sequence search for LIVIA Atlas: shared 5-mers between the query and every sequence of a dataset (seqs.fasta),
// encoded as integers so a query of any length scans ~11M residues in well under a second. Exact and
// substring (fragment) matches are reported first; the rest are ranked by shared 5-mers.
const AA = 'ACDEFGHIKLMNPQRSTVWY', K = 5, BASE = 21, SPACE = BASE ** K;
const CODE = new Uint8Array(128).fill(20);
for (let i = 0; i < AA.length; i++) { CODE[AA.charCodeAt(i)] = i; CODE[AA.toLowerCase().charCodeAt(i)] = i; }
let IDS = null, SEQS = null, ENC = null;

function encode(s) { const a = new Uint8Array(s.length); for (let i = 0; i < s.length; i++) a[i] = CODE[s.charCodeAt(i) & 127]; return a; }

async function init(url) {
  const text = await (await fetch(url)).text();
  IDS = []; SEQS = [];
  let id = null, buf = [];
  for (const line of text.split('\n')) {
    if (line.startsWith('>')) { if (id) { IDS.push(id); SEQS.push(buf.join('')); } id = line.slice(1).trim(); buf = []; }
    else if (line) buf.push(line.trim());
  }
  if (id) { IDS.push(id); SEQS.push(buf.join('')); }
  ENC = SEQS.map(encode);
}

function search(q) {
  q = q.toUpperCase().replace(/[^A-Z]/g, '');
  const qe = encode(q), n = qe.length - K + 1;
  if (n < 1) return [];
  const mark = new Uint8Array(SPACE);
  let nq = 0;
  for (let i = 0; i < n; i++) {
    let c = 0, bad = false;
    for (let j = 0; j < K; j++) { const v = qe[i + j]; if (v === 20) { bad = true; break; } c = c * BASE + v; }
    if (!bad && !mark[c]) { mark[c] = 1; nq++; }
  }
  const hits = [], stamp = new Int32Array(SPACE);          // stamp[c] = s + 1 once k-mer c was counted for sequence s
  for (let s = 0; s < ENC.length; s++) {
    const e = ENC[s], m = e.length - K + 1;
    if (m < 1) continue;
    let shared = 0;
    for (let i = 0; i < m; i++) {
      let c = 0, bad = false;
      for (let j = 0; j < K; j++) { const v = e[i + j]; if (v === 20) { bad = true; break; } c = c * BASE + v; }
      if (!bad && mark[c] && stamp[c] !== s + 1) { stamp[c] = s + 1; shared++; }
    }
    if (shared) hits.push([s, shared]);
  }
  hits.sort((a, b) => b[1] - a[1]);
  return hits.slice(0, 12).map(([s, shared]) => {
    const seq = SEQS[s], exact = seq === q, within = !exact && seq.includes(q), contains = !exact && !within && q.includes(seq);
    return { id: IDS[s], length: seq.length, shared, queryKmers: nq, frac: nq ? shared / nq : 0,
      kind: exact ? 'identical' : within ? 'fragment of this protein' : contains ? 'contains this protein' : 'similar' };
  });
}

self.onmessage = async (ev) => {
  const { type, url, seq, id } = ev.data;
  try {
    if (type === 'init') { if (!ENC) await init(url); self.postMessage({ type: 'ready', n: IDS.length }); }
    else if (type === 'query') { if (!ENC) await init(url); self.postMessage({ type: 'result', id, hits: search(seq) }); }
  } catch (e) { self.postMessage({ type: 'error', id, message: String(e && e.message || e) }); }
};
