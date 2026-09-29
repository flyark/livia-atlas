#!/usr/bin/env python3
"""Every pair folded in a species' screens, for the network views: data/species/<sp>/tested.bin.

Input: the pair-level tables of tested pairs (one row per unordered pair per screen, keyed like the Atlas; built on O2
from the warehouse's clean tables) and the species index (data/species/<sp>/proteins.json, whose row order the file uses).
Output, little-endian:
  bytes 0-5   b'LVTP2\\0'
  bytes 6-9   uint32 N, the rows of the species index
  bytes 10-13 uint32 P, the pairs (each once, i < j; homodimers i == j included)
  bytes 14-25 ASCII, the first 12 hex digits of SHA-1 over the index keys joined by newlines (as biogrid.tsv carries)
  bytes 26-31 zero
  then uint32 offsets[N + 1], uint16 partners[P] (row i's tested partners j >= i, sorted, at partners[offsets[i]:offsets[i + 1]])
  and uint8 scores[P]: the pair's best iLIS over every model of every screen, as round(iLIS * 255)
Keys that are not rows of the index (constructs the index does not list) are counted and left out.
Usage: tested_index.py <species dir> <pairs.parquet> [<pairs.parquet> ...]"""
import sys, json, hashlib, struct
import numpy as np, pandas as pd

sp_dir, files = sys.argv[1], sys.argv[2:]
prot = json.load(open(f'{sp_dir}/proteins.json'))
kc = prot['columns'].index('key'); keys = [r[kc] for r in prot['rows']]; N = len(keys)
assert N < 65536, 'row numbers must fit in 16 bits'
row = {k: i for i, k in enumerate(keys)}
sha = hashlib.sha1('\n'.join(keys).encode()).hexdigest()[:12]
parts, lost = [], 0
for f in files:
    cols = pd.read_parquet(f).columns; t = pd.read_parquet(f, columns=['protein_a', 'protein_b'] + (['iLIS_best'] if 'iLIS_best' in cols else ['iLIS']))
    t['s'] = t[t.columns[2]].fillna(0).clip(0, 1)
    a, b = t.protein_a.map(row), t.protein_b.map(row); ok = a.notna() & b.notna(); lost += int((~ok).sum())
    a, b = a[ok].astype(np.int64).values, b[ok].astype(np.int64).values
    sc = t.s[ok].values
    parts.append(pd.DataFrame({'i': np.minimum(a, b), 'j': np.maximum(a, b), 's': sc}))
    print(f.split('/')[-1], len(t), 'pairs,', int((~ok).sum()), 'with a key outside the index', flush=True)
df = pd.concat(parts).groupby(['i', 'j'], as_index=False)['s'].max().sort_values(['i', 'j'])   # a pair folded in two screens counts once, its best iLIS
pr = df[['i', 'j']].values.astype(np.int64); scores = np.round(df.s.values * 255).astype(np.uint8)
off = np.zeros(N + 1, dtype=np.uint32); np.add.at(off, pr[:, 0] + 1, 1); off = np.cumsum(off, dtype=np.uint32)
nbr = pr[:, 1].astype(np.uint16)
with open(f'{sp_dir}/tested.bin', 'wb') as out:
    out.write(b'LVTP2\0' + struct.pack('<II', N, len(pr)) + sha.encode() + b'\0' * 6); out.write(off.astype('<u4').tobytes()); out.write(nbr.astype('<u2').tobytes()); out.write(scores.tobytes())
man = json.load(open(f'{sp_dir}/manifest.json')); man.setdefault('files', {})['tested'] = 'tested.bin'
man['counts']['pairsTested'] = int(len(pr))
json.dump(man, open(f'{sp_dir}/manifest.json', 'w'), indent=1)   # the manifests' own format: ASCII escapes, no final newline
print(f'{sp_dir}/tested.bin: N {N}, {len(pr)} pairs ({int((pr[:, 0] == pr[:, 1]).sum())} homodimers), {lost} rows left out, index {sha}')
