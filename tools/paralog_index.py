#!/usr/bin/env python3
"""Paralog families of a species, for the nested network: data/species/<sp>/paralogs.json.

Input: an MMseqs2 easy-cluster table (representative <tab> member, keys as the Atlas uses) and the species index
(data/species/<sp>/proteins.json). Output: {"index": first 12 hex of SHA-1 over the index keys, "identity": ...,
"families": [[row, row, ...], ...]}: every family of two or more rows of the index. A protein in no family is its own.
Usage: paralog_index.py <species dir> <cluster.tsv> "<identity note>" """
import sys, json, hashlib, collections
sp_dir, tsv, note = sys.argv[1], sys.argv[2], sys.argv[3]
prot = json.load(open(f'{sp_dir}/proteins.json')); kc = prot['columns'].index('key'); keys = [r[kc] for r in prot['rows']]
row = {k: i for i, k in enumerate(keys)}
fam, lost = collections.defaultdict(set), 0
for line in open(tsv):
    rep, mem = line.rstrip('\n').split('\t')
    if rep in row and mem in row: fam[rep].add(row[mem]); fam[rep].add(row[rep])
    else: lost += 1
fams = sorted((sorted(v) for v in fam.values() if len(v) >= 2), key=lambda f: f[0])
out = {'index': hashlib.sha1('\n'.join(keys).encode()).hexdigest()[:12], 'identity': note, 'families': fams}
json.dump(out, open(f'{sp_dir}/paralogs.json', 'w'), separators=(',', ':'))
man = json.load(open(f'{sp_dir}/manifest.json')); man.setdefault('files', {})['paralogs'] = 'paralogs.json'; json.dump(man, open(f'{sp_dir}/manifest.json', 'w'), indent=1)
print(f'{sp_dir}/paralogs.json: {len(fams)} families of 2+ covering {sum(map(len, fams))} rows; {lost} table rows outside the index')
