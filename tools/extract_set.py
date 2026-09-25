"""Extract one thematic set of a LIVIA Atlas dataset zip as a table of predictions (CSV on stdout).

  python3 extract_set.py livia-atlas-flypredictome.zip kinase-tf > kinase-tf.csv
  python3 extract_set.py livia-atlas-flypredictome.zip            # lists the sets

A set is a list of prediction runs in sets.json. Each gene bundle b/<FBgn>.zip holds every run of the gene, so the set
is the rows whose `batch` is one of the set's runs; only the set's genes are read. A prediction sits in the bundles of
both of its proteins and is written once. Columns are LIVIA lis.py's, plus set (source category) and batch (run).
"""
import csv, io, json, sys, zipfile


def main():
    z = zipfile.ZipFile(sys.argv[1])
    sets = json.loads(z.read('sets.json'))['sets']
    if len(sys.argv) < 3:
        for s in sets:
            print(f"{s['id']:26} {s['type']:9} {s['counts']['predictions']:>10,} predictions  {s['title']}")
        return
    S = next((s for s in sets if s['id'] == sys.argv[2]), None)
    if not S:
        sys.exit(f"No set {sys.argv[2]!r}; the sets are: {', '.join(s['id'] for s in sets)}")
    runs = {str(r) for r in S['runs']}
    genes = [r[0] for r in json.loads(z.read(S['files']['proteins']))['rows']]
    out, seen, header, n = csv.writer(sys.stdout, lineterminator='\n'), set(), None, 0
    for k in genes:
        b = zipfile.ZipFile(io.BytesIO(z.read(f'b/{k}.zip')))
        rows = csv.reader(io.StringIO(b.read(next(f for f in b.namelist() if f.endswith('.csv'))).decode()))
        h = next(rows)
        if header is None:
            header = h; out.writerow(h)
        bi, ni, ri = h.index('batch'), h.index('name'), h.index('rank')
        for row in rows:
            if row[bi] in runs and (row[ni], row[ri], row[bi]) not in seen:
                seen.add((row[ni], row[ri], row[bi])); out.writerow(row); n += 1
    print(f"{n:,} predictions of {S['title']} ({len(genes):,} proteins)", file=sys.stderr)


if __name__ == '__main__':
    main()
