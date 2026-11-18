---
name: livia-atlas
description: >
  Query the proteome-scale LIVIA Atlas (https://flyark.github.io/livia-atlas/) and LIVIA
  (https://flyark.github.io/LIVIA/) for predicted protein-protein interactions, iLIS false-positive-rate
  scores, cLIP binding-site clusters, contact residues (cLIR), local interaction residues (LIR),
  AlphaFold DB pLDDT context, and ortholog conservation across 50+ species.
---

# LIVIA Atlas & LIVIA Agent Skill

Use **LIVIA Atlas** (`https://flyark.github.io/livia-atlas/`) to inspect proteome-scale AlphaFold-Multimer screens resolved to individual residues across human, mouse (`mus-musculus`), rat (`rattus-norvegicus`), *Drosophila* (`fly`), zebrafish, *C. elegans* (`worm`), yeast, *Arabidopsis*, *E. coli*, and viral proteomes.

## Metric & Cutoff Reference

- **LIR (local interaction residues)**: Cross-chain residue pairs with `PAE ≤ 12 Å`.
- **cLIR (contact residues)**: Subset of LIR that also satisfies `Cβ ≤ 8 Å` (`Cα` for Gly).
- **iLIS**: `sqrt(LIS * cLIS)` — primary interaction confidence score.
  - `10% FPR`: `iLIS ≥ 0.223`
  - `5% FPR`: `iLIS ≥ 0.339`
  - `1% FPR`: `iLIS ≥ 0.551`
- **cLIP footprint**: Residues contacted (`cLIR`) by `≥ 30%` of the passing models in a cluster.

## 1. Zero-Dependency CLI / Python Sandbox Queries (`atlas_query.py`)

When running in a terminal or Python code-execution environment without a browser, run `scripts/atlas_query.py` (uses only Python 3 standard library; fetches small byte-range slices directly from `https://flyark.github.io/livia-atlas/`):

```bash
# Resolve a gene symbol, UniProt accession, FlyBase ID, or alias (e.g., UBC6e -> UBE2J1)
python3 scripts/atlas_query.py resolve human UBC6e

# Protein cLIP binding sites, top partners, Pfam domains, and region evaluation (Markdown or --json)
python3 scripts/atlas_query.py protein human UBE2J1 --res 171-184 --markdown
python3 scripts/atlas_query.py protein human CXCR2 --res 10-19,188-202 --cut 10 --json

# Binary pair ranked models and exact cLIR / LIR spans + sequence stretches
python3 scripts/atlas_query.py pair human CXCR2 CXCL8 --res-a 10-19 --markdown
```

If you do not have a local clone of `livia-atlas`, download the single script first:
`curl -fsSL https://flyark.github.io/livia-atlas/scripts/atlas_query.py -o atlas_query.py`

## 2. Direct Hash URL Deep-Links (Browser Navigation)

- **Protein page**: `https://flyark.github.io/livia-atlas/#/<species>/<gene_or_accession>`
  - `?res=<query>`: Marks one or more regions (`?res=10-19`, `?res=10-19,188-202`, `?res=R175,R248,R273`, peptide `?res=SFEDFWKGED`, or motif `?res=LXXLL`, `?res=PX(2,3)P`) and populates `#res-look` in plain text with sequence, domain context, AlphaFold DB `pLDDT` mean/min–max, `cLIR` contacting partners/clusters, and `LIR`-only partners.
  - `?cut=10|5|1`: Sets the active FPR cutoff (`10%`, `5%`, or `1%`).
- **Pair page**: `https://flyark.github.io/livia-atlas/#/<species>/<proteinA>/<proteinB>`
  - `?rank=<n>`: Selects ranked model `n` (1-based by `iLIS`).
  - `?resA=<query>&resB=<query>`: Marks and evaluates specific residues/ranges/motifs on Protein A and/or Protein B in `#pair-res-look`.
- **Network builder**: `https://flyark.github.io/livia-atlas/#/network/<species>?ids=<g1,g2,...>&cut=5`
  - Add `&rc=<g1,g2,...>&focus=rc` to open the circular **Residue chord** diagram directly.

## 3. Browser JavaScript APIs (`window.LIVIA_ATLAS` & `window.LIVIA`)

When driving a live browser tab via Playwright, DevTools, or a browser subagent, call the built-in read-only JS APIs instead of scraping `<canvas>` elements:

```js
// On LIVIA Atlas (https://flyark.github.io/livia-atlas/):
const snapshot = await window.LIVIA_ATLAS.waitReady();
const md = window.LIVIA_ATLAS.toMarkdown();
const region = window.LIVIA_ATLAS.queryRegion('10-19,188-202', { mark: true });

// On LIVIA (universal.html, clip.html, flypredictome.html, dimer.html, monomer.html, network.html):
const pred = await window.LIVIA.waitReady();
const reportMd = window.LIVIA.toMarkdown();
```

Users can also click **Copy Markdown** (`#copy-md` in LIVIA Atlas, `#livia-copy-md` in LIVIA) on protein, pair, prediction, cLIP, and network pages to copy the structured report directly into an LLM chat.
