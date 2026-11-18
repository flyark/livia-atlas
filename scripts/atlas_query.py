#!/usr/bin/env python3
"""Zero-dependency CLI and Python sandbox client for LIVIA Atlas.

Queries https://flyark.github.io/livia-atlas/ (or a local checkout via --base)
using only the Python 3 standard library. Fetches compact species indexes and
byte-range per-protein bundle slices (never full screen archives), resolves
UniProt synonyms, computes cLIP contact-fingerprint clusters, and evaluates
residue/range/motif queries in Markdown or JSON.

Examples:
  python3 scripts/atlas_query.py species
  python3 scripts/atlas_query.py resolve human UBC6e
  python3 scripts/atlas_query.py protein human UBE2J1 --res 171-184 --markdown
  python3 scripts/atlas_query.py protein human CXCR2 --res 10-19,188-202 --json
  python3 scripts/atlas_query.py pair human CXCR2 CXCL8 --res-a 10-19 --markdown
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import pathlib
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from typing import Any

DEFAULT_BASE = "https://flyark.github.io/livia-atlas/"
CUT = {10: 0.223, 5: 0.339, 1: 0.551}
KIND_DROP = {"phosphosite", "mutant", "variant", "construct", "engineered"}
UA = "LIVIA-Atlas-Query/1.1 (Python-urllib)"


class AtlasSource:
    """Reads JSON, text, and byte-range zip slices from a remote URL or local directory."""

    def __init__(self, base: str = DEFAULT_BASE, timeout: float = 35.0):
        self.timeout = timeout
        self.is_local = not re.match(r"^https?://", base, re.I)
        if self.is_local:
            self.local_root = pathlib.Path(base).expanduser().resolve()
            self.remote_root = DEFAULT_BASE
        else:
            self.local_root = None
            self.remote_root = base if base.endswith("/") else base + "/"
        self._json_cache: dict[str, Any] = {}

    def _http_bytes(self, url: str, byte_range: tuple[int, int] | None = None) -> bytes:
        headers = {"User-Agent": UA}
        if byte_range is not None:
            s, length = byte_range
            headers["Range"] = f"bytes={s}-{s + length - 1}"
        req = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(req, timeout=self.timeout) as resp:
            status = getattr(resp, "status", 200)
            if byte_range is not None and status != 206:
                raise RuntimeError(f"Expected HTTP 206 for Range read on {url}, got {status}")
            return resp.read()

    def read_bytes(self, rel_or_url: str, base_rel: str = "", byte_range: tuple[int, int] | None = None) -> bytes:
        if re.match(r"^https?://", rel_or_url, re.I):
            return self._http_bytes(rel_or_url, byte_range=byte_range)
        joined = urllib.parse.urljoin(base_rel, rel_or_url) if base_rel else rel_or_url
        clean = joined.split("?", 1)[0].lstrip("/")
        if self.is_local and self.local_root is not None:
            fp = (self.local_root / clean).resolve()
            if fp.exists() and fp.is_file():
                if byte_range is None:
                    return fp.read_bytes()
                with fp.open("rb") as f:
                    f.seek(byte_range[0])
                    return f.read(byte_range[1])
        url = urllib.parse.urljoin(self.remote_root, joined)
        return self._http_bytes(url, byte_range=byte_range)

    def read_json(self, rel_or_url: str, base_rel: str = "") -> Any:
        key = f"{base_rel}|{rel_or_url}"
        if key not in self._json_cache:
            raw = self.read_bytes(rel_or_url, base_rel=base_rel)
            self._json_cache[key] = json.loads(raw.decode("utf-8"))
        return self._json_cache[key]

    def read_text(self, rel_or_url: str, base_rel: str = "") -> str:
        return self.read_bytes(rel_or_url, base_rel=base_rel).decode("utf-8", errors="replace")


def expand_indices(val: str | None) -> list[int]:
    if not val:
        return []
    s = re.sub(r"[\[\]\"'\s]", "", str(val))
    if not s:
        return []
    out: list[int] = []
    for tok in re.split(r"[;,]", s):
        if not tok:
            continue
        parts = tok.split("-")
        try:
            a = int(parts[0])
            b = int(parts[1]) if len(parts) > 1 else a
        except ValueError:
            continue
        if b >= a and a >= 1:
            out.extend(range(a, b + 1))
        elif a >= 1:
            out.append(a)
    return out


def runs_of(nums: list[int] | set[int]) -> list[tuple[int, int]]:
    s = sorted({int(x) for x in nums if int(x) >= 1})
    if not s:
        return []
    out: list[tuple[int, int]] = []
    a = p = s[0]
    for x in s[1:]:
        if x == p + 1:
            p = x
        else:
            out.append((a, p))
            a = p = x
    out.append((a, p))
    return out


def stretches_of(rs: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for a, b in rs:
        if out and a - out[-1][1] <= 4:
            out[-1] = (out[-1][0], b)
        else:
            out.append((a, b))
    return out


def fmt_runs(rs: list[tuple[int, int]]) -> str:
    return ";".join(str(a) if a == b else f"{a}-{b}" for a, b in rs)


def seq_snippets(seq: str, rs: list[tuple[int, int]], limit: int = 6) -> str:
    if not seq or not rs:
        return ""
    st = stretches_of(rs)
    parts: list[str] = []
    for a, b in st[:limit]:
        sub = seq[max(0, a - 1) : min(len(seq), b)]
        if not sub:
            continue
        disp = sub if len(sub) <= 28 else f"{sub[:12]}...{sub[-8:]}"
        parts.append(f"{a if a == b else f'{a}-{b}'}:{disp}")
    if len(st) > limit:
        parts.append(f"+{len(st) - limit} more")
    return "; ".join(parts)


def parse_fasta(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    cur_id: str | None = None
    buf: list[str] = []
    for line in (text or "").splitlines():
        if line.startswith(">"):
            if cur_id is not None:
                out[cur_id] = "".join(buf)
            cur_id = line[1:].strip()
            buf = []
        elif line.strip():
            buf.append(line.strip())
    if cur_id is not None:
        out[cur_id] = "".join(buf)
    return out


def load_species(src: AtlasSource, sp_id: str) -> dict[str, Any]:
    reg = src.read_json("datasets.json")
    sp_reg = next((s for s in reg.get("species", []) if s["id"] == sp_id), None)
    if not sp_reg:
        known = ", ".join(s["id"] for s in reg.get("species", [])[:15])
        raise ValueError(f"Unknown species '{sp_id}'. Examples: {known}")
    base = sp_reg["base"]
    manifest = src.read_json("manifest.json", base_rel=base)
    prot = src.read_json(sp_reg.get("index", "proteins.json"), base_rel=base)
    cmp_fmt = prot.get("format") == "livia-atlas-index/1"
    cols = {k: i for i, k in enumerate(prot["columns"])}
    ce = len(prot["columns"]) - 1 if cmp_fmt else -1
    rows: list[dict[str, Any]] = []
    for idx, raw_r in enumerate(prot["rows"]):
        r = list(raw_r)
        if cmp_fmt:
            if r[cols["acc"]] == 0:
                r[cols["acc"]] = r[cols["key"]]
            if r[cols["clen"]] == 0:
                r[cols["clen"]] = r[cols["len"]]
        occ: list[dict[str, Any]] = []
        for tok in str(r[cols["occ"]] or "").split():
            if ":" not in tok:
                dot = tok.endswith(".")
                di = int(tok[:-1] if dot else tok)
                occ.append({"di": di, "name": r[ce] if dot else r[cols["key"]]})
            else:
                di_s, nm = tok.split(":", 1)
                occ.append({"di": int(di_s), "name": nm})
        gene_raw = r[cols["gene"]]
        gene = (gene_raw if gene_raw and gene_raw != "?" else "") or r[cols["acc"]] or r[cols["key"]]
        rows.append(
            {
                "i": idx,
                "key": r[cols["key"]],
                "gene": gene,
                "acc": r[cols["acc"]] or "",
                "name": r[cols["name"]] or "",
                "syn": r[cols["syn"]] or "",
                "len": r[cols["len"]] or 0,
                "clen": r[cols["clen"]] or 0,
                "status": r[cols["status"]] or "",
                "occ": occ,
                "id": occ[0]["name"] if occ else r[cols["key"]],
                "partners": r[cols["partners"]] or 0,
                "pos10": r[cols["pos10"]] or 0,
                "pos5": r[cols["pos5"]] or 0,
                "pos1": r[cols["pos1"]] or 0,
                "best": r[cols["bestIlis"]],
            }
        )
    by_key = {r["key"]: r for r in rows}
    by_name: dict[str, dict[str, Any]] = {}
    by_gene: dict[str, dict[str, Any]] = {}
    for r in rows:
        for o in r["occ"]:
            by_name[o["name"]] = r
        by_gene.setdefault(r["gene"], r)
    ds_ids = [d["id"] for d in manifest.get("datasets", [])]
    ds_short = [d["short"] for d in manifest.get("datasets", [])]
    ds_regs = {d["id"]: d for d in reg.get("datasets", [])}
    return {
        "id": sp_id,
        "reg": sp_reg,
        "base": base,
        "manifest": manifest,
        "rows": rows,
        "by_key": by_key,
        "by_name": by_name,
        "by_gene": by_gene,
        "ds_ids": ds_ids,
        "ds_short": ds_short,
        "ds_regs": ds_regs,
    }


def uniprot_aliases(raw: str, taxon: int | None) -> list[dict[str, str]]:
    t = re.sub(r"[^\w\s.-]", " ", str(raw or "")).strip()
    if len(t) < 2 or len(t) > 80:
        return []
    q_esc = t.replace('"', "")
    alnum = re.sub(r"[^A-Za-z0-9]", "", t)
    clauses = [f'gene:"{q_esc}"', f'protein_name:"{q_esc}"']
    if len(alnum) >= 3 and alnum.lower() != q_esc.lower():
        clauses.append(f'gene:"{alnum}"')
        clauses.append(f'protein_name:"{alnum}"')
    if len(alnum) >= 4 and " " not in t:
        clauses.append(f"protein_name:*{alnum}*")
    field_q = "(" + " OR ".join(clauses) + ")"
    tax_part = f" AND (organism_id:{taxon})" if taxon else ""

    def run_search(extra: str) -> list[dict[str, Any]]:
        q = f"{field_q}{tax_part}{extra}"
        params = urllib.parse.urlencode(
            {"query": q, "fields": "accession,gene_names,organism_id,protein_name", "size": "10", "format": "json"}
        )
        url = f"https://rest.uniprot.org/uniprotkb/search?{params}"
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        try:
            with urllib.request.urlopen(req, timeout=8.0) as resp:
                return json.loads(resp.read().decode("utf-8")).get("results", [])
        except Exception:
            return []

    hits = run_search(" AND (reviewed:true)") or run_search("")
    out: list[dict[str, str]] = []
    for e in hits:
        acc = e.get("primaryAccession") or ""
        g0 = (e.get("genes") or [{}])[0]
        gene = (g0.get("geneName") or {}).get("value") or ""
        syns = [s.get("value", "") for s in (g0.get("synonyms") or []) if s.get("value")]
        if acc:
            out.append({"acc": acc, "gene": gene, "syns": syns})
    return out


def resolve_protein(sp: dict[str, Any], raw: str, allow_uniprot: bool = True) -> tuple[dict[str, Any] | None, str]:
    q = str(raw or "").strip()
    if not q:
        return None, ""
    if q in sp["by_key"]:
        return sp["by_key"][q], "key"
    if q in sp["by_gene"]:
        return sp["by_gene"][q], "gene"
    if q in sp["by_name"]:
        return sp["by_name"][q], "screen_name"
    base_iso = re.sub(r"-\d+$", "", q)
    if base_iso != q and base_iso in sp["by_key"]:
        return sp["by_key"][base_iso], "isoform_accession"
    ql = q.lower()
    for r in sp["rows"]:
        if r["gene"].lower() == ql or r["key"].lower() == ql:
            return r, "gene_case_insensitive"
        if any(tok.lower() == ql for tok in str(r["acc"]).split()):
            return r, "accession"
        if any(o["name"].lower() == ql for o in r["occ"]):
            return r, "screen_name"
    for r in sp["rows"]:
        if any(s.lower() == ql for s in str(r["syn"]).split()):
            return r, "synonym"
    if allow_uniprot:
        taxon = (sp["manifest"].get("species") or {}).get("taxon") or sp["reg"].get("taxon")
        for u in uniprot_aliases(q, taxon):
            for cand in [u["acc"], u["gene"], *u["syns"]]:
                if cand:
                    hit, _ = resolve_protein(sp, cand, allow_uniprot=False)
                    if hit is not None:
                        return hit, f"UniProt synonym ({q} -> {hit['gene']} / {hit['key']})"
    return None, ""


def offset_for_bundle(src: AtlasSource, ds_manifest: dict[str, Any], offsets_path: str, rel: str) -> list[int] | None:
    if (ds_manifest.get("files") or {}).get("offsetShards"):
        nm0 = re.sub(r"\.[^.]+$", "", rel.rsplit("/", 1)[-1])
        m = re.match(r"^([A-Z0-9]+)_[A-Z0-9]+$", nm0)
        nm = m.group(1) if m else nm0
        xx = re.sub(r"[^a-z0-9_-]", "_", nm[-2:].lower() if len(nm) >= 2 else "_")
        shard_path = re.sub(r"\.json$", f"/{xx}.json", offsets_path)
        table = src.read_json(shard_path)
    else:
        table = src.read_json(offsets_path)
    return table.get(rel)


def read_bundle_zip_bytes(src: AtlasSource, ds_reg: dict[str, Any], ds_manifest: dict[str, Any], name: str) -> bytes:
    rel = (ds_manifest.get("files") or {}).get("bundle", "b/{id}.zip").replace("{id}", name)
    z = ds_reg.get("zip")
    if not z:
        base_rel = ds_manifest.get("bundleBase") or ds_reg["base"]
        return src.read_bytes(rel, base_rel=base_rel)
    at = offset_for_bundle(src, ds_manifest, z["offsets"], rel)
    if at:
        try:
            return src.read_bytes(z["url"], byte_range=(int(at[0]), int(at[1])))
        except Exception:
            fb = z.get("fallback") or ({"url": z["zenodo"]} if z.get("zenodo") else None)
            if fb:
                at_fb = offset_for_bundle(src, ds_manifest, fb.get("offsets") or z["offsets"], rel)
                if at_fb:
                    return src.read_bytes(fb["url"], byte_range=(int(at_fb[0]), int(at_fb[1])))
            raise
    raise FileNotFoundError(f"No bundle offset for {rel} in {ds_reg['id']}")


def parse_cons_tsv(text: str, seqs: dict[str, str]) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for line in text.splitlines()[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        while len(parts) < 10:
            parts.append("")
        name, key, kind, label, ln, off, exact, mut, seg, seq = parts[:10]
        seg_parsed = [[int(x) for x in t.split(":")] for t in seg.split(",")] if seg else None
        out[name] = {
            "key": key,
            "kind": kind,
            "label": label,
            "len": int(ln or 0),
            "off": None if off == "" else int(off),
            "exact": exact == "y",
            "mut": mut,
            "seg": seg_parsed,
            "seq": seq or seqs.get(name, ""),
        }
    return out


def parse_bundle_bytes(raw_bytes: bytes, bundle_name: str) -> dict[str, Any]:
    with zipfile.ZipFile(io.BytesIO(raw_bytes)) as zf:
        names = zf.namelist()
        csv_name = next((f for f in names if f.lower().endswith(".csv") and "identity" not in f.lower()), None)
        fa_name = next((f for f in names if re.search(r"\.(fa|fasta)$", f, re.I)), None)
        con_name = next((f for f in names if f.endswith("constructs.tsv")), None)
        iso_name = "isoforms.json" if "isoforms.json" in names else None
        seqs = parse_fasta(zf.read(fa_name).decode("utf-8", errors="replace")) if fa_name else {}
        cons = parse_cons_tsv(zf.read(con_name).decode("utf-8", errors="replace"), seqs) if con_name else None
        if cons:
            for nm in list(seqs.keys()):
                c = cons.get(nm)
                if c and not c["seq"]:
                    c["seq"] = seqs[nm]
                if nm != bundle_name and (not c or c["key"] != bundle_name):
                    del seqs[nm]
        rows: list[list[str]] = []
        header: list[str] = []
        if csv_name:
            txt = zf.read(csv_name).decode("utf-8", errors="replace")
            reader = list(csv.reader(io.StringIO(txt)))
            if reader:
                header = reader[0]
                last_k = None
                for r in reader[1:]:
                    if len(r) <= 10:
                        continue
                    k = "\x01".join(r)
                    if k == last_k:
                        continue
                    last_k = k
                    rows.append(r)
        isoforms = json.loads(zf.read(iso_name).decode("utf-8")) if iso_name else None
        return {"header": header, "H": {k: i for i, k in enumerate(header)}, "rows": rows, "seqs": seqs, "cons": cons, "isoforms": isoforms}


def to_float(row: list[str], H: dict[str, int], col: str) -> float:
    idx = H.get(col)
    if idx is None or idx >= len(row) or row[idx] == "":
        return math.nan
    try:
        return float(row[idx])
    except ValueError:
        return math.nan


def assemble_protein(src: AtlasSource, sp: dict[str, Any], P: dict[str, Any], whole_isoforms: bool = False) -> dict[str, Any]:
    parts: list[dict[str, Any]] = []
    for o in P["occ"]:
        ds_id = sp["ds_ids"][o["di"]]
        ds_reg = sp["ds_regs"][ds_id]
        ds_man = src.read_json("manifest.json", base_rel=ds_reg["base"])
        try:
            b_bytes = read_bundle_zip_bytes(src, ds_reg, ds_man, o["name"])
            raw = parse_bundle_bytes(b_bytes, o["name"])
            parts.append({"di": o["di"], "name": o["name"], "ds_id": ds_id, "ds_reg": ds_reg, "ds_man": ds_man, "raw": raw})
            if whole_isoforms and raw.get("isoforms"):
                for ch in raw["isoforms"].get("choices", []):
                    sub_b = read_bundle_zip_bytes(src, ds_reg, ds_man, ch["file"])
                    parts.append({"di": o["di"], "name": ch["file"], "ds_id": ds_id, "ds_reg": ds_reg, "ds_man": ds_man, "raw": parse_bundle_bytes(sub_b, ch["file"])})
        except Exception:
            continue
    if not parts:
        raise RuntimeError(f"Could not read interaction bundles for {P['gene']} ({P['key']}).")

    preds: list[dict[str, Any]] = []
    cons: dict[str, dict[str, Any]] = {}
    seqs: dict[str, str] = {}
    runs_map: dict[str, dict[str, Any]] = {}
    seq_no: dict[str, int] = {}

    def sno(s: str) -> int:
        if s not in seq_no:
            seq_no[s] = len(seq_no)
        return seq_no[s]

    def seq_id(C: dict[str, Any] | None, nm: str, ln: float) -> str:
        c = C.get(nm) if C else None
        if c and c["exact"] and c["off"] is not None and not c["mut"] and c["kind"] not in ("mutant", "variant"):
            return f"{c['key']}@{c['off']}+{c['len'] or int(ln or 0)}"
        return f"{nm}:{int(ln or 0)}"

    for part in parts:
        C = part["raw"]["cons"]
        if C:
            cons.update(C)

        def key_of(nm: str) -> str:
            if C and nm in C:
                return C[nm]["key"]
            r = sp["by_name"].get(nm)
            return r["key"] if r else nm

        for nm, sq in part["raw"]["seqs"].items():
            seqs[nm if C else key_of(nm)] = sq
        H = part["raw"]["H"]
        for r in part["raw"]["rows"]:
            nm = r[H["name"]]
            if "___" not in nm:
                continue
            a, b = nm.split("___", 1)
            qi = (key_of(a) == P["key"]) if C else (a == part["name"])
            if not qi and ((key_of(b) != P["key"]) if C else (b != part["name"])):
                continue
            qc = a if qi else b
            pc = b if qi else a
            pkey = key_of(pc)
            if part["ds_id"].startswith("afdb-het-") and pkey not in sp["by_key"]:
                continue
            batch_s = r[H["batch"]] if "batch" in H and H["batch"] < len(r) else ""
            rid = f"{part['di']}|{nm}" + (f"|{batch_s}" if "batch" in H else "")
            if rid not in runs_map:
                runs_map[rid] = {"id": rid, "di": part["di"], "ds": part["ds_id"], "short": sp["ds_short"][part["di"]], "name": nm, "qi": qi, "key": pkey, "qc": qc, "pc": pc}
            pick = lambda x, y: x if qi else y
            q_len = to_float(r, H, pick("len_i", "len_j"))
            p_len = to_float(r, H, pick("len_j", "len_i"))
            ilis = to_float(r, H, "iLIS")
            ilia = to_float(r, H, "iLIA")
            ilisa = to_float(r, H, "iLISA")
            if math.isnan(ilisa) and not math.isnan(ilis):
                ilisa = ilis * (0.0 if math.isnan(ilia) else ilia)
            sx = part["raw"]["seqs"].get(qc)
            sy = part["raw"]["seqs"].get(pc)
            sp_pair = f"{sno(sx)}|{sno(sy)}" if sx and sy else f"{seq_id(C, qc, q_len)}|{seq_id(C, pc, p_len)}"
            preds.append(
                {
                    "partner": pkey,
                    "run": rid,
                    "di": part["di"],
                    "qi": qi,
                    "qc": qc,
                    "pc": pc,
                    "rank": to_float(r, H, "rank"),
                    "model": to_float(r, H, "model"),
                    "iLIS": ilis,
                    "iLIA": ilia,
                    "iLISA": ilisa,
                    "ipTM": to_float(r, H, "ipTM"),
                    "pTM": to_float(r, H, "pTM"),
                    "LIS": to_float(r, H, "LIS"),
                    "cLIS": to_float(r, H, "cLIS"),
                    "ipSAE": to_float(r, H, "ipSAE"),
                    "actifpTM": to_float(r, H, "actifpTM"),
                    "pDockQ": to_float(r, H, "pDockQ"),
                    "pDockQ2": to_float(r, H, "pDockQ2"),
                    "qLen": int(q_len) if not math.isnan(q_len) else 0,
                    "pLen": int(p_len) if not math.isnan(p_len) else 0,
                    "qLIR": int(to_float(r, H, pick("LIR_i", "LIR_j")) or 0) if not math.isnan(to_float(r, H, pick("LIR_i", "LIR_j"))) else 0,
                    "pLIR": int(to_float(r, H, pick("LIR_j", "LIR_i")) or 0) if not math.isnan(to_float(r, H, pick("LIR_j", "LIR_i"))) else 0,
                    "qcLIR": int(to_float(r, H, pick("cLIR_i", "cLIR_j")) or 0) if not math.isnan(to_float(r, H, pick("cLIR_i", "cLIR_j"))) else 0,
                    "pcLIR": int(to_float(r, H, pick("cLIR_j", "cLIR_i")) or 0) if not math.isnan(to_float(r, H, pick("cLIR_j", "cLIR_i"))) else 0,
                    "qL": r[H[pick("LIR_indices_i", "LIR_indices_j")]] if pick("LIR_indices_i", "LIR_indices_j") in H else "",
                    "pL": r[H[pick("LIR_indices_j", "LIR_indices_i")]] if pick("LIR_indices_j", "LIR_indices_i") in H else "",
                    "qC": r[H[pick("cLIR_indices_i", "cLIR_indices_j")]] if pick("cLIR_indices_i", "cLIR_indices_j") in H else "",
                    "pC": r[H[pick("cLIR_indices_j", "cLIR_indices_i")]] if pick("cLIR_indices_j", "cLIR_indices_i") in H else "",
                    "sp": sp_pair,
                }
            )

    dropped: list[dict[str, Any]] = []
    if cons:
        kept_preds: list[dict[str, Any]] = []
        for p in preds:
            if (p["qc"] in cons and cons[p["qc"]]["kind"] in KIND_DROP) or (p["pc"] in cons and cons[p["pc"]]["kind"] in KIND_DROP):
                dropped.append(p)
            else:
                kept_preds.append(p)
        preds = kept_preds

    top_run: dict[str, tuple[float, str]] = {}
    for p in preds:
        v = 0.0 if math.isnan(p["iLIS"]) else p["iLIS"]
        if p["run"] not in top_run or v > top_run[p["run"]][0]:
            top_run[p["run"]] = (v, p["sp"])

    def better_run(rx_id: str, ry_id: str) -> bool:
        vx, _ = top_run[rx_id]
        vy, _ = top_run[ry_id]
        if vx != vy:
            return vx > vy
        rx, ry = runs_map[rx_id], runs_map[ry_id]
        if rx["di"] != ry["di"]:
            return rx["di"] < ry["di"]
        if rx["qi"] != ry["qi"]:
            return rx["qi"]
        return rx_id < ry_id

    kept_by_sp: dict[str, str] = {}
    for rid, (_, sp_pair) in top_run.items():
        if sp_pair not in kept_by_sp or better_run(rid, kept_by_sp[sp_pair]):
            kept_by_sp[sp_pair] = rid
    for p in preds:
        p["rep"] = kept_by_sp.get(top_run[p["run"]][1]) != p["run"]
    for p in dropped:
        p["rep"] = False

    def aggregate(plist: list[dict[str, Any]]) -> list[dict[str, Any]]:
        by_p: dict[str, list[dict[str, Any]]] = {}
        for p in plist:
            by_p.setdefault(p["partner"], []).append(p)
        out_p: list[dict[str, Any]] = []
        for key, ps in by_p.items():
            cs = [p for p in ps if not p["rep"]] or ps
            rep_only = not any(not p["rep"] for p in ps)
            bm = max(cs, key=lambda x: -1.0 if math.isnan(x["iLIS"]) else x["iLIS"])
            best_ilis = max((0.0 if math.isnan(x["iLIS"]) else x["iLIS"]) for x in cs)
            best_ilisa = max((0.0 if math.isnan(x["iLISA"]) else x["iLISA"]) for x in cs)
            iptms = [x["ipTM"] for x in cs if not math.isnan(x["ipTM"])]
            best_iptm = max(iptms) if iptms else math.nan
            out_p.append(
                {
                    "id": key,
                    "row": sp["by_key"].get(key),
                    "preds": ps,
                    "counted": cs,
                    "rep": rep_only,
                    "best": best_ilis,
                    "ilisaBest": best_ilisa,
                    "iptmBest": best_iptm,
                    "bm": bm,
                    "pass10": sum(1 for x in cs if not math.isnan(x["iLIS"]) and x["iLIS"] >= CUT[10]),
                }
            )
        return out_p

    q_seq = seqs.get(P["key"], "")
    if not q_seq:
        for o in P["occ"]:
            if o["name"] in seqs:
                q_seq = seqs[o["name"]]
                break
            ds_id = sp["ds_ids"][o["di"]]
            ds_reg = sp["ds_regs"][ds_id]
            ds_man = src.read_json("manifest.json", base_rel=ds_reg["base"])
            seq_tmpl = (ds_man.get("files") or {}).get("sequence")
            if seq_tmpl:
                try:
                    rel_fa = seq_tmpl.replace("{id}", o["name"])
                    z = ds_reg.get("zip")
                    if z:
                        at = offset_for_bundle(src, ds_man, z["offsets"], rel_fa)
                        fa_txt = src.read_bytes(z["url"], byte_range=(int(at[0]), int(at[1]))).decode("utf-8", errors="replace") if at else ""
                    else:
                        fa_txt = src.read_text(rel_fa, base_rel=ds_reg["base"])
                    vals = list(parse_fasta(fa_txt).values())
                    if vals:
                        q_seq = vals[0]
                        break
                except Exception:
                    pass

    return {
        "preds": preds,
        "partners": aggregate(preds),
        "dropped": dropped,
        "droppedPartners": aggregate(dropped) if dropped else [],
        "runs": runs_map,
        "seqs": seqs,
        "cons": cons,
        "qSeq": q_seq,
    }


def to_ref_coord(c: dict[str, Any] | None, r: int) -> int | None:
    if not c:
        return r
    if c["exact"] and c["off"] is not None:
        return r + c["off"]
    if c["seg"]:
        for a, b, n in c["seg"]:
            if a <= r < a + n:
                return b + r - a
    return None


def map_pred_residues(p: dict[str, Any], cons: dict[str, Any], P: dict[str, Any], side: str = "q") -> tuple[list[int], list[int]]:
    cname = p["qc"] if side == "q" else p["pc"]
    clen = p["qLen"] if side == "q" else p["pLen"]
    l_raw = expand_indices(p["qL"] if side == "q" else p["pL"])
    c_raw = expand_indices(p["qC"] if side == "q" else p["pC"])
    c = cons.get(cname) if cons else None
    ref_len = P.get("len") or 0
    if c and c["len"] == clen and c["kind"] not in ("mutant", "variant"):
        if c["exact"] and c["off"] is not None and (not ref_len or c["off"] + clen <= ref_len):
            return [r + c["off"] for r in l_raw], [r + c["off"] for r in c_raw]
        if c["seg"]:
            lm = sorted({x for x in (to_ref_coord(c, r) for r in l_raw) if x is not None})
            cm = sorted({x for x in (to_ref_coord(c, r) for r in c_raw) if x is not None})
            return lm, cm
    return l_raw, c_raw


def cluster_fingerprints(fps: list[list[int]], max_k: int = 30, drop_threshold: float = 0.05) -> tuple[int, list[int]]:
    n = len(fps)
    if n == 0:
        return 0, []
    if n == 1:
        return 1, [1]
    sets = [set(f) for f in fps]
    norms = [math.sqrt(len(s)) for s in sets]
    D = [[0.0] * n for _ in range(n)]
    for i in range(n):
        si, ni = sets[i], norms[i]
        for j in range(i + 1, n):
            inter = len(si.intersection(sets[j]))
            d = 1.0 if (ni == 0.0 or norms[j] == 0.0) else max(0.0, 1.0 - inter / (ni * norms[j]))
            D[i][j] = D[j][i] = d

    size = {i: 1 for i in range(n)}
    dist = {i: {j: D[i][j] for j in range(n) if j != i} for i in range(n)}
    active = list(range(n))
    Z: list[tuple[int, int, float, int]] = []
    next_id = n
    while len(active) > 1:
        active.sort()
        best = math.inf
        bi = bj = -1
        for x in range(len(active)):
            i = active[x]
            di = dist[i]
            for y in range(x + 1, len(active)):
                j = active[y]
                d = di[j]
                if d < best - 1e-12:
                    best, bi, bj = d, i, j
        ni, nj = size[bi], size[bj]
        nid = next_id
        next_id += 1
        size[nid] = ni + nj
        dist[nid] = {}
        active.remove(bi)
        active.remove(bj)
        for k in active:
            dn = (ni * dist[bi][k] + nj * dist[bj][k]) / (ni + nj)
            dist[nid][k] = dn
            dist[k][nid] = dn
        active.append(nid)
        Z.append((bi, bj, best, ni + nj))

    def fcluster(k_val: int) -> list[int]:
        uf: dict[int, int] = {}

        def find(x: int) -> int:
            while x in uf:
                x = uf[x]
            return x

        for m in range(n - k_val):
            a, b, _, _ = Z[m]
            ra, rb = find(a), find(b)
            if ra != rb:
                uf[ra] = rb
            uf[n + m] = find(b)
        root_lab: dict[int, int] = {}
        nxt = 1
        labs: list[int] = []
        for i in range(n):
            r = find(i)
            if r not in root_lab:
                root_lab[r] = nxt
                nxt += 1
            labs.append(root_lab[r])
        return labs

    def silhouette(labs: list[int]) -> float:
        groups: dict[int, list[int]] = {}
        for idx, lb in enumerate(labs):
            groups.setdefault(lb, []).append(idx)
        if len(groups) < 2:
            return -1.0
        tot = 0.0
        for i in range(n):
            own = groups[labs[i]]
            if len(own) == 1:
                continue
            a = sum(D[i][j] for j in own if j != i) / (len(own) - 1)
            b = min(sum(D[i][j] for j in grp) / len(grp) for cid, grp in groups.items() if cid != labs[i])
            denom = max(a, b)
            tot += 0.0 if denom == 0.0 else (b - a) / denom
        return tot / n

    k_cap = min(max_k, n - 1)
    if k_cap < 2:
        return 2 if n >= 2 else 1, fcluster(min(2, n))
    k_scores: list[tuple[int, float]] = []
    best_s = -1.0
    opt_k = 2
    for kv in range(2, k_cap + 1):
        lb = fcluster(kv)
        if len(set(lb)) < 2:
            continue
        sc = silhouette(lb)
        k_scores.append((kv, sc))
        if sc > best_s:
            best_s = sc
            opt_k = kv
    for idx in range(1, len(k_scores)):
        if k_scores[idx - 1][1] - k_scores[idx][1] > drop_threshold:
            opt_k = k_scores[idx - 1][0]
            break
    raw_labels = fcluster(opt_k)
    counts: dict[int, int] = {}
    for lb in raw_labels:
        counts[lb] = counts.get(lb, 0) + 1
    ordered_old = sorted(counts.keys(), key=lambda c: (-counts[c], c))
    remap = {old: idx + 1 for idx, old in enumerate(ordered_old)}
    return len(ordered_old), [remap[lb] for lb in raw_labels]


def fetch_domains(src: AtlasSource, acc: str) -> list[dict[str, Any]]:
    if not acc:
        return []
    primary = acc.split()[0].split("-")[0]
    if len(primary) < 2:
        return []
    try:
        shard = src.read_json(f"data/dom/{primary[-2:].lower()}.json")
        entry = shard.get(primary)
        if not entry:
            return []
        nums = [int(entry[i]) for i in range(3, len(entry), 3)]
        name_tables = {k: src.read_json(f"data/dom/names/{k}.json") for k in {n % 128 for n in nums}}
        doms: list[dict[str, Any]] = []
        for i in range(1, len(entry) - 2, 3):
            st, en, pf_n = int(entry[i]), int(entry[i + 1]), int(entry[i + 2])
            pf_id = f"PF{pf_n:05d}"
            nm = (name_tables.get(pf_n % 128) or {}).get(str(pf_n)) or pf_id
            doms.append({"name": nm, "start": st, "end": en, "pfam": pf_id, "source": "Pfam"})
        return sorted(doms, key=lambda d: d["start"])
    except Exception:
        return []


def fetch_afdb_plddt(acc: str, seq_len: int) -> list[float] | None:
    if not acc:
        return None
    primary = acc.split()[0].split("-")[0]
    try:
        req = urllib.request.Request(f"https://alphafold.ebi.ac.uk/api/prediction/{urllib.parse.quote(primary)}", headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=8.0) as resp:
            items = json.loads(resp.read().decode("utf-8"))
        if isinstance(items, dict):
            items = [items]
        entry = next((x for x in items if x.get("uniprotAccession") == primary and (x.get("uniprotStart") or 1) == 1), items[0] if items else None)
        if not entry:
            return None
        if entry.get("plddtDocUrl"):
            r2 = urllib.request.Request(entry["plddtDocUrl"], headers={"User-Agent": UA})
            with urllib.request.urlopen(r2, timeout=8.0) as resp2:
                doc = json.loads(resp2.read().decode("utf-8"))
            scores = doc.get("confidenceScore")
            if isinstance(scores, list) and (not seq_len or len(scores) == seq_len):
                return [float(x) for x in scores]
        cif_url = entry.get("cifUrl")
        if cif_url:
            r3 = urllib.request.Request(cif_url, headers={"User-Agent": UA})
            with urllib.request.urlopen(r3, timeout=10.0) as resp3:
                cif_txt = resp3.read().decode("utf-8", errors="replace")
            by_res: dict[int, float] = {}
            for line in cif_txt.splitlines():
                if line.startswith("ATOM "):
                    toks = line.split()
                    if len(toks) >= 15 and toks[3] == "CA":
                        try:
                            by_res[int(toks[8])] = float(toks[14])
                        except ValueError:
                            pass
            if by_res and (not seq_len or max(by_res) == seq_len):
                L = seq_len or max(by_res)
                return [by_res.get(i, math.nan) for i in range(1, L + 1)]
    except Exception:
        return None
    return None


def split_res_parts(v: str) -> list[str]:
    out: list[str] = []
    cur: list[str] = []
    depth = 0
    for ch in str(v or ""):
        if ch in "({[":
            depth += 1
            cur.append(ch)
        elif ch in ")}]":
            depth = max(0, depth - 1)
            cur.append(ch)
        elif ch in ",;" and depth == 0:
            s = "".join(cur).strip()
            if s:
                out.append(s)
            cur = []
        else:
            cur.append(ch)
    tail = "".join(cur).strip()
    if tail:
        out.append(tail)
    return out


def parse_motif_regex(pat_raw: str) -> re.Pattern[str] | None:
    s = re.sub(r"\s+", "", pat_raw.upper())
    prosite = "-" in s
    s = s.replace("-", "")
    if prosite and s.endswith("."):
        s = s[:-1]
    aa = set("ACDEFGHIKLMNPQRSTVWYUO")
    amb = {"B": "DN", "Z": "EQ", "J": "IL"}
    out = ""
    i = 0
    fixed = False
    min_len = 0
    while i < len(s):
        c = s[i]
        if c in "^<" and i == 0:
            out += "^"
            i += 1
            continue
        if c in "$>" and i == len(s) - 1:
            out += "$"
            i += 1
            continue
        if c in "X.":
            el = "."
            i += 1
        elif c in amb:
            el = f"[{amb[c]}]"
            fixed = True
            i += 1
        elif c in aa:
            el = c
            fixed = True
            i += 1
        elif c in "[{":
            end_c = "]" if c == "[" else "}"
            j = s.find(end_c, i + 1)
            if j < 0:
                return None
            body = s[i + 1 : j]
            neg = c == "{"
            if body.startswith("^"):
                neg = True
                body = body[1:]
            if not body or re.search(r"[^A-Z]|X", body):
                return None
            expanded = "".join(amb.get(ch, ch) for ch in body)
            el = f"[{'^' if neg else ''}{expanded}]"
            fixed = True
            i = j + 1
        else:
            return None
        m = re.match(r"^[({](\d+)(?:,(\d+))?[)}]", s[i:])
        if m:
            a = int(m.group(1))
            b = int(m.group(2)) if m.group(2) is not None else a
            if b < a or b > 99:
                return None
            out += f"{el}{{{a},{b}}}"
            min_len += a
            i += len(m.group(0))
        else:
            out += el
            min_len += 1
    return re.compile(out) if fixed and min_len >= 2 else None


def resolve_res_spans(query: str, seq: str, length: int) -> list[dict[str, Any]]:
    spans: list[dict[str, Any]] = []
    for part in split_res_parts(query):
        t = re.sub(r"^P\.(?=[A-Z]{0,3}\d)", "", part.strip().upper())
        mr = re.match(r"^[A-Z]{0,3}(\d+)\s*(?:[-–—:]|\.\.)\s*[A-Z]{0,3}(\d+)$", t)
        if mr:
            a, b = sorted((int(mr.group(1)), int(mr.group(2))))
            if a >= 1 and a <= length:
                e = min(length, b)
                spans.append({"start": a, "end": e, "mismatches": 0, "sequence": seq[a - 1 : e] if seq else ""})
            continue
        mn = re.match(r"^[A-Z]{0,3}(\d+)[A-Z*]{0,3}$", t)
        if mn:
            n = int(mn.group(1))
            if 1 <= n <= length:
                spans.append({"start": n, "end": n, "mismatches": 0, "sequence": seq[n - 1 : n] if seq else ""})
            continue
        if not seq:
            continue
        if re.search(r"[\[\]{}()^$<>.]|[XBZJ]", t):
            rx = parse_motif_regex(t)
            if rx:
                for i in range(len(seq)):
                    m = rx.match(seq, i)
                    if m and m.end() > i:
                        st, en = i + 1, m.end()
                        if spans and st <= spans[-1]["end"]:
                            spans[-1]["end"] = max(spans[-1]["end"], en)
                            spans[-1]["sequence"] = seq[spans[-1]["start"] - 1 : spans[-1]["end"]]
                        else:
                            spans.append({"start": st, "end": en, "mismatches": 0, "sequence": seq[st - 1 : en]})
            continue
        q_aa = re.sub(r"[^A-Z]", "", t)
        if len(q_aa) >= 3 and len(q_aa) <= len(seq):
            allow = len(q_aa) // 5
            best: tuple[int, int, int] | None = None
            for i in range(len(seq) - len(q_aa) + 1):
                mm = sum(1 for j, ch in enumerate(q_aa) if ch != "X" and ch != seq[i + j])
                if mm <= allow and (best is None or mm < best[2]):
                    best = (i + 1, i + len(q_aa), mm)
                    if mm == 0:
                        break
            if best:
                spans.append({"start": best[0], "end": best[1], "mismatches": best[2], "sequence": seq[best[0] - 1 : best[1]]})
    spans.sort(key=lambda x: (x["start"], x["end"]))
    return spans


def run_protein_command(src: AtlasSource, args: argparse.Namespace) -> dict[str, Any]:
    sp = load_species(src, args.species)
    P, via = resolve_protein(sp, args.query)
    if not P:
        return {"ok": False, "error": f"No protein matching '{args.query}' in {sp['reg']['label']}."}
    cut_fpr = args.cut if args.cut in CUT else 10
    cut_ilis = CUT[cut_fpr]
    B = assemble_protein(src, sp, P, whole_isoforms=False)
    q_seq = B["qSeq"]
    q_len = P["len"] or P["clen"] or len(q_seq) or max((p["qLen"] for p in B["preds"]), default=0)

    others = [pt for pt in B["partners"] if pt["id"] != P["key"] and not pt["rep"]]
    counts = {
        "folded": len(others),
        "fpr10": sum(1 for x in others if x["best"] >= CUT[10]),
        "fpr5": sum(1 for x in others if x["best"] >= CUT[5]),
        "fpr1": sum(1 for x in others if x["best"] >= CUT[1]),
        "activeCutoffFpr": cut_fpr,
        "activeCutoffILIS": cut_ilis,
    }

    passing_models: list[dict[str, Any]] = []
    fps: list[list[int]] = []
    for p in B["preds"]:
        if p["rep"] or math.isnan(p["iLIS"]) or p["iLIS"] < cut_ilis:
            continue
        q_lir, q_clir = map_pred_residues(p, B["cons"], P, side="q")
        q_clir = [r for r in q_clir if 1 <= r <= q_len]
        q_lir = [r for r in q_lir if 1 <= r <= q_len]
        if q_clir:
            passing_models.append({**p, "_qCLIR": q_clir, "_qLIR": q_lir})
            fps.append(q_clir)

    if len(fps) > 4000:
        by_pair: dict[str, tuple[dict[str, Any], list[int]]] = {}
        for pm, fp in zip(passing_models, fps):
            cur = by_pair.get(pm["partner"])
            if cur is None or pm["iLIS"] > cur[0]["iLIS"]:
                by_pair[pm["partner"]] = (pm, fp)
        top_items = sorted(by_pair.values(), key=lambda x: -x[0]["iLIS"])[:4000]
        passing_models = [x[0] for x in top_items]
        fps = [x[1] for x in top_items]

    k_opt, labels = cluster_fingerprints(fps) if len(fps) >= 2 else (len(fps), [1] * len(fps))
    partner_cluster: dict[str, tuple[float, int]] = {}
    by_cl: dict[int, dict[str, Any]] = {}
    for pm, fp, lb in zip(passing_models, fps, labels):
        pk = pm["partner"]
        if pk != P["key"] and (pk not in partner_cluster or pm["iLIS"] > partner_cluster[pk][0]):
            partner_cluster[pk] = (pm["iLIS"], lb)
        g = by_cl.setdefault(lb, {"c": lb, "n": 0, "hits": {}, "best": {}})
        g["n"] += 1
        for r in fp:
            g["hits"][r] = g["hits"].get(r, 0) + 1
        if pk != P["key"] and pm["iLIS"] > g["best"].get(pk, -1.0):
            g["best"][pk] = pm["iLIS"]

    doms = fetch_domains(src, P["acc"]) if q_len and (not P["len"] or P["len"] == q_len) else []
    plddt_arr = None if args.no_afdb else fetch_afdb_plddt(P["acc"], q_len)
    plddt_vals = [v for v in (plddt_arr or []) if not math.isnan(v)]
    plddt_summary = (
        {
            "mapped": True,
            "mean": round(sum(plddt_vals) / len(plddt_vals), 1),
            "min": round(min(plddt_vals), 1),
            "max": round(max(plddt_vals), 1),
        }
        if plddt_vals
        else {"mapped": False}
    )

    def gname(k: str) -> str:
        r = sp["by_key"].get(k)
        return r["gene"] if r else k

    clusters_out: list[dict[str, Any]] = []
    for g in sorted(by_cl.values(), key=lambda x: (-len(x["best"]), -x["n"], x["c"])):
        foot = sorted(r for r, h in g["hits"].items() if h / g["n"] >= 0.3)
        st = stretches_of(runs_of(foot))
        c_doms = list(dict.fromkeys(d["name"] for d in doms if any(d["start"] <= r <= d["end"] for r in foot)))
        clusters_out.append(
            {
                "cluster": g["c"],
                "label": f"C{g['c']}",
                "models": g["n"],
                "partnerCount": len(g["best"]),
                "footprintStretches": [str(a) if a == b else f"{a}-{b}" for a, b in st],
                "footprintSequence": seq_snippets(q_seq, st, limit=5),
                "domains": c_doms[:3],
                "partners": [{"key": k, "gene": gname(k), "bestILIS": round(v, 3)} for k, v in sorted(g["best"].items(), key=lambda kv: -kv[1])],
            }
        )

    region_eval: dict[str, Any] | None = None
    if args.res:
        spans = resolve_res_spans(args.res, q_seq, q_len)
        if not spans:
            region_eval = {"ok": False, "query": args.res, "reason": "no_match", "length": q_len}
        else:
            res_set = {r for spn in spans for r in range(spn["start"], spn["end"] + 1)}
            dom_hits = [d for d in doms if any(d["start"] <= r <= d["end"] for r in res_set)]
            if dom_hits:
                dom_ctx = "in " + ", ".join(f"{d['name']} ({d['start']}-{d['end']})" for d in dom_hits)
            elif doms and len(spans) == 1:
                s0, e0 = spans[0]["start"], spans[0]["end"]
                if e0 < doms[0]["start"]:
                    dom_ctx = f"outside annotated domains (N-terminal to {doms[0]['name']} {doms[0]['start']}-{doms[0]['end']})"
                elif s0 > doms[-1]["end"]:
                    dom_ctx = f"outside annotated domains (C-terminal to {doms[-1]['name']} {doms[-1]['start']}-{doms[-1]['end']})"
                else:
                    dom_ctx = "outside annotated domains"
            else:
                dom_ctx = ""
            reg_pvs = [plddt_arr[r - 1] for r in sorted(res_set) if plddt_arr and 1 <= r <= len(plddt_arr) and not math.isnan(plddt_arr[r - 1])]
            reg_plddt = (
                {"mean": round(sum(reg_pvs) / len(reg_pvs), 1), "min": round(min(reg_pvs), 1), "max": round(max(reg_pvs), 1), "count": len(reg_pvs)}
                if reg_pvs
                else None
            )
            clir_models = 0
            clir_map: dict[str, dict[str, Any]] = {}
            for pm, fp in zip(passing_models, fps):
                hit_r = [r for r in fp if r in res_set]
                if not hit_r:
                    continue
                clir_models += 1
                pk = pm["partner"]
                if pk == P["key"]:
                    continue
                cur = clir_map.setdefault(pk, {"key": pk, "gene": gname(pk), "best": -1.0, "res": set()})
                if pm["iLIS"] > cur["best"]:
                    cur["best"] = pm["iLIS"]
                cur["res"].update(hit_r)
            lir_models = 0
            lir_map: dict[str, dict[str, Any]] = {}
            for p in B["preds"]:
                if p["rep"] or math.isnan(p["iLIS"]) or p["iLIS"] < cut_ilis:
                    continue
                q_lir, _ = map_pred_residues(p, B["cons"], P, side="q")
                hit_l = [r for r in q_lir if r in res_set]
                if not hit_l:
                    continue
                lir_models += 1
                pk = p["partner"]
                if pk == P["key"]:
                    continue
                cur = lir_map.setdefault(pk, {"key": pk, "gene": gname(pk), "best": -1.0, "res": set()})
                if p["iLIS"] > cur["best"]:
                    cur["best"] = p["iLIS"]
                cur["res"].update(hit_l)
            clir_partners = [
                {"key": x["key"], "gene": x["gene"], "bestILIS": round(x["best"], 3), "contactedResidues": sorted(x["res"])}
                for x in sorted(clir_map.values(), key=lambda z: -z["best"])
            ]
            lir_only_partners = [
                {"key": x["key"], "gene": x["gene"], "bestILIS": round(x["best"], 3), "interactionResidues": sorted(x["res"])}
                for x in sorted(lir_map.values(), key=lambda z: -z["best"])
                if x["key"] not in clir_map
            ]
            region_eval = {
                "ok": True,
                "query": args.res,
                "cutoffFpr": cut_fpr,
                "cutoffILIS": cut_ilis,
                "spans": spans,
                "domainContext": dom_ctx,
                "plddt": reg_plddt,
                "totalPassingModels": len(passing_models),
                "clirModels": clir_models,
                "clirPartners": clir_partners,
                "lirModels": lir_models,
                "lirOnlyPartners": lir_only_partners,
            }

    top_partners: list[dict[str, Any]] = []
    passing_others = [pt for pt in others if pt["best"] >= cut_ilis]
    for pt in sorted(passing_others, key=lambda x: -x["best"])[: args.top]:
        bm = pt["bm"]
        q_l, q_c = map_pred_residues(bm, B["cons"], P, side="q")
        p_row = pt["row"] or {}
        p_l, p_c = map_pred_residues(bm, B["cons"], p_row, side="p")
        cl_id = partner_cluster.get(pt["id"], (0.0, 0))[1]
        top_partners.append(
            {
                "key": pt["id"],
                "gene": p_row.get("gene") or pt["id"],
                "acc": p_row.get("acc") or "",
                "cluster": f"C{cl_id}" if cl_id else None,
                "bestILIS": round(pt["best"], 3),
                "iLISA": None if math.isnan(pt["ilisaBest"]) else round(pt["ilisaBest"], 1),
                "ipTM": None if math.isnan(pt["iptmBest"]) else round(pt["iptmBest"], 2),
                "queryCLIRRanges": fmt_runs(runs_of(q_c)),
                "queryLIRRanges": fmt_runs(runs_of(q_l)),
                "partnerCLIRRanges": fmt_runs(runs_of(p_c)),
                "partnerLIRRanges": fmt_runs(runs_of(p_l)),
                "modelsPast10": pt["pass10"],
                "modelsTotal": len(pt["counted"]),
                "pairUrl": f"https://flyark.github.io/livia-atlas/#/{sp['id']}/{P['key']}/{pt['id']}",
            }
        )

    url_qs = []
    if args.res:
        url_qs.append(f"res={urllib.parse.quote(args.res)}")
    if cut_fpr != 10:
        url_qs.append(f"cut={cut_fpr}")
    atlas_url = f"https://flyark.github.io/livia-atlas/#/{sp['id']}/{P['key']}" + (f"?{'&'.join(url_qs)}" if url_qs else "")

    return {
        "ok": True,
        "view": "protein",
        "resolvedVia": via,
        "url": atlas_url,
        "species": {"id": sp["id"], "label": sp["reg"]["label"], "name": sp["reg"]["name"]},
        "protein": {"key": P["key"], "gene": P["gene"], "acc": P["acc"], "name": P["name"], "length": q_len, "sequence": q_seq},
        "counts": counts,
        "domains": doms,
        "plddt": plddt_summary,
        "activeRegionQuery": region_eval,
        "clusters": clusters_out,
        "partners": top_partners,
    }


def run_pair_command(src: AtlasSource, args: argparse.Namespace) -> dict[str, Any]:
    sp = load_species(src, args.species)
    P, via_a = resolve_protein(sp, args.query_a)
    if not P:
        return {"ok": False, "error": f"No protein matching '{args.query_a}' in {sp['reg']['label']}."}
    O0, via_b = resolve_protein(sp, args.query_b)
    o_key = O0["key"] if O0 else args.query_b
    O = O0 or {"key": o_key, "gene": o_key, "acc": "", "name": "", "len": 0, "clen": 0, "occ": []}

    B = assemble_protein(src, sp, P, whole_isoforms=True)
    part0 = next((pt for pt in B["partners"] if pt["id"] == O["key"]), None)
    dpart = next((pt for pt in B["droppedPartners"] if pt["id"] == O["key"]), None)
    part = part0 or dpart
    if not part:
        return {"ok": False, "error": f"{P['gene']} and {O['gene']} were not predicted together in {sp['reg']['label']} screens."}

    seq_a = B["qSeq"]
    seq_b = B["seqs"].get(O["key"], "")
    if not seq_b and O0:
        try:
            B2 = assemble_protein(src, sp, O0, whole_isoforms=False)
            seq_b = B2["qSeq"]
        except Exception:
            seq_b = ""

    sorted_preds = sorted(part["preds"], key=lambda x: -1.0 if math.isnan(x["iLIS"]) else -x["iLIS"])
    models_out: list[dict[str, Any]] = []
    for idx, p in enumerate(sorted_preds, 1):
        q_l, q_c = map_pred_residues(p, B["cons"], P, side="q")
        p_l, p_c = map_pred_residues(p, B["cons"], O, side="p")
        q_cr, p_cr = runs_of(q_c), runs_of(p_c)
        models_out.append(
            {
                "rank": idx,
                "source": B["runs"][p["run"]]["short"],
                "predictorRank": None if math.isnan(p["rank"]) else int(p["rank"]),
                "iLIS": None if math.isnan(p["iLIS"]) else round(p["iLIS"], 3),
                "iLISA": None if math.isnan(p["iLISA"]) else round(p["iLISA"], 1),
                "ipTM": None if math.isnan(p["ipTM"]) else round(p["ipTM"], 2),
                "LIS": None if math.isnan(p["LIS"]) else round(p["LIS"], 3),
                "cLIS": None if math.isnan(p["cLIS"]) else round(p["cLIS"], 3),
                "ipSAE": None if math.isnan(p["ipSAE"]) else round(p["ipSAE"], 3),
                "actifpTM": None if math.isnan(p["actifpTM"]) else round(p["actifpTM"], 3),
                "pDockQ2": None if math.isnan(p["pDockQ2"]) else round(p["pDockQ2"], 3),
                "queryCLIRRanges": fmt_runs(q_cr),
                "queryLIRRanges": fmt_runs(runs_of(q_l)),
                "partnerCLIRRanges": fmt_runs(p_cr),
                "partnerLIRRanges": fmt_runs(runs_of(p_l)),
                "queryCLIRSequence": seq_snippets(seq_a, q_cr, limit=6),
                "partnerCLIRSequence": seq_snippets(seq_b, p_cr, limit=6),
                "_qC": set(q_c),
                "_qL": set(q_l),
                "_pC": set(p_c),
                "_pL": set(p_l),
            }
        )

    sel_idx = max(0, min(len(models_out) - 1, (args.rank or 1) - 1))
    sel_model = {k: v for k, v in models_out[sel_idx].items() if not k.startswith("_")}
    len_a = P["clen"] or P["len"] or len(seq_a) or (sorted_preds[0]["qLen"] if sorted_preds else 0)
    len_b = O["clen"] or O["len"] or len(seq_b) or (sorted_preds[0]["pLen"] if sorted_preds else 0)

    def eval_pair_query(q_str: str | None, side: str) -> dict[str, Any] | None:
        if not q_str:
            return None
        sq = seq_a if side == "A" else seq_b
        ln = len_a if side == "A" else len_b
        gene = P["gene"] if side == "A" else O["gene"]
        spans = resolve_res_spans(q_str, sq, ln)
        if not spans:
            return {"ok": False, "gene": gene, "query": q_str, "reason": "no_match"}
        rset = {r for s in spans for r in range(s["start"], s["end"] + 1)}
        ck = "_qC" if side == "A" else "_pC"
        lk = "_qL" if side == "A" else "_pL"
        sel_m = models_out[sel_idx]
        return {
            "ok": True,
            "gene": gene,
            "query": q_str,
            "spans": spans,
            "selectedModelCLIR": sorted(sel_m[ck].intersection(rset)),
            "selectedModelLIR": sorted(sel_m[lk].intersection(rset)),
            "clirModels": sum(1 for m in models_out if m[ck].intersection(rset)),
            "lirModels": sum(1 for m in models_out if m[lk].intersection(rset)),
            "totalModels": len(models_out),
        }

    reg_a = eval_pair_query(args.res_a, "A")
    reg_b = eval_pair_query(args.res_b, "B")
    clean_models = [{k: v for k, v in m.items() if not k.startswith("_")} for m in models_out]

    qs = []
    if args.rank and args.rank > 1:
        qs.append(f"rank={args.rank}")
    if args.res_a:
        qs.append(f"resA={urllib.parse.quote(args.res_a)}")
    if args.res_b:
        qs.append(f"resB={urllib.parse.quote(args.res_b)}")
    pair_url = f"https://flyark.github.io/livia-atlas/#/{sp['id']}/{P['key']}/{O['key']}" + (f"?{'&'.join(qs)}" if qs else "")

    return {
        "ok": True,
        "view": "pair",
        "resolvedVia": {"A": via_a, "B": via_b or "raw"},
        "url": pair_url,
        "species": {"id": sp["id"], "label": sp["reg"]["label"], "name": sp["reg"]["name"]},
        "proteinA": {"key": P["key"], "gene": P["gene"], "acc": P["acc"], "length": len_a, "sequence": seq_a},
        "proteinB": {"key": O["key"], "gene": O["gene"], "acc": O["acc"], "length": len_b, "sequence": seq_b},
        "bestILIS": round(part["best"], 3),
        "bestIpTM": None if math.isnan(part["iptmBest"]) else round(part["iptmBest"], 2),
        "selectedRank": sel_model["rank"],
        "selectedModel": sel_model,
        "activeRegionQueries": {"resA": reg_a, "resB": reg_b} if (reg_a or reg_b) else None,
        "models": clean_models,
    }


def format_markdown(d: dict[str, Any]) -> str:
    if not d.get("ok"):
        return f"# LIVIA Atlas Error\n\n- {d.get('error', 'Unknown error')}"
    view = d.get("view")
    lines: list[str] = []
    if view == "protein":
        p, sp, c = d["protein"], d["species"], d["counts"]
        lines.append(f"# LIVIA Atlas Protein: {p['gene']}" + (f" ({p['acc']})" if p["acc"] else "") + f" — {sp['label']}")
        if d.get("resolvedVia") and "synonym" in d["resolvedVia"].lower():
            lines.append(f"- **Resolved via**: {d['resolvedVia']}")
        if p.get("name"):
            lines.append(f"- **Protein**: {p['name']} · **Length**: {p['length']} aa")
        lines.append(
            f"- **Partners**: {c['folded']} folded · **≥ 10% FPR** (iLIS ≥ {CUT[10]}): {c['fpr10']} · "
            f"**≥ 5% FPR** (iLIS ≥ {CUT[5]}): {c['fpr5']} · **≥ 1% FPR** (iLIS ≥ {CUT[1]}): {c['fpr1']}"
        )
        if d.get("domains"):
            lines.append("- **Domains**: " + ", ".join(f"{x['name']} ({x['start']}-{x['end']})" for x in d["domains"]))
        if d.get("plddt", {}).get("mapped"):
            pl = d["plddt"]
            lines.append(f"- **AlphaFold DB pLDDT**: mean {pl['mean']:.1f} (min {pl['min']:.1f}–max {pl['max']:.1f})")
        lines.append(f"- **Atlas URL**: {d['url']}")
        rq = d.get("activeRegionQuery")
        if rq:
            lines.append(f"\n## Queried Region (`{rq['query']}` at {c['activeCutoffFpr']}% FPR, iLIS ≥ {c['activeCutoffILIS']})")
            if not rq.get("ok"):
                lines.append(f"- **Result**: {rq.get('reason', 'no_match')}")
            else:
                sp_txt = ", ".join(
                    f"{s['start'] if s['start'] == s['end'] else f'{s['start']}-{s['end']}'}"
                    + (f" ({s['sequence']})" if s.get("sequence") and len(s["sequence"]) <= 32 else "")
                    for s in rq["spans"]
                )
                meta = [sp_txt]
                if rq.get("domainContext"):
                    meta.append(rq["domainContext"])
                if rq.get("plddt"):
                    meta.append(f"AlphaFold DB pLDDT mean {rq['plddt']['mean']:.1f} ({rq['plddt']['min']:.1f}–{rq['plddt']['max']:.1f})")
                lines.append("- **Spans**: " + " · ".join(meta))
                c_p = ", ".join(f"{x['gene']} ({x['bestILIS']:.3f})" for x in rq["clirPartners"][:12])
                lines.append(
                    f"- **Contact models (cLIR, PAE ≤ 12 Å & Cβ ≤ 8 Å)**: {rq['clirModels']} of {rq['totalPassingModels']} passing models "
                    f"({len(rq['clirPartners'])} partners)" + (f": {c_p}" if c_p else "")
                )
                if rq["lirOnlyPartners"]:
                    l_p = ", ".join(f"{x['gene']} ({x['bestILIS']:.3f})" for x in rq["lirOnlyPartners"][:10])
                    lines.append(f"- **Local interaction residues only (LIR, PAE ≤ 12 Å)**: {len(rq['lirOnlyPartners'])} partners: {l_p}")
        if d.get("clusters"):
            lines.append(f"\n## Predicted Binding-Site Clusters ({c['activeCutoffFpr']}% FPR, footprint = ≥ 30% of cluster models)")
            for cl in d["clusters"]:
                fp = ", ".join(cl["footprintStretches"]) if cl["footprintStretches"] else "none ≥ 30%"
                tops = ", ".join(f"{x['gene']} ({x['bestILIS']:.3f})" for x in cl["partners"][:8])
                dom_s = f" ({'; '.join(cl['domains'])})" if cl.get("domains") else ""
                seq_s = f" (`{cl['footprintSequence']}`)" if cl.get("footprintSequence") else ""
                lines.append(f"- **{cl['label']}** ({cl['partnerCount']} partners, {cl['models']} models): footprint `{fp}`{seq_s}{dom_s}" + (f" · top partners: {tops}" if tops else ""))
        if d.get("partners"):
            pts = d["partners"]
            lines.append(f"\n## Top Predicted Interactors (showing {len(pts)})")
            lines.append(f"| Partner | UniProt | Cluster | Best iLIS | iLISA | ipTM | {p['gene']} cLIR | Partner cLIR | Models ≥ 10% FPR |")
            lines.append("|---|---|---|---:|---:|---:|---|---|---|")
            for pt in pts:
                ilis_s = f"{pt['bestILIS']:.3f}"
                ilisa_s = f"{pt['iLISA']:.1f}" if pt["iLISA"] is not None else "–"
                iptm_s = f"{pt['ipTM']:.2f}" if pt["ipTM"] is not None else "–"
                lines.append(
                    f"| {pt['gene']} | {pt['acc'] or '–'} | {pt['cluster'] or '–'} | {ilis_s} | {ilisa_s} | {iptm_s} | "
                    f"{pt['queryCLIRRanges'] or 'none'} | {pt['partnerCLIRRanges'] or 'none'} | {pt['modelsPast10']}/{pt['modelsTotal']} |"
                )
        return "\n".join(lines)

    if view == "pair":
        a, b, sp = d["proteinA"], d["proteinB"], d["species"]
        lines.append(f"# LIVIA Atlas Pair: {a['gene']} × {b['gene']} — {sp['label']}")
        lines.append(f"- **{a['gene']}**: {a['acc'] or a['key']} ({a['length']} aa) · **{b['gene']}**: {b['acc'] or b['key']} ({b['length']} aa)")
        iptm_b = f"{d['bestIpTM']:.2f}" if d["bestIpTM"] is not None else "–"
        lines.append(f"- **Best iLIS**: {d['bestILIS']:.3f} · **Best ipTM**: {iptm_b} · **Selected model**: rank {d['selectedRank']}")
        lines.append(f"- **Atlas URL**: {d['url']}")
        sm = d["selectedModel"]
        lines.append(f"\n## Selected Model Residue Spans (Rank {sm['rank']}: {sm['source']})")
        lines.append(
            f"- **{a['gene']} cLIR** (PAE ≤ 12 Å & Cβ ≤ 8 Å): `{sm['queryCLIRRanges'] or 'none'}`"
            + (f" — `{sm['queryCLIRSequence']}`" if sm.get("queryCLIRSequence") else "")
        )
        lines.append(f"- **{a['gene']} LIR** (PAE ≤ 12 Å): `{sm['queryLIRRanges'] or 'none'}`")
        lines.append(
            f"- **{b['gene']} cLIR** (PAE ≤ 12 Å & Cβ ≤ 8 Å): `{sm['partnerCLIRRanges'] or 'none'}`"
            + (f" — `{sm['partnerCLIRSequence']}`" if sm.get("partnerCLIRSequence") else "")
        )
        lines.append(f"- **{b['gene']} LIR** (PAE ≤ 12 Å): `{sm['partnerLIRRanges'] or 'none'}`")
        if d.get("activeRegionQueries"):
            lines.append("\n## Queried Regions")
            for k in ("resA", "resB"):
                rq = d["activeRegionQueries"].get(k)
                if not rq:
                    continue
                if not rq.get("ok"):
                    lines.append(f"- **{rq['gene']} `{rq['query']}`**: {rq.get('reason', 'no_match')}")
                else:
                    lines.append(
                        f"- **{rq['gene']} `{rq['query']}`**: selected model cLIR {len(rq['selectedModelCLIR'])} aa, "
                        f"LIR {len(rq['selectedModelLIR'])} aa · across {rq['totalModels']} models: "
                        f"{rq['clirModels']} with cLIR, {rq['lirModels']} with LIR"
                    )
        lines.append("\n## Ranked Models")
        lines.append(f"| Rank | Source | iLIS | iLISA | ipTM | LIS | cLIS | ipSAE | {a['gene']} cLIR | {b['gene']} cLIR |")
        lines.append("|---:|---|---:|---:|---:|---:|---:|---:|---|---|")
        for m in d["models"]:
            fn = lambda v, dec=3: f"{v:.{dec}f}" if v is not None else "–"
            lines.append(
                f"| {m['rank']} | {m['source']} | {fn(m['iLIS'], 3)} | {fn(m['iLISA'], 1)} | {fn(m['ipTM'], 2)} | "
                f"{fn(m['LIS'], 3)} | {fn(m['cLIS'], 3)} | {fn(m['ipSAE'], 3)} | "
                f"{m['queryCLIRRanges'] or 'none'} | {m['partnerCLIRRanges'] or 'none'} |"
            )
        return "\n".join(lines)
    return json.dumps(d, indent=2)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Zero-dependency CLI client for LIVIA Atlas.")
    ap.add_argument("--base", default=DEFAULT_BASE, help="Atlas base URL or local repo path (default: %(default)s)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    sp_cmd = sub.add_parser("species", help="List available species in LIVIA Atlas")
    sp_cmd.add_argument("--json", action="store_true", help="Output JSON")

    rs_cmd = sub.add_parser("resolve", help="Resolve a gene symbol, accession, or synonym in a species")
    rs_cmd.add_argument("species", help="Species ID (e.g., human, fly, mouse / mus-musculus, yeast, worm)")
    rs_cmd.add_argument("query", help="Gene symbol, UniProt accession, FlyBase ID, or alias")
    rs_cmd.add_argument("--json", action="store_true", help="Output JSON")

    pr_cmd = sub.add_parser("protein", help="Fetch a protein's cLIP binding sites, interactors, and optional residue lookup")
    pr_cmd.add_argument("species", help="Species ID (e.g., human, fly, zebrafish, worm, yeast)")
    pr_cmd.add_argument("query", help="Gene symbol, UniProt accession, FlyBase ID, or alias")
    pr_cmd.add_argument("--res", default="", help="Residue number, range (10-19), comma list (10-19,188-202), sequence, or motif")
    pr_cmd.add_argument("--cut", type=int, choices=[10, 5, 1], default=10, help="FPR cutoff percentage (10, 5, or 1)")
    pr_cmd.add_argument("--top", type=int, default=25, help="Max top interactors to return (default: 25)")
    pr_cmd.add_argument("--no-afdb", action="store_true", help="Skip live AlphaFold DB pLDDT lookup")
    pr_cmd.add_argument("--json", action="store_true", help="Output structured JSON")
    pr_cmd.add_argument("--markdown", action="store_true", help="Output Markdown (default)")

    pa_cmd = sub.add_parser("pair", help="Fetch all ranked models and exact residue spans for a protein pair")
    pa_cmd.add_argument("species", help="Species ID (e.g., human, fly)")
    pa_cmd.add_argument("query_a", help="First protein (gene, accession, or alias)")
    pa_cmd.add_argument("query_b", help="Second protein (gene, accession, or alias)")
    pa_cmd.add_argument("--rank", type=int, default=1, help="1-based ranked model to highlight (default: 1)")
    pa_cmd.add_argument("--res-a", default="", help="Residue/range/motif query on protein A")
    pa_cmd.add_argument("--res-b", default="", help="Residue/range/motif query on protein B")
    pa_cmd.add_argument("--json", action="store_true", help="Output structured JSON")
    pa_cmd.add_argument("--markdown", action="store_true", help="Output Markdown (default)")

    args = ap.parse_args(argv)
    alias_sp = {"mouse": "mus-musculus", "rat": "rattus-norvegicus", "arabidopsis": "arabidopsis-thaliana", "ecoli": "escherichia-coli-83333"}
    if hasattr(args, "species") and args.species.lower() in alias_sp:
        args.species = alias_sp[args.species.lower()]

    src = AtlasSource(base=args.base)
    try:
        if args.cmd == "species":
            reg = src.read_json("datasets.json")
            items = [
                {
                    "id": s["id"],
                    "label": s["label"],
                    "name": s.get("name", ""),
                    "taxon": s.get("taxon"),
                    "proteins": (s.get("counts") or {}).get("proteins"),
                    "pairsFpr10": (s.get("counts") or {}).get("pairsFpr10"),
                }
                for s in reg.get("species", [])
            ]
            if args.json:
                print(json.dumps(items, indent=2))
            else:
                print("# LIVIA Atlas Species\n\n| ID | Label | Scientific Name | Taxon | Proteins | Pairs ≥ 10% FPR |")
                print("|---|---|---|---:|---:|---:|")
                for s in items:
                    print(f"| {s['id']} | {s['label']} | {s['name']} | {s['taxon'] or '–'} | {s['proteins'] or '–'} | {s['pairsFpr10'] or '–'} |")
            return 0

        if args.cmd == "resolve":
            sp = load_species(src, args.species)
            row, via = resolve_protein(sp, args.query)
            if not row:
                out = {"ok": False, "species": sp["id"], "query": args.query}
                if args.json:
                    print(json.dumps(out, indent=2))
                else:
                    print(f"No match for '{args.query}' in {sp['reg']['label']}.")
                return 1
            out = {
                "ok": True,
                "species": sp["id"],
                "query": args.query,
                "via": via,
                "key": row["key"],
                "gene": row["gene"],
                "acc": row["acc"],
                "name": row["name"],
                "length": row["clen"] or row["len"],
                "partnersFolded": row["partners"],
                "partnersFpr10": row["pos10"],
                "url": f"https://flyark.github.io/livia-atlas/#/{sp['id']}/{row['key']}",
            }
            if args.json:
                print(json.dumps(out, indent=2))
            else:
                print(f"- **{out['gene']}** (`{out['key']}`, UniProt `{out['acc'] or '–'}`) — {out['name']}")
                print(f"- Resolved via `{via}` · {out['length']} aa · {out['partnersFpr10']} partners past 10% FPR (of {out['partnersFolded']} folded)")
                print(f"- URL: {out['url']}")
            return 0

        if args.cmd == "protein":
            d = run_protein_command(src, args)
            if not args.json:
                for p_obj in [d.get("protein")]:
                    if p_obj and "sequence" in p_obj and not args.res:
                        pass
            if args.json:
                print(json.dumps(d, indent=2))
            else:
                print(format_markdown(d))
            return 0 if d.get("ok") else 1

        if args.cmd == "pair":
            d = run_pair_command(src, args)
            if args.json:
                print(json.dumps(d, indent=2))
            else:
                print(format_markdown(d))
            return 0 if d.get("ok") else 1
    except Exception as exc:
        if getattr(args, "json", False):
            print(json.dumps({"ok": False, "error": str(exc)}, indent=2))
        else:
            print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
