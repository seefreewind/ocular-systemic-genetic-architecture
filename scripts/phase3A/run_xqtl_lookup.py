#!/usr/bin/env python3
"""Perform descriptive GTEx v8 variant-to-gene xQTL lookups.

The records are lookup support only.  They are not colocalization evidence and
are never used to redefine a locus or assign causality.
"""
from __future__ import annotations

import csv
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[2]
FREEZE = ROOT / "config/phase3A_LOCUS_FREEZE.tsv"
OUT = ROOT / "results/phase3A/XQTL_LOOKUP.tsv"
CACHE = ROOT / "results/phase3A/gtex_v8_lookup_cache.json"
BASE = "https://gtexportal.org/api/v2"
RELEVANT = ("Brain", "Artery", "Heart", "Kidney", "Liver", "Adipose", "Whole_Blood", "Nerve", "Muscle", "Spleen", "Pancreas", "Small_Intestine", "Skin")


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def get_json(url):
    req = Request(url, headers={"Accept": "application/json", "User-Agent": "ocular-systemic-antagonistic-pleiotropy/phase3A"})
    with urlopen(req, timeout=60) as response:
        return json.loads(response.read().decode())


def fetch_variant(rsid):
    url = BASE + "/dataset/variant?" + urlencode({"snpId": rsid, "format": "json"})
    try:
        data = get_json(url).get("data", [])
        v8 = [x for x in data if x.get("datasetId") == "gtex_v8"]
        return rsid, {"data": v8[:1], "status": "OK"}
    except (HTTPError, URLError, TimeoutError, OSError, json.JSONDecodeError) as error:
        return rsid, {"data": [], "status": f"API_ERROR: {type(error).__name__}: {error}"}


def fetch_eqtl(internal_id):
    query = urlencode({"variantId": internal_id, "datasetId": "gtex_v8", "format": "json"})
    try:
        data = get_json(BASE + "/association/singleTissueEqtl?" + query)
        return internal_id, {"data": data.get("data", []), "status": "OK"}
    except (HTTPError, URLError, TimeoutError, OSError, json.JSONDecodeError) as error:
        return internal_id, {"data": [], "status": f"API_ERROR: {type(error).__name__}: {error}"}


def main():
    loci = read_tsv(FREEZE)
    assert len(loci) == 58
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {"variant": {}, "eqtl": {}}
    cache.setdefault("variant", {})
    cache.setdefault("eqtl", {})
    rsids = sorted({l["lead_SNP"] for l in loci})

    missing = [x for x in rsids if x not in cache["variant"]]
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(fetch_variant, x) for x in missing]
        for future in as_completed(futures):
            key, value = future.result()
            cache["variant"][key] = value
    CACHE.write_text(json.dumps(cache, indent=2, sort_keys=True))

    internal_ids = sorted({
        item["data"][0]["variantId"]
        for item in cache["variant"].values()
        if item.get("data")
    })
    missing = [x for x in internal_ids if x not in cache["eqtl"]]
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(fetch_eqtl, x) for x in missing]
        for future in as_completed(futures):
            key, value = future.result()
            cache["eqtl"][key] = value
    CACHE.write_text(json.dumps(cache, indent=2, sort_keys=True))

    loci_by_rs = {l["lead_SNP"]: l for l in loci}
    fields = [
        "lead_SNP", "pair", "direction_class", "source", "dataset", "variant_id", "tissue",
        "gene_symbol", "gencode_id", "p_value", "nes", "status", "query_date", "build_note", "lookup_scope",
    ]
    rows = []
    for rsid in rsids:
        locus = loci_by_rs[rsid]
        variant = cache["variant"].get(rsid, {"data": [], "status": "MISSING_CACHE"})
        if not variant.get("data"):
            rows.append({
                "lead_SNP": rsid, "pair": locus["pair"], "direction_class": locus["direction_class"],
                "source": "GTEx_v8", "dataset": "gtex_v8", "status": variant.get("status", "NO_VARIANT_RECORD"),
                "query_date": str(date.today()), "build_note": "GTEx API variant identifier lookup", "lookup_scope": "relevant_tissues",
            })
            continue
        internal = variant["data"][0]["variantId"]
        eqtl = cache["eqtl"].get(internal, {"data": [], "status": "MISSING_CACHE"})
        relevant = [x for x in eqtl.get("data", []) if any(str(x.get("tissueSiteDetailId", "")).startswith(prefix) for prefix in RELEVANT)]
        if not relevant:
            rows.append({
                "lead_SNP": rsid, "pair": locus["pair"], "direction_class": locus["direction_class"],
                "source": "GTEx_v8", "dataset": "gtex_v8", "variant_id": internal,
                "status": "NO_RELEVANT_TISSUE_RECORD" if eqtl.get("status") == "OK" else eqtl.get("status", "NO_RECORD"),
                "query_date": str(date.today()), "build_note": "GTEx v8 API; relevant tissues filtered descriptively", "lookup_scope": "relevant_tissues",
            })
        else:
            for item in relevant:
                rows.append({
                    "lead_SNP": rsid, "pair": locus["pair"], "direction_class": locus["direction_class"],
                    "source": "GTEx_v8", "dataset": item.get("datasetId", "gtex_v8"), "variant_id": internal,
                    "tissue": item.get("tissueSiteDetailId", ""), "gene_symbol": item.get("geneSymbol", ""),
                    "gencode_id": item.get("gencodeId", ""), "p_value": item.get("pValue", ""), "nes": item.get("nes", ""),
                    "status": "LOOKUP_RETURNED", "query_date": str(date.today()),
                    "build_note": "GTEx v8 API; raw p-values are descriptive and not used as colocalization evidence",
                    "lookup_scope": "relevant_tissues",
                })
        # EyeGEx/retina is recorded as an explicit optional source status so
        # absence is not silently interpreted as a negative biological result.
        rows.append({
            "lead_SNP": rsid, "pair": locus["pair"], "direction_class": locus["direction_class"],
            "source": "EyeGEx/retina", "dataset": "", "variant_id": "",
            "status": "NOT_ASSESSED_STABLE_PUBLIC_API_NOT_CONFIGURED", "query_date": str(date.today()),
            "build_note": "Optional retina xQTL lookup; no stable endpoint was used in this reproducible workflow",
            "lookup_scope": "optional_source_status",
        })
    write_tsv(OUT, rows, fields)
    print(f"xQTL lookup complete: {len(rsids)} lead variants, {len(rows)} rows")


if __name__ == "__main__":
    main()
