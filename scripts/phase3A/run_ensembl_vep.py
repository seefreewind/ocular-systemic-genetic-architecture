#!/usr/bin/env python3
"""Annotate the 58 frozen lead variants with Ensembl GRCh37 VEP.

The output is an annotation lookup, not a causal-variant or causal-gene call.
"""
from __future__ import annotations

import csv
import json
import time
from datetime import date
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

ROOT = Path(__file__).resolve().parents[2]
FREEZE = ROOT / "config/phase3A_LOCUS_FREEZE.tsv"
OUT = ROOT / "results/phase3A/VEP_LEAD_VARIANT_ANNOTATION.tsv"
CACHE = ROOT / "results/phase3A/ensembl_vep_cache.json"
URL = "https://grch37.rest.ensembl.org/vep/human/id"


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def request_batch(ids):
    body = json.dumps({"ids": ids}).encode()
    request = Request(URL, data=body, headers={
        "Content-Type": "application/json",
        "Accept": "application/json",
        "User-Agent": "ocular-systemic-antagonistic-pleiotropy/phase3A",
    }, method="POST")
    with urlopen(request, timeout=90) as response:
        return json.loads(response.read().decode())


def flatten_annotation(item, meta):
    base = {
        "lead_SNP": item.get("id", item.get("input", "")),
        "input": item.get("input", ""),
        "assembly": "GRCh37",
        "vep_release": "116",
        "annotation_date": str(date.today()),
        "most_severe_consequence": item.get("most_severe_consequence", ""),
        "seq_region_name": item.get("seq_region_name", ""),
        "start": item.get("start", ""),
        "end": item.get("end", ""),
        "allele_string": item.get("allele_string", ""),
        "consequence": "",
        "impact": "",
        "gene_id": "",
        "gene_symbol": "",
        "transcript": "",
        "biotype": "",
        "protein_position": "",
        "amino_acids": "",
        "codons": "",
        "canonical": "",
        "distance": "",
        "regulatory_feature": "",
        "source": "Ensembl VEP GRCh37 REST",
        "status": "ANNOTATED",
    }
    transcripts = item.get("transcript_consequences") or []
    regulatory = item.get("regulatory_feature_consequences") or []
    intergenic = item.get("intergenic_consequences") or []
    rows = []
    for consequence in transcripts:
        row = dict(base)
        row.update({
            "consequence": ";".join(consequence.get("consequence_terms", []) or []),
            "impact": consequence.get("impact", ""),
            "gene_id": consequence.get("gene_id", ""),
            "gene_symbol": consequence.get("gene_symbol", ""),
            "transcript": consequence.get("transcript_id", ""),
            "biotype": consequence.get("biotype", ""),
            "protein_position": consequence.get("protein_start", consequence.get("protein_end", "")),
            "amino_acids": consequence.get("amino_acids", ""),
            "codons": consequence.get("codons", ""),
            "canonical": consequence.get("canonical", ""),
            "distance": ";".join(str(x) for x in consequence.get("distance", []) or []) if isinstance(consequence.get("distance"), list) else consequence.get("distance", ""),
            "regulatory_feature": ";".join(consequence.get("regulatory_feature", []) or []) if isinstance(consequence.get("regulatory_feature"), list) else consequence.get("regulatory_feature", ""),
        })
        rows.append(row)
    for consequence in regulatory:
        row = dict(base)
        row.update({
            "consequence": ";".join(consequence.get("consequence_terms", []) or []),
            "impact": consequence.get("impact", ""),
            "regulatory_feature": consequence.get("regulatory_feature_id", ""),
            "status": "REGULATORY_ANNOTATION",
        })
        rows.append(row)
    if not rows and intergenic:
        for consequence in intergenic:
            row = dict(base)
            row.update({
                "consequence": ";".join(consequence.get("consequence_terms", []) or ["intergenic_variant"]),
                "status": "INTERGENIC_ANNOTATION",
            })
            rows.append(row)
    if not rows:
        base["status"] = "NO_TRANSCRIPT_OR_REGULATORY_CONSEQUENCE_RETURNED"
        rows.append(base)
    return rows


def main():
    loci = read_tsv(FREEZE)
    assert len(loci) == 58, len(loci)
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    ids = [r["lead_SNP"] for r in loci]
    for start in range(0, len(ids), 50):
        batch = [x for x in ids[start:start + 50] if x not in cache]
        if not batch:
            continue
        try:
            result = request_batch(batch)
            for item in result:
                key = item.get("id", item.get("input", ""))
                cache[key] = item
        except (HTTPError, URLError, TimeoutError, OSError, json.JSONDecodeError) as error:
            for key in batch:
                cache[key] = {"id": key, "_error": f"{type(error).__name__}: {error}"}
        CACHE.write_text(json.dumps(cache, indent=2, sort_keys=True))
        time.sleep(0.2)

    fields = [
        "lead_SNP", "input", "assembly", "vep_release", "annotation_date", "most_severe_consequence",
        "seq_region_name", "start", "end", "allele_string", "consequence", "impact", "gene_id",
        "gene_symbol", "transcript", "biotype", "protein_position", "amino_acids", "codons",
        "canonical", "distance", "regulatory_feature", "source", "status",
    ]
    out = []
    for locus in loci:
        item = cache.get(locus["lead_SNP"], {"id": locus["lead_SNP"], "_error": "MISSING_RESPONSE"})
        if "_error" in item:
            out.append({
                "lead_SNP": locus["lead_SNP"], "input": locus["lead_SNP"], "assembly": "GRCh37",
                "vep_release": "116", "annotation_date": str(date.today()), "source": "Ensembl VEP GRCh37 REST",
                "status": "API_ERROR: " + item["_error"],
            })
        else:
            out.extend(flatten_annotation(item, locus))
    write_tsv(OUT, out, fields)
    print(f"VEP annotation complete: {len(loci)} lead variants, {len(out)} rows")


if __name__ == "__main__":
    main()
