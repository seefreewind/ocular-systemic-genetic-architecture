#!/usr/bin/env python3
"""Map all 2,096 valid Phase 2A-R lead loci to descriptive VEP genes.

The background is restricted to genes observed among the valid clumped loci;
it is never replaced by the full human gene universe.  Variants without a
protein-coding VEP transcript are retained as unmapped loci and are reported.
"""
from __future__ import annotations

import csv
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[2]
CLUMP_DIR = ROOT / "results/phase2A/ld_clump"
OUT = ROOT / "results/phase3A/FUNCTIONAL_BACKGROUND_GENES.tsv"
CACHE = ROOT / "results/phase3A/background_vep_cache.json"
URL = "https://grch37.rest.ensembl.org/vep/human/id"


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def request_batch(ids):
    body = json.dumps({"ids": ids}).encode()
    request = Request(URL, data=body, headers={
        "Content-Type": "application/json", "Accept": "application/json",
        "User-Agent": "ocular-systemic-antagonistic-pleiotropy/phase3A",
    }, method="POST")
    with urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode())


def main():
    records = []
    for path in sorted(CLUMP_DIR.glob("*.clumps")):
        if path.name.startswith("._"):
            continue
        pair = path.name.split(".chr", 1)[0]
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                records.append({
                    "pair": pair, "chr": row.get("#CHROM", ""), "lead_BP": row.get("POS", ""),
                    "lead_SNP": row.get("ID", ""), "PLACO_clump_P": row.get("P", ""), "source_clump_file": path.name,
                })
    assert len(records) == 2096, f"Expected 2096 valid clump lead loci, found {len(records)}"
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    ids = sorted({r["lead_SNP"] for r in records if r["lead_SNP"]})
    batches = [[x for x in ids[start:start + 50] if x not in cache] for start in range(0, len(ids), 50)]
    batches = [b for b in batches if b]
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(request_batch, batch): batch for batch in batches}
        for future in as_completed(futures):
            batch = futures[future]
            try:
                response = future.result()
                for item in response:
                    cache[item.get("id", item.get("input", ""))] = item
            except (HTTPError, URLError, TimeoutError, OSError, json.JSONDecodeError) as error:
                for key in batch:
                    cache[key] = {"id": key, "_error": f"{type(error).__name__}: {error}"}
            CACHE.write_text(json.dumps(cache, indent=2, sort_keys=True)
                            )

    fields = [
        "pair", "chr", "lead_BP", "lead_SNP", "PLACO_clump_P", "source_clump_file", "gene_symbol",
        "gene_id", "transcript", "biotype", "consequence", "mapping_level", "mapping_status",
        "annotation_source", "assembly", "vep_release", "annotation_date",
    ]
    out = []
    for record in records:
        item = cache.get(record["lead_SNP"], {"_error": "MISSING_RESPONSE"})
        transcripts = [] if "_error" in item else item.get("transcript_consequences", []) or []
        coding = []
        for tx in transcripts:
            if tx.get("biotype") != "protein_coding" or not tx.get("gene_symbol"):
                continue
            terms = tx.get("consequence_terms", []) or []
            high = any(x in terms for x in ("splice", "stop_gained", "frameshift", "missense", "inframe", "start_lost", "stop_lost", "transcript_ablation", "coding_sequence_variant"))
            coding.append((tx, high, terms))
        if coding:
            # Preserve all protein-coding transcript mappings but collapse
            # duplicate gene/transcript/consequence rows deterministically.
            seen = set()
            for tx, high, terms in coding:
                key = (tx.get("gene_symbol"), tx.get("gene_id"), tx.get("transcript_id"), ";".join(terms))
                if key in seen:
                    continue
                seen.add(key)
                out.append({
                    **record, "gene_symbol": tx.get("gene_symbol", ""), "gene_id": tx.get("gene_id", ""),
                    "transcript": tx.get("transcript_id", ""), "biotype": tx.get("biotype", ""),
                    "consequence": ";".join(terms), "mapping_level": "LEVEL1" if high else "LEVEL2",
                    "mapping_status": "VEP_PROTEIN_CODING_TRANSCRIPT", "annotation_source": "Ensembl VEP GRCh37 REST",
                    "assembly": "GRCh37", "vep_release": "116", "annotation_date": str(date.today()),
                })
        else:
            status = "UNMAPPED_NO_PROTEIN_CODING_VEP_TRANSCRIPT"
            if "_error" in item:
                status = "VEP_API_ERROR: " + item["_error"]
            out.append({
                **record, "mapping_level": "UNMAPPED", "mapping_status": status,
                "annotation_source": "Ensembl VEP GRCh37 REST", "assembly": "GRCh37", "vep_release": "116",
                "annotation_date": str(date.today()),
            })
    write_tsv(OUT, out, fields)
    mapped_loci = len({(r["pair"], r["chr"], r["lead_SNP"]) for r in out if r.get("gene_symbol")})
    print(f"Functional background complete: {len(records)} loci; {mapped_loci} loci with protein-coding VEP mappings; {len(out)} rows")


if __name__ == "__main__":
    main()
