#!/usr/bin/env python3
"""Create a descriptive protein-coding gene inventory for 41 supported regions."""
from __future__ import annotations

import csv
import json
import time
from datetime import date
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[2]
REGIONS = ROOT / "results/phase2AR/REGION_LEVEL_DIRECTIONAL_CONVERGENCE.tsv"
FREEZE = ROOT / "config/phase3A_LOCUS_FREEZE.tsv"
VEP = ROOT / "results/phase3A/VEP_LEAD_VARIANT_ANNOTATION.tsv"
OUT = ROOT / "results/phase3A/REGION_GENE_INVENTORY.tsv"
CACHE = ROOT / "results/phase3A/ensembl_region_gene_cache.json"
BASE = "https://grch37.rest.ensembl.org"


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def api_json(path):
    request = Request(BASE + path, headers={
        "Content-Type": "application/json",
        "Accept": "application/json",
        "User-Agent": "ocular-systemic-antagonistic-pleiotropy/phase3A",
    })
    with urlopen(request, timeout=90) as response:
        return json.loads(response.read().decode())


def gene_symbol(item):
    return item.get("external_name") or item.get("display_name") or item.get("id", "")


def distance_to_gene(bp, gene):
    start, end = int(gene["start"]), int(gene["end"])
    if start <= bp <= end:
        return 0
    return min(abs(bp - start), abs(bp - end))


def main():
    regions = [r for r in read_tsv(REGIONS) if int(float(r.get("N_independent_PLACO_loci", 0) or 0)) > 0]
    loci = read_tsv(FREEZE)
    vep = read_tsv(VEP) if VEP.exists() else []
    assert len(regions) == 41
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}

    def cached(key, path):
        if key in cache:
            return cache[key]
        try:
            value = api_json(path)
            cache[key] = value
        except (HTTPError, URLError, TimeoutError, OSError, json.JSONDecodeError) as error:
            cache[key] = {"_error": f"{type(error).__name__}: {error}"}
        CACHE.write_text(json.dumps(cache, indent=2, sort_keys=True))
        time.sleep(0.12)
        return cache[key]

    fields = [
        "pair", "genomic_region", "chr", "region_start", "region_end", "SUPERGNOVA_rho",
        "direction_class", "lead_SNP", "lead_BP", "gene", "gene_id", "biotype",
        "gene_start", "gene_end", "distance_to_lead", "mapping_type", "mapping_level",
        "source", "build", "annotation_date", "status",
    ]
    rows = []
    for region in regions:
        chrom = region["chr"]
        start = int(float(region["region_start"]))
        end = int(float(region["region_end"]))
        region_name = region["SUPERGNOVA_region"]
        members = [l for l in loci if l["SUPERGNOVA_region"] == region_name]
        lead_names = {l["lead_SNP"] for l in members}
        direction_class = "NEGATIVE_DIRECTION_CONCORDANT" if float(region["SUPERGNOVA_rho"]) < 0 else "POSITIVE_DIRECTION_CONCORDANT"
        key = f"region:{chrom}:{start}-{end}"
        data = cached(key, f"/overlap/region/human/{quote(str(chrom))}:{start}-{end}?feature=gene")
        if isinstance(data, dict) and "_error" in data:
            rows.append({
                "pair": region["pair"], "genomic_region": region_name, "chr": chrom,
                "region_start": start, "region_end": end, "SUPERGNOVA_rho": region["SUPERGNOVA_rho"],
                "direction_class": direction_class, "gene": "", "mapping_type": "GENE_BODY_OVERLAP",
                "mapping_level": "LEVEL2", "source": "Ensembl overlap/region GRCh37",
                "build": "GRCh37", "annotation_date": str(date.today()), "status": "API_ERROR: " + data["_error"],
            })
            genes = []
        else:
            genes = [g for g in (data or []) if g.get("feature_type") == "gene" and g.get("biotype") == "protein_coding"]
            for gene in genes:
                rows.append({
                    "pair": region["pair"], "genomic_region": region_name, "chr": chrom,
                    "region_start": start, "region_end": end, "SUPERGNOVA_rho": region["SUPERGNOVA_rho"],
                    "direction_class": direction_class, "lead_SNP": ";".join(sorted(lead_names)),
                    "gene": gene_symbol(gene), "gene_id": gene.get("id", ""), "biotype": gene.get("biotype", ""),
                    "gene_start": gene.get("start", ""), "gene_end": gene.get("end", ""),
                    "mapping_type": "GENE_BODY_OVERLAP", "mapping_level": "LEVEL2",
                    "source": "Ensembl overlap/region GRCh37", "build": "GRCh37",
                    "annotation_date": str(date.today()), "status": "ANNOTATED",
                })

        for locus in members:
            bp = int(float(locus["lead_BP"]))
            window_start = max(1, bp - 1_000_000)
            window_end = bp + 1_000_000
            key = f"nearest:{chrom}:{bp}"
            near = cached(key, f"/overlap/region/human/{quote(str(chrom))}:{window_start}-{window_end}?feature=gene")
            if isinstance(near, dict) and "_error" in near:
                rows.append({
                    "pair": locus["pair"], "genomic_region": region_name, "chr": chrom,
                    "region_start": start, "region_end": end, "SUPERGNOVA_rho": locus["SUPERGNOVA_rho"],
                    "direction_class": locus["direction_class"], "lead_SNP": locus["lead_SNP"], "lead_BP": bp,
                    "mapping_type": "NEAREST_PROTEIN_CODING_GENE", "mapping_level": "LEVEL3",
                    "source": "Ensembl overlap/region GRCh37", "build": "GRCh37",
                    "annotation_date": str(date.today()), "status": "API_ERROR: " + near["_error"],
                })
                near_genes = []
            else:
                near_genes = [g for g in (near or []) if g.get("feature_type") == "gene" and g.get("biotype") == "protein_coding"]
            if near_genes:
                distances = [(distance_to_gene(bp, g), g) for g in near_genes]
                best_distance = min(d for d, _ in distances)
                best = sorted({(gene_symbol(g), g.get("id", ""), d, g.get("start", ""), g.get("end", "")) for d, g in distances if d == best_distance})
                for symbol, gid, dist, gstart, gend in best:
                    rows.append({
                        "pair": locus["pair"], "genomic_region": region_name, "chr": chrom,
                        "region_start": start, "region_end": end, "SUPERGNOVA_rho": locus["SUPERGNOVA_rho"],
                        "direction_class": locus["direction_class"], "lead_SNP": locus["lead_SNP"], "lead_BP": bp,
                        "gene": symbol, "gene_id": gid, "biotype": "protein_coding", "gene_start": gstart,
                        "gene_end": gend, "distance_to_lead": dist,
                        "mapping_type": "NEAREST_PROTEIN_CODING_GENE", "mapping_level": "LEVEL3",
                        "source": "Ensembl overlap/region GRCh37", "build": "GRCh37",
                        "annotation_date": str(date.today()), "status": "ANNOTATED",
                    })

    # Add VEP-supported protein-coding transcript annotations as Level 1 rows.
    for ann in vep:
        if not ann.get("gene_symbol") or ann.get("status", "").startswith("API_ERROR"):
            continue
        locus = next((x for x in loci if x["lead_SNP"] == ann["lead_SNP"]), None)
        if not locus:
            continue
        consequences = ann.get("consequence", "")
        high = any(x in consequences for x in ("splice", "stop_gained", "frameshift", "missense", "inframe", "start_lost", "stop_lost", "transcript_ablation", "coding_sequence_variant"))
        level = "LEVEL1" if high else "LEVEL4"
        if ann.get("biotype") == "protein_coding" and high:
            rows.append({
                "pair": locus["pair"], "genomic_region": locus["SUPERGNOVA_region"], "chr": locus["chr"],
                "region_start": locus["SUPERGNOVA_region"].split(":", 1)[1].split("-", 1)[0],
                "region_end": locus["SUPERGNOVA_region"].rsplit("-", 1)[1], "SUPERGNOVA_rho": locus["SUPERGNOVA_rho"],
                "direction_class": locus["direction_class"], "lead_SNP": locus["lead_SNP"], "lead_BP": locus["lead_BP"],
                "gene": ann["gene_symbol"], "gene_id": ann.get("gene_id", ""), "biotype": ann.get("biotype", ""),
                "distance_to_lead": 0, "mapping_type": "VEP_TRANSCRIPT_CONSEQUENCE", "mapping_level": level,
                "source": "Ensembl VEP GRCh37 REST", "build": "GRCh37", "annotation_date": ann.get("annotation_date", str(date.today())),
                "status": ann.get("status", "ANNOTATED"),
            })
    write_tsv(OUT, rows, fields)
    print(f"Region gene inventory complete: {len(rows)} rows for {len(regions)} regions")


if __name__ == "__main__":
    main()
