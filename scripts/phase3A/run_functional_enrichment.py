#!/usr/bin/env python3
"""Build the descriptive locus matrix and direction-specific pathway tables."""
from __future__ import annotations

import csv
import json
import math
import urllib.parse
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FREEZE = ROOT / "config/phase3A_LOCUS_FREEZE.tsv"
INV = ROOT / "results/phase3A/REGION_GENE_INVENTORY.tsv"
BG = ROOT / "results/phase3A/FUNCTIONAL_BACKGROUND_GENES.tsv"
XQTL = ROOT / "results/phase3A/XQTL_LOOKUP.tsv"
OUT = ROOT / "results/phase3A"

DOMAINS = "RETINAL/OCULAR;VASCULAR/ENDOTHELIAL;IMMUNE/INFLAMMATORY;LIPID/METABOLIC;ECM;NEURONAL/NEURODEGENERATIVE;RENAL;OXIDATIVE/MITOCHONDRIAL"


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def panther(genes, background):
    if not genes:
        return {"status": "NO_INPUT_GENES", "result": []}
    params = {
        "geneInputList": ",".join(sorted(genes)),
        "refInputList": ",".join(sorted(background)),
        "organism": "9606",
        "refOrganism": "9606",
        "annotDataSet": "GO:0008150",
        "enrichmentTestType": "FISHER",
        "correction": "FDR",
    }
    url = "https://pantherdb.org/services/oai/pantherdb/enrich/overrep"
    try:
        req = urllib.request.Request(
            url, data=urllib.parse.urlencode(params).encode(), method="POST",
            headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json", "Content-Type": "application/x-www-form-urlencoded"},
        )
        with urllib.request.urlopen(req, timeout=120) as response:
            payload = json.loads(response.read().decode())
        result_block = payload.get("results", {})
        if "result" not in result_block:
            error = result_block.get("search", {}).get("error", "NO_RESULT_FIELD") if isinstance(result_block, dict) else "NO_RESULTS_BLOCK"
            return {"status": "API_ERROR: " + str(error), "result": [], "meta": result_block}
        return {"status": "CUSTOM_BACKGROUND", "result": result_block.get("result", []), "meta": result_block}
    except Exception as error:
        return {"status": f"API_ERROR: {type(error).__name__}: {error}", "result": []}


def reactome(genes):
    if not genes:
        return {"status": "NO_INPUT_GENES", "pathways": []}
    body = ("\n".join(sorted(genes)) + "\n").encode()
    request = urllib.request.Request(
        "https://reactome.org/AnalysisService/identifiers/", data=body,
        headers={"Content-Type": "text/plain", "Accept": "application/json", "User-Agent": "ocular-systemic-antagonistic-pleiotropy/phase3A"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            payload = json.loads(response.read().decode())
        return {"status": "DESCRIPTIVE_DEFAULT_BACKGROUND", "pathways": payload.get("pathways", [])}
    except Exception as error:
        return {"status": f"API_ERROR: {type(error).__name__}: {error}", "pathways": []}


def main():
    loci = read_tsv(FREEZE)
    inv = read_tsv(INV)
    bg = read_tsv(BG)
    xqtl = read_tsv(XQTL)
    assert len(loci) == 58
    locus_keys = {(r["pair"], r["chr"], r["lead_SNP"]): r for r in loci}
    by_locus = defaultdict(lambda: {"LEVEL1": set(), "LEVEL2": set(), "LEVEL3": set()})
    for row in inv:
        key = (row.get("pair", ""), row.get("chr", ""), row.get("lead_SNP", ""))
        if key not in locus_keys:
            # Region-level rows without a single lead are used below only as
            # region inventory and not as locus-level mappings.
            continue
        if row.get("gene") and row.get("mapping_level") in {"LEVEL1", "LEVEL2", "LEVEL3"}:
            by_locus[key][row["mapping_level"]].add(row["gene"])
    qtl_by_locus = defaultdict(set)
    for row in xqtl:
        if row.get("gene_symbol") and row.get("status") == "LOOKUP_RETURNED":
            qtl_by_locus[row.get("lead_SNP", "")].add(row["gene_symbol"])

    matrix_fields = [
        "PLACO_locus_id", "pair", "hidden_mixed_pair", "chr", "lead_SNP", "lead_BP", "SUPERGNOVA_region",
        "direction_class", "SUPERGNOVA_q", "PLACO_P", "positional_priority", "level1_genes", "level2_genes",
        "level3_genes", "primary_mapping_genes", "one_locus_one_gene", "qtl_lookup_genes",
        "qtl_lookup_status", "predefined_functional_domains", "domain_assignment_status", "interpretation_status",
    ]
    matrix = []
    one_gene = {}
    all_gene_sets = {"NEGATIVE_ALL_LEVELS": set(), "POSITIVE_ALL_LEVELS": set(), "NEGATIVE_ONE_LOCUS_ONE_GENE": set(), "POSITIVE_ONE_LOCUS_ONE_GENE": set()}
    for locus in loci:
        key = (locus["pair"], locus["chr"], locus["lead_SNP"])
        levels = by_locus[key]
        if levels["LEVEL1"]:
            primary_level, primary = "LEVEL1", sorted(levels["LEVEL1"])
        elif levels["LEVEL2"]:
            primary_level, primary = "LEVEL2", sorted(levels["LEVEL2"])
        elif levels["LEVEL3"]:
            primary_level, primary = "LEVEL3", sorted(levels["LEVEL3"])
        else:
            primary_level, primary = "UNMAPPED", []
        one = primary[0] if primary else ""
        one_gene[key] = one
        direction = "NEGATIVE" if "NEGATIVE" in locus["direction_class"] else "POSITIVE"
        all_gene_sets[f"{direction}_ALL_LEVELS"].update(set().union(*levels.values()))
        if one:
            all_gene_sets[f"{direction}_ONE_LOCUS_ONE_GENE"].add(one)
        qtls = sorted(qtl_by_locus.get(locus["lead_SNP"], set()))
        matrix.append({
            "PLACO_locus_id": locus["PLACO_locus_id"], "pair": locus["pair"], "hidden_mixed_pair": locus["hidden_mixed_pair"],
            "chr": locus["chr"], "lead_SNP": locus["lead_SNP"], "lead_BP": locus["lead_BP"],
            "SUPERGNOVA_region": locus["SUPERGNOVA_region"], "direction_class": locus["direction_class"],
            "SUPERGNOVA_q": locus["SUPERGNOVA_q"], "PLACO_P": locus["PLACO_P"],
            "positional_priority": primary_level, "level1_genes": ";".join(sorted(levels["LEVEL1"])),
            "level2_genes": ";".join(sorted(levels["LEVEL2"])), "level3_genes": ";".join(sorted(levels["LEVEL3"])),
            "primary_mapping_genes": ";".join(primary), "one_locus_one_gene": one,
            "qtl_lookup_genes": ";".join(qtls), "qtl_lookup_status": "LOOKUP_RETURNED" if qtls else "NO_GTEx_V8_GENE_RECORD",
            "predefined_functional_domains": DOMAINS, "domain_assignment_status": "NOT_ASSIGNED_FROM_POSITIONAL_OR_XQTL_LOOKUP",
            "interpretation_status": "DESCRIPTIVE_SUPPORT_ONLY_NO_CAUSAL_GENE_CLAIM",
        })
    write_tsv(OUT / "LOCUS_BIOLOGICAL_MATRIX.tsv", matrix, matrix_fields)

    background_genes = {r["gene_symbol"] for r in bg if r.get("gene_symbol")}
    # Rank recurrent regions by independent locus count, then pair count.
    region_members = defaultdict(list)
    for locus in loci:
        region_members[locus["SUPERGNOVA_region"]].append(locus)
    recurrence = []
    for region, members in region_members.items():
        recurrence.append((region, len({r["pair"] for r in members}), len(members)))
    recurrence.sort(key=lambda x: (-x[1], -x[2], x[0]))
    top1 = {x[0] for x in recurrence[:1]}
    top3 = {x[0] for x in recurrence[:3]}

    # Sets used for the primary analysis and predefined major-locus sensitivity.
    analysis_sets = {}
    for direction in ("NEGATIVE", "POSITIVE"):
        for mode in ("ALL_LEVELS", "ONE_LOCUS_ONE_GENE"):
            analysis_sets[f"{direction}_{mode}"] = set(all_gene_sets[f"{direction}_{mode}"])
        for label, excluded in (("EXCLUDE_TOP1", top1), ("EXCLUDE_TOP3", top3)):
            genes = set()
            for locus in loci:
                if direction in locus["direction_class"] and locus["SUPERGNOVA_region"] not in excluded:
                    key = (locus["pair"], locus["chr"], locus["lead_SNP"])
                    genes.add(one_gene.get(key, ""))
            analysis_sets[f"{direction}_{label}"] = {x for x in genes if x}

    pathway_fields = [
        "analysis_set", "direction", "annotation_source", "background_scope", "n_input_genes", "n_background_genes",
        "term_id", "term_name", "n_in_list", "n_in_reference", "expected", "fold_enrichment", "p_value", "fdr",
        "status", "query_date",
    ]
    pathway_rows = []
    sensitivity_rows = []
    results_cache = {}
    with ThreadPoolExecutor(max_workers=4) as pool:
        panther_futures = {pool.submit(panther, genes, background_genes): name for name, genes in analysis_sets.items()}
        reactome_futures = {pool.submit(reactome, genes): name for name, genes in analysis_sets.items()}
        panther_results = {}
        reactome_results = {}
        for future in as_completed(panther_futures):
            panther_results[panther_futures[future]] = future.result()
        for future in as_completed(reactome_futures):
            reactome_results[reactome_futures[future]] = future.result()

    for name, genes in analysis_sets.items():
        direction = "NEGATIVE" if name.startswith("NEGATIVE") else "POSITIVE"
        p = panther_results.get(name, {"status": "MISSING_RESULT", "result": []})
        results_cache[name] = {"genes": genes, "panther": p}
        returned = 0
        best = None
        for term in p.get("result", []):
            returned += 1
            fdr = term.get("fdr", "")
            if best is None or (isinstance(fdr, (float, int)) and fdr < best[0]):
                best = (fdr, term)
            pathway_rows.append({
                "analysis_set": name, "direction": direction, "annotation_source": "PANTHER_GO_BIOLOGICAL_PROCESS",
                "background_scope": "FULL_2096_VALID_LOCI_MAPPED_PROTEIN_CODING_GENES",
                "n_input_genes": len(genes), "n_background_genes": len(background_genes),
                "term_id": term.get("term", {}).get("id", ""), "term_name": term.get("term", {}).get("label", ""),
                "n_in_list": term.get("number_in_list", ""), "n_in_reference": term.get("number_in_reference", ""),
                "expected": term.get("expected", ""), "fold_enrichment": term.get("fold_enrichment", ""),
                "p_value": term.get("pValue", ""), "fdr": fdr, "status": p.get("status", ""), "query_date": str(date.today()),
            })
        sensitivity_rows.append({
            "analysis_set": name, "direction": direction, "excluded_regions": "TOP1" if "EXCLUDE_TOP1" in name else "TOP3" if "EXCLUDE_TOP3" in name else "NONE",
            "n_loci": sum(1 for l in loci if direction in l["direction_class"] and (l["SUPERGNOVA_region"] not in (top1 if "EXCLUDE_TOP1" in name else top3 if "EXCLUDE_TOP3" in name else set()))),
            "n_genes": len(genes), "n_background_genes": len(background_genes),
            "n_terms_returned": returned, "n_terms_fdr_lt_0_05": sum(1 for x in p.get("result", []) if isinstance(x.get("fdr"), (float, int)) and x["fdr"] < 0.05),
            "top_term": best[1].get("term", {}).get("label", "") if best else "",
            "top_term_fdr": best[0] if best else "", "status": p.get("status", ""),
        })

        # Reactome is retained as a separate descriptive lookup.  Its API does
        # not expose the same custom background in this workflow, so it is not
        # used for the primary enrichment conclusion.
        r = reactome_results.get(name, {"status": "MISSING_RESULT", "pathways": []})
        for term in r.get("pathways", [])[:500]:
            entities = term.get("entities", {})
            pathway_rows.append({
                "analysis_set": name, "direction": direction, "annotation_source": "REACTOME_ANALYSIS_SERVICE",
                "background_scope": "REACTOME_DEFAULT_BACKGROUND_NOT_PRIMARY", "n_input_genes": len(genes),
                "n_background_genes": "", "term_id": term.get("stId", ""), "term_name": term.get("name", ""),
                "n_in_list": entities.get("found", ""), "n_in_reference": entities.get("total", ""),
                "expected": "", "fold_enrichment": entities.get("ratio", ""), "p_value": entities.get("pValue", ""),
                "fdr": entities.get("fdr", ""), "status": r.get("status", ""), "query_date": str(date.today()),
            })
    write_tsv(OUT / "PATHWAY_ENRICHMENT.tsv", pathway_rows, pathway_fields)
    sens_fields = [
        "analysis_set", "direction", "excluded_regions", "n_loci", "n_genes", "n_background_genes",
        "n_terms_returned", "n_terms_fdr_lt_0_05", "top_term", "top_term_fdr", "status",
    ]
    write_tsv(OUT / "MAJOR_LOCUS_ENRICHMENT_SENSITIVITY.tsv", sensitivity_rows, sens_fields)

    # Update recurrence table with locus-level primary mapping genes.
    recurrence_rows = []
    for region, n_pairs, n_loci in recurrence:
        members = [l for l in loci if l["SUPERGNOVA_region"] == region]
        mapped = []
        for locus in members:
            key = (locus["pair"], locus["chr"], locus["lead_SNP"])
            if one_gene.get(key):
                mapped.append(one_gene[key])
        recurrence_rows.append({
            "genomic_region": region, "chr": members[0]["chr"], "region_start": region.split(":", 1)[1].split("-", 1)[0],
            "region_end": region.rsplit("-", 1)[1], "n_loci": n_loci, "n_pairs": n_pairs,
            "pairs": ";".join(sorted({l["pair"] for l in members})),
            "n_negative_loci": sum("NEGATIVE" in l["direction_class"] for l in members),
            "n_positive_loci": sum("POSITIVE" in l["direction_class"] for l in members),
            "lead_variants": ";".join(sorted({l["lead_SNP"] for l in members})),
            "mapped_genes": ";".join(sorted(set(mapped))),
            "biological_support_status": "POSITIONAL_OR_XQTL_SUPPORT_ONLY" if mapped else "NO_PRIMARY_GENE_MAPPING",
        })
    recurrence_fields = [
        "genomic_region", "chr", "region_start", "region_end", "n_loci", "n_pairs", "pairs", "n_negative_loci",
        "n_positive_loci", "lead_variants", "mapped_genes", "biological_support_status",
    ]
    write_tsv(OUT / "RECURRENT_CROSS_PAIR_LOCI.tsv", recurrence_rows, recurrence_fields)

    (OUT / "FUNCTIONAL_ANALYSIS_AUDIT.txt").write_text(
        "primary_gene_background=all_unique_protein_coding_genes_from_2096_valid_loci\n"
        f"background_gene_count={len(background_genes)}\n"
        "primary_enrichment=PANTHER_GO_BP_with_custom_reference_list\n"
        "reactome=descriptive_default_background_only\n"
        "causal_gene_claims=NOT_MADE\n"
    )
    print(f"Functional matrix and enrichment complete: matrix={len(matrix)} loci, background_genes={len(background_genes)}, pathway_rows={len(pathway_rows)}")


if __name__ == "__main__":
    main()
