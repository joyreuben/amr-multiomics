#!/usr/bin/env python3
"""
Merge Pathogenwatch's summary, metrics, kleborate, and plasmidfinder exports
into one clean per-genome table for the West Africa comparison dataset.

Drops:
  - Cameroon genomes (out of scope: Central/Middle Africa, not West Africa)
  - Genomes that failed Pathogenwatch's own QC check

Usage (run from the folder containing the 8 pathogenwatch-*.csv files):
    python clean_pathogenwatch_west_africa.py
"""
import csv
from collections import defaultdict, Counter


def load(path):
    with open(path, encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def main():
    summary = load("pathogenwatch-summary.csv")
    metrics = {r["Genome ID"]: r for r in load("pathogenwatch-metrics.csv")}
    kleborate = {r["Genome ID"]: r for r in load("pathogenwatch-kleborate.csv")}

    plasmid_rows = load("pathogenwatch-plasmidfinder.csv")
    plasmids_by_genome = defaultdict(set)
    for row in plasmid_rows:
        plasmids_by_genome[row["Genome ID"]].add(row["Inc Match"])

    out_rows = []
    excluded_cameroon = 0
    excluded_qc_failed = 0

    for row in summary:
        country_raw = row["Country"]
        if "Cameroon" in country_raw:
            excluded_cameroon += 1
            continue
        if row["QC"] != "Passed":
            excluded_qc_failed += 1
            continue

        gid = row["Genome ID"]
        country = country_raw.split(" - ", 1)[-1] if " - " in country_raw else country_raw
        kb = kleborate.get(gid, {})
        mt = metrics.get(gid, {})

        out_rows.append({
            "genome_id": gid,
            "genome_name": row["Genome Name"],
            "country": country,
            "qc_status": row["QC"],
            "species_prediction": row["Species Prediction"],
            "contig_count": mt.get("No. Contigs", ""),
            "genome_length": mt.get("Genome Length", ""),
            "n50": mt.get("N50", ""),
            "st": kb.get("ST", ""),
            "resistance_score": kb.get("resistance_score", ""),
            "resistance_gene_count": kb.get("resistance_gene_count", ""),
            "bla_acquired": kb.get("Bla_acquired", ""),
            "bla_esbl_acquired": kb.get("Bla_ESBL_acquired", ""),
            "bla_carb_acquired": kb.get("Bla_Carb_acquired", ""),
            "virulence_score": kb.get("virulence_score", ""),
            "k_type": kb.get("K_type", ""),
            "o_type": kb.get("O_type", ""),
            "plasmid_replicons": ";".join(sorted(plasmids_by_genome.get(gid, []))),
        })

    fields = list(out_rows[0].keys())
    with open("pathogenwatch_west_africa_summary.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(out_rows)

    print(f"Excluded {excluded_cameroon} Cameroon genomes")
    print(f"Excluded {excluded_qc_failed} QC-failed genomes")
    print(f"Final row count: {len(out_rows)}")
    print(Counter(r["country"] for r in out_rows))


if __name__ == "__main__":
    main()
