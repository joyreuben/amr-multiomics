#!/usr/bin/env python3
"""
Combine AMRFinderPlus, RGI, PlasmidFinder, and MOB-suite output for one sample
into a single summary row, appended to a shared CSV.

Usage:
    python resistome/scripts/summarize_resistome.py --sample kpn_test \
        --amrfinder resistome/results/kpn_amrfinder.tsv \
        --rgi resistome/results/kpn_rgi.txt \
        --plasmidfinder resistome/results/kpn_plasmidfinder/results_tab.tsv \
        --mob_contig_report resistome/results/kpn_mob_recon/contig_report.txt \
        --out resistome/results/resistome_summary.csv
"""
import argparse
import os
import pandas as pd


def summarize_amrfinder(path):
    df = pd.read_csv(path, sep="\t")
    genes = sorted(df["Element symbol"].dropna().unique().tolist())
    classes = sorted(df["Class"].dropna().unique().tolist())
    return len(df), genes, classes


def summarize_rgi(path):
    df = pd.read_csv(path, sep="\t")
    genes = sorted(df["Best_Hit_ARO"].dropna().unique().tolist())
    drug_classes = sorted({
        c.strip() for entry in df["Drug Class"].dropna() for c in entry.split(";")
    })
    return len(df), genes, drug_classes


def summarize_plasmidfinder(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return []
    return sorted(df["Plasmid"].dropna().unique().tolist())


def summarize_mob_recon(path):
    df = pd.read_csv(path, sep="\t")
    plasmid_rows = df[df["molecule_type"] == "plasmid"]
    mobility = sorted(plasmid_rows["predicted_mobility"].dropna().unique().tolist())
    return len(plasmid_rows), mobility


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    p.add_argument("--amrfinder", required=True)
    p.add_argument("--rgi", required=True)
    p.add_argument("--plasmidfinder", required=True)
    p.add_argument("--mob_contig_report", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    amr_count, amr_genes, amr_classes = summarize_amrfinder(args.amrfinder)
    rgi_count, rgi_genes, rgi_drug_classes = summarize_rgi(args.rgi)
    plasmid_replicons = summarize_plasmidfinder(args.plasmidfinder)
    num_plasmid_contigs, mobility = summarize_mob_recon(args.mob_contig_report)

    row = {
        "sample_id": args.sample,
        "amrfinder_gene_count": amr_count,
        "amrfinder_genes": ";".join(amr_genes),
        "amrfinder_classes": ";".join(amr_classes),
        "rgi_gene_count": rgi_count,
        "rgi_genes": ";".join(rgi_genes),
        "rgi_drug_classes": ";".join(rgi_drug_classes),
        "plasmidfinder_replicons": ";".join(plasmid_replicons),
        "mob_recon_plasmid_contigs": num_plasmid_contigs,
        "mob_recon_predicted_mobility": ";".join(mobility),
    }

    out_df = pd.DataFrame([row])
    write_header = not os.path.exists(args.out)
    out_df.to_csv(args.out, mode="a", header=write_header, index=False)
    print(f"Wrote summary for sample '{args.sample}' to {args.out}")


if __name__ == "__main__":
    main()
