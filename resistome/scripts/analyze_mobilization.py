#!/usr/bin/env python3
"""
Objective 2: is the carbapenemase on a plasmid, and is it the same plasmid in
lineages too distant for clonal descent to explain?

Joins AMRFinderPlus gene coordinates to mob_recon's verdict on each contig, so
every carbapenemase call gets a location (chromosome or plasmid), a mob_suite
cluster, a replicon type where one resolves, and a predicted mobility.

The cluster is what carries the argument. Replicon typing returned "-" on the
carbapenemase contigs of the verified test genome, so replicon names alone
cannot be relied on; mob_suite primary cluster IDs are assigned by mash distance
to a reference plasmid database and are comparable across genomes, which is what
a shared-vehicle claim needs. Both are reported and the script says how often
each resolves.

Three comparisons decide the objective:

  ST17 vs ST234    OXA-181 in a 28-genome clonal outbreak and in an unrelated
                   lineage. Same cluster means one vehicle crossing lineages;
                   different clusters mean OXA-181 arrived twice by separate
                   routes.
  ST340            two genomes 4 alleles apart carrying different families.
                   Different clusters confirm two separate arrivals into what is
                   effectively the same host background.
  ST395            a positive and a negative 0 alleles apart. If the negative
                   lacks the cluster entirely, the plasmid left (or never
                   arrived) rather than the gene being lost from it.

Sources: resistome/results/pw_*_amrfinder.tsv and pw_*_mob_recon/contig_report.txt
         produced by run_mobilization_tools.sh, plus mobilization_subset.tsv
         for ST and family.

Usage:
    python resistome/scripts/analyze_mobilization.py
"""
import argparse
import glob
import os
import re

import pandas as pd

# The collection contains only NDM and OXA-48-like, but the pattern covers the
# other acquired carbapenemases so a future genome is not silently missed.
CARBAPENEMASE = re.compile(
    r"bla(?:NDM|KPC|VIM|IMP|GES)|blaOXA-(?:48|181|232|204|244|245)", re.IGNORECASE)


def first_token(value):
    return str(value).split()[0] if str(value).strip() else ""


def family_of(gene):
    gene = gene.replace("bla", "")
    return gene.split("-")[0].upper()


def load_sample(sample, results_dir):
    """Join one genome's AMR calls to its contig assignments."""
    amr_path = os.path.join(results_dir, f"{sample}_amrfinder.tsv")
    mob_path = os.path.join(results_dir, f"{sample}_mob_recon", "contig_report.txt")
    if not (os.path.exists(amr_path) and os.path.exists(mob_path)):
        return None

    amr = pd.read_csv(amr_path, sep="\t", dtype=str)
    mob = pd.read_csv(mob_path, sep="\t", dtype=str)

    gene_column = next((c for c in amr.columns
                        if c.lower() in ("element symbol", "gene symbol")), None)
    contig_column = next((c for c in amr.columns if "contig" in c.lower()), None)
    if gene_column is None or contig_column is None:
        return None

    amr["contig"] = amr[contig_column].map(first_token)
    mob["contig"] = mob["contig_id"].map(first_token)

    keep = ["contig", "molecule_type", "primary_cluster_id", "rep_type(s)",
            "predicted_mobility"]
    keep = [c for c in keep if c in mob.columns]
    joined = amr.merge(mob[keep], on="contig", how="left")
    joined["sample"] = sample
    joined["gene"] = joined[gene_column]
    return joined


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results_dir", default="resistome/results")
    parser.add_argument("--subset", default="resistome/results/mobilization_subset.tsv")
    parser.add_argument("--out", default="resistome/results/mobilization_placement.tsv")
    args = parser.parse_args()

    samples = sorted(
        os.path.basename(p).replace("_amrfinder.tsv", "")
        for p in glob.glob(os.path.join(args.results_dir, "pw_*_amrfinder.tsv")))
    frames = [f for f in (load_sample(s, args.results_dir) for s in samples) if f is not None]
    if not frames:
        raise SystemExit("No processed genomes found. Run run_mobilization_tools.sh first.")
    everything = pd.concat(frames, ignore_index=True)
    print(f"{len(frames)} genomes processed, {len(everything)} AMR gene calls")

    subset = pd.read_csv(args.subset, sep="\t", dtype=str)
    subset["sample"] = "pw_" + subset["run_accession"]
    meta = subset.set_index("sample")[["ST", "family", "carbapenemase", "rationale"]]
    everything = everything.join(meta, on="sample")

    carb = everything[everything["gene"].str.contains(CARBAPENEMASE, na=False)].copy()
    carb["gene_family"] = carb["gene"].map(family_of)

    columns = ["sample", "ST", "gene", "gene_family", "contig", "molecule_type",
               "primary_cluster_id", "rep_type(s)", "predicted_mobility"]
    columns = [c for c in columns if c in carb.columns]
    carb[columns].sort_values(["gene_family", "ST", "sample"]).to_csv(
        args.out, sep="\t", index=False)

    print(f"\n--- where each carbapenemase sits ---")
    print(f"{'sample':<18}{'ST':<9}{'gene':<12}{'location':<12}{'cluster':<10}"
          f"{'replicon':<22}{'mobility':<12}")
    for _, r in carb[columns].sort_values(["gene_family", "ST"]).iterrows():
        print(f"{r['sample']:<18}{str(r.get('ST','?')):<9}{r['gene']:<12}"
              f"{str(r.get('molecule_type','?')):<12}{str(r.get('primary_cluster_id','-')):<10}"
              f"{str(r.get('rep_type(s)','-'))[:20]:<22}{str(r.get('predicted_mobility','-')):<12}")

    on_plasmid = (carb["molecule_type"] == "plasmid").sum()
    print(f"\n{on_plasmid}/{len(carb)} carbapenemase calls sit on a plasmid contig")
    resolved_rep = (~carb["rep_type(s)"].isin(["-", ""]) & carb["rep_type(s)"].notna()).sum() \
        if "rep_type(s)" in carb.columns else 0
    print(f"replicon type resolves for {resolved_rep}/{len(carb)}; "
          f"mob_suite cluster for "
          f"{int(carb['primary_cluster_id'].notna().sum())}/{len(carb)}")

    # --- is one vehicle crossing lineages? ---
    print("\n--- clusters carrying the same gene family across STs ---")
    plasmid = carb[carb["molecule_type"] == "plasmid"]
    for fam, group in plasmid.groupby("gene_family"):
        for cluster, rows in group.groupby("primary_cluster_id"):
            sts = sorted(set(rows["ST"].dropna()))
            flag = "  <-- crosses lineages" if len(sts) > 1 else ""
            print(f"  {fam:<5} cluster {cluster:<8} in {len(rows)} genome(s), "
                  f"ST(s): {', '.join(sts)}{flag}")

    # --- what else rides on the carbapenemase plasmid ---
    # Aggregated over the whole mob_suite cluster, not the single contig. These
    # are short-read assemblies of 167-215 contigs, so one plasmid is routinely
    # split across several of them; asking what shares a contig understates the
    # cargo, while the cluster is mob_recon's reconstruction of the replicon.
    print("\n--- cargo on carbapenemase-carrying plasmid clusters ---")
    for (sample, cluster), rows in carb[carb["molecule_type"] == "plasmid"].groupby(
            ["sample", "primary_cluster_id"]):
        mates = everything[(everything["sample"] == sample) &
                           (everything["primary_cluster_id"] == cluster)]["gene"]
        others = sorted(set(mates) - set(rows["gene"]))
        st = rows.iloc[0].get("ST", "?")
        print(f"  {sample} ({st}) cluster {cluster} {', '.join(sorted(set(rows['gene'])))}: "
              f"{len(others)} co-located — {', '.join(others[:10]) if others else 'none'}")

    print("\n--- is the objective 3 cargo actually on the carbapenemase plasmid? ---")
    # Objective 3 found armA, rmtB, mphE/msrE and aph(3')-VI tracking the
    # carbapenemase families at genome level. Genome level cannot say whether
    # they ride the same replicon, and this can.
    for marker in ["armA", "rmtB", "aph(3')-VI", "mphE", "msrE", "blaCTX-M-15", "qnrS1"]:
        hits = everything[everything["gene"] == marker]
        if hits.empty:
            continue
        same = 0
        for _, hit in hits.iterrows():
            carb_clusters = set(carb[(carb["sample"] == hit["sample"]) &
                                     (carb["molecule_type"] == "plasmid")]["primary_cluster_id"])
            if hit["primary_cluster_id"] in carb_clusters:
                same += 1
        print(f"  {marker:<14} present in {len(hits):>2} genome(s); "
              f"on the carbapenemase cluster in {same}")

    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
