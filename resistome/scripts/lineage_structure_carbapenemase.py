#!/usr/bin/env python3
"""
Objective 1: lineage structure of carbapenemase carriage.

For each sequence type that carries an acquired carbapenemase, ask whether the
positives are one clone that expanded or several independent acquisitions of
the gene. The discriminating evidence is cgMLST allele distance: positives that
sit within a few alleles of each other are one clone; positives scattered
across the ST's diversity, interleaved with negatives, are separate
acquisitions of the same gene into the same background.

Distance is the number of loci where two genomes have different alleles,
counted only over loci called in both (loci per genome range 590-628 of 629),
then rescaled to the full scheme so genomes with different amounts of missing
data stay comparable.

Carbapenemase is resolved to gene family, never to a binary. NDM and OXA-48-like
are two independent gene flows and an ST carrying both families has by
definition acquired the phenotype twice. See claude.md.

Sources (Pathogenwatch West Africa export, downloaded 2026-09):
    pathogenwatch-cgmlst.csv     629-locus allele calls, long format
    pathogenwatch-kleborate.csv  Kleborate v3.2.4, ST and Bla_Carb_acquired

Usage:
    python resistome/scripts/lineage_structure_carbapenemase.py \
        --pathogenwatch_dir resistome/results/pathogenwatch_west_africa \
        --out resistome/results/lineage_structure.tsv
"""
import argparse
import os

import numpy as np
import pandas as pd

# Single-linkage threshold, in allele differences over the 629-locus scheme,
# below which two genomes are treated as the same clone. 10 is the threshold in
# common use for K. pneumoniae cgMLST; --threshold reports the sensitivity.
DEFAULT_THRESHOLD = 10
MISSING = -1


def has_acquired_carbapenemase(value):
    return str(value or "").strip() not in ("", "-", "nan")


def gene_families(value):
    families = sorted({g.strip().split("-")[0] for g in str(value).split(";") if g.strip()})
    return "+".join(families)


def load_allele_matrix(pathogenwatch_dir):
    """Return (genome_ids, integer allele matrix with MISSING for uncalled loci)."""
    cgmlst = pd.read_csv(
        os.path.join(pathogenwatch_dir, "pathogenwatch-cgmlst.csv"),
        usecols=["Genome ID", "Gene", "Allele ID"], dtype=str, low_memory=False,
    )
    wide = cgmlst.pivot_table(
        index="Genome ID", columns="Gene", values="Allele ID", aggfunc="first",
    )
    genomes = wide.index.to_numpy()
    # Factorize per locus: allele identity only has to be comparable within a
    # locus, so per-column integer codes are enough and keep the matrix small.
    codes = np.full(wide.shape, MISSING, dtype=np.int32)
    for j, locus in enumerate(wide.columns):
        column = wide[locus]
        called = column.notna().to_numpy()
        codes[called, j] = pd.factorize(column[called])[0]
    return genomes, codes


def pairwise_distances(codes):
    """Scaled allele distance and shared-locus count for every pair of genomes."""
    n, n_loci = codes.shape
    distance = np.zeros((n, n), dtype=np.float64)
    shared = np.zeros((n, n), dtype=np.int32)
    for i in range(n - 1):
        rest = codes[i + 1:]
        both_called = (codes[i] != MISSING) & (rest != MISSING)
        mismatches = ((codes[i] != rest) & both_called).sum(axis=1)
        comparable = both_called.sum(axis=1)
        # Rescale to the full scheme so a pair sharing fewer called loci is not
        # credited with artificially few differences.
        scaled = np.where(comparable > 0, mismatches * (n_loci / np.maximum(comparable, 1)), np.nan)
        distance[i, i + 1:] = distance[i + 1:, i] = scaled
        shared[i, i + 1:] = shared[i + 1:, i] = comparable
    return distance, shared


def single_linkage(distance_block, threshold):
    """Cluster labels under single linkage at `threshold` allele differences."""
    n = len(distance_block)
    label = list(range(n))

    def find(x):
        while label[x] != x:
            label[x] = label[label[x]]
            x = label[x]
        return x

    for i in range(n):
        for j in range(i + 1, n):
            if distance_block[i, j] <= threshold:
                ri, rj = find(i), find(j)
                if ri != rj:
                    label[rj] = ri
    return [find(i) for i in range(n)]


def load_typing(pathogenwatch_dir, index):
    kleborate = pd.read_csv(
        os.path.join(pathogenwatch_dir, "pathogenwatch-kleborate.csv"),
        encoding="utf-8-sig", dtype=str,
    )
    kleborate["positive"] = kleborate["Bla_Carb_acquired"].map(has_acquired_carbapenemase)
    kleborate["family"] = kleborate["Bla_Carb_acquired"].where(kleborate["positive"], "").map(
        lambda v: gene_families(v) if str(v).strip() else "")

    # Contributing study travels with every genome: a tight cluster confined to
    # one BioProject is an outbreak that one study sampled, not evidence that
    # the clone is widespread. See DATA_PROVENANCE.md.
    metadata = pd.read_csv(
        os.path.join(pathogenwatch_dir, "pathogenwatch-metadata.csv"),
        encoding="utf-8-sig", dtype=str,
    )
    kleborate = kleborate.merge(
        metadata[["Genome ID", "INSDC Study Accession"]], on="Genome ID", how="left")
    kleborate["study"] = kleborate["INSDC Study Accession"].fillna("(unrecorded)").str.strip()

    kleborate = kleborate[kleborate["Genome ID"].isin(index)].copy()
    kleborate["row"] = kleborate["Genome ID"].map(index)
    return kleborate


def summarize(distance, kleborate, threshold):
    """Per-ST clonality summary at one single-linkage threshold."""
    rows = []
    for st, group in kleborate.groupby("ST"):
        positives = group[group["positive"]]
        if len(positives) == 0:
            continue
        pos_rows = positives["row"].to_numpy()
        block = distance[np.ix_(pos_rows, pos_rows)]
        off_diagonal = block[np.triu_indices(len(pos_rows), k=1)]
        families = sorted(set(positives["family"]))

        if len(pos_rows) > 1:
            labels = single_linkage(block, threshold)
            n_clusters = len(set(labels))
            sizes = pd.Series(labels).value_counts()
            largest = int(sizes.iloc[0])
        else:
            labels, n_clusters, largest = [0], 1, 1

        # Minimum independent acquisitions. Host clonality alone understates it:
        # two near-identical genomes carrying different carbapenemase families
        # acquired the phenotype twice, and one family appearing in two distant
        # clusters is also two acquisitions. Counting distinct (clone, family)
        # pairs captures both.
        acquisitions = len(set(zip(labels, positives["family"])))

        studies = positives["study"].value_counts()

        # How far the positives sit from same-ST negatives: if positives are a
        # distinct subclade this is large relative to the within-positive spread.
        negatives = group[~group["positive"]]
        if len(negatives) and len(pos_rows):
            cross = distance[np.ix_(pos_rows, negatives["row"].to_numpy())]
            median_to_negative = float(np.nanmedian(cross))
        else:
            median_to_negative = np.nan

        rows.append({
            "st": st,
            "n_genomes": len(group),
            "n_positive": len(positives),
            "families": "|".join(f for f in families if f),
            "alleles": "|".join(sorted(set(positives["Bla_Carb_acquired"].str.strip()))),
            "median_within_positive": round(float(np.nanmedian(off_diagonal)), 1)
                                      if len(off_diagonal) else np.nan,
            "min_within_positive": round(float(np.nanmin(off_diagonal)), 1)
                                   if len(off_diagonal) else np.nan,
            "max_within_positive": round(float(np.nanmax(off_diagonal)), 1)
                                   if len(off_diagonal) else np.nan,
            "median_positive_to_negative": round(median_to_negative, 1)
                                           if not np.isnan(median_to_negative) else "",
            f"clusters_at_{threshold}": n_clusters,
            "largest_cluster": largest,
            "min_acquisitions": acquisitions,
            "n_studies": len(studies),
            "studies": "; ".join(f"{s}={n}" for s, n in studies.items()),
            "single_study_clone": bool(n_clusters == 1 and len(studies) == 1 and len(pos_rows) > 1),
        })

    return pd.DataFrame(rows).sort_values("n_positive", ascending=False).reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pathogenwatch_dir",
                        default="resistome/results/pathogenwatch_west_africa")
    parser.add_argument("--out", default="resistome/results/lineage_structure.tsv")
    parser.add_argument("--threshold", type=int, default=DEFAULT_THRESHOLD,
                        help="allele differences below which two genomes are one clone")
    args = parser.parse_args()

    genomes, codes = load_allele_matrix(args.pathogenwatch_dir)
    print(f"allele matrix: {codes.shape[0]} genomes x {codes.shape[1]} loci")
    distance, shared = pairwise_distances(codes)
    print(f"pairwise distances computed; median shared loci "
          f"{int(np.median(shared[np.triu_indices_from(shared, k=1)]))}")

    index = {g: i for i, g in enumerate(genomes)}
    kleborate = load_typing(args.pathogenwatch_dir, index)

    profile = summarize(distance, kleborate, args.threshold)
    profile.to_csv(args.out, sep="\t", index=False)

    cluster_col = f"clusters_at_{args.threshold}"
    multi = profile[profile["n_positive"] >= 2]
    print(f"\n--- STs with >=2 carbapenemase-positive genomes "
          f"(threshold {args.threshold} alleles) ---")
    print(f"{'ST':<10}{'n+':>4}{'/n':>5}  {'families':<9}{'within+':>9}{'to neg':>8}"
          f"{'clust':>7}{'acq':>5}{'studies':>9}  note")
    for _, r in multi.iterrows():
        note = "single-study clone" if r["single_study_clone"] else ""
        if r["min_acquisitions"] > r[cluster_col]:
            note = (note + " " if note else "") + "2 families in one clone"
        print(f"{r['st']:<10}{r['n_positive']:>4}{r['n_genomes']:>5}  {r['families']:<9}"
              f"{r['median_within_positive']:>9}{str(r['median_positive_to_negative']):>8}"
              f"{r[cluster_col]:>7}{r['min_acquisitions']:>5}{r['n_studies']:>9}  {note}")

    print(f"\n{int((multi['min_acquisitions'] == 1).sum())} STs consistent with a single "
          f"acquisition, {int((multi['min_acquisitions'] > 1).sum())} requiring more than one.")
    print(f"total minimum independent acquisitions across all STs: "
          f"{int(profile['min_acquisitions'].sum())} for {int(profile['n_positive'].sum())} "
          f"positive genomes")
    confined = multi[multi["single_study_clone"]]
    print(f"\n{len(confined)} of the single-clone STs are confined to one BioProject "
          f"({', '.join(confined['st'])}) —")
    print("  an outbreak that one study sampled cannot be distinguished from a widespread clone.")
    print(f"singleton-positive STs (uninformative on clonality): "
          f"{int((profile['n_positive'] == 1).sum())}")

    print("\n--- threshold sensitivity ---")
    for t in sorted({5, 10, 20, 50, args.threshold}):
        p = profile if t == args.threshold else summarize(distance, kleborate, t)
        m = p[p["n_positive"] >= 2]
        print(f"  <= {t:>3} alleles: {int((m[f'clusters_at_{t}'] == 1).sum())}/{len(m)} STs "
              f"single-clone, {int(p['min_acquisitions'].sum())} total acquisitions")

    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
