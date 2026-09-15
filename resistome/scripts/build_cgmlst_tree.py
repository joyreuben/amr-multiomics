#!/usr/bin/env python3
"""
Build a neighbour-joining tree from the cgMLST allele distances, plus the
annotation table needed to colour it by carbapenemase family.

This is the figure for objective 1. The numbers already say that ST17's
positives are one clone while ST147's are scattered; the tree is what makes
that visible — a tight monophyletic OXA-181 block against NDM-1 tips
interleaved with negatives in several places.

Neighbour joining is implemented here directly rather than via Bio.Phylo: its
DistanceTreeConstructor is pure Python and O(n^3) over 557 taxa is far too slow
on a loaded machine. Each iteration here is vectorised, so the whole tree takes
seconds.

Distances are the same as lineage_structure_carbapenemase.py: allele
differences over loci called in both genomes, rescaled to the 629-locus scheme.
The matrix is cached as .npy so the tree can be rebuilt without recomputing it.

Outputs:
    results/cgmlst_tree.nwk          Newick, tip labels are Genome IDs
    results/cgmlst_tree_annot.tsv    per-tip ST, carbapenemase family and
                                     allele, country, study — load alongside
                                     the tree in iTOL or FigTree

Usage:
    python resistome/scripts/build_cgmlst_tree.py
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lineage_structure_carbapenemase import (  # noqa: E402
    load_allele_matrix, pairwise_distances, gene_families, has_acquired_carbapenemase,
)


def neighbour_joining(distance, labels):
    """Standard neighbour joining, vectorised per iteration. Returns Newick."""
    d = distance.astype(np.float64).copy()
    np.fill_diagonal(d, 0.0)
    # Any pair that shared no loci cannot be placed; fall back to the matrix max
    # so the join is possible but maximally distant.
    if np.isnan(d).any():
        d[np.isnan(d)] = np.nanmax(d)

    n = len(labels)
    nodes = list(labels)
    active = np.ones(n, dtype=bool)

    while active.sum() > 2:
        idx = np.flatnonzero(active)
        m = len(idx)
        sub = d[np.ix_(idx, idx)]
        totals = sub.sum(axis=1)

        # Q_ij = (m-2) d_ij - r_i - r_j, minimised off the diagonal.
        q = (m - 2) * sub - totals[:, None] - totals[None, :]
        np.fill_diagonal(q, np.inf)
        flat = np.argmin(q)
        a, b = divmod(flat, m)
        i, j = idx[a], idx[b]

        # Branch lengths to the new internal node.
        delta = (totals[a] - totals[b]) / (m - 2)
        limb_i = max(0.0, 0.5 * (sub[a, b] + delta))
        limb_j = max(0.0, 0.5 * (sub[a, b] - delta))

        # Distances from the new node to every other active taxon.
        new_row = 0.5 * (d[i, :] + d[j, :] - d[i, j])
        new_row[i] = new_row[j] = 0.0

        nodes[i] = f"({nodes[i]}:{limb_i:.4f},{nodes[j]}:{limb_j:.4f})"
        d[i, :] = new_row
        d[:, i] = new_row
        d[i, i] = 0.0
        active[j] = False

    left, right = np.flatnonzero(active)
    return f"({nodes[left]}:{d[left, right] / 2:.4f},{nodes[right]}:{d[left, right] / 2:.4f});"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pathogenwatch_dir",
                        default="resistome/results/pathogenwatch_west_africa")
    parser.add_argument("--tree", default="resistome/results/cgmlst_tree.nwk")
    parser.add_argument("--annot", default="resistome/results/cgmlst_tree_annot.tsv")
    parser.add_argument("--cache", default="resistome/results/cgmlst_distance.npy")
    args = parser.parse_args()

    if os.path.exists(args.cache):
        cached = np.load(args.cache, allow_pickle=True).item()
        genomes, distance = cached["genomes"], cached["distance"]
        print(f"loaded cached distance matrix: {len(genomes)} genomes")
    else:
        genomes, codes = load_allele_matrix(args.pathogenwatch_dir)
        print(f"allele matrix: {codes.shape[0]} genomes x {codes.shape[1]} loci")
        distance, _ = pairwise_distances(codes)
        np.save(args.cache, {"genomes": genomes, "distance": distance})
        print(f"cached distance matrix to {args.cache}")

    print("building neighbour-joining tree...")
    newick = neighbour_joining(distance, list(genomes))
    with open(args.tree, "w") as handle:
        handle.write(newick + "\n")
    print(f"wrote {args.tree} ({len(newick):,} chars)")

    kleborate = pd.read_csv(
        os.path.join(args.pathogenwatch_dir, "pathogenwatch-kleborate.csv"),
        dtype=str, encoding="utf-8-sig")
    metadata = pd.read_csv(
        os.path.join(args.pathogenwatch_dir, "pathogenwatch-metadata.csv"),
        dtype=str, encoding="utf-8-sig")
    annot = kleborate[["Genome ID", "ST", "Bla_Carb_acquired"]].merge(
        metadata[["Genome ID", "Country", "INSDC Study Accession"]], on="Genome ID", how="left")
    annot["carbapenemase"] = annot["Bla_Carb_acquired"].where(
        annot["Bla_Carb_acquired"].map(has_acquired_carbapenemase), "").str.strip()
    annot["family"] = annot["carbapenemase"].map(
        lambda v: gene_families(v) if v else "negative")
    annot = annot[annot["Genome ID"].isin(set(genomes))]
    annot[["Genome ID", "ST", "family", "carbapenemase", "Country",
           "INSDC Study Accession"]].rename(
        columns={"Genome ID": "genome_id", "ST": "st", "Country": "country",
                 "INSDC Study Accession": "study"}
    ).to_csv(args.annot, sep="\t", index=False)
    print(f"wrote {args.annot}")
    print("\nfamily counts on the tree:")
    print(annot["family"].value_counts().to_string())


if __name__ == "__main__":
    main()
