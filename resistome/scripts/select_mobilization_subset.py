#!/usr/bin/env python3
"""
Choose the ~20 genomes worth assembling for objective 2, instead of all 99.

Assembling the full cohort costs about three days and mostly buys prevalence,
which is not what objective 2 needs. The claim that needs physical evidence is
narrower: that the same vehicle carries a carbapenemase into lineages too
distant to be explained by clonal descent. That is answered by genomes chosen
at the points where the argument turns, not by every positive genome.

Each pick has a stated reason, in the `rationale` column of the output:

  ST17   3  the monophyletic 28-genome OXA-181 outbreak — is one plasmid behind it?
  ST234  2  OXA-181 in an unrelated lineage — same plasmid as ST17, or convergent?
  ST147  4  NDM-1 arriving three separate times in one ST — same plasmid each time?
  ST464  2  NDM-7 rather than NDM-1, so a different allele on possibly the same vehicle
  ST392  2  one ST carrying both families
  ST307  2  both families again, in a globally important lineage
  ST340  2  two genomes 4 alleles apart carrying different families — the single
            most informative pair in the collection
  ST395  2  a positive and the negative sitting 0 alleles from it — did the
            plasmid leave, or was the gene never there?

Where an ST contributes several picks they are drawn from different cgMLST
clusters where clusters exist, so the picks span the ST's diversity rather than
resampling one clone.

The output is written in the same columns as carbapenemase_cohort.tsv, so it
drops straight into the assembler:

    MANIFEST=resistome/results/mobilization_subset.tsv \
        bash resistome/scripts/assemble_carbapenemase_cohort.sh

Usage:
    python resistome/scripts/select_mobilization_subset.py
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lineage_structure_carbapenemase import (  # noqa: E402
    load_allele_matrix, pairwise_distances, load_typing, single_linkage,
)
from build_carbapenemase_cohort import fetch_study_runs  # noqa: E402

CLUSTER_THRESHOLD = 10

# (ST, family filter or None for any, how many, why)
TARGETS = [
    ("ST17",  "OXA", 3, "monophyletic 28-genome OXA-181 outbreak: one plasmid behind it?"),
    ("ST234", "OXA", 2, "OXA-181 in an unrelated lineage: same plasmid as ST17?"),
    ("ST147", "NDM", 4, "NDM-1 acquired 3 separate times in one ST: same plasmid each time?"),
    ("ST464", "NDM", 2, "NDM-7 not NDM-1: different allele, possibly same vehicle"),
    ("ST392", None,  2, "one ST carrying both carbapenemase families"),
    ("ST307", None,  2, "both families in a globally important lineage"),
    ("ST340", None,  2, "two genomes 4 alleles apart with different families: the key pair"),
]
# Handled separately: needs the negative, which is not in the positives-only cohort.
PAIRED_NEGATIVE_ST = "ST395"


def pick_spanning_clusters(candidates, distance, count):
    """Take `count` genomes from `candidates`, spreading across cgMLST clusters."""
    rows = candidates["row"].to_numpy()
    if len(rows) <= count:
        return list(candidates["Genome ID"])
    block = distance[np.ix_(rows, rows)]
    labels = single_linkage(block, CLUSTER_THRESHOLD)
    candidates = candidates.assign(cluster=labels)

    picked, seen = [], set()
    # One from each distinct cluster first, largest clusters first.
    for cluster in candidates["cluster"].value_counts().index:
        if len(picked) >= count:
            break
        member = candidates[candidates["cluster"] == cluster].iloc[0]
        picked.append(member["Genome ID"])
        seen.add(cluster)
    # Then fill from the largest clusters.
    for _, row in candidates.iterrows():
        if len(picked) >= count:
            break
        if row["Genome ID"] not in picked:
            picked.append(row["Genome ID"])
    return picked[:count]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pathogenwatch_dir",
                        default="resistome/results/pathogenwatch_west_africa")
    parser.add_argument("--cohort", default="resistome/results/carbapenemase_cohort.tsv")
    parser.add_argument("--out", default="resistome/results/mobilization_subset.tsv")
    parser.add_argument("--cache", default="resistome/results/cgmlst_distance.npy")
    args = parser.parse_args()

    if os.path.exists(args.cache):
        cached = np.load(args.cache, allow_pickle=True).item()
        genomes, distance = cached["genomes"], cached["distance"]
    else:
        genomes, codes = load_allele_matrix(args.pathogenwatch_dir)
        distance, _ = pairwise_distances(codes)
    index = {g: i for i, g in enumerate(genomes)}
    typing = load_typing(args.pathogenwatch_dir, index)

    selected = []
    for st, family, count, rationale in TARGETS:
        pool = typing[(typing["ST"] == st) & typing["positive"]]
        if family:
            pool = pool[pool["family"] == family]
        if pool.empty:
            print(f"!!! {st}: no candidates", file=sys.stderr)
            continue

        # Where the point of the ST is that it carries both families, take one
        # of each before spreading across clusters — otherwise cluster spanning
        # can quietly draw two of the same family and miss the comparison.
        families = [f for f in pool["family"].unique() if f]
        if family is None and len(families) > 1:
            chosen, per_family = [], max(1, count // len(families))
            for fam in families:
                subset = pool[pool["family"] == fam]
                chosen += pick_spanning_clusters(subset, distance, per_family)
            for genome in pick_spanning_clusters(pool, distance, count):
                if len(chosen) >= count:
                    break
                if genome not in chosen:
                    chosen.append(genome)
            picks_for_st = chosen[:count]
        else:
            picks_for_st = pick_spanning_clusters(pool, distance, count)

        for genome in picks_for_st:
            selected.append({"Genome ID": genome, "rationale": f"{st}: {rationale}"})

    # ST395: one positive and the negative that sits 0 alleles from it.
    st395 = typing[typing["ST"] == PAIRED_NEGATIVE_ST]
    positive = st395[st395["positive"]]
    negative = st395[~st395["positive"]]
    for frame, kind in ((positive, "positive"), (negative, "negative")):
        if len(frame):
            selected.append({
                "Genome ID": frame.iloc[0]["Genome ID"],
                "rationale": f"{PAIRED_NEGATIVE_ST}: {kind} of a pair 0 alleles apart — "
                             "did the plasmid leave, or was the gene never there?",
            })

    picks = pd.DataFrame(selected).drop_duplicates("Genome ID")
    print(f"selected {len(picks)} genomes")

    cohort = pd.read_csv(args.cohort, sep="\t", dtype=str)
    merged = picks.merge(cohort, left_on="Genome ID", right_on="genome_id", how="left")

    # Anything not in the positives-only cohort (the ST395 negative) needs its
    # own ENA lookup.
    missing = merged[merged["run_accession"].isna()]
    if len(missing):
        metadata = pd.read_csv(
            os.path.join(args.pathogenwatch_dir, "pathogenwatch-metadata.csv"),
            dtype=str, encoding="utf-8-sig")
        extra = metadata[metadata["Genome ID"].isin(missing["Genome ID"])]
        print(f"looking up ENA reads for {len(extra)} genome(s) outside the cohort")
        runs = pd.concat(
            [fetch_study_runs(s) for s in extra["INSDC Study Accession"].dropna().unique()],
            ignore_index=True)
        extra = extra.merge(runs, left_on="INSDC Run Accession",
                            right_on="run_accession", how="left")
        typed = typing.set_index("Genome ID")
        for _, row in extra.iterrows():
            urls = str(row.get("fastq_ftp") or "").split(";")
            mask = merged["Genome ID"] == row["Genome ID"]
            merged.loc[mask, "genome_id"] = row["Genome ID"]
            merged.loc[mask, "run_accession"] = row["INSDC Run Accession"]
            merged.loc[mask, "sample_accession"] = row["INSDC Sample Accession"]
            merged.loc[mask, "study"] = row["INSDC Study Accession"]
            merged.loc[mask, "ST"] = typed.loc[row["Genome ID"], "ST"]
            merged.loc[mask, "carbapenemase"] = ""
            merged.loc[mask, "family"] = "negative"
            merged.loc[mask, "country"] = row["Country"]
            merged.loc[mask, "fastq_1"] = urls[0] if urls else ""
            merged.loc[mask, "fastq_2"] = urls[-1] if len(urls) > 1 else ""
            merged.loc[mask, "paired"] = bool(len(urls) > 1 and urls[0] != urls[-1])

    columns = ["genome_id", "run_accession", "sample_accession", "study", "ST",
               "carbapenemase", "family", "country", "date", "fastq_1", "fastq_2",
               "fastq_bytes", "read_count", "paired", "rationale"]
    for column in columns:
        if column not in merged.columns:
            merged[column] = ""
    out = merged[columns]
    out.to_csv(args.out, sep="\t", index=False)

    resolved = int(out["paired"].astype(str).str.lower().eq("true").sum())
    print(f"paired FASTQ resolved for {resolved}/{len(out)}")
    print("\nby ST:")
    print(out.groupby(["ST", "family"]).size().to_string())
    gb = sum(int(b) for row in out["fastq_bytes"].dropna()
             for b in str(row).split(";") if b.isdigit()) / 1e9
    print(f"\ndownload {gb:.1f} GB, roughly {len(out) * 45 / 60:.0f} h of assembly "
          f"at the measured 45 min/genome")
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
