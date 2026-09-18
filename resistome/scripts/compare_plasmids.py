#!/usr/bin/env python3
"""
Turn "same mob_suite cluster" into "same plasmid", or show that it is not.

analyze_mobilization.py found one cluster carrying a carbapenemase across
unrelated lineages — AA038 with OXA-181 in ST17, ST234 and ST340, AA405 with
NDM-1 in four STs. A cluster is a mash-similarity group against a reference
database, not proof that two genomes carry the same element, and AA038 also
carries NDM-7 in ST464, which is what a loose grouping looks like. The claim
needs direct sequence comparison.

Each pair of reconstructed plasmids is compared two ways:

  mash distance  whole-sequence similarity, cheap, gives an ANI estimate
  blastn         aligned fraction of each sequence, which mash cannot give

Both are needed. Two plasmids sharing a transposon score well on identity over
the aligned part while sharing very little of their length, so identity without
coverage would call any pair of plasmids with a common mobile element "the
same". The convention applied here is >=95% identity AND >=80% reciprocal
coverage; the thresholds are reported alongside the raw numbers so a reader can
apply their own.

The comparison that matters is within-ST versus cross-ST. Plasmids from one
clonal outbreak should be near-identical simply by descent; the question is
whether cross-lineage pairs reach the same similarity, because that is what
distinguishes one plasmid moving between lineages from two related plasmids
arriving separately.

Usage:
    python resistome/scripts/compare_plasmids.py --cluster AA038
    python resistome/scripts/compare_plasmids.py          # every shared cluster
"""
import argparse
import glob
import itertools
import os
import re
import subprocess
import tempfile

import pandas as pd

IDENTITY_THRESHOLD = 95.0
COVERAGE_THRESHOLD = 80.0


def sequence_lengths(path):
    lengths, name, total = {}, None, 0
    for line in open(path):
        if line.startswith(">"):
            if name:
                lengths[name] = total
            name, total = line[1:].split()[0], 0
        else:
            total += len(line.strip())
    if name:
        lengths[name] = total
    return lengths


def mash_distance(path_a, path_b):
    """Whole-sequence distance; ANI estimate is (1 - distance)."""
    try:
        out = subprocess.run(["mash", "dist", path_a, path_b],
                             capture_output=True, text=True, timeout=300, check=True)
        fields = out.stdout.split("\t")
        return float(fields[2]), float(fields[3])
    except Exception:  # noqa: BLE001 - a failed pair should not stop the sweep
        return float("nan"), float("nan")


def merge_intervals(intervals):
    if not intervals:
        return 0
    intervals.sort()
    total, current_start, current_end = 0, *intervals[0]
    for start, end in intervals[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start + 1
            current_start, current_end = start, end
    return total + current_end - current_start + 1


def blast_coverage(query, subject):
    """Aligned fraction of `query` against `subject`, and weighted identity."""
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "db")
        subprocess.run(["makeblastdb", "-in", subject, "-dbtype", "nucl", "-out", db],
                       capture_output=True, check=True)
        out = subprocess.run(
            ["blastn", "-query", query, "-db", db, "-outfmt",
             "6 qseqid qstart qend pident length", "-evalue", "1e-10",
             "-perc_identity", "80"],
            capture_output=True, text=True, timeout=900, check=True).stdout

    per_contig, weighted, aligned_total = {}, 0.0, 0
    for line in out.strip().splitlines():
        if not line:
            continue
        qid, qstart, qend, pident, length = line.split("\t")
        start, end = sorted((int(qstart), int(qend)))
        per_contig.setdefault(qid, []).append((start, end))
        weighted += float(pident) * int(length)
        aligned_total += int(length)

    lengths = sequence_lengths(query)
    covered = sum(merge_intervals(v) for v in per_contig.values())
    total = sum(lengths.values()) or 1
    identity = weighted / aligned_total if aligned_total else float("nan")
    return 100.0 * covered / total, identity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results_dir", default="resistome/results")
    parser.add_argument("--placement",
                        default="resistome/results/mobilization_placement.tsv")
    parser.add_argument("--cluster", help="restrict to one cluster, e.g. AA038")
    parser.add_argument("--out", default="resistome/results/plasmid_comparison.tsv")
    args = parser.parse_args()

    placement = pd.read_csv(args.placement, sep="\t", dtype=str)
    placement = placement[placement["molecule_type"] == "plasmid"]
    carrying = placement.groupby("primary_cluster_id")["sample"].apply(set).to_dict()
    st_of = dict(zip(placement["sample"], placement["ST"]))
    gene_of = dict(zip(placement["sample"], placement["gene"]))

    clusters = [args.cluster] if args.cluster else sorted(carrying)
    rows = []
    for cluster in clusters:
        files = {}
        for path in sorted(glob.glob(
                os.path.join(args.results_dir, f"pw_*_mob_recon/plasmid_{cluster}.fasta"))):
            sample = re.search(r"(pw_[^/]+)_mob_recon", path).group(1)
            if sample in carrying.get(cluster, set()):
                files[sample] = path
        if len(files) < 2:
            continue

        print(f"\n=== cluster {cluster}: {len(files)} carbapenemase-carrying genomes ===")
        for a, b in itertools.combinations(sorted(files), 2):
            distance, pvalue = mash_distance(files[a], files[b])
            cov_ab, id_ab = blast_coverage(files[a], files[b])
            cov_ba, id_ba = blast_coverage(files[b], files[a])
            same_st = st_of.get(a) == st_of.get(b)
            identity = (id_ab + id_ba) / 2
            min_cov = min(cov_ab, cov_ba)
            same_plasmid = (identity >= IDENTITY_THRESHOLD and min_cov >= COVERAGE_THRESHOLD)
            rows.append({
                "cluster": cluster, "sample_a": a, "sample_b": b,
                "st_a": st_of.get(a), "st_b": st_of.get(b),
                "gene_a": gene_of.get(a), "gene_b": gene_of.get(b),
                "comparison": "within-ST" if same_st else "cross-ST",
                "mash_distance": round(distance, 5),
                "mash_ani_estimate": round(100 * (1 - distance), 2),
                "blast_identity": round(identity, 2),
                "coverage_a_in_b": round(cov_ab, 1),
                "coverage_b_in_a": round(cov_ba, 1),
                "min_coverage": round(min_cov, 1),
                "same_plasmid": same_plasmid,
            })
            tag = "SAME" if same_plasmid else ""
            print(f"  {a.replace('pw_',''):<14}{st_of.get(a):<8} vs "
                  f"{b.replace('pw_',''):<14}{st_of.get(b):<8} "
                  f"{'within' if same_st else 'CROSS ':<7} "
                  f"ANI~{100*(1-distance):5.1f}  id {identity:5.1f}%  "
                  f"cov {min_cov:5.1f}%  {tag}")

    if not rows:
        raise SystemExit("No cluster had two or more carbapenemase-carrying genomes.")

    table = pd.DataFrame(rows)
    table.to_csv(args.out, sep="\t", index=False)

    print("\n--- within-ST versus cross-ST ---")
    for (cluster, comparison), group in table.groupby(["cluster", "comparison"]):
        print(f"  {cluster} {comparison:<10} n={len(group):<3} "
              f"median identity {group['blast_identity'].median():.1f}%  "
              f"median min-coverage {group['min_coverage'].median():.1f}%  "
              f"same-plasmid {int(group['same_plasmid'].sum())}/{len(group)}")

    cross = table[table["comparison"] == "cross-ST"]
    print(f"\n{int(cross['same_plasmid'].sum())}/{len(cross)} cross-lineage pairs meet "
          f">={IDENTITY_THRESHOLD}% identity and >={COVERAGE_THRESHOLD}% reciprocal coverage.")
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
