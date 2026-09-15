#!/usr/bin/env python3
"""
Objective 3: resistome architecture — which acquired resistance genes travel
with a carbapenemase, and which plasmid replicons they travel with.

Resistance genes cluster on the same mobile element, so one transfer can hand
over several drug classes at once. This asks which genes are carried
disproportionately by carbapenemase-positive genomes, and whether the answer
differs between the two carbapenemase families, which are expected on different
vehicles.

Two guards against reading study composition as biology:

1. Every association is reported twice — a crude collection-wide odds ratio,
   and a Cochran-Mantel-Haenszel odds ratio stratified by contributing
   BioProject. The 557 are 16 studies stacked together with carbapenemase rates
   from 0% to 94% (DATA_PROVENANCE.md), so a gene common in one high-carbapenem
   study would look associated collection-wide without being associated within
   any study. Where the crude and CMH estimates disagree, the CMH one stands.
2. Benjamini-Hochberg correction across the genes tested.

Co-occurrence within a genome is genome-intrinsic and so survives the sampling
frame; how often a pair is *seen* does not, which is what the stratification is
for.

Sources (Pathogenwatch West Africa export, downloaded 2026-09):
    pathogenwatch-kleborate.csv      Kleborate v3.2.4 acquired determinants
    pathogenwatch-plasmidfinder.csv  replicon hits per genome
    pathogenwatch-metadata.csv       INSDC study accession

Usage:
    python resistome/scripts/resistome_architecture.py \
        --pathogenwatch_dir resistome/results/pathogenwatch_west_africa \
        --out resistome/results/resistome_architecture.tsv
"""
import argparse
import os
import re

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

# Kleborate acquired-determinant columns. Chromosomal genes (Bla_chr) and
# mutation columns are excluded: they are intrinsic or point mutations, not
# cargo that can ride a plasmid, so they cannot co-travel with a carbapenemase.
ACQUIRED_COLUMNS = [
    "AGly_acquired", "Col_acquired", "Fcyn_acquired", "Flq_acquired",
    "Gly_acquired", "MLS_acquired", "Phe_acquired", "Rif_acquired",
    "Sul_acquired", "Tet_acquired", "Tgc_acquired", "Tmt_acquired",
    "Bla_acquired", "Bla_inhR_acquired", "Bla_ESBL_acquired",
    "Bla_ESBL_inhR_acquired",
]
CARB_COLUMN = "Bla_Carb_acquired"
MIN_GENOMES = 10        # don't test genes too rare to say anything about
MIN_STUDY_SIZE = 10     # strata smaller than this contribute little to CMH


def clean_gene(token):
    """Strip Kleborate's allele-version and confidence marks: aac(3)-IIa.v1*^ -> aac(3)-IIa."""
    token = token.strip().rstrip("*?^")
    token = re.sub(r"\.v\d+$", "", token)
    return token.rstrip("*?^")


def parse_genes(value):
    if str(value).strip() in ("", "-", "nan"):
        return set()
    return {clean_gene(t) for t in str(value).split(";") if clean_gene(t)}


def gene_families(value):
    return "+".join(sorted({g.split("-")[0] for g in parse_genes(value)}))


def benjamini_hochberg(pvalues):
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    order = np.argsort(p)
    adjusted = np.empty(n, dtype=float)
    previous = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        value = min(previous, p[idx] * n / (n - rank + 1))
        adjusted[idx] = previous = value
    return adjusted


def cmh_odds_ratio(table_by_stratum):
    """Cochran-Mantel-Haenszel common odds ratio across 2x2 strata."""
    numerator = denominator = 0.0
    for a, b, c, d in table_by_stratum:
        n = a + b + c + d
        if n == 0:
            continue
        numerator += a * d / n
        denominator += b * c / n
    if denominator == 0:
        return np.inf if numerator > 0 else np.nan
    return numerator / denominator


def cmh_test(table_by_stratum):
    """CMH chi-square statistic and its p-value, testing a common OR of 1."""
    from scipy.stats import chi2
    observed = expected = variance = 0.0
    for a, b, c, d in table_by_stratum:
        n = a + b + c + d
        if n < 2:
            continue
        row1, row2 = a + b, c + d
        col1, col2 = a + c, b + d
        observed += a
        expected += row1 * col1 / n
        if n > 1:
            variance += row1 * row2 * col1 * col2 / (n * n * (n - 1))
    if variance == 0:
        return np.nan, np.nan
    statistic = (abs(observed - expected) - 0.5) ** 2 / variance
    return statistic, float(chi2.sf(statistic, 1))


def load(pathogenwatch_dir):
    kleborate = pd.read_csv(
        os.path.join(pathogenwatch_dir, "pathogenwatch-kleborate.csv"),
        dtype=str, encoding="utf-8-sig",
    )
    metadata = pd.read_csv(
        os.path.join(pathogenwatch_dir, "pathogenwatch-metadata.csv"),
        dtype=str, encoding="utf-8-sig",
    )
    df = kleborate.merge(
        metadata[["Genome ID", "INSDC Study Accession"]], on="Genome ID", how="left")
    df["study"] = df["INSDC Study Accession"].fillna("(unrecorded)").str.strip()

    df["carb_family"] = df[CARB_COLUMN].map(gene_families)
    df["carb_positive"] = df["carb_family"] != ""

    present = [c for c in ACQUIRED_COLUMNS if c in df.columns]
    df["genes"] = df[present].apply(
        lambda row: set().union(*(parse_genes(v) for v in row)), axis=1)
    return df


def associations(df, label, mask_positive, mask_reference):
    """Test every gene for enrichment in `mask_positive` versus `mask_reference`."""
    universe = mask_positive | mask_reference
    subset = df[universe]
    positive = mask_positive[universe]
    counts = {}
    for genes in subset["genes"]:
        for gene in genes:
            counts[gene] = counts.get(gene, 0) + 1
    testable = sorted(g for g, n in counts.items() if n >= MIN_GENOMES)

    rows = []
    for gene in testable:
        has_gene = subset["genes"].map(lambda s, g=gene: g in s)
        a = int((has_gene & positive).sum())
        b = int((has_gene & ~positive).sum())
        c = int((~has_gene & positive).sum())
        d = int((~has_gene & ~positive).sum())
        if a + c == 0 or b + d == 0:
            continue

        odds, p_crude = fisher_exact([[a, b], [c, d]])

        strata = []
        for _, group in subset.groupby("study"):
            if len(group) < MIN_STUDY_SIZE:
                continue
            gh = group["genes"].map(lambda s, g=gene: g in s)
            gp = positive.loc[group.index]
            strata.append((int((gh & gp).sum()), int((gh & ~gp).sum()),
                           int((~gh & gp).sum()), int((~gh & ~gp).sum())))
        cmh_or = cmh_odds_ratio(strata) if strata else np.nan
        _, p_cmh = cmh_test(strata) if strata else (np.nan, np.nan)

        rows.append({
            "comparison": label,
            "gene": gene,
            "n_positive_with": a, "n_positive": a + c,
            "n_reference_with": b, "n_reference": b + d,
            "pct_positive": round(100 * a / (a + c), 1),
            "pct_reference": round(100 * b / (b + d), 1),
            "crude_or": round(float(odds), 2) if np.isfinite(odds) else np.inf,
            "crude_p": p_crude,
            "cmh_or": round(float(cmh_or), 2) if np.isfinite(cmh_or) else cmh_or,
            "cmh_p": p_cmh,
            "n_strata": len(strata),
        })

    result = pd.DataFrame(rows)
    if len(result):
        result["crude_q"] = benjamini_hochberg(result["crude_p"])
        valid = result["cmh_p"].notna()
        result.loc[valid, "cmh_q"] = benjamini_hochberg(result.loc[valid, "cmh_p"])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pathogenwatch_dir",
                        default="resistome/results/pathogenwatch_west_africa")
    parser.add_argument("--out", default="resistome/results/resistome_architecture.tsv")
    parser.add_argument("--alpha", type=float, default=0.05)
    args = parser.parse_args()

    df = load(args.pathogenwatch_dir)
    negative = ~df["carb_positive"]
    print(f"{len(df)} genomes: {int(df['carb_positive'].sum())} carbapenemase-positive "
          f"({', '.join(f'{k}={v}' for k, v in df[df.carb_positive].carb_family.value_counts().items())})")
    print(f"acquired resistance genes seen in >= {MIN_GENOMES} genomes are tested\n")

    blocks = [
        ("any", df["carb_positive"], negative),
        ("NDM", df["carb_family"].eq("NDM"), negative),
        ("OXA", df["carb_family"].eq("OXA"), negative),
    ]
    tables = []
    for label, positive, reference in blocks:
        table = associations(df, label, positive, reference)
        tables.append(table)

        significant = table[(table["cmh_q"] < args.alpha) & (table["cmh_or"] > 1)]
        significant = significant.sort_values("cmh_or", ascending=False)
        print(f"--- genes enriched with {label} carbapenemase "
              f"(study-stratified, q<{args.alpha}) ---")
        if not len(significant):
            print("  none survive stratification\n")
            continue
        print(f"{'gene':<18}{'in carb+':>10}{'in carb-':>10}{'crude OR':>10}{'CMH OR':>9}{'CMH q':>10}")
        for _, r in significant.head(15).iterrows():
            print(f"{r['gene']:<18}{r['pct_positive']:>9.1f}%{r['pct_reference']:>9.1f}%"
                  f"{r['crude_or']:>10}{r['cmh_or']:>9}{r['cmh_q']:>10.2g}")
        print()

    combined = pd.concat(tables, ignore_index=True)
    combined.to_csv(args.out, sep="\t", index=False)

    # How often does the crude estimate overstate the association?
    tested = combined[combined["cmh_or"].notna() & np.isfinite(combined["cmh_or"])]
    inflated = tested[(tested["crude_q"] < args.alpha) & (tested["cmh_q"] >= args.alpha)]
    print(f"--- sampling-frame check ---")
    print(f"{len(inflated)} of {len(tested)} gene tests are significant collection-wide but "
          f"not after stratifying by study.")
    if len(inflated):
        print("  " + ", ".join(sorted(set(inflated['gene']))[:15]))
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
