#!/usr/bin/env python3
"""
Quantify how far the West Africa Pathogenwatch collection is confounded by its
contributing studies, before any country- or year-stratified claim is made from it.

For each contributing BioProject this reports n, acquired-carbapenemase rate,
the years it spans and the countries it covers. If the between-study spread in
carbapenemase rate is wide, a naive "carbapenemase by country" or "by year"
figure is reporting study composition, not epidemiology.

Sources (Pathogenwatch West Africa export, downloaded 2026-09):
    pathogenwatch-metadata.csv   country, collection date, isolation source, study accession
    pathogenwatch-kleborate.csv  Kleborate v3.2.4 (wrapper 6) typing and AMR determinants

Usage:
    python resistome/scripts/profile_west_africa_sampling.py \
        --pathogenwatch_dir resistome/results/pathogenwatch_west_africa \
        --out resistome/results/west_africa_sampling_frame.tsv
"""
import argparse
import os

import pandas as pd

# Pathogenwatch free-texts the country field; these are the same three countries.
COUNTRY_ALIASES = {
    "NG - Nigeria": "Nigeria",
    "NG": "Nigeria",
    "SN": "Senegal",
    "GH": "Ghana",
}


def normalize_country(value):
    value = str(value or "").strip()
    return COUNTRY_ALIASES.get(value, value) or "(unrecorded)"


def collection_year(row):
    """Prefer the curated collection date, fall back to the Pathogenwatch Date field.

    `Collection date` is filled for 222/557 genomes but `Date` for 510/557, so
    falling back roughly doubles the temporally usable set.
    """
    for field in ("Collection date", "Date"):
        value = str(row.get(field) or "").strip()
        if len(value) >= 4 and value[:4].isdigit():
            return int(value[:4])
    return None


def has_acquired_carbapenemase(value):
    value = str(value or "").strip()
    return value not in ("", "-", "nan")


def build_profile(pathogenwatch_dir):
    metadata = pd.read_csv(
        os.path.join(pathogenwatch_dir, "pathogenwatch-metadata.csv"),
        encoding="utf-8-sig", dtype=str,
    )
    kleborate = pd.read_csv(
        os.path.join(pathogenwatch_dir, "pathogenwatch-kleborate.csv"),
        encoding="utf-8-sig", dtype=str,
    )

    df = metadata.merge(
        kleborate[["Genome ID", "ST", "Bla_Carb_acquired", "Bla_ESBL_acquired"]],
        on="Genome ID", how="inner",
    )
    df["study"] = df["INSDC Study Accession"].fillna("(unrecorded)").str.strip()
    df["country"] = df["Country"].map(normalize_country)
    df["year"] = df.apply(collection_year, axis=1)
    df["carbapenemase"] = df["Bla_Carb_acquired"].map(has_acquired_carbapenemase)

    rows = []
    for study, group in df.groupby("study"):
        years = group["year"].dropna().astype(int)
        countries = group["country"].value_counts()
        rows.append({
            "study": study,
            "n": len(group),
            "carbapenemase_n": int(group["carbapenemase"].sum()),
            "carbapenemase_pct": round(100 * group["carbapenemase"].mean(), 1),
            "year_min": int(years.min()) if len(years) else "",
            "year_max": int(years.max()) if len(years) else "",
            "countries": "; ".join(f"{c}={n}" for c, n in countries.items()),
        })

    profile = pd.DataFrame(rows).sort_values("n", ascending=False).reset_index(drop=True)
    return df, profile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pathogenwatch_dir",
                        default="resistome/results/pathogenwatch_west_africa")
    parser.add_argument("--out",
                        default="resistome/results/west_africa_sampling_frame.tsv")
    parser.add_argument("--min_n", type=int, default=15,
                        help="only print studies at least this large (all are written to --out)")
    args = parser.parse_args()

    df, profile = build_profile(args.pathogenwatch_dir)
    profile.to_csv(args.out, sep="\t", index=False)

    print(f"{len(df)} genomes across {len(profile)} contributing studies\n")

    shown = profile[profile["n"] >= args.min_n]
    print(f"--- studies with n >= {args.min_n} ---")
    for _, r in shown.iterrows():
        span = f"{r['year_min']}-{r['year_max']}" if r["year_min"] != "" else "no dates"
        print(f"  {r['study']:<14} n={r['n']:<4} carb={r['carbapenemase_n']:>3} "
              f"({r['carbapenemase_pct']:>5.1f}%)  {span:<10} {r['countries']}")

    rates = shown["carbapenemase_pct"]
    print(f"\nBetween-study carbapenemase rate spans {rates.min():.1f}%-{rates.max():.1f}% "
          f"across studies with n >= {args.min_n}.")

    print("\n--- the same rate, aggregated the naive way ---")
    for key in ("country", "year"):
        print(f"  by {key}:")
        agg = df.dropna(subset=[key]).groupby(key)["carbapenemase"].agg(["sum", "count"])
        for value, r in agg.iterrows():
            label = int(value) if key == "year" else value
            print(f"    {str(label):<12} {int(r['sum']):>3}/{int(r['count']):<4} "
                  f"{100 * r['sum'] / r['count']:>5.1f}%")

    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
