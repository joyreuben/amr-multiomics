#!/usr/bin/env python3
"""
Build the download manifest for objective 2 (mobilization): the carbapenemase-
positive subset of the West Africa Pathogenwatch collection, with the ENA read
files needed to re-assemble each one.

Objective 2 asks whether each acquired carbapenemase sits on a plasmid contig
and under which replicon. The Pathogenwatch export cannot answer that: it gives
replicon-to-contig (pathogenwatch-plasmidfinder.csv) but Kleborate reports no
contig coordinates for AMR genes, so there is nothing to join on. That requires
the assemblies themselves, which requires the reads.

Only carbapenemase-positive genomes are in scope, so this is 99 genomes rather
than all 557.

Sources (Pathogenwatch West Africa export, downloaded 2026-09):
    pathogenwatch-metadata.csv   country, date, INSDC run/sample/study accessions
    pathogenwatch-kleborate.csv  Kleborate v3.2.4, ST and Bla_Carb_acquired
    ENA portal API               run_accession -> FASTQ URLs, one query per study

Usage:
    python resistome/scripts/build_carbapenemase_cohort.py \
        --pathogenwatch_dir resistome/results/pathogenwatch_west_africa \
        --out resistome/results/carbapenemase_cohort.tsv
"""
import argparse
import io
import os
import sys
import time
import urllib.request

import pandas as pd

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FIELDS = "run_accession,sample_accession,fastq_ftp,fastq_bytes,read_count"


def has_acquired_carbapenemase(value):
    return str(value or "").strip() not in ("", "-", "nan")


def gene_families(value):
    """Collapse allele calls to families: NDM-1;OXA-48 -> NDM+OXA.

    Objective 2 is asked per family, not per allele: NDM and OXA-48-like are
    two independent gene flows expected on different vehicles (OXA-181 on
    ColKP3/IncX3, NDM on IncF/IncC). See claude.md.
    """
    families = sorted({g.strip().split("-")[0] for g in str(value).split(";") if g.strip()})
    return "+".join(families)


def fetch_study_runs(study, retries=3):
    """Return the ENA read_run table for one study accession."""
    url = f"{ENA_PORTAL}?accession={study}&result=read_run&fields={ENA_FIELDS}&format=tsv"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                text = response.read().decode("utf-8")
            if not text.strip():
                return pd.DataFrame()
            return pd.read_csv(io.StringIO(text), sep="\t", dtype=str)
        except Exception as exc:  # noqa: BLE001 - report and retry, ENA rate-limits
            if attempt == retries - 1:
                print(f"  !!! {study}: {exc}", file=sys.stderr)
                return pd.DataFrame()
            time.sleep(3 * (attempt + 1))
    return pd.DataFrame()


def build_cohort(pathogenwatch_dir):
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
    positives = df[df["Bla_Carb_acquired"].map(has_acquired_carbapenemase)].copy()
    positives["carbapenemase"] = positives["Bla_Carb_acquired"].str.strip()
    positives["family"] = positives["carbapenemase"].map(gene_families)
    positives["study"] = positives["INSDC Study Accession"].fillna("(unrecorded)").str.strip()

    studies = sorted(positives["study"].unique())
    print(f"{len(positives)} carbapenemase-positive genomes across {len(studies)} studies")

    runs = []
    for study in studies:
        table = fetch_study_runs(study)
        print(f"  {study:<14} ENA returned {len(table):>4} runs")
        if len(table):
            runs.append(table)
        time.sleep(0.5)

    if not runs:
        raise SystemExit("ENA returned nothing for any study; cannot build manifest.")

    ena = pd.concat(runs, ignore_index=True).drop_duplicates("run_accession")
    cohort = positives.merge(
        ena, left_on="INSDC Run Accession", right_on="run_accession", how="left",
    )

    cohort["fastq_1"] = cohort["fastq_ftp"].fillna("").str.split(";").str[0]
    cohort["fastq_2"] = cohort["fastq_ftp"].fillna("").str.split(";").str[-1]
    cohort["paired"] = cohort["fastq_1"].ne("") & cohort["fastq_2"].ne("") \
        & cohort["fastq_1"].ne(cohort["fastq_2"])

    columns = [
        "Genome ID", "INSDC Run Accession", "INSDC Sample Accession", "study",
        "ST", "carbapenemase", "family", "Country", "Date",
        "fastq_1", "fastq_2", "fastq_bytes", "read_count", "paired",
    ]
    return cohort[columns].rename(columns={
        "Genome ID": "genome_id",
        "INSDC Run Accession": "run_accession",
        "INSDC Sample Accession": "sample_accession",
        "Country": "country",
        "Date": "date",
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pathogenwatch_dir",
                        default="resistome/results/pathogenwatch_west_africa")
    parser.add_argument("--out",
                        default="resistome/results/carbapenemase_cohort.tsv")
    args = parser.parse_args()

    cohort = build_cohort(args.pathogenwatch_dir)
    cohort.to_csv(args.out, sep="\t", index=False)

    resolved = int(cohort["paired"].sum())
    print(f"\npaired FASTQ resolved for {resolved}/{len(cohort)} genomes")
    if resolved < len(cohort):
        missing = cohort[~cohort["paired"]][["genome_id", "run_accession", "study"]]
        print("unresolved:")
        print(missing.to_string(index=False))

    total_gb = sum(
        int(b) for row in cohort["fastq_bytes"].dropna()
        for b in str(row).split(";") if b.isdigit()
    ) / 1e9
    print(f"\ntotal read download: {total_gb:.1f} GB")
    print("\nby carbapenemase family:")
    print(cohort["family"].value_counts().to_string())
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
