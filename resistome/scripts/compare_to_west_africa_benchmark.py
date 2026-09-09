#!/usr/bin/env python3
"""
Compare this study's assembled genomes against the public West Africa
Pathogenwatch clinical benchmark (Nigeria + Ghana + Benin, n=430).

Two things make the comparison non-trivial and are handled explicitly here:

1. Deduplication. resistome/raw/ holds several Ghana genomes under BOTH their
   GenBank (GCA_) and RefSeq (GCF_) accessions. These are the same biological
   assembly: all 11 duplicated pairs give identical MLST allele profiles and
   identical AMRFinder gene calls. Counting both would inflate n and
   double-weight those isolates, so only the GCA representative is kept.

2. ESBL / carbapenemase calling. Status is taken from AMRFinderPlus's
   structured Class/Subclass fields (Type=AMR, Class=BETA-LACTAM, then
   Subclass CEPHALOSPORIN => ESBL, CARBAPENEM => carbapenemase). Do NOT
   substitute a text search over RGI output: "NDM" occurs by chance in raw
   protein-sequence columns, and RGI's carbapenem drug-class tag also fires on
   core genes every K. pneumoniae carries (LptD, OmpK porins), neither of
   which is an acquired carbapenemase.

Tools: AMRFinderPlus 4.2.7 (db 2026-05-15.1), mlst 2.11 (PubMLST snapshot,
kpneumoniae scheme). Benchmark: Pathogenwatch export, Kleborate calls.

Usage (from project root):
    python resistome/scripts/compare_to_west_africa_benchmark.py
"""
import csv
import glob
import os
from collections import Counter

RESULTS = "resistome/results"
MLST_CALLS = f"{RESULTS}/mlst_st_calls.tsv"
BENCHMARK = f"{RESULTS}/pathogenwatch_west_africa/pathogenwatch_west_africa_clinical_only.csv"
OUT = f"{RESULTS}/west_africa_benchmark_comparison.csv"

# Genomes excluded from the primary figures on quality/identity grounds; see
# the notes emitted at the bottom of this script's output.
SUSPECT = {
    # Bimodal GC: ~1.4 Mbp of contigs at 50-52% GC alongside ~2.6 Mbp at
    # 56-58% GC, i.e. two organisms co-assembled. mlst finds gapA/pgi/phoE at
    # ~100% but rpoB at only 95.2% and infB/mdh/tonB absent entirely; mob_recon
    # mash assigns some contigs to K. pneumoniae and others to E. coli.
    # NB: mlst prints every locus as "-" here, which looks like "no loci
    # found" but is really mlst blanking a partial result (SCORE=47,
    # 4/-/-/52/1/~146/-) that fell below its acceptance threshold.
    "ghana_GCA_025660315.1": "contaminated: bimodal GC (50-52% + 56-58%), 3/7 MLST loci absent, rpoB 95.2%",
    # GC smeared flat across 45-62% with no dominant peak, 2594 contigs,
    # 6.08 Mbp (K. pneumoniae is ~5.3 Mbp): multi-organism contamination.
    "ghana_GCA_022968745.1": "contaminated: no dominant GC peak, 2594 contigs, 6.08 Mbp (oversized)",
}


def beta_lactam_status(path):
    """Return (esbl_genes, carbapenemase_genes) from one AMRFinder TSV."""
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    def genes(subclass):
        return sorted({
            r["Element symbol"] for r in rows
            if r["Type"] == "AMR" and r["Class"] == "BETA-LACTAM"
            and subclass in r["Subclass"]
        })
    return genes("CEPHALOSPORIN"), genes("CARBAPENEM")


def main():
    amr = {
        os.path.basename(p).replace("_amrfinder.tsv", ""): beta_lactam_status(p)
        for p in glob.glob(f"{RESULTS}/*_amrfinder.tsv")
    }
    st = {
        r["sample_id"]: r["ST"]
        for r in csv.DictReader(open(MLST_CALLS), delimiter="\t")
    }

    # Drop RefSeq duplicates of GenBank assemblies, and the public reference
    # genome (kpn), which is not part of this study's isolate collection.
    distinct = sorted(
        s for s in amr if not s.startswith("ghana_GCF_") and s != "kpn"
    )
    clean = [s for s in distinct if s not in SUSPECT]

    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample_id", "st", "esbl_genes", "carbapenemase_genes",
                    "esbl_positive", "carbapenemase_positive", "quality_flag"])
        for s in distinct:
            esbl, carb = amr[s]
            w.writerow([s, st.get(s, ""), ";".join(esbl), ";".join(carb),
                        bool(esbl), bool(carb), SUSPECT.get(s, "")])

    bench = list(csv.DictReader(open(BENCHMARK)))
    bn = len(bench)
    b_esbl = sum(1 for r in bench if r["bla_esbl_acquired"].strip() not in ("", "-"))
    b_carb = sum(1 for r in bench if r["bla_carb_acquired"].strip() not in ("", "-"))

    print(f"Benchmark (Pathogenwatch clinical): n={bn}  "
          f"ESBL+ {b_esbl/bn:.1%}  carbapenemase+ {b_carb/bn:.1%}\n")

    for label, gs in [("This study (all distinct)", distinct),
                      ("This study (quality-filtered)", clean)]:
        n = len(gs)
        e = sum(1 for s in gs if amr[s][0])
        c = sum(1 for s in gs if amr[s][1])
        print(f"{label:<32} n={n:<4} ESBL+ {e/n:6.1%}   carbapenemase+ {c/n:.1%}")

    print("\nST distribution (this study, all distinct):")
    for s_, cnt in Counter(st.get(s, "?") for s in distinct).most_common():
        print(f"  {'untypeable' if s_ == '-' else 'ST' + s_:<12} {cnt}")

    # Pathogenwatch/Kleborate writes STs as "ST17"; mlst writes bare "17".
    # Normalise to bare digits on both sides before comparing.
    def bare(x):
        return x[2:] if x.startswith("ST") else x

    print("\nBenchmark top-8 clones vs this study:")
    b_st = Counter(bare(r["st"]) for r in bench)
    mine = Counter(st.get(s) for s in distinct)
    for s_, cnt in b_st.most_common(8):
        print(f"  ST{s_:<6} benchmark n={cnt:<4} this study n={mine.get(s_, 0)}")

    print("\nQuality-suspect genomes excluded from the filtered figures:")
    for s, why in SUSPECT.items():
        print(f"  {s}: {why}")
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
