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

2. ESBL / carbapenemase calling. Headline figures come from Kleborate, which
   is the tool the Pathogenwatch benchmark itself used, so both sides of the
   comparison apply the same definitions. AMRFinder is carried alongside as a
   cross-check only. The two agree on 14/15 genomes; they part company on
   ghana_GCA_022968315.1, whose only cephalosporinase is DHA-1. AMRFinder tags
   DHA-1 Subclass=CEPHALOSPORIN, but DHA-1 is an AmpC, not an ESBL, and
   Kleborate correctly keeps it out of Bla_ESBL_acquired. Kleborate likewise
   excludes narrow-spectrum OXA-1. Prefer the Kleborate column when they
   disagree.

   Whichever tool is used, do NOT substitute a text search over RGI output:
   "NDM" occurs by chance in raw protein-sequence columns, and RGI's
   carbapenem drug-class tag also fires on core genes every K. pneumoniae
   carries (LptD, OmpK porins), neither of which is an acquired carbapenemase.

The zero-carbapenemase result is supported by three independent lines of
evidence: Kleborate, AMRFinder, and read-level mapping of the raw trimmed
reads against blaNDM-1/KPC-2/OXA-48/VIM-1/IMP-1 references (which returned
zero hits while the blaCTX-M-15 internal control returned 142 and 136 reads in
the two isolates known to carry it, confirming nothing was lost in assembly).

Tools: Kleborate 3.1.3 (kpsc preset), AMRFinderPlus 4.2.7 (db 2026-05-15.1),
mlst 2.11 (PubMLST snapshot). Benchmark: Pathogenwatch export, Kleborate calls.

Usage (from project root):
    python resistome/scripts/compare_to_west_africa_benchmark.py
"""
import csv
import glob
import os
from collections import Counter

RESULTS = "resistome/results"
MLST_CALLS = f"{RESULTS}/mlst_st_calls.tsv"
KLEBORATE = f"{RESULTS}/kleborate/klebsiella_pneumo_complex_output.txt"
BENCHMARK = f"{RESULTS}/pathogenwatch_west_africa/pathogenwatch_west_africa_clinical_only.csv"
OUT = f"{RESULTS}/west_africa_benchmark_comparison.csv"

# Genomes excluded from the primary figures on quality/identity grounds; see
# the notes emitted at the bottom of this script's output.
SUSPECT = {
    # Missing 3 of 7 core housekeeping loci (infB, mdh, tonB absent; rpoB only
    # an inexact match). mlst and Kleborate independently agree on the same
    # partial profile 4/-/-/52/1/146*/-, so this is real, not a tool artifact.
    # infB is essential, so a complete genome cannot lack it. GC also splits in
    # two: ~1 Mbp of chromosome-assigned sequence sits at 50-52% against the
    # ~57% expected of K. pneumoniae, and mob_recon accounts for only 0.39 Mbp
    # of plasmid, so plasmid content does not explain the low-GC mass.
    # NB: mlst prints every locus as "-" here, which reads as "nothing found"
    # but is really mlst blanking a partial call below its score threshold
    # (visible only under --debug). Kleborate reports the same data honestly
    # as ST307-4LV.
    "ghana_GCA_025660315.1": "incomplete/likely mixed: 3/7 core MLST loci absent, ~1 Mbp of chromosome at 50-52% GC",
    # Severely fragmented: N50 3,641 bp across 2594 contigs, 6.08 Mbp against
    # an expected ~5.3 Mbp. Kleborate raises its own N50 QC warning here.
    # NOTE: an earlier revision of this file called this genome "contaminated"
    # on the basis of a smeared GC histogram. That was wrong: only 2.18 Mbp of
    # the 6.08 Mbp sits on contigs >=5 kb, so the histogram covered a third of
    # the assembly and short-contig GC is noisy. Fragmentation is established;
    # contamination is not.
    "ghana_GCA_022968745.1": "severely fragmented: N50 3641 bp, 2594 contigs, 6.08 Mbp (oversized)",
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


def load_kleborate():
    """Kleborate calls, keyed by strain. This is the benchmark's own tool, so
    its Bla_ESBL_acquired / Bla_Carb_acquired / ST columns are what make the
    comparison like-for-like."""
    rows = list(csv.DictReader(open(KLEBORATE), delimiter="\t"))
    def col(row, suffix):
        return row[[c for c in row if c.endswith(suffix)][0]]
    return {
        r["strain"]: {
            "st": col(r, "mlst__ST"),
            "esbl": col(r, "Bla_ESBL_acquired"),
            "carb": col(r, "Bla_Carb_acquired"),
            "species": col(r, "species__species"),
            "qc": col(r, "QC_warnings"),
        }
        for r in rows
    }


def positive(value):
    """Kleborate and Pathogenwatch both write '-' for absent, not ''."""
    return value.strip() not in ("", "-")


def main():
    amr = {
        os.path.basename(p).replace("_amrfinder.tsv", ""): beta_lactam_status(p)
        for p in glob.glob(f"{RESULTS}/*_amrfinder.tsv")
    }
    kb = load_kleborate()
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
        w.writerow(["sample_id", "st", "kleborate_esbl", "kleborate_carbapenemase",
                    "esbl_positive", "carbapenemase_positive",
                    "amrfinder_esbl_genes", "amrfinder_carbapenemase_genes",
                    "quality_flag"])
        for s_ in distinct:
            k = kb.get(s_, {})
            a_esbl, a_carb = amr[s_]
            w.writerow([s_, k.get("st", ""), k.get("esbl", ""), k.get("carb", ""),
                        positive(k.get("esbl", "-")), positive(k.get("carb", "-")),
                        ";".join(a_esbl), ";".join(a_carb), SUSPECT.get(s_, "")])

    bench = list(csv.DictReader(open(BENCHMARK)))
    bn = len(bench)
    b_esbl = sum(1 for r in bench if positive(r["bla_esbl_acquired"]))
    b_carb = sum(1 for r in bench if positive(r["bla_carb_acquired"]))

    print(f"Benchmark (Pathogenwatch clinical, Kleborate calls): n={bn}  "
          f"ESBL+ {b_esbl/bn:.1%}  carbapenemase+ {b_carb/bn:.1%}\n")

    for label, gs in [("This study (all distinct)", distinct),
                      ("This study (quality-filtered)", clean)]:
        n = len(gs)
        e = sum(1 for x in gs if positive(kb.get(x, {}).get("esbl", "-")))
        c = sum(1 for x in gs if positive(kb.get(x, {}).get("carb", "-")))
        print(f"{label:<32} n={n:<4} ESBL+ {e/n:6.1%}   carbapenemase+ {c/n:.1%}   [Kleborate]")

    # Cross-check: AMRFinder's CEPHALOSPORIN subclass is broader than
    # Kleborate's ESBL definition (it also picks up AmpC such as DHA-1 and
    # narrow-spectrum OXA-1), so the two disagree on a small number of calls.
    disagree = [x for x in distinct
                if positive(kb.get(x, {}).get("esbl", "-")) != bool(amr[x][0])]
    print(f"\nAMRFinder vs Kleborate ESBL agreement: "
          f"{len(distinct) - len(disagree)}/{len(distinct)}"
          + (f"  (differs: {', '.join(disagree)})" if disagree else ""))

    print("\nST distribution (this study, all distinct):")
    for s_, cnt in Counter(kb.get(x, {}).get("st", "?") for x in distinct).most_common():
        print(f"  {s_:<12} {cnt}")

    def bare(x):
        return x[2:] if x.startswith("ST") else x

    print("\nBenchmark top-8 clones vs this study:")
    b_st = Counter(bare(r["st"]) for r in bench)
    mine = Counter(bare(kb.get(x, {}).get("st", "")) for x in distinct)
    for s_, cnt in b_st.most_common(8):
        print(f"  ST{s_:<6} benchmark n={cnt:<4} this study n={mine.get(s_, 0)}")

    print("\nQuality-suspect genomes excluded from the filtered figures:")
    for s, why in SUSPECT.items():
        print(f"  {s}: {why}")
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
