#!/usr/bin/env bash
set -euo pipefail

# Run the 4-tool resistome pipeline (AMRFinderPlus, RGI, PlasmidFinder,
# MOB-recon) on every assembled genome in resistome/raw/, then fold each
# result into resistome/results/resistome_summary.csv via
# summarize_resistome.py.
#
# Skips samples already present in resistome_summary.csv (so it's safe to
# re-run after adding new genomes to raw/).
#
# Must be run from the project root (RGI's --local flag looks for ./localDB
# relative to cwd).

cd "$(dirname "$0")/../.."

PLASMIDFINDER_DB=/home/joy/miniconda3/envs/amr-multiomics/share/plasmidfinder-2.1.6/database
SUMMARY=resistome/results/resistome_summary.csv
THREADS=8

already_done() {
  [ -f "$SUMMARY" ] && cut -d, -f1 "$SUMMARY" | tail -n +2 | grep -qx "$1"
}

for fasta in resistome/raw/ghana_*.fna resistome/raw/nigeria_*.fna; do
  sample=$(basename "$fasta" .fna)

  if already_done "$sample"; then
    echo "=== $sample: already in $SUMMARY, skipping ==="
    continue
  fi

  echo "=== $sample: AMRFinderPlus ==="
  amrfinder -n "$fasta" -O Klebsiella_pneumoniae --threads "$THREADS" \
    -o "resistome/results/${sample}_amrfinder.tsv"

  echo "=== $sample: RGI ==="
  rgi main -i "$fasta" -o "resistome/results/${sample}_rgi" \
    -t contig -n "$THREADS" --local --clean

  echo "=== $sample: PlasmidFinder ==="
  mkdir -p "resistome/results/${sample}_plasmidfinder"
  plasmidfinder.py -i "$fasta" -o "resistome/results/${sample}_plasmidfinder" \
    -p "$PLASMIDFINDER_DB" -mp blastn -x

  echo "=== $sample: mob_recon ==="
  rm -rf "resistome/results/${sample}_mob_recon"
  mob_recon -i "$fasta" -o "resistome/results/${sample}_mob_recon" \
    -s "$sample" -n "$THREADS" -f

  echo "=== $sample: summarizing ==="
  python resistome/scripts/summarize_resistome.py --sample "$sample" \
    --amrfinder "resistome/results/${sample}_amrfinder.tsv" \
    --rgi "resistome/results/${sample}_rgi.txt" \
    --plasmidfinder "resistome/results/${sample}_plasmidfinder/results_tab.tsv" \
    --mob_contig_report "resistome/results/${sample}_mob_recon/contig_report.txt" \
    --out "$SUMMARY"

  echo "=== $sample: done ==="
done

echo "=== All samples processed. $SUMMARY now has: ==="
cut -d, -f1 "$SUMMARY"
