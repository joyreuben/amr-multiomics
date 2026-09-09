#!/usr/bin/env bash
set -euo pipefail

# Download, QC, and assemble Nigeria K. pneumoniae runs from the filtered
# PRJEB29739 accession list.
#
# N_GENOMES controls how many accessions to process (default 30). With only 3
# isolates, a zero-carbapenemase observation carries no weight: if the true
# rate matched the West Africa clinical benchmark (21.4%), seeing zero in 3
# genomes would happen 48.6% of the time by chance. At n=30 that falls to 0.1%.
#
# Safe to re-run: any accession whose final assembly already exists is
# skipped, and intermediate reads are removed after a successful assembly to
# keep the raw read footprint bounded.
#
#   N_GENOMES=50 bash resistome/scripts/nigeria_qc_assemble.sh

cd "$(dirname "$0")/../.."

N_GENOMES="${N_GENOMES:-30}"
ACCESSION_LIST=resistome/raw/PRJEB29739_klebsiella_only.tsv
# Trimmed reads are kept by default: they are what makes read-level checks
# possible (e.g. mapping reads against carbapenemase references to prove an
# absent gene is really absent rather than lost during assembly). Set to 0 to
# discard them if disk ever gets tight.
KEEP_READS="${KEEP_READS:-1}"

mkdir -p resistome/raw/nigeria_reads resistome/processed/nigeria_qc \
         resistome/processed/nigeria_assembly

mapfile -t ACCESSIONS < <(cut -f1 "$ACCESSION_LIST" | head -"$N_GENOMES")
echo "=== ${#ACCESSIONS[@]} accessions selected (N_GENOMES=$N_GENOMES) ==="

done_count=0
for acc in "${ACCESSIONS[@]}"; do
  if [ -s "resistome/raw/nigeria_${acc}.fna" ]; then
    echo "=== $acc: assembly already present, skipping ==="
    done_count=$((done_count + 1))
    continue
  fi

  echo "=== $acc: downloading reads ==="
  urls=$(awk -F'\t' -v acc="$acc" '$1==acc {print $6}' "$ACCESSION_LIST")
  IFS=';' read -r url1 url2 <<< "$urls"
  if [ -z "${url1:-}" ] || [ -z "${url2:-}" ]; then
    echo "!!! $acc: no FASTQ URLs in $ACCESSION_LIST, skipping" >&2
    continue
  fi
  wget -q "https://$url1" -O "resistome/raw/nigeria_reads/${acc}_1.fastq.gz"
  wget -q "https://$url2" -O "resistome/raw/nigeria_reads/${acc}_2.fastq.gz"

  echo "=== $acc: quality control ==="
  fastp -i "resistome/raw/nigeria_reads/${acc}_1.fastq.gz" -I "resistome/raw/nigeria_reads/${acc}_2.fastq.gz" \
        -o "resistome/processed/nigeria_qc/${acc}_1.trim.fastq.gz" -O "resistome/processed/nigeria_qc/${acc}_2.trim.fastq.gz" \
        -j "resistome/processed/nigeria_qc/${acc}_fastp.json" -h "resistome/processed/nigeria_qc/${acc}_fastp.html"

  echo "=== $acc: assembly (this is the slow step, be patient) ==="
  # shovill rejects --ram below 8 regardless of how much the machine has.
  shovill --outdir "resistome/processed/nigeria_assembly/${acc}" \
          --R1 "resistome/processed/nigeria_qc/${acc}_1.trim.fastq.gz" \
          --R2 "resistome/processed/nigeria_qc/${acc}_2.trim.fastq.gz" \
          --gsize 5.5M --ram 8 --force

  cp "resistome/processed/nigeria_assembly/${acc}/contigs.fa" "resistome/raw/nigeria_${acc}.fna"

  # Untrimmed reads are re-downloadable and dominate disk use; drop them once
  # the assembly they produced is safely in place.
  rm -f "resistome/raw/nigeria_reads/${acc}_1.fastq.gz" \
        "resistome/raw/nigeria_reads/${acc}_2.fastq.gz"
  if [ "$KEEP_READS" != "1" ]; then
    rm -f "resistome/processed/nigeria_qc/${acc}_1.trim.fastq.gz" \
          "resistome/processed/nigeria_qc/${acc}_2.trim.fastq.gz"
  fi

  done_count=$((done_count + 1))
  echo "=== $acc: done ($done_count/${#ACCESSIONS[@]}), saved as resistome/raw/nigeria_${acc}.fna ==="
done

echo "=== finished: $done_count/${#ACCESSIONS[@]} assemblies present ==="
ls resistome/raw/nigeria_*.fna | wc -l
