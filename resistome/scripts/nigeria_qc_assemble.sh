#!/usr/bin/env bash
set -euo pipefail

# Take the first 3 real Klebsiella run accessions from the filtered Nigeria list
mapfile -t ACCESSIONS < <(cut -f1 resistome/raw/PRJEB29739_klebsiella_only.tsv | head -3)

for acc in "${ACCESSIONS[@]}"; do
  echo "=== $acc: downloading reads ==="
  urls=$(awk -F'\t' -v acc="$acc" '$1==acc {print $6}' resistome/raw/PRJEB29739_klebsiella_only.tsv)
  IFS=';' read -r url1 url2 <<< "$urls"
  wget -q "https://$url1" -O "resistome/raw/nigeria_reads/${acc}_1.fastq.gz"
  wget -q "https://$url2" -O "resistome/raw/nigeria_reads/${acc}_2.fastq.gz"

  echo "=== $acc: quality control ==="
  fastp -i "resistome/raw/nigeria_reads/${acc}_1.fastq.gz" -I "resistome/raw/nigeria_reads/${acc}_2.fastq.gz" \
        -o "resistome/processed/nigeria_qc/${acc}_1.trim.fastq.gz" -O "resistome/processed/nigeria_qc/${acc}_2.trim.fastq.gz" \
        -j "resistome/processed/nigeria_qc/${acc}_fastp.json" -h "resistome/processed/nigeria_qc/${acc}_fastp.html"

  echo "=== $acc: assembly (this is the slow step, be patient) ==="
  shovill --outdir "resistome/processed/nigeria_assembly/${acc}" \
          --R1 "resistome/processed/nigeria_qc/${acc}_1.trim.fastq.gz" \
          --R2 "resistome/processed/nigeria_qc/${acc}_2.trim.fastq.gz" \
          --gsize 5.5M --ram 8 --force

  cp "resistome/processed/nigeria_assembly/${acc}/contigs.fa" "resistome/raw/nigeria_${acc}.fna"
  echo "=== $acc: done, saved as resistome/raw/nigeria_${acc}.fna ==="
done
