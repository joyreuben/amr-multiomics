#!/usr/bin/env bash
set -euo pipefail

# Download, QC and assemble the carbapenemase-positive cohort listed in
# resistome/results/carbapenemase_cohort.tsv (built by
# build_carbapenemase_cohort.py), so that objective 2 can ask whether each
# carbapenemase sits on a plasmid contig.
#
# Every genome is assembled the same way on purpose. 40 of the 99 have a public
# assembly in NCBI, but mixing those in would mean comparing contig structure
# across assemblers, and plasmid-vs-chromosome assignment is sensitive to
# exactly that. Uniform shovill assembly from the same reads keeps the
# mobilization call comparable across the cohort.
#
# Safe to re-run: any run whose final assembly already exists is skipped, and
# untrimmed reads are deleted once the assembly they produced is in place.
#
#   N_GENOMES=5 bash resistome/scripts/assemble_carbapenemase_cohort.sh   # smoke test
#   THREADS=4 bash resistome/scripts/assemble_carbapenemase_cohort.sh     # full run

cd "$(dirname "$0")/../.."

MANIFEST="${MANIFEST:-resistome/results/carbapenemase_cohort.tsv}"
N_GENOMES="${N_GENOMES:-0}"          # 0 = all
THREADS="${THREADS:-4}"              # leave headroom; this box has 8 cores
RAM="${RAM:-8}"                      # shovill refuses anything below 8
KEEP_READS="${KEEP_READS:-0}"        # trimmed reads are large; drop by default

if [ ! -s "$MANIFEST" ]; then
  echo "!!! manifest not found: $MANIFEST" >&2
  echo "    build it with resistome/scripts/build_carbapenemase_cohort.py" >&2
  exit 1
fi

mkdir -p resistome/raw/pw_reads resistome/processed/pw_qc resistome/processed/pw_assembly

# Columns: genome_id run_accession sample_accession study ST carbapenemase
#          family country date fastq_1 fastq_2 ...
mapfile -t ROWS < <(awk -F'\t' 'NR>1 && $14=="True" {print $2"\t"$10"\t"$11}' "$MANIFEST")
if [ "$N_GENOMES" -gt 0 ]; then
  ROWS=("${ROWS[@]:0:$N_GENOMES}")
fi

echo "=== ${#ROWS[@]} genomes queued (THREADS=$THREADS RAM=${RAM}G) ==="
started=$(date +%s)
done_count=0
failed=()

for row in "${ROWS[@]}"; do
  IFS=$'\t' read -r acc url1 url2 <<< "$row"
  target="resistome/raw/pw_${acc}.fna"

  if [ -s "$target" ]; then
    echo "=== $acc: assembly already present, skipping ==="
    done_count=$((done_count + 1))
    continue
  fi

  # wget exits 0 on a truncated file, so the gzip stream is verified before use:
  # an unverified short read reached fastp as "invalid gzip header" and, being
  # unguarded, took the whole run down with it under set -e.
  echo "=== $acc: downloading reads ==="
  download_ok=0
  for attempt in 1 2; do
    if wget -q --tries=3 --timeout=120 "https://$url1" -O "resistome/raw/pw_reads/${acc}_1.fastq.gz" &&
       wget -q --tries=3 --timeout=120 "https://$url2" -O "resistome/raw/pw_reads/${acc}_2.fastq.gz" &&
       gzip -t "resistome/raw/pw_reads/${acc}_1.fastq.gz" 2>/dev/null &&
       gzip -t "resistome/raw/pw_reads/${acc}_2.fastq.gz" 2>/dev/null; then
      download_ok=1
      break
    fi
    echo "!!! $acc: download incomplete or corrupt (attempt $attempt)" >&2
    rm -f "resistome/raw/pw_reads/${acc}_1.fastq.gz" "resistome/raw/pw_reads/${acc}_2.fastq.gz"
  done
  if [ "$download_ok" -ne 1 ]; then
    echo "!!! $acc: download failed twice, skipping" >&2
    failed+=("$acc:download")
    continue
  fi

  echo "=== $acc: quality control ==="
  if ! fastp -i "resistome/raw/pw_reads/${acc}_1.fastq.gz" -I "resistome/raw/pw_reads/${acc}_2.fastq.gz" \
             -o "resistome/processed/pw_qc/${acc}_1.trim.fastq.gz" -O "resistome/processed/pw_qc/${acc}_2.trim.fastq.gz" \
             -j "resistome/processed/pw_qc/${acc}_fastp.json" -h "resistome/processed/pw_qc/${acc}_fastp.html" \
             --thread "$THREADS" > "resistome/processed/pw_qc/${acc}_fastp.log" 2>&1; then
    echo "!!! $acc: fastp failed, skipping (see ${acc}_fastp.log)" >&2
    failed+=("$acc:qc")
    rm -f "resistome/raw/pw_reads/${acc}"_*.fastq.gz
    continue
  fi
  tail -3 "resistome/processed/pw_qc/${acc}_fastp.log"

  echo "=== $acc: assembly (slow step) ==="
  if ! shovill --outdir "resistome/processed/pw_assembly/${acc}" \
               --R1 "resistome/processed/pw_qc/${acc}_1.trim.fastq.gz" \
               --R2 "resistome/processed/pw_qc/${acc}_2.trim.fastq.gz" \
               --gsize 5.5M --cpus "$THREADS" --ram "$RAM" --force; then
    echo "!!! $acc: shovill failed, skipping" >&2
    failed+=("$acc:assembly")
    rm -f "resistome/raw/pw_reads/${acc}"_*.fastq.gz
    continue
  fi

  cp "resistome/processed/pw_assembly/${acc}/contigs.fa" "$target"

  rm -f "resistome/raw/pw_reads/${acc}_1.fastq.gz" "resistome/raw/pw_reads/${acc}_2.fastq.gz"
  if [ "$KEEP_READS" != "1" ]; then
    rm -f "resistome/processed/pw_qc/${acc}_1.trim.fastq.gz" \
          "resistome/processed/pw_qc/${acc}_2.trim.fastq.gz"
  fi

  done_count=$((done_count + 1))
  elapsed=$(( ($(date +%s) - started) / 60 ))
  echo "=== $acc: done ($done_count/${#ROWS[@]}, ${elapsed}m elapsed) -> $target ==="
done

echo "=== finished: $done_count/${#ROWS[@]} assemblies present ==="
if [ ${#failed[@]} -gt 0 ]; then
  echo "=== ${#failed[@]} failed: ${failed[*]} ==="
fi
ls resistome/raw/pw_*.fna 2>/dev/null | wc -l
