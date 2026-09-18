#!/usr/bin/env bash
set -euo pipefail

# Run the two tools objective 2 needs over the assembled mobilization subset:
# AMRFinderPlus for gene-to-contig, and mob_recon to say whether that contig is
# chromosome or plasmid and which cluster it belongs to.
#
# RGI and PlasmidFinder are deliberately not run here. run_resistome_pipeline.sh
# runs all four for the study genomes, but objective 2 only needs the gene's
# coordinates and the contig's identity, and the other two roughly double the
# runtime without contributing to that join.
#
# Safe to re-run: a genome whose mob_recon contig report already exists is
# skipped, and a failure on one genome does not stop the rest.
#
#   bash resistome/scripts/run_mobilization_tools.sh
#   THREADS=4 bash resistome/scripts/run_mobilization_tools.sh

cd "$(dirname "$0")/../.."

THREADS="${THREADS:-6}"
RESULTS=resistome/results

shopt -s nullglob
assemblies=(resistome/raw/pw_*.fna)
if [ ${#assemblies[@]} -eq 0 ]; then
  echo "!!! no resistome/raw/pw_*.fna assemblies found" >&2
  exit 1
fi

echo "=== ${#assemblies[@]} assemblies found (THREADS=$THREADS) ==="
failed=()
done_count=0

for fasta in "${assemblies[@]}"; do
  sample=$(basename "$fasta" .fna)

  if [ -s "${RESULTS}/${sample}_mob_recon/contig_report.txt" ] &&
     [ -s "${RESULTS}/${sample}_amrfinder.tsv" ]; then
    echo "=== $sample: already processed, skipping ==="
    done_count=$((done_count + 1))
    continue
  fi

  if [ ! -s "${RESULTS}/${sample}_amrfinder.tsv" ]; then
    echo "=== $sample: AMRFinderPlus ==="
    if ! amrfinder -n "$fasta" -O Klebsiella_pneumoniae --threads "$THREADS" \
                   -o "${RESULTS}/${sample}_amrfinder.tsv" \
                   > "${RESULTS}/${sample}_amrfinder.log" 2>&1; then
      echo "!!! $sample: AMRFinder failed, skipping" >&2
      failed+=("$sample:amrfinder")
      continue
    fi
  fi

  echo "=== $sample: mob_recon ==="
  rm -rf "${RESULTS}/${sample}_mob_recon"
  if ! mob_recon -i "$fasta" -o "${RESULTS}/${sample}_mob_recon" \
                 -s "$sample" -n "$THREADS" -f \
                 > "${RESULTS}/${sample}_mob_recon.log" 2>&1; then
    echo "!!! $sample: mob_recon failed, skipping" >&2
    failed+=("$sample:mob_recon")
    continue
  fi

  done_count=$((done_count + 1))
  echo "=== $sample: done ($done_count/${#assemblies[@]}) ==="
done

echo "=== finished: $done_count/${#assemblies[@]} processed ==="
if [ ${#failed[@]} -gt 0 ]; then
  echo "=== ${#failed[@]} failed: ${failed[*]} ==="
fi
