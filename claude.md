# AMR multi-omics project

## What this project is
Genomic determinants and mobilization of carbapenem resistance in *Klebsiella
pneumoniae* in West Africa (Nigeria, Ghana, Senegal).

The question is **how** carbapenem resistance spreads — which lineages carry
which carbapenemases, whether those genes sit on mobile elements, and whether
their distribution reflects clonal expansion or repeated horizontal acquisition.

### Objectives
1. **Lineage structure of carbapenemase carriage.** Which STs carry acquired
   carbapenemases, and are positives within an ST a single clone or independent
   acquisitions? Uses Kleborate ST calls plus the cgMLST allele profiles in the
   Pathogenwatch export.
2. **Mobilization.** For each carbapenemase-positive genome, is the gene on a
   plasmid contig, and which replicon type? Uses MOB-suite and PlasmidFinder.
   Shared replicons across distant lineages are evidence for horizontal spread.
3. **Resistome architecture.** Which AMR genes co-occur, and on which replicon
   types. Genome-intrinsic, so robust to how the collection was assembled.
4. **Sampling-frame heterogeneity as a result.** Quantify how much of the
   apparent country and year variation in carbapenemase prevalence is explained
   by contributing study rather than epidemiology, and report it as a limit on
   what public genome collections can support for AMR surveillance in the region.

### Objectives that were dropped, and why
The original objective was multi-omics (host genomics + microbiome + resistome)
modeling of AMR emergence in **CAP-causing pathogens**. Changed 2026-09-14 with
approval. Both halves failed on data availability, not on interest:

- **CAP / respiratory.** No genome in this repository is a documented
  respiratory isolate, and the entire West African Pathogenwatch export contains
  ~16 sputum isolates, none carbapenemase-positive. See
  `resistome/DATA_PROVENANCE.md`.
- **Host genomics + microbiome.** No data was ever obtained for either track.
  `host_genomics/` and `microbiome/` hold only placeholder files. "Multi-omics"
  now means layers of pathogen genomic evidence — resistome, plasmidome,
  virulome, core genome — not multiple organisms. The folders are kept in case
  those tracks are revived; nothing currently depends on them.
- **Emergence over time.** Ruled out as a headline claim: year is nearly nested
  within contributing study. Retained only as objective 4, where the confounding
  is the subject rather than a flaw.

Also ruled out: hypervirulence–MDR convergence. Only 2 of 557 genomes have
Kleborate virulence score ≥3 together with resistance score ≥1.

## Environment
- conda env: amr-multiomics (python 3.11)
- Activate with: conda activate amr-multiomics

## Rules
- Never fabricate or simulate data to fill gaps. If real data isn't available yet,
  say so explicitly.
- Cite the source/tool version for any bioinformatics result.
- Keep each omics track's code in its own folder: host_genomics/, microbiome/, resistome/
- Shared feature tables go in integration/
- Before stratifying the Pathogenwatch collection by country, year or isolation
  source, check `resistome/DATA_PROVENANCE.md`. Both the sample compartment and
  the contributing study confound those comparisons.

## Data
- **Pathogenwatch West Africa export, n=557** — the primary dataset. All
  *K. pneumoniae*. Kleborate v3.2.4 typing, cgMLST (~627 loci), PlasmidFinder
  replicons, MLST. Country 557/557, isolation source 409/557, year 510/557
  (2009–2022). In `resistome/results/pathogenwatch_west_africa/`.
- **27 study genomes** (12 distinct Ghana + Nigeria assemblies) — processed
  through the 4-tool pipeline. Ghana isolates are hospital-environmental and
  rectal carriage; Nigeria isolates have no isolation metadata. Useful for
  gene-presence work, not for stratified comparison.

## Current status
Resistome track only. 4-tool pipeline (AMRFinderPlus, RGI, PlasmidFinder,
MOB-suite) run end-to-end on all 27 genomes; MLST and Kleborate typing done;
sampling-frame profiling done. Gene detection is independently verified three
ways (AMRFinderPlus, Kleborate, read-level mapping with a blaCTX-M-15 positive
control).

Next: objectives 1 and 2 — join carbapenemase calls to cgMLST distances and to
plasmid replicon assignments across the 557.

Superseded: `resistome/scripts/compare_to_west_africa_benchmark.py` compares
study genomes against a clinical-only benchmark across mismatched compartments.
Do not extend it; its framing is retired with the CAP objective.
