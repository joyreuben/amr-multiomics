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

### Resolve the carbapenemase to gene family, not a binary
99/557 genomes carry an acquired carbapenemase, in two families that barely
overlap: **NDM** (NDM-1 46, NDM-7 12, NDM-5 3) and **OXA-48-like** (OXA-181 34,
OXA-48 5). One ST152 genome carries both. No KPC, VIM, IMP or GES anywhere in
the collection.

The two families track different lineages. ST17 is 28/28 OXA-181 with zero NDM;
the next four largest positive lineages — ST464 (NDM-7), ST147, ST395, ST442
(NDM-1) — are pure NDM. Carriage is spread across 28 distinct STs.

So objectives 1–3 must stratify by family. These are two independent gene flows
in different clonal backgrounds, on vehicles expected to differ (OXA-181 on
ColKP3/IncX3, NDM on IncF/IncC); a `carbapenemase+` flag averages across exactly
the contrast those objectives exist to measure, and would make a shared-replicon
signal unreadable. The binary in `profile_west_africa_sampling.py` is correct
**only** for objective 4, where gene identity does not bear on the question.

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
Resistome track only. Objectives 1, 3 and 4 are done; objective 2 is prepared
and method-verified but not run.

**Objective 1 — done.** `lineage_structure_carbapenemase.py` →
`results/lineage_structure.tsv`. The 99 positives need a minimum of 37
independent acquisitions across 28 STs, counted as distinct (clone, family)
pairs at 10 cgMLST alleles. Robust to threshold: 40 acquisitions at 5 alleles,
31 at 50. ST17 is the one true clone (28 × OXA-181, median 1 allele apart, 117
from its own ST17 negatives); ST147/ST392/ST15 are repeated acquisition into
one lineage; ST395 positives sit 0 alleles from same-ST negatives, so the gene
varies within a clone. Caveat carried in the output: 5 single-clone STs, ST17
included, come from one BioProject each, so an outbreak one study sampled
cannot be told from a widespread clone.

**Objective 3 — done.** `resistome_architecture.py` →
`results/resistome_architecture.tsv`. NDM and OXA-48-like carry disjoint cargo,
and the aminoglycoside markers are perfectly exclusive: armA 25/60 NDM vs 0/38
OXA, rmtB 0/60 vs 28/38, aph(3')-VI 27/60 vs 0/38, tet(G) 0/60 vs 28/38. armA
and rmtB are interchangeable 16S methyltransferases and no genome has both —
two vehicles, argued without any assembly. Associations are reported crude and
CMH-stratified by study; 20 of 99 tests are significant collection-wide but not
within study, so always read the CMH column.

**Objective 4 — done.** See `DATA_PROVENANCE.md`.

**Objective 2 — blocked on compute, not on method.** The export has no contig
coordinates for AMR genes, so the assemblies are required.
`carbapenemase_cohort.tsv` pins the 99 genomes (all resolve to paired ENA
reads, 27.2 GB) and `assemble_carbapenemase_cohort.sh` pins the pipeline. All
99 are assembled from reads even though 40 have public NCBI assemblies, because
mixing assemblers would confound the plasmid-vs-chromosome call that objective
2 rests on. Verified end-to-end on ERR4783440 (ST392, NDM-1): 167 contigs,
N50 215 kb, and blaNDM-1 lands on a plasmid contig in mob_suite cluster AA405
alongside CTX-M-15, OXA-1 and TEM-1 — 21 AMR genes on plasmid contigs against 7
on the chromosome. Budget ~45 min/genome, so ~3 days for the cohort. Note
`rep_type(s)` came back `-` on those contigs, so replicon naming may have to
fall back on mob_suite cluster IDs; measure how often before relying on it.

Cohort assemblies (`raw/pw_*.fna`) are gitignored — 99 × ~5.8 MB against a
`.git` already at 357 MB, and they regenerate from the manifest plus the
script. Derived evidence is committed.

Superseded: `resistome/scripts/compare_to_west_africa_benchmark.py` compares
study genomes against a clinical-only benchmark across mismatched compartments.
Do not extend it; its framing is retired with the CAP objective.
