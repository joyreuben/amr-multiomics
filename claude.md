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
Resistome track only. All four objectives are done.

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

**Objective 2 — done, on 19 genomes not 99.** The export has no contig
coordinates for AMR genes, so assemblies were required.
`select_mobilization_subset.py` picked 19 genomes at the points where the
argument turns (rationale column in `mobilization_subset.tsv`),
`assemble_carbapenemase_cohort.sh` built them uniformly from reads, and
`analyze_mobilization.py` → `results/mobilization_placement.tsv` did the join.

16/18 carbapenemase calls sit on a plasmid contig. **One vehicle does cross
lineages**: cluster AA038 carries OXA-181 in ST17, ST234 and ST340; AA405
carries NDM-1 in ST147, ST307, ST340 and ST392; AA002 carries OXA-48 in ST307
and ST392. ST340 settles the mechanism — its two genomes are 4 alleles apart
and carry NDM-1 on AA405 and OXA-181 on AA038, two different vehicles into one
host background.

The two families ride very different plasmids. NDM clusters carry 10–20
co-located AMR genes; the OXA-181 cluster carries 1–4 (usually just qnrS1) and
the OXA-48 cluster carries none. Objective 3's specific markers are confirmed
physically co-located: armA on the carbapenemase cluster in 4/4 genomes that
have it, aph(3')-VI in 5/5 — whereas blaCTX-M-15, common everywhere, is only
9/19. Genome-level co-occurrence and same-replicon carriage are not the same
claim, and only this analysis separates them.

**Same cluster confirmed as the same plasmid.** `compare_plasmids.py` →
`results/plasmid_comparison.tsv` aligns every pair of reconstructed plasmids
(mash for ANI, blastn for reciprocal coverage; "same plasmid" = ≥95% identity
and ≥80% reciprocal coverage — identity alone would call any two plasmids
sharing a transposon identical).

The OXA-181 plasmid is the same element in ST17, ST234 and ST340: 11
cross-lineage pairs at **100% identity** and 79–96% coverage, 10/11 meeting the
threshold, and **statistically indistinguishable from within-ST pairs** (median
coverage 89.3% cross vs 91.7% within). Clonal descent cannot explain that — it
is one plasmid moving between lineages. AA405 behaves the same for NDM-1 across
ST147/ST307/ST340/ST392 (5/6 pairs, 97–99% identity) and AA002 for OXA-48
across ST307/ST392.

The plasmid-family caution was right and now has a number. AA038's ST464
NDM-7 member is **0/6** against the OXA-181 members at 68–80% coverage, so the
cluster splits exactly along gene lines: the OXA-181 members are one plasmid,
the NDM-7 member is a relative. Also note AA019 within ST147 is only 1/3, so a
shared cluster inside a single ST is not automatically a shared plasmid either.

Remaining caution: `rep_type(s)` resolves for only 6/18 carbapenemase contigs
while mob_suite cluster resolves 18/18, so clusters are the unit throughout.
Cargo is aggregated per cluster, never per contig: these assemblies run 167–215
contigs and one plasmid is routinely split across several.

Cohort assemblies (`raw/pw_*.fna`) are gitignored — 99 × ~5.8 MB against a
`.git` already at 357 MB, and they regenerate from the manifest plus the
script. Derived evidence is committed.

Superseded: `resistome/scripts/compare_to_west_africa_benchmark.py` compares
study genomes against a clinical-only benchmark across mismatched compartments.
Do not extend it; its framing is retired with the CAP objective.
