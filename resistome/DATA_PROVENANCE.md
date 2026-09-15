# Data provenance — resistome track

Verified 2026-09-11 against NCBI BioProject/BioSample and ENA records. Read this
before interpreting any comparison in `resistome/results/`.

## The short version

None of the genomes in this repository are documented respiratory or
pneumonia isolates. The Ghana genomes are hospital *environmental* and *rectal
carriage* isolates; the Nigeria genomes carry no isolation metadata at all.
Any comparison of these genomes against a *clinical* benchmark is comparing
different sample compartments and must say so explicitly.

## Ghana genomes (n=12 distinct, 23 files incl. GCA/GCF duplicates)

**BioProject PRJNA823741** — "Molecular characterization of multi-drug
resistant (MDR) Gram-negative pathogens from hospital [environment]".
Submitter: University of KwaZulu Natal. Location: Ghana, Kumasi.
Collected 2021. Released 2022-04-13.

Isolation sources, from BioSample attributes:

| source | n | compartment |
|---|---|---|
| rectal | 5 | gut carriage (not infection) |
| dripstand | 2 | hospital surface |
| bed | 1 | hospital surface |
| door | 1 | hospital surface |
| tap | 1 | hospital water |
| hand | 2 | hand swab (`host=Homo sapiens`) |

No blood, sputum, urine or other clinical infection isolates. No respiratory
specimens. This collection was designed to study the hospital environment as
an AMR reservoir, which is a different question from clinical infection.

Note: this data is likely already published — the Kumasi location matches
existing literature on MDR *K. pneumoniae* from that setting. Check for a
primary publication before framing any analysis as novel.

## Nigeria genomes (PRJEB29739, 746 K. pneumoniae read sets)

**NIHR Global Health Research Unit on Genomic Surveillance of AMR**
(Wellcome Sanger Institute). General AMR surveillance across 2,331 runs and
12+ species — including *Salmonella* Typhi and *E. faecalis*, which are not
respiratory pathogens. This is not a pneumonia study.

BioSample records (e.g. SAMEA9844719 for ERR7226174) contain organism name,
strain ID and deposition dates only. Empty across all 746:

- `isolation_source`
- `host`, `host_body_site`
- `collection_date`
- `country`

So genomes assembled from this list can support gene-presence analysis, but
not stratification by source, time or host. 139 of the 746 also appear in the
Pathogenwatch export below, where 112 of them *do* have isolation source —
those 139 are the analytically useful subset.

## Pathogenwatch West Africa export (n=557)

The only collection here with usable metadata. After normalising the free-text
country field (`NG`, `NG - Nigeria` → Nigeria; `SN` → Senegal; `GH` → Ghana):
Nigeria 335, Senegal 122, Ghana 97, Cameroon 2, Benin 1.

| field | filled |
|---|---|
| Country | 557/557 |
| Isolation source | 409/557 |
| Collection date | 222/557 |
| Collection date *or* `Date` | 510/557 (2009–2022) |
| Host status | 87/557 |
| Host sex | 0/557 |
| Age | not a field |

Note the two date fields. `Collection date` is filled for 222 genomes, but the
Pathogenwatch `Date` field is filled for 510 and agrees with `Collection date`
wherever both exist. Use the fallback — it roughly doubles the temporally
usable set. All 557 are *K. pneumoniae* by Kleborate v3.2.4.

Isolation sources: Blood 162, Faecal 101, Urine 46, **sputum 16**, Carriage 18,
ENRH 17, swabs 11, blank 148.

**Only ~16 sputum isolates exist in the entire West African set.** A
CAP-specific analysis is not achievable at meaningful sample size from this
data.

## The 557 are 16 studies stacked together, not a survey

Verified 2026-09-14. Reproduce with:

```
python resistome/scripts/profile_west_africa_sampling.py
```

Output: `resistome/results/west_africa_sampling_frame.tsv`.

This export is an aggregation of 16 BioProjects with different designs,
catchments and inclusion criteria. Acquired-carbapenemase rate per contributing
study (Kleborate v3.2.4 `Bla_Carb_acquired`, studies n≥15):

| study | n | carbapenemase | years | countries |
|---|---|---|---|---|
| PRJEB29739 | 139 | 12 (8.6%) | 2016–2018 | Nigeria 139 |
| PRJEB29143 | 101 | 0 (0.0%) | 2016 | Senegal 101 |
| PRJNA351846 | 79 | 1 (1.3%) | 2009–2016 | Nigeria 79 |
| PRJEB33565 | 56 | 24 (42.9%) | 2016 | Nigeria 56 |
| PRJEB56918 | 48 | 5 (10.4%) | 2020–2021 | Nigeria 30; Ghana 17; Benin 1 |
| PRJEB37523 | 39 | 29 (74.4%) | 2017–2019 | Ghana 39 |
| PRJNA473419 | 22 | 0 (0.0%) | 2015–2016 | Ghana 22 |
| PRJEB78995 | 21 | 5 (23.8%) | 2018–2021 | Senegal 21 |
| PRJEB27707 | 16 | 15 (93.8%) | 2016–2018 | Nigeria 16 |

**Between-study carbapenemase rate spans 0.0% to 93.8%.** That range is far
wider than any country or year contrast computed from the same genomes, so
aggregate stratifications are reporting study composition:

- *By country* — Ghana 35.1%, Nigeria 17.9%, Senegal 4.1%. But Ghana's rate is
  largely PRJEB37523 (74.4%), and Senegal's is largely PRJEB29143 (0.0%).
  Country is nearly nested within study.
- *By year* — 0% before 2014 rising to 36.0% in 2017. But every genome
  collected 2009–2014 comes from a single BioProject (PRJNA351846, itself
  1.3% carbapenemase), and 2016 is 200 genomes dominated by PRJEB29143 and
  PRJEB33565. Year is nearly nested within study too.

Only PRJNA351846 (n=79, 2009–2016) spans enough years for a within-study
temporal series, and it contains one carbapenemase-positive genome. There is no
usable within-study time trend here.

**Consequence.** Do not present carbapenemase prevalence by country or by year
from this collection as an epidemiological finding. Either condition on study —
e.g. `carbapenemase ~ country + year + (1|study)` — and report how much variance
study absorbs, or restrict claims to genome-intrinsic questions (which lineage
carries which gene, whether it sits on a plasmid, which genes co-occur) that do
not depend on the sampling frame being representative.

This is the same error as the compartment mismatch below, one level up: the
earlier problem was comparing across sample types, this one is comparing across
study designs.

## Consequence for `compare_to_west_africa_benchmark.py`

That script compares this study's genomes against
`pathogenwatch_west_africa_clinical_only.csv`, a benchmark deliberately
filtered to clinical isolates (Senegal carriage samples were excluded to make
it clinical). The study-side genomes are environmental/carriage (Ghana) or of
unknown source (Nigeria).

Environmental and carriage isolates are expected to carry fewer acquired
carbapenemases than invasive clinical isolates, so the reported gap
(0% vs 21.4%) is confounded by sample compartment and must not be presented as
a regional or temporal finding.

The gene-detection result itself is sound and independently verified three ways
(AMRFinderPlus, Kleborate, and read-level mapping with a blaCTX-M-15 positive
control). It is the choice of reference population that requires requalifying,
not the measurement.
