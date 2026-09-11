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

The only collection here with usable metadata. Nigeria 333, Senegal 121,
Ghana 96, Benin 1.

| field | filled |
|---|---|
| Country | 557/557 |
| Isolation source | 409/557 |
| Collection date | 222/557 (2009–2022) |
| Host status | 87/557 |
| Host sex | 0/557 |
| Age | not a field |

Isolation sources: Blood 162, Faecal 101, Urine 46, **sputum 16**, Carriage 18,
ENRH 17, swabs 11, blank 148.

**Only ~16 sputum isolates exist in the entire West African set.** A
CAP-specific analysis is not achievable at meaningful sample size from this
data.

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
