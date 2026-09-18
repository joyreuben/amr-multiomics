# Genomic determinants and mobilization of carbapenem resistance in *Klebsiella pneumoniae* in West Africa

---

## Methods

### Genome collection and metadata

The primary dataset was a Pathogenwatch export of 557 *Klebsiella pneumoniae*
genomes from West Africa, comprising isolates from Nigeria, Senegal, Ghana,
Cameroon and Benin. Species identity was confirmed for all 557 by Pathogenwatch
Speciator. The export supplied multi-locus sequence types (MLST), core-genome
MLST (cgMLST) allele profiles over a 629-locus scheme, PlasmidFinder replicon
calls, and antimicrobial resistance determinants called by Kleborate v3.2.4.

Country was recorded for 557/557 isolates, isolation source for 409/557, and
collection year for 510/557 (2009–2022). Two date fields were present in the
export: `Collection date`, populated for 222 isolates, and `Date`, populated for
510. The two agreed wherever both were present, and the latter was used with the
former as a fallback, approximately doubling the temporally usable set. The
free-text country field required normalisation (`NG`, `NG - Nigeria` → Nigeria;
`SN` → Senegal; `GH` → Ghana), after which the composition was Nigeria 335,
Senegal 122, Ghana 97, Cameroon 2 and Benin 1.

Acquired carbapenemases were taken from the Kleborate `Bla_Carb_acquired` field.
Throughout, carbapenemases were resolved to gene family (NDM or OXA-48-like)
rather than collapsed to a presence/absence indicator, except where the analysis
question was explicitly independent of gene identity.

### Assessment of the sampling frame

Because the collection aggregates isolates contributed by multiple independent
studies, the extent to which apparent geographic and temporal variation reflects
study composition rather than epidemiology was quantified before any stratified
comparison was made. Isolates were grouped by INSDC study accession, and
acquired-carbapenemase prevalence, year span and country composition were
tabulated per contributing BioProject. The between-study range in prevalence was
then compared against the country- and year-stratified estimates computed from
the same genomes.

### Lineage structure and independent acquisition counting

cgMLST allele profiles were reshaped from the long-format export into a
genome-by-locus matrix. Loci called per genome ranged from 590 to 628 of 629
(median 628). Pairwise distance between two genomes was defined as the number of
loci at which they carried different alleles, counted only over loci called in
both, and rescaled to the full 629-locus scheme so that genomes differing in
missing-data burden remained comparable.

Within each sequence type carrying at least one carbapenemase, positive genomes
were clustered by single linkage at a threshold of 10 allele differences, a
threshold in common use for *K. pneumoniae* cgMLST. The minimum number of
independent acquisition events was then counted as the number of distinct
(cluster, gene family) pairs rather than the number of clusters. This distinction
matters: two genomes within a single cgMLST cluster that carry different
carbapenemase families must have acquired the phenotype on two separate
occasions, and a cluster count alone would score such a case as one event.
Sensitivity to the clustering threshold was assessed at 5, 10, 20 and 50 allele
differences.

Distance from carbapenemase-positive to carbapenemase-negative genomes of the
same sequence type was computed to distinguish positives forming a distinct
subclade from positives interleaved with negatives. The contributing study of
each positive genome was carried through the analysis so that clusters confined
to a single BioProject could be identified.

### Phylogenetic reconstruction

A neighbour-joining tree was built from the same rescaled cgMLST distance matrix
over all 557 genomes. Neighbour joining was implemented directly with vectorised
per-iteration operations, as the pure-Python implementation in Bio.Phylo is
prohibitively slow at O(n³) for 557 taxa. Tree topology was validated against the
independently derived single-linkage clustering by testing whether defined groups
of interest were recovered as monophyletic clades. Tip annotations comprised
sequence type, carbapenemase family and allele, country and contributing study.

### Resistome architecture and co-occurrence testing

Acquired resistance determinants were extracted from the Kleborate acquired-gene
columns spanning aminoglycoside, colistin, fosfomycin, fluoroquinolone,
glycopeptide, macrolide, phenicol, rifampicin, sulphonamide, tetracycline,
tigecycline, trimethoprim and beta-lactam classes. Chromosomal determinants
(`Bla_chr`) and point-mutation columns were excluded, as intrinsic genes and
resistance mutations cannot be carried on a plasmid and therefore cannot
co-transfer with a carbapenemase. Allele version and confidence marks were
stripped so that alleles of the same gene were pooled.

Every gene present in at least 10 genomes was tested for association with
carbapenemase carriage, separately for each family, using Fisher's exact test.
Because the collection aggregates studies with markedly different carbapenemase
prevalence, each association was additionally estimated as a
Cochran–Mantel–Haenszel (CMH) common odds ratio stratified by contributing
BioProject, restricted to strata of at least 10 genomes. Where crude and
stratified estimates disagreed, the stratified estimate was taken as
authoritative. p-values were corrected for multiple testing by the
Benjamini–Hochberg procedure, and significance was assessed on the
study-stratified q-value.

Genome-level association between carbapenemase family and PlasmidFinder replicon
type was tabulated from the export as a preliminary, assembly-free assessment of
plasmid vehicle.

### Selection of a subset for mobilization analysis

The Pathogenwatch export reports replicon calls per contig but reports no contig
coordinates for resistance genes, so the export alone cannot establish whether a
carbapenemase is physically located on a plasmid. Resolving this required
assemblies.

Rather than assembling all 99 carbapenemase-positive genomes, 19 were selected at
the points where the inferential argument turns, each with a documented
rationale: three from the largest clonal carbapenemase-positive lineage; two
carrying the same allele in an unrelated sequence type; four spanning the
separate cgMLST clusters of a sequence type showing repeated acquisition; two
carrying a different allele of the same family; one of each family from each of
three sequence types carrying both; and a positive–negative pair separated by
zero cgMLST alleles. Where a sequence type contributed several genomes, these
were drawn from different cgMLST clusters so as to span the diversity of that
lineage rather than resample a single clone.

### Sequencing data retrieval and genome assembly

Paired-end short reads were retrieved from the European Nucleotide Archive for
all selected genomes. Downloaded files were verified for gzip stream integrity
before use, as the archive's generated FASTQ path was found to be broken for one
run, returning an HTML directory listing under a successful HTTP status; for that
run the originally submitted read files were retrieved instead.

Reads were quality- and adapter-trimmed with fastp v1.3.6 and assembled with
Shovill v1.1.0 (SPAdes v3.15.5) using an expected genome size of 5.5 Mb. All
genomes were assembled from reads under identical parameters, including those for
which public assemblies already existed, because plasmid-versus-chromosome
assignment is sensitive to assembler and to contig fragmentation, and mixing
assembly provenance would confound precisely the measurement of interest. The
resulting assemblies comprised 126–329 contigs (median 169), totalling
5.35–5.97 Mb (median 5.73 Mb), with N50 78.6–374.0 kb (median 214.7 kb).

### Resistance gene localisation and plasmid reconstruction

Resistance genes were called with AMRFinderPlus v4.2.7 (database 2026-05-15.1)
under the *Klebsiella pneumoniae* organism model, which reports the contig on
which each gene lies. Contigs were assigned to chromosome or plasmid, grouped
into putative replicons, and typed with MOB-suite `mob_recon` v3.1.8. Gene calls
were joined to contig assignments on contig identifier.

Analyses were conducted at the level of the MOB-suite primary cluster rather than
the individual contig or the replicon type, for two reasons. First, replicon type
resolved for only 6 of 18 carbapenemase-bearing contigs, whereas a primary
cluster was assigned in 18 of 18. Second, short-read assemblies of this
fragmentation routinely distribute a single plasmid across several contigs, so
co-location assessed per contig systematically understates plasmid cargo.

### Plasmid sequence comparison

Assignment to a shared MOB-suite cluster indicates similarity to a common
reference group and does not by itself demonstrate that two isolates carry the
same element. Reconstructed plasmid sequences were therefore compared directly,
pairwise, for every pair of genomes carrying a carbapenemase on the same cluster.

Two complementary measures were computed: Mash v2.3 distance, giving a
whole-sequence average nucleotide identity estimate, and BLASTN v2.16.0+ derived
reciprocal alignment coverage with length-weighted identity. Both are necessary:
two plasmids sharing only a common transposon achieve near-perfect identity over
the aligned region while sharing little of their total length, so an
identity-based criterion alone would classify any two plasmids carrying a shared
mobile element as identical. Aligned intervals were merged on the query before
coverage was summed, to avoid double-counting overlapping high-scoring pairs.

A pair was classified as the same plasmid at ≥95% length-weighted identity and
≥80% reciprocal coverage. The decisive comparison was between cross-lineage and
within-sequence-type pairs: plasmids from genomes of a single clonal lineage are
expected to be near-identical through vertical inheritance alone, so the question
is whether cross-lineage pairs attain comparable similarity.

### Software and reproducibility

Analyses were performed in Python 3.11.15 with pandas, NumPy and SciPy. All
analysis scripts, intermediate tables and derived outputs are version-controlled.
Assemblies and the cgMLST distance cache are excluded from version control but
regenerate deterministically from the committed accession manifests and pipeline
scripts.

---

## Results

### Composition of the collection and its sampling frame

The 557 genomes derive from 16 contributing BioProjects, nine of which contributed
15 or more isolates. Acquired-carbapenemase prevalence varied from 0.0% to 93.8%
between these nine studies (Table 1). This range substantially exceeds any
contrast obtainable by stratifying the same genomes geographically or temporally.

Stratified by country, carbapenemase prevalence was 35.1% in Ghana, 17.9% in
Nigeria and 4.1% in Senegal. However, the Ghanaian estimate is dominated by a
single study with 74.4% prevalence, and the Senegalese estimate by a single study
with 0.0%; country is therefore close to nested within study. The temporal
pattern behaves similarly: prevalence rose from 0% before 2014 to 36.0% in 2017,
but every genome collected between 2009 and 2014 originates from one BioProject,
itself only 1.3% carbapenemase-positive, and the 200 genomes from 2016 are
dominated by two studies at opposite extremes of prevalence. Only one study
(n = 79) spans sufficient years to support a within-study temporal series, and it
contains a single carbapenemase-positive genome.

Accordingly, carbapenemase prevalence by country or by year is not interpretable
as an epidemiological quantity in this collection, and subsequent analyses were
restricted to genome-intrinsic questions or were conditioned on contributing
study.

### Carbapenemase diversity

Ninety-nine of 557 genomes (17.8%) carried an acquired carbapenemase, distributed
across two families. NDM was present in 60 genomes (NDM-1 n = 46, NDM-7 n = 12,
NDM-5 n = 3) and OXA-48-like in 38 (OXA-181 n = 34, OXA-48 n = 5). A single ST152
genome carried both NDM-1 and OXA-48. No KPC, VIM, IMP or GES determinant was
detected in any genome.

Carbapenemase carriage was distributed across 28 distinct sequence types, of
which 13 contained two or more positive genomes and 15 contained a single
positive.

### Carbapenem resistance arises predominantly by repeated acquisition

The 99 carbapenemase-positive genomes require a minimum of 37 independent
acquisition events. This estimate is robust to the clustering threshold: 40
events at 5 allele differences, 37 at 10, 32 at 20 and 31 at 50. A tenfold
relaxation in the definition of a clone therefore alters the total by less than a
quarter and at no point approaches a single origin.

Three distinct patterns were observed among lineages with multiple positives.

**Clonal expansion.** ST17 contained 28 carbapenemase-positive genomes, all
carrying OXA-181, separated by a median of 1 allele (range 0–7) and lying a
median of 116.6 alleles from the 20 carbapenemase-negative ST17 genomes in the
collection. The neighbour-joining tree recovered these 28 genomes as a perfectly
monophyletic clade (100% purity) nested within the broader ST17 background,
independently confirming the distance-based result.

**Repeated acquisition within a lineage.** ST147 contained nine NDM-1-positive
genomes resolving into three separate cgMLST clusters, drawn from three
contributing studies and two countries. The smallest tree clade containing all
nine also contained three carbapenemase-negative ST147 genomes, so positives are
interleaved with negatives rather than forming a distinct subclade. ST392 and
ST15 showed the same pattern.

**Gene variation within a clone.** ST395 positives lay a median of 0 alleles from
same-sequence-type negatives, indicating that carbapenemase presence varies
between genomes that are otherwise indistinguishable at the core genome.

The clearest single instance of independent acquisition is ST340, in which two
genomes separated by 4 cgMLST alleles carry different carbapenemase families. On
host genotype alone these would be scored as one clone; the gene calls establish
two separate events.

Of the lineages consistent with a single clone, five — including ST17 — derive
entirely from one BioProject each. For these, an outbreak captured by a single
study cannot be distinguished from a genuinely widespread clone.

### The two carbapenemase families carry disjoint accessory resistomes

Testing every acquired resistance gene present in at least 10 genomes revealed
that the two families travel with almost non-overlapping accessory gene content.
The aminoglycoside determinants are mutually exclusive (Table 2): *armA* was
present in 25 of 60 NDM genomes and 0 of 38 OXA-48-like genomes, while *rmtB* was
present in 0 of 60 and 28 of 38 respectively. No genome in the collection carried
both. Both encode 16S rRNA methyltransferases conferring essentially
pan-aminoglycoside resistance and are therefore functionally interchangeable;
their strict partition between families is not explicable by selection for
phenotype.

The same partition held for *aph(3′)-VI* and *mphE*/*msrE* (NDM-associated) and
*tet(G)* and *aac(3)-IId* (OXA-48-like-associated).

Study stratification materially changed the conclusions. Of 99 gene-level tests,
20 were significant across the collection but not after conditioning on
contributing BioProject, including *sul2*, *strA*, *strB*, *aac(6′)-Ib-cr*,
*qnrB2* and *arr-3*. These reflect study composition rather than co-transfer and
are not reported as associations.

### Genome-level replicon association

Genomes carried a median of 4 PlasmidFinder replicon types (range 0–10).
Replicon distribution differed sharply by carbapenemase family: ColKP3 was
present in 82% of OXA-48-like genomes, 0% of NDM genomes and 0% of the 458
carbapenemase-negative genomes; IncX3 in 89%, 25% and 0.2%; and IncQ1 in 74%, 0%
and 1%. Conversely IncFIB(pQil) and ColpVC were confined to NDM genomes (43% and
30%) and absent from OXA-48-like genomes. Because genomes carry several plasmids,
these associations are suggestive of vehicle but cannot establish physical
linkage.

### Carbapenemases are predominantly plasmid-borne

Across the 19 assembled genomes, 16 of 18 carbapenemase calls were located on
plasmid contigs; two (an NDM-1 in ST395 and an NDM-7 in ST464) were assigned to
the chromosome.

Four MOB-suite clusters carried carbapenemases, three of which spanned multiple
sequence types (Table 3): AA038 carried OXA-181 in ST17, ST234 and ST340; AA405
carried NDM-1 in ST147, ST307, ST340 and ST392; and AA002 carried OXA-48 in ST307
and ST392. AA019 carried NDM-1 within ST147 only.

ST340, in which two near-identical genomes carry different families, carried
NDM-1 on AA405 and OXA-181 on AA038 — two distinct vehicles entering effectively
the same host background.

Plasmid cargo differed markedly between families. NDM-carrying clusters carried
10–20 additional resistance genes, including *armA*, *aph(3′)-VI*, *blaCTX-M-15*,
*blaOXA-1*, *blaTEM-1*, *aadA2*, *dfrA12* and *mph(E)*. The OXA-181 cluster
carried 1–4 additional genes, most often only *qnrS1*, and the OXA-48 cluster
carried none.

Contig-level localisation refined the co-occurrence results. Genes showing a
family-specific genome-level association were confirmed to be physically
co-resident on the carbapenemase-bearing cluster: *armA* in 4 of 4 genomes
carrying it, and *aph(3′)-VI* in 5 of 5. By contrast *blaCTX-M-15*, which is
common throughout the collection, was co-resident in only 9 of 19. Genome-level
co-occurrence and same-replicon carriage are therefore distinguishable, and only
the latter supports an inference of co-transfer.

### A single plasmid carries OXA-181 across unrelated lineages

Direct pairwise alignment established that shared cluster assignment reflects a
genuinely shared element for the OXA-181 plasmid. Eleven cross-lineage pairs
spanning ST17, ST234 and ST340 showed 100% length-weighted identity and 79.2–95.8%
reciprocal coverage, with 10 of 11 meeting the same-plasmid criterion.

Critically, these cross-lineage pairs were not distinguishable from
within-sequence-type pairs, which are near-identical by vertical inheritance:
median reciprocal coverage was 89.3% across lineages versus 91.7% within (both at
100% median identity). Since clonal descent cannot account for the cross-lineage
similarity, the parsimonious interpretation is that a single plasmid has moved
between these lineages.

AA405 behaved comparably for NDM-1, with 5 of 6 cross-lineage pairs meeting the
criterion at 97.1–99.2% identity, as did AA002 for OXA-48 across ST307 and ST392.

The comparison also delimited the shared-vehicle claim. Within AA038, the ST464
genome carrying NDM-7 scored 0 of 6 against the OXA-181-carrying members, at
67.9–79.6% reciprocal coverage. The cluster therefore partitions exactly along
gene lines: its OXA-181 members constitute one plasmid, while the NDM-7 member is
a related element sharing a backbone. Similarly, AA019 within ST147 yielded only
1 of 3 same-plasmid pairs, demonstrating that shared cluster assignment within a
single sequence type does not guarantee a shared plasmid either.

---

## Discussion

### Principal findings

This study addressed whether carbapenem resistance in West African
*K. pneumoniae* propagates predominantly through clonal expansion or through
horizontal transfer of resistance genes, and reaches three principal conclusions.

First, repeated acquisition dominates. The 99 carbapenemase-positive genomes
require a minimum of 37 independent acquisition events across 28 sequence types,
and this figure is stable across a tenfold variation in clustering stringency.
Carbapenem resistance in this collection is therefore not the dissemination of a
single successful clone.

Second, two independent gene flows are in circulation. NDM and OXA-48-like
determinants occupy different lineages, ride different plasmids and carry
disjoint accessory resistomes. The most striking evidence is the strict mutual
exclusivity of *armA* and *rmtB*. These genes are functionally equivalent, so
their partition cannot be explained by phenotypic selection and is most
parsimoniously attributed to physical linkage on two separately circulating
plasmid lineages — an inference reached from gene content alone and subsequently
confirmed by direct plasmid alignment.

Third, a single plasmid demonstrably crosses lineage boundaries. The OXA-181
plasmid is indistinguishable, at 100% identity and 79–96% reciprocal coverage,
between ST17, ST234 and ST340, and is no less similar across lineages than within
a single clonal lineage. As vertical inheritance cannot produce this pattern in
unrelated hosts, horizontal transfer of an intact element is the most economical
explanation.

These findings are mutually reinforcing and were obtained by methods that fail in
different ways — allele-distance clustering, tree reconstruction, statistical
gene association, and direct sequence alignment — which strengthens confidence in
the shared conclusion.

### Relationship to existing knowledge

The individual genetic associations reported here are largely consistent with the
international literature: OXA-181 is well documented on ColKP3/IncX3 backbones,
NDM determinants on IncF and IncHI1B replicons, and ST17, ST147 and ST307 are
recognised high-risk clones. The contribution of this study is therefore
principally regional characterisation — establishing which of these globally
described elements are circulating in Nigeria, Ghana and Senegal, in which
lineages, and with what accessory content — together with the methodological
finding described below.

### Implications for control

If carbapenem resistance were spreading principally as a clone, conventional
infection-control measures directed at that clone would be expected to contain
it. The predominance of independent acquisition argues otherwise: the resistance
genes are circulating within a mobile element pool accessible to many lineages,
and interventions aimed at individual strains will not interrupt that reservoir.

The asymmetry in plasmid cargo carries direct therapeutic significance. The
NDM-associated plasmids carry 10–20 additional resistance determinants, so a
single conjugative transfer can confer simultaneous resistance to carbapenems,
extended-spectrum cephalosporins, aminoglycosides (via *armA*) and macrolides,
leaving severely constrained treatment options. OXA-48-like plasmids, carrying
little beyond the carbapenemase, present a comparatively narrower resistance
burden. Surveillance and stewardship that treat "carbapenem resistance" as a
single entity will obscure this distinction, and the two families warrant
separate monitoring.

### Sampling-frame heterogeneity as a substantive finding

The collection analysed here is an aggregation of 16 independent studies rather
than a survey, with between-study carbapenemase prevalence spanning 0.0% to
93.8%. Both country and collection year are close to nested within contributing
study, so prevalence estimates stratified on either quantity principally reflect
study composition.

The practical consequence was quantified: 20 of 99 gene-level association tests
were significant across the pooled collection but not after conditioning on
study. An analysis that pooled these studies without stratification would have
reported those associations as biological findings.

This has implications beyond the present work. Public genome repositories are
increasingly used as a basis for regional AMR surveillance in settings where
systematic surveillance is unavailable, and such collections are typically
convenience aggregations with the structure described here. The results argue
that prevalence estimates derived from them should be regarded as uninterpretable
unless conditioned on contributing study, while genome-intrinsic questions —
which lineage carries which gene, on which element, with what co-resident content
— remain answerable because they do not require the sample to be representative.

### Limitations

Several limitations qualify these conclusions.

The analysis is entirely secondary, based on publicly deposited genomes and their
associated metadata; no isolates were sequenced for this study and no phenotypic
susceptibility testing was performed, so all resistance is inferred genotypically.

The sampling frame is not representative, as documented above. Most consequentially,
the largest single signal — the 28-genome ST17 OXA-181 clonal expansion — derives
entirely from one BioProject, so it cannot be determined whether this clone is
disseminated across the region or was simply well sampled in one setting.

Mobilization analysis rests on 19 genomes selected to span the key comparisons
rather than on the full set of 99 positives. This design supports conclusions
about the identity and lineage distribution of vehicles but does not support
prevalence estimates for particular plasmid types.

Short-read assembly cannot fully resolve plasmid structure. Assemblies comprised
126–329 contigs, so plasmids are fragmented and reconstructed by inference;
analyses were conducted at cluster level to mitigate this, but long-read
sequencing would be required to establish complete plasmid structures, to confirm
circularity, and to characterise the immediate genetic context of the
carbapenemase genes. The two chromosomal assignments should be regarded with
particular caution, as failure to recognise a plasmid contig would produce
exactly this observation.

MOB-suite clusters are similarity groupings against a reference database rather
than definitive replicon identities. This limitation was addressed by direct
pairwise alignment, which showed that cluster AA038 in fact comprises at least
two distinct elements.

Finally, the analysis has no temporal or global dimension. Acquisition events
cannot be dated, and no comparison was made to international collections, so
whether these plasmids were introduced from elsewhere or have circulated
regionally remains unresolved.

### Future work

The most direct extensions are long-read sequencing of representative isolates to
resolve complete plasmid structures and the genetic context of the carbapenemase
genes; comparison of the identified plasmids against global collections to
establish their international relationships; and conjugation assays to confirm
that the elements inferred to be mobile are in fact transmissible.

### Conclusion

Carbapenem resistance in this West African *K. pneumoniae* collection is
principally driven by repeated horizontal acquisition rather than clonal
dissemination, mediated by two independently circulating plasmid lineages with
markedly different accessory resistance content, one of which is demonstrably the
same element in unrelated host lineages. Interventions directed at individual
clones are unlikely to interrupt this process. Separately, the collection's
structure as an aggregation of heterogeneous studies materially affects what can
be inferred from it, and conditioning on contributing study should be regarded as
a prerequisite for analyses of comparable public genome collections.

---

## Tables

**Table 1.** Acquired-carbapenemase prevalence by contributing BioProject
(studies with n ≥ 15).

| Study | n | Carbapenemase-positive | Years | Countries |
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

**Table 2.** Accessory resistance genes distinguishing the two carbapenemase
families.

| Gene | NDM (n = 60) | OXA-48-like (n = 38) | Carbapenemase-negative (n = 458) |
|---|---|---|---|
| *armA* | 25 | 0 | 3 |
| *rmtB* | 0 | 28 | 0 |
| *aph(3′)-VI* | 27 | 0 | 0 |
| *mphE* / *msrE* | 30 | 0 | 3 |
| *tet(G)* | 0 | 28 | 0 |
| *aac(3)-IId* | 6 | 24 | 16 |

**Table 3.** MOB-suite clusters carrying carbapenemases across the 19-genome
mobilization subset.

| Cluster | Carbapenemase | Sequence types | Genomes | Cross-lineage |
|---|---|---|---|---|
| AA038 | OXA-181 | ST17, ST234, ST340 | 6 | Yes |
| AA405 | NDM-1 | ST147, ST307, ST340, ST392 | 4 | Yes |
| AA002 | OXA-48 | ST307, ST392 | 2 | Yes |
| AA019 | NDM-1 | ST147 | 3 | No |
