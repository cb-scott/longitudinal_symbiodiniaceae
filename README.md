# Longitudinal Coral-Symbiont Flexibility

This repository contains a curated longitudinal database of coral-Symbiodiniaceae
associations and the reproducible analysis pipeline used to evaluate symbiont
turnover across environmental change, coral phylogeny, and coral functional
traits. The analyses operate at two complementary resolutions:

- dominant symbiont genus states, used to classify genus-level turnover and
  estimate turnover, reversion, and pairwise transition probabilities; and
- SymPortal defining intragenomic variants (DIVs), used to quantify continuous
  community change with Bray-Curtis and weighted UniFrac distances.

Contact: Carly Scott, cabscott@unc.edu

> [!NOTE]
> Do you have a publication or colony-level dataset with repeated symbiont
> observations through time? Please contact me at cabscott@unc.edu or use this
> [data contribution form](https://forms.gle/nG9z4qEJfJ7kMbBUA). Published data
> will be cited, and contributors of unpublished data are welcome to discuss
> participation in future analyses.

## Analysis overview

1. Published supplemental data, figure-derived records, and raw ITS2
   datasets are standardized into the longitudinal database.
2. SymPortal DIV abundance tables are harmonized, linked to colony metadata,
   and converted into sequential Bray-Curtis and weighted UniFrac distances.
3. Satellite sea-surface temperature, heat-stress, ocean-color, and
   precipitation data are summarized around each sampling location and date.
   Environmental variables are reduced to three principal components.
4. Sequential genus observations are classified as no change, emergent
   co-dominance, or complete genus replacement.
5. Bayesian hierarchical models estimate environmental effects and partition
   host-species variance among phylogenetic, trait-associated, and residual
   components. The main community-distance model uses DIV-level Bray-Curtis
   distance, with weighted UniFrac as a sensitivity analysis.
6. Multistate Markov models estimate annual turnover and reversion
   probabilities and pairwise transition rates among dominant symbiont genera.

## Repository structure

```text
db/                       dated releases of the longitudinal database
data/                     working database, source data, intermediate objects,
                          environmental extractions, and saved model fits
Figures/                  final figures published with the paper
remote_run_resources/     resources used for remote preprocessing
results/                  figures and tables grouped by producing script
scripts/                  numbered pipeline scripts, shared helpers, and Stan code
supplemental_litsearchlogs/ search strategies and records from the literature review
```

## The longitudinal genus-level database

The database has one row per coral-colony observation at a sampling time. The
current release contains 10,398 observations from 2,995 colonies, 67 coral
species, and 42 study labels, spanning August 1998 through April 2024.

Two equivalent copies are provided for different purposes:

- `db/LongitudinalSymbiontDB_Nov172026.xlsx` is the dated version of the database which was used in this
  analysis (current release).
- `data/Long_Symbiont_Manual.xlsx` is the working copy read by the analysis
  pipeline.

### Genus-level database columns

| Column | Type | Description |
|---|---|---|
| `SPECIES` | character | Host coral species in uppercase `GENUS_SPECIES` format, for example `PORITES_LOBATA`. |
| `COLONY_ID` | character | Unique colony identifier, generally prefixed by a study code. |
| `STUDY` | character | Source-study citation or study identifier. |
| `DATE` | Excel date | Sampling date. Exact dates are retained when reported; month-only observations use the first of the month, and approximate mid- or end-of-month observations use the 15th or 30th. |
| `GENUS` | character | Dominant Symbiodiniaceae genus. Co-dominant genera are joined with underscores and sorted alphabetically. |
| `SUBTYPE` | character | Legacy dominant species/type label retained for provenance and backward compatibility; see the note below. |
| `BLEACHING_OBS` | character | Observed bleaching status (`Y`, `N`, or blank/not reported). |
| `DEAD` | character | Colony mortality status (`Y`, `N`, or blank/not reported). |
| `DEPTH_M` | numeric | Sampling depth in meters. |
| `LATITUDE` | numeric | Latitude in decimal degrees. |
| `LONGITUDE` | numeric | Longitude in decimal degrees. |
| `NOTES` | character | Study-design, sampling, or data-extraction notes. |
| `NOTES2` | character | Additional notes. |

### Legacy `SUBTYPE` field

`SUBTYPE` records the dominant symbiont species or type labels available from
the source studies. When relative-abundance data were available, a lineage was
treated as dominant when it exceeded 90% of the community; samples were treated
as co-dominant when two or more lineages each exceeded 10%. For figure-only
datasets, the same 10% threshold was estimated from plotted proportions. Where
possible, strain-level suffixes were removed to improve consistency among
studies, because older ITS2 labels and newer SymPortal type-profile identifiers
are not always directly comparable. Multiple retained labels were joined with
underscores.

This field is now considered legacy. It may remain in intermediate objects for
consistency, but it is not used as the within-genus response in the current
models. Within-genus community change is instead quantified from the underlying
DIV relative-abundance data using Bray-Curtis and weighted UniFrac distances.

## DIV-level database

`data/ColonyDistances.Rdata` is a second data product for
analyses on finer-scale Symbiodiniaceae community turnover. Its
primary object, `all.colonydists`, contains 3,127 sequential colony intervals
and includes the following columns:

### DIV-level database columns

| Column | Type | Description |
|---|---|---|
| `STUDY` | character | Short study identifier used to link the interval to its source dataset. |
| `COLONY_ID` | character | Unique coral-colony identifier. |
| `SPECIES` | character | Host coral species in uppercase `GENUS_SPECIES` format. |
| `FROM_DATE` | date | Starting date of the sequential sampling interval. |
| `TO_DATE` | date | Ending date of the sequential sampling interval. |
| `FROM_DOM_DIV` | character | Dominant DIV at the start of the interval. |
| `FROM_DOM_DIV_RA` | numeric | Relative abundance of the starting dominant DIV at `FROM_DATE`. |
| `TO_DOM_DIV` | character | Dominant DIV at the end of the interval. |
| `TO_DOM_DIV_RA` | numeric | Relative abundance of the ending dominant DIV at `TO_DATE`. |
| `TO_DIV` | character | Ending dominant DIV retained as the lookup key for its abundance at the start of the interval. |
| `TO_DIV_FROM_RA` | numeric | Relative abundance at `FROM_DATE` of the DIV that is dominant at `TO_DATE`; zero indicates that it was not detected initially. |
| `DIST` | numeric | Bray-Curtis distance between the sequential `FROM_DATE` and `TO_DATE` communities. |
| `DIST_INITIAL` | numeric | Bray-Curtis distance between the colony's initial community and the community at `TO_DATE`. |
| `UNIFRAC` | numeric | Weighted UniFrac distance between the sequential `FROM_DATE` and `TO_DATE` communities. |
| `UNIFRAC_INITIAL` | numeric | Weighted UniFrac distance between the colony's initial community and the community at `TO_DATE`. |
| `UNIFRAC_METHOD` | character | Method used for the sequential UniFrac value: sequence-based calculation or calibrated hybrid proxy. |
| `UNIFRAC_INITIAL_METHOD` | character | Method used for the initial-reference UniFrac value: sequence-based calculation or calibrated hybrid proxy. |
| `last_genus` | character | Dominant or co-dominant genus state at the start of the sequential interval. |
| `initial_genus` | character | Dominant or co-dominant genus state at the colony's first observation. |
| `LATITUDE` | numeric | Sampling latitude in decimal degrees. |
| `LONGITUDE` | numeric | Sampling longitude in decimal degrees. |

The starting and ending DIV identities are finer-resolution observations that
are not converted into discrete DIV-switch calls by this pipeline. The current
models use the continuous Bray-Curtis and weighted UniFrac distance columns,
while the retained DIV identities allow other researchers to ask taxon-specific
questions that were outside the scope of this analysis.

The file also contains long-format taxon-change, turnover, turnover-partition,
and method-summary objects generated during DIV harmonization.

## Pipeline order

Scripts are numbered by analysis stage. Letter suffixes indicate execution
order within a stage.

| Prefix | Stage |
|---|---|
| `00` | Raw-read and database preprocessing |
| `01` | Environmental-data extraction and PCA |
| `02` | Analysis-dataset cleaning and genus-level transition calling |
| `03` | Bayesian hierarchical models |
| `04` | Multistate Markov models |

The `00` scripts document how the databases were assembled and generally do not
need to be rerun to reproduce the analysis. Scripts `01` through `04` use
checked-in intermediate data products and should be run from the repository
root.

| Script | Main inputs | Main outputs |
|---|---|---|
| `00a_RawReads_SymPortal_Preprocessing.sh` | NCBI/SRA records and raw ITS2 reads | DADA2 taxonomy objects, SymPortal outputs, and DIV phylogeny inputs |
| `00b_RawData_Cleaning.R` | Study supplemental data and metadata | Study-specific cleaned tables used to assemble `data/Long_Symbiont_Manual.xlsx` |
| `00c_SymPortal_BrayCurtis_UniFrac_Distances.R` | SymPortal abundance tables, sample metadata, and DIV phylogeny | `data/ColonyDistances.Rdata` |
| `01a_EnvironmentalData_Extraction_PCA.R` | Longitudinal database, DIV intervals, and satellite products | `data/EnvData_PCs.Rdata` and `results/01a_EnvironmentalData_Extraction_PCA/` |
| `02a_AnalysisDataset_Cleaning.R` | Longitudinal database and environmental PCs | `data/LongZooxInputData_wEnv.Rdata` and `results/02a_AnalysisDataset_Cleaning/` |
| `03a_Bayesian_GenusOrdinal.R` | Genus-level intervals, host phylogeny, traits, and environmental PCs | Main genus-level ordinal model results in `results/03a_Bayesian_GenusOrdinal/` |
| `03b_Bayesian_GenusOrdinal_ReducedTimeframe.R` | Same inputs as `03a`, restricted to the DIV sampling period | Reduced-timeframe sensitivity results in `results/03b_Bayesian_GenusOrdinal_ReducedTimeframe/` |
| `03c_Bayesian_FirstTurnoverBernoulli.R` | Genus-level intervals, host phylogeny, traits, and environmental PCs | First-turnover Bernoulli sensitivity results in `results/03c_Bayesian_FirstTurnoverBernoulli/` |
| `03d_Bayesian_DIVBrayCurtis.R` | DIV intervals, environmental PCs, host phylogeny, and traits | Main DIV-distance model results in `results/03d_Bayesian_DIVBrayCurtis/` |
| `03e_Bayesian_DIVUniFrac.R` | DIV intervals, environmental PCs, host phylogeny, and traits | Weighted-UniFrac sensitivity results in `results/03e_Bayesian_DIVUniFrac/` |
| `03f_Bayesian_GenusOrdinal_InitialSingleGenus.R` | Same inputs as `03a`, excluding colonies whose first genus state was mixed | Initial-single-genus sensitivity results in `results/03f_Bayesian_GenusOrdinal_InitialSingleGenus/` |
| `03g_Bayesian_GenusOrdinal_NoLastGenus.R` | Same inputs as `03a`, omitting prior genus as a predictor | No-prior-genus sensitivity results in `results/03g_Bayesian_GenusOrdinal_NoLastGenus/` |
| `03h_Bayesian_GenusOrdinal_TraitPCInteractions.R` | Genus-level intervals, TraitPC1-TraitPC3, environmental PCs, and host phylogeny | Trait-by-environment interaction results in `results/03h_Bayesian_GenusOrdinal_TraitPCInteractions/` |
| `04a_MSM_Modeling.R` | Genus-level intervals with environmental PCs | turnover/reversion and pairwise-transition results in `results/04a_MSM_Modeling/` |

## Script details

### `00a_RawReads_SymPortal_Preprocessing.sh`

Documents the remote and local steps used to retrieve raw reads, process ITS2
amplicons with DADA2, rerun selected studies through SymPortal, and construct
the DIV phylogeny used for weighted UniFrac. This is a reproducibility record,
not a single portable command: its sections were run in different local and
cluster environments and contain paths that must be adapted before reuse.

External requirements include Entrez Direct, SRA Toolkit, DADA2, SymPortal,
MAFFT, and FastTree.

### `00b_RawData_Cleaning.R`

Cleans study-specific supplemental datasets into the 13-column database schema.
Only colonies observed more than once are retained for longitudinal analysis.
Some studies were extracted manually from figures and therefore do not have a
corresponding automated cleaning block.

R dependencies: `tidyverse`, `readxl`, `writexl`, `phyloseq`,
`janitor`, `vegan`, and `pheatmap`.

### `00c_SymPortal_BrayCurtis_UniFrac_Distances.R`

Combines DIV-level SymPortal outputs across studies, assigns DIVs to genera,
checks phylogenetic coverage, identifies sequential samples from the same
colony, and computes Bray-Curtis and weighted UniFrac community distances. It
also retains the underlying start/end DIV identities and relative abundances in
`data/ColonyDistances.Rdata`.

R dependencies: `tidyverse`, `readxl`, `phyloseq`, `ape`, `phangorn`,
`vegan`, and `janitor`.

### `01a_EnvironmentalData_Extraction_PCA.R`

Queries environmental data for each unique sampling location and date. The
script summarizes nearby valid ocean cells consistently through time, caches
the extraction results, imputes the remaining missing environmental values,
and reduces the standardized variables to three environmental PCs. It also
produces maps and summaries of the geographic, temporal, taxonomic, and study
coverage of the database.

Environmental sources include NOAA Coral Reef Watch temperature, anomaly, and
degree-heating-week products; ESA Ocean Colour Climate Change Initiative
chlorophyll-a and light-attenuation products; and IMERG precipitation.

R dependencies: `tidyverse`, `lubridate`, `missForest`, `future`,
`furrr`, `rerddap`, `sf`, `geosphere`, `pheatmap`, and `writexl`.

### `02a_AnalysisDataset_Cleaning.R`

Joins environmental PCs to the longitudinal observations, standardizes selected
host-species names, summarizes study coverage, and constructs sequential
genus-level intervals. The current turnover outcome distinguishes no genus
change, emergent co-dominance, and complete genus replacement; legacy subtype
differences are not treated as modeled turnover events.

R dependencies: `tidyverse`, `readxl`, and `writexl`.

### Bayesian hierarchical models (`03a`-`03h`)

All Bayesian scripts use `brms` to prepare model data and custom Stan programs
for fitting. Shared dependencies include `cmdstanr`, `posterior`, `ape`,
`phangorn`, `phytools`, `Matrix`, `ggdist`, `ggtree`, `patchwork`, `Ternary`, `viridisLite`,
and `writexl`. Saved model-fit objects in `data/` allow post-processing to be
reproduced without refitting every model.

#### `03a_Bayesian_GenusOrdinal.R`

Fits the genus-level cumulative ordinal model. Fixed effects include elapsed
time between observations, changes in environmental PC1-PC3, and the previous
dominant symbiont genus. Species-level intercepts and environmental slopes are
partitioned among host phylogeny, trait similarity, and residual species
variation, with additional grouping structure for study, site, and colony.

#### `03b_Bayesian_GenusOrdinal_ReducedTimeframe.R`

Repeats the genus-level ordinal analysis over the shorter period represented by
the DIV-level dataset. This provides a time-window sensitivity analysis for
comparison with the community-distance models.

#### `03c_Bayesian_FirstTurnoverBernoulli.R`

Restricts each colony to observations through its first modeled genus-level
turnover and fits a Bernoulli sensitivity model for comparison with annual
transition probabilities from the Markov models. It uses the same fixed effects,
hierarchical grouping structure, and phylogenetic, trait-associated, and residual
species-variance partition as the ordinal model. The script loads the saved fit
from `data/stan_dirichlet_explicit_bernoulli_first.RDS` when available, so the
posterior summaries and figures can be reproduced without refitting the model.

#### `03d_Bayesian_DIVBrayCurtis.R`

Fits the primary continuous community-change model using sequential DIV-level
Bray-Curtis distance. It uses the same environmental predictors and the same
phylogenetic, trait-associated, and residual species-variance partition as the
genus-level model. Outputs include fixed and random effects, variance
partitions, expected species-level distances, and species-specific
environmental-response summaries.

#### `03e_Bayesian_DIVUniFrac.R`

Repeats the DIV-distance analysis using weighted UniFrac, incorporating DIV
phylogenetic relatedness as a sensitivity analysis.

#### `03f_Bayesian_GenusOrdinal_InitialSingleGenus.R`

Repeats the genus-level ordinal analysis after excluding colonies whose first
observed genus state was mixed. This sensitivity analysis tests whether
initially co-dominant communities drive the main genus-level patterns.

#### `03g_Bayesian_GenusOrdinal_NoLastGenus.R`

Repeats the genus-level ordinal analysis while omitting previous dominant genus
as a predictor. This tests whether prior-genus effects mask associations between
environmental change and turnover.

#### `03h_Bayesian_GenusOrdinal_TraitPCInteractions.R`

Tests whether TraitPC1-TraitPC3 predict species differences in genus-level
responses to environmental PC1-PC3. The trait PCs and all nine
trait-by-environment interactions enter as fixed effects. Because traits are
tested explicitly, this model omits the trait-covariance random-effect block
used by `03a`; remaining species variation is partitioned between phylogenetic
and idiosyncratic components.

### `04a_MSM_Modeling.R`

Fits two multistate Markov analyses. The first estimates annual probabilities
of leaving and returning to a colony's original dominant-genus state as a
function of environmental PC changes. The second estimates pairwise transition
rates among the dominant Symbiodiniaceae genera while treating mixed states as
censored. Outputs include hazard-ratio plots, annual transition-probability
summaries, transition networks, and pairwise rate comparisons.

R dependencies: `msm`, `tidyverse`, and `igraph`.

## Shared code and Stan models

`scripts/helpers.R` contains only the two shared functions called by the active
numbered analyses: `assign_ocean_basin()` for geographic classification and
`extract_re()` for extracting species-level random-effect draws from the custom
Stan models.

The active Stan programs are:

- `scripts/coral_dirichlet_partition_site.stan` for the genus-level ordinal
  model;
- `scripts/coral_dirichlet_partition_bernoulli_site.stan` for the first-turnover
  Bernoulli model;
- `scripts/coral_dirichlet_partition_distance_beta.stan` for the Bray-Curtis
  and weighted-UniFrac distance models;
- `scripts/coral_dirichlet_partition_site_no_last_genus.stan` for the
  no-prior-genus sensitivity model; and
- `scripts/coral_traitpc_interactions_ordinal.stan` for the explicit TraitPC by
  environmental-PC interaction model.

## Reproducing the analyses

1. Clone the repository and open it at the repository root.
2. Confirm that the required R packages, CmdStan installation, and external
   preprocessing tools are available for the stage you intend to run.
3. Use the checked-in database and intermediate objects unless you specifically
   intend to rebuild them from source studies.
4. Run scripts in numeric and letter order. Each analysis script creates its
   own subdirectory under `results/`.
5. For Bayesian scripts, either fit the models or load the corresponding saved
   `.RDS` objects before running posterior summaries and figures.

The analysis scripts use project-relative paths and are intended to be run from
the repository root.

