#Need each file and it's metatdata
library(tidyverse)
library(readxl)
library(vegan)
#the plan: for each study, want to calculate the bray-curtis community dissimilarity between timepoints
#we will also need to join the enviornmental data, which we should have for every study if we use the overall metadata
#Some studies are new; we will need to rerun the enviornmental data generation and add them to the main database
#we'll then model community dissimilarity as a function of environmental disimmilarity w the same bayesian methods

#------------------------- Lee et al. 2026 ----------------------------------
#New inclusion, needs to be added to overall database too

lee.divs <- read.table("data/symportal_profiles/Lee2026_DIVs_corr.txt", sep = "\t", header = T)
lee.meta <- read_xlsx("data/symportal_profiles/Lee_Symportal_Metadata.xlsx", skip=1)
head(lee.meta)
relevant.meta = lee.meta %>% select(sample_name, colonyID, collection_date) %>%
  mutate(collection_date = as.Date(collection_date, format = "%Y%m%d"))

#5.611236805538839, 103.06104205143536 coordinates for study

#first do colony level calls to include in database
lee.dists <- lee.divs %>%
  select(sample_name, A6b:G3x) %>% data.frame() %>% column_to_rownames(var = "sample_name")

# All within-colony pairwise distances
all.dists <- as.matrix(vegdist(lee.dists, method = "bray")) %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(!from, names_to = "to", values_to = "dist") %>%
  left_join(relevant.meta, by = c("from" = "sample_name")) %>%
  left_join(relevant.meta, by = c("to" = "sample_name")) %>%
  filter(colonyID.x == colonyID.y) %>%
  group_by(colonyID.x) %>%
  mutate(
    date.rank.x = dense_rank(collection_date.x),
    date.rank.y = dense_rank(collection_date.y)
  ) %>%
  ungroup()

# Distance between sequential observations
sequential <- all.dists %>%
  filter(date.rank.y == date.rank.x + 1) %>%
  select(colonyID.x, collection_date.x, collection_date.y, dist)

# Distance from initial community to each subsequent observation
initial <- all.dists %>%
  filter(date.rank.x == 1, date.rank.y > 1) %>%
  select(colonyID.x, collection_date.y, DIST_INITIAL = dist)

# Combine and format
lee.colonydists <- sequential %>%
  left_join(initial, by = c("colonyID.x", "collection_date.y")) %>%
  arrange(colonyID.x, collection_date.x) %>%
  mutate(colonyID.x = paste("LEE", "ADIG", colonyID.x, sep = "_")) %>%
  rename(
    COLONYID = colonyID.x,
    FROM_DATE = collection_date.x,
    TO_DATE = collection_date.y,
    DIST = dist
  ) %>%
  mutate(
    LATITUDE = 5.611236805538839,
    LONGITUDE = 103.06104205143536
  ) %>%
  select(COLONYID, FROM_DATE, TO_DATE, DIST, DIST_INITIAL,
         LATITUDE, LONGITUDE)

lee.colonydists$SPECIES = "ACROPORA_DIGITIFERA"

#----------- Terraneo et al., 2026 --------------
terr.divs <- read.table("data/symportal_profiles/Terraneo_DIVs.txt", sep = "\t", header = T)
terr.meta <- read_xlsx("data/symportal_profiles/Terraneo_SymPortal_Metadata.xlsx", skip=1)
terr.meta2 <- read_excel("data/Tarraneo_Meta.xlsx")
terr.meta2 = terr.meta2[-1,]
terr.meta2 = janitor::clean_names(terr.meta2)

terr.meta.full = terr.meta %>% select(colonyID, collection_date, sample_name) %>% 
  left_join(terr.meta2, by= c("colonyID" = "tagged_colony")) %>% 
  select(colonyID, collection_date, sample_name, genus, x4, sampling_locality) %>% distinct() %>%
  mutate(SPECIES = toupper(paste(genus, x4, sep = "_"))) %>%
  mutate(DATE = as.Date(paste(collection_date, "01", sep = ""), format = "%Y%m%d")) %>%
  mutate(COLONY_ID = paste("TERR", toupper(substr(genus, 1, 1)),
                           toupper(substr(x4, 1, 3)), colonyID, sep = "_")) 
terr.meta.full
lonlat.key = data.frame(sampling_locality = sort(unique(terr.meta.full$sampling_locality)), 
                        LATITUDE = c(-22.477683, -22.2295633, -22.285833, -22.285833, -22.312217),
                        LONGITUDE = c(166.463033, 166.435, 166.183167, 166.183167,166.633833))
terr.meta.full = terr.meta.full %>% left_join(lonlat.key)
relevant.meta = terr.meta.full %>% select(sample_name, COLONY_ID, DATE) %>% distinct()
colony.meta = terr.meta.full %>% select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>% distinct()


terr.dists <- terr.divs %>%
  select(sample_name, A1:X28138_D) %>% data.frame() %>% column_to_rownames(var = "sample_name")

# All within-colony pairwise distances
all.dists <- as.matrix(vegdist(terr.dists, method = "bray")) %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(!from, names_to = "to", values_to = "dist") %>%
  left_join(relevant.meta, by = c("from" = "sample_name")) %>%
  left_join(relevant.meta, by = c("to" = "sample_name")) %>%
  filter(COLONY_ID.x == COLONY_ID.y) %>%
  group_by(COLONY_ID.x) %>%
  mutate(
    date.rank.x = dense_rank(DATE.x),
    date.rank.y = dense_rank(DATE.y)
  ) %>%
  ungroup() %>% distinct()

# Distance between sequential observations
sequential <- all.dists %>%
  filter(date.rank.y == date.rank.x + 1) %>%
  select(COLONY_ID.x, DATE.x, DATE.y, dist) %>% distinct()

# Distance from initial community to each subsequent observation
initial <- all.dists %>%
  filter(date.rank.x == 1, date.rank.y > 1) %>%
  select(COLONY_ID.x, DATE.y, DIST_INITIAL = dist) %>% distinct()

# Combine and format
terrcolonydists <- sequential %>%
  left_join(initial, by = c("COLONY_ID.x", "DATE.y")) %>%
  arrange(COLONY_ID.x, DATE.x) %>%
  rename(
    COLONY_ID = COLONY_ID.x,
    FROM_DATE = DATE.x,
    TO_DATE = DATE.y,
    DIST = dist
  ) %>% left_join(colony.meta)

terrcolonydists

#----------- Claar et al., 2020 ---------
#should have been many more samples in this raw rerun; not sure what the deal is with that.


claar.divs <- read.table(
  "data/symportal_profiles/Claar_DIVs.txt",
  sep = "\t", header = TRUE, check.names = FALSE
)

claar.meta <- read_xlsx(
  "data/symportal_profiles/Claar_Symportal_BioPRJ_Meta.xlsx",
  skip = 1
)

# Prepare metadata
claar.meta.full <- claar.meta %>%
  select(
    sample_name, colonyID, collection_date,
    host_genus, host_species,
    collection_latitude, collection_longitude
  ) %>%
  distinct() %>%
  mutate(
    SPECIES = toupper(paste(host_genus, host_species, sep = "_")),
    DATE = as.Date(
      paste(collection_date, "01", sep = ""),
      format = "%Y%m%d"
    ),
    COLONY_ID = paste(
      "CLAAR",
      toupper(substr(host_genus, 1, 1)),
      toupper(substr(host_species, 1, 3)),
      colonyID, sep = "_"
    ),
    LATITUDE = collection_latitude,
    LONGITUDE = collection_longitude
  )

relevant.meta <- claar.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- claar.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare community composition matrix
claar.dists <- claar.divs %>%
  filter(!is.na(sample_name)) %>%
  semi_join(relevant.meta, by = "sample_name") %>%
  select(sample_name, A5b:D17c) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# All within-colony pairwise distances
all.dists <- as.matrix(vegdist(claar.dists, method = "bray")) %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(!from, names_to = "to", values_to = "dist") %>%
  left_join(relevant.meta, by = c("from" = "sample_name")) %>%
  left_join(relevant.meta, by = c("to" = "sample_name")) %>%
  filter(COLONY_ID.x == COLONY_ID.y) %>%
  group_by(COLONY_ID.x) %>%
  mutate(
    date.rank.x = dense_rank(DATE.x),
    date.rank.y = dense_rank(DATE.y)
  ) %>%
  ungroup() %>%
  distinct()

# Distance between sequential observations
sequential <- all.dists %>%
  filter(date.rank.y == date.rank.x + 1) %>%
  select(COLONY_ID.x, DATE.x, DATE.y, dist) %>%
  distinct()

# Distance from initial community to each subsequent observation
initial <- all.dists %>%
  filter(date.rank.x == 1, date.rank.y > 1) %>%
  select(COLONY_ID.x, DATE.y, DIST_INITIAL = dist) %>%
  distinct()

# Combine and format
claarcolonydists <- sequential %>%
  left_join(initial, by = c("COLONY_ID.x", "DATE.y")) %>%
  arrange(COLONY_ID.x, DATE.x) %>%
  rename(
    COLONY_ID = COLONY_ID.x,
    FROM_DATE = DATE.x,
    TO_DATE = DATE.y,
    DIST = dist
  ) %>%
  left_join(colony.meta, by = "COLONY_ID") %>%
  select(
    COLONY_ID, FROM_DATE, TO_DATE,
    DIST, DIST_INITIAL,
    LATITUDE, LONGITUDE, SPECIES
  )

claarcolonydists

#-------- Buzzoni et al., 2026 --------
#New, add to DB as well


buzz.divs <- read.csv("data/symportal_profiles/Buzzoni2026_DIVs.csv")
buzz.meta <- read.csv("data/symportal_profiles/Buzzoni2026_metadata.csv")

# Only retain samples with both metadata and DIVs
shared.samps <- intersect(buzz.divs$sample_name, buzz.meta$sample_name)

# Prepare metadata
buzz.meta.full <- buzz.meta %>%
  filter(sample_name %in% shared.samps) %>%
  mutate(
    SPECIES = "PLATYGYRA_RYUKYUENSIS",
    DATE = as.Date(as.character(collection_date), format = "%Y%m%d"),
    COLONY_ID = paste("BUZZ", "PRYUK", coral_tag, sep = "_"),
    LATITUDE = as.numeric(collection_latitude),
    LONGITUDE = as.numeric(collection_longitude)
  ) %>%
  distinct()

relevant.meta <- buzz.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- buzz.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare community composition matrix
buzz.dists <- buzz.divs %>%
  filter(sample_name %in% shared.samps) %>%
  select(sample_name, A5b:`X1269165_G`) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Convert counts to relative abundances
buzz.dists <- sweep(buzz.dists, 1, rowSums(buzz.dists), "/")

# All within-colony pairwise distances
all.dists <- as.matrix(vegdist(buzz.dists, method = "bray")) %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(!from, names_to = "to", values_to = "dist") %>%
  left_join(relevant.meta, by = c("from" = "sample_name")) %>%
  left_join(relevant.meta, by = c("to" = "sample_name")) %>%
  filter(COLONY_ID.x == COLONY_ID.y) %>%
  group_by(COLONY_ID.x) %>%
  mutate(
    date.rank.x = dense_rank(DATE.x),
    date.rank.y = dense_rank(DATE.y)
  ) %>%
  ungroup() %>%
  distinct()

# Distance between sequential observations
sequential <- all.dists %>%
  filter(date.rank.y == date.rank.x + 1) %>%
  select(COLONY_ID.x, DATE.x, DATE.y, dist) %>%
  distinct()

# Distance from initial community to each subsequent observation
initial <- all.dists %>%
  filter(date.rank.x == 1, date.rank.y > 1) %>%
  select(COLONY_ID.x, DATE.y, DIST_INITIAL = dist) %>%
  distinct()

# Combine and format
buzzcolonydists <- sequential %>%
  left_join(initial, by = c("COLONY_ID.x", "DATE.y")) %>%
  arrange(COLONY_ID.x, DATE.x) %>%
  rename(
    COLONY_ID = COLONY_ID.x,
    FROM_DATE = DATE.x,
    TO_DATE = DATE.y,
    DIST = dist
  ) %>%
  left_join(colony.meta, by = "COLONY_ID") %>%
  select(
    COLONY_ID, FROM_DATE, TO_DATE,
    DIST, DIST_INITIAL,
    LATITUDE, LONGITUDE, SPECIES
  )

buzzcolonydists

#------------------ Claar et al., 2020 ------------------
#some overlap with Buzzoni; keep only Favites, worried about double counting platygyra

claar.divs <- read.table(
  "data/symportal_profiles/Claar_DIVs.txt",
  sep = "\t", header = TRUE
)

claar.meta <- read_xlsx(
  "data/symportal_profiles/Claar_Symportal_BioPRJ_Meta.xlsx",
  skip = 1
)

# Prepare metadata, excluding Platygyra
claar.meta.full <- claar.meta %>%
  filter(tolower(host_genus) != "platygyra") %>%
  mutate(
    SPECIES = toupper(paste(host_genus, host_species, sep = "_")),
    DATE = as.Date(
      paste(collection_date, "01", sep = ""),
      format = "%Y%m%d"
    ),
    COLONY_ID = paste(
      "CLAAR",
      toupper(substr(host_genus, 1, 1)),
      toupper(substr(host_species, 1, 3)),
      colonyID, sep = "_"
    ),
    LATITUDE = as.numeric(collection_latitude),
    LONGITUDE = as.numeric(collection_longitude)
  ) %>%
  distinct()

relevant.meta <- claar.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- claar.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare community composition matrix
claar.dists <- claar.divs %>%
  filter(!is.na(sample_name)) %>%
  semi_join(relevant.meta, by = "sample_name") %>%
  select(sample_name, A5b:D17c) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# All within-colony pairwise distances
all.dists <- as.matrix(vegdist(claar.dists, method = "bray")) %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(!from, names_to = "to", values_to = "dist") %>%
  left_join(relevant.meta, by = c("from" = "sample_name")) %>%
  left_join(relevant.meta, by = c("to" = "sample_name")) %>%
  filter(COLONY_ID.x == COLONY_ID.y) %>%
  group_by(COLONY_ID.x) %>%
  mutate(
    date.rank.x = dense_rank(DATE.x),
    date.rank.y = dense_rank(DATE.y)
  ) %>%
  ungroup() %>%
  distinct()

# Distance between sequential observations
sequential <- all.dists %>%
  filter(date.rank.y == date.rank.x + 1) %>%
  select(COLONY_ID.x, DATE.x, DATE.y, dist) %>%
  distinct()

# Distance from initial community to each subsequent observation
initial <- all.dists %>%
  filter(date.rank.x == 1, date.rank.y > 1) %>%
  select(COLONY_ID.x, DATE.y, DIST_INITIAL = dist) %>%
  distinct()

# Combine and format
claarcolonydists <- sequential %>%
  left_join(initial, by = c("COLONY_ID.x", "DATE.y")) %>%
  arrange(COLONY_ID.x, DATE.x) %>%
  rename(
    COLONY_ID = COLONY_ID.x,
    FROM_DATE = DATE.x,
    TO_DATE = DATE.y,
    DIST = dist
  ) %>%
  left_join(colony.meta, by = "COLONY_ID") %>%
  select(
    COLONY_ID, FROM_DATE, TO_DATE,
    DIST, DIST_INITIAL,
    LATITUDE, LONGITUDE, SPECIES
  )

claarcolonydists


#------------- Van Nynatten et al., 2025 ----------- 

vn.divs <- read.table(
  "data/symportal_profiles/van_nynatten_DIVs_corr.txt",
  sep = "\t", header = TRUE, check.names = FALSE
)

vn.seqid <- read.csv("data/van_nynatten_sample_metadata.csv")
vn.sites <- read.csv("data/van_nynatten_survey_metadata.csv")

vn.exps <- read.csv("data/van_nynatten_expeditions.csv") %>%
  select(Expedition.code, Mean.date) %>%
  rename(
    expedition = Expedition.code,
    Mean.date = Mean.date
  ) %>%
  mutate(
    expedition = ifelse(expedition == "2023", "2023b", expedition),
    DATE = as.Date(Mean.date, format = "%m/%d/%Y")
  )

# Only retain samples with both metadata and DIVs
shared.samps <- intersect(
  vn.divs$sample_name,
  vn.seqid$library_ID
)

# Prepare metadata
vn.meta.full <- vn.seqid %>%
  filter(library_ID %in% shared.samps) %>%
  left_join(vn.exps, by = "expedition") %>%
  left_join(vn.sites, by = "site_name") %>%
  mutate(
    sample_name = library_ID,
    SPECIES = toupper(paste(host_genus, host_species, sep = "_")),
    COLONY_ID = paste0(
      "VN_",
      substr(toupper(host_genus), 1, 1),
      substr(toupper(host_species), 1, 3),
      "_", coral_tag
    ),
    LATITUDE = lat,
    LONGITUDE = lon
  ) %>%
  distinct()

relevant.meta <- vn.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- vn.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare community composition matrix
vn.dists <- vn.divs %>%
  filter(sample_name %in% shared.samps) %>%
  select(sample_name, A5z:G3aq) %>%
  filter(rowSums(across(-sample_name)) > 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Convert counts to relative abundances
vn.dists <- sweep(vn.dists, 1, rowSums(vn.dists), "/")

# All within-colony pairwise distances
all.dists <- as.matrix(vegdist(vn.dists, method = "bray")) %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(!from, names_to = "to", values_to = "dist") %>%
  left_join(relevant.meta, by = c("from" = "sample_name")) %>%
  left_join(relevant.meta, by = c("to" = "sample_name")) %>%
  filter(COLONY_ID.x == COLONY_ID.y) %>%
  group_by(COLONY_ID.x) %>%
  mutate(
    date.rank.x = dense_rank(DATE.x),
    date.rank.y = dense_rank(DATE.y)
  ) %>%
  ungroup() %>%
  distinct()

# Distance between sequential observations
sequential <- all.dists %>%
  filter(date.rank.y == date.rank.x + 1) %>%
  select(COLONY_ID.x, DATE.x, DATE.y, dist) %>%
  distinct()

# Distance from initial community to each subsequent observation
initial <- all.dists %>%
  filter(date.rank.x == 1, date.rank.y > 1) %>%
  select(COLONY_ID.x, DATE.y, DIST_INITIAL = dist) %>%
  distinct()

# Combine and format
vncolonydists <- sequential %>%
  left_join(initial, by = c("COLONY_ID.x", "DATE.y")) %>%
  arrange(COLONY_ID.x, DATE.x) %>%
  rename(
    COLONY_ID = COLONY_ID.x,
    FROM_DATE = DATE.x,
    TO_DATE = DATE.y,
    DIST = dist
  ) %>%
  left_join(colony.meta, by = "COLONY_ID") %>%
  select(
    COLONY_ID, FROM_DATE, TO_DATE,
    DIST, DIST_INITIAL,
    LATITUDE, LONGITUDE, SPECIES
  )

vncolonydists

## ------------- Leinbach et al., 2023 ------------------------


lbach.meta <- read.csv("data/Leinbach_profs.csv")

lbach.divs <- read.table(
  "data/symportal_profiles/Leinbach_DIVs.txt",
  sep = "\t", header = TRUE, check.names = FALSE
)

# Prepare metadata
lbach.meta.full <- lbach.meta %>%
  select(Sample, Individual, month, reef_zone) %>%
  distinct() %>%
  mutate(
    DATE = as.Date(
      paste(month, "01", "2019"),
      format = "%B %d %Y"
    ),
    COLONY_ID = paste("LBACH_AHYC", Individual, sep = "_"),
    SPECIES = "ACROPORA_HYACINTHUS",
    LATITUDE = case_when(
      reef_zone == "crest" ~ -17.4741,
      reef_zone == "forereef" ~ -17.4731
    ),
    LONGITUDE = case_when(
      reef_zone == "crest" ~ -149.8173,
      reef_zone == "forereef" ~ -149.8176
    )
  )

relevant.meta <- lbach.meta.full %>%
  select(Sample, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- lbach.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare DIV community composition matrix
lbach.dists <- lbach.divs %>%
  semi_join(relevant.meta, by = c("sample_name" = "Sample")) %>%
  select(sample_name, A1:`1767177_G`) %>%
  filter(rowSums(across(-sample_name)) > 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# All within-colony pairwise distances
all.dists <- as.matrix(vegdist(lbach.dists, method = "bray")) %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(!from, names_to = "to", values_to = "dist") %>%
  left_join(relevant.meta, by = c("from" = "Sample")) %>%
  left_join(relevant.meta, by = c("to" = "Sample")) %>%
  filter(COLONY_ID.x == COLONY_ID.y) %>%
  group_by(COLONY_ID.x) %>%
  mutate(
    date.rank.x = dense_rank(DATE.x),
    date.rank.y = dense_rank(DATE.y)
  ) %>%
  ungroup() %>%
  distinct()

# Distance between sequential observations
sequential <- all.dists %>%
  filter(date.rank.y == date.rank.x + 1) %>%
  select(COLONY_ID.x, DATE.x, DATE.y, dist) %>%
  distinct()

# Distance from initial community to each subsequent observation
initial <- all.dists %>%
  filter(date.rank.x == 1, date.rank.y > 1) %>%
  select(COLONY_ID.x, DATE.y, DIST_INITIAL = dist) %>%
  distinct()

# Combine and format
lbachcolonydists <- sequential %>%
  left_join(initial, by = c("COLONY_ID.x", "DATE.y")) %>%
  arrange(COLONY_ID.x, DATE.x) %>%
  rename(
    COLONY_ID = COLONY_ID.x,
    FROM_DATE = DATE.x,
    TO_DATE = DATE.y,
    DIST = dist
  ) %>%
  left_join(colony.meta, by = "COLONY_ID") %>%
  select(
    COLONY_ID, FROM_DATE, TO_DATE,
    DIST, DIST_INITIAL,
    LATITUDE, LONGITUDE, SPECIES
  )

lbachcolonydists

#------------- HOWELLS ET AL 2020 --------------------

how.divs <- read_xlsx(
  "data/symportal_profiles/Howells_DIVs_Meta.xlsx",
  skip = 1
)

# Prepare metadata

# Published sampling coordinates: Howells et al. (2020)
how.sites <- tibble(
  SITE = c("Delma", "Saadiyat", "Ras Ghanada", "Dibba", "Al Aqah"),
  LATITUDE = c(24.5208, 24.5986, 24.8482, 25.6031, 25.4929),
  LONGITUDE = c(52.2781, 54.4201, 54.6903, 56.3485, 56.3635)
)

# Prepare metadata
how.meta.full <- how.divs %>%
  transmute(
    SITE = trimws(Site),
    SPECIES = toupper(gsub(" ", "_", trimws(Species))),
    ORIGINAL_ID = trimws(`Colony ID`),
    DATE = as.Date(
      paste0("01-", sub("^DNA\\s+", "", Timepoint)),
      format = "%d-%b-%Y"
    )
  ) %>%
  filter(
    !is.na(SITE), SITE != "",
    !is.na(ORIGINAL_ID), ORIGINAL_ID != "",
    !is.na(DATE)
  ) %>%
  mutate(
    SITE = sub("^.*,\\s*", "", SITE),
    COLONY_ID = paste("HOW", ORIGINAL_ID, sep = "_"),
    sample_name = paste(COLONY_ID, format(DATE, "%Y%m"), sep = "_")
  ) %>%
  left_join(how.sites, by = "SITE") %>%
  distinct()

relevant.meta <- how.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- how.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare community composition matrix
# Columns 5 onward contain ITS2 DIV abundances
how.dists <- how.divs %>%
  mutate(
    ORIGINAL_ID = trimws(`Colony ID`),
    COLONY_ID = paste("HOW", ORIGINAL_ID, sep = "_"),
    DATE = as.Date(
      paste0("01-", sub("^DNA\\s+", "", Timepoint)),
      format = "%d-%b-%Y"
    ),
    sample_name = paste(COLONY_ID, format(DATE, "%Y%m"), sep = "_")
  ) %>%
  semi_join(relevant.meta, by = "sample_name") %>%
  select(sample_name, all_of(names(how.divs)[5:ncol(how.divs)])) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Normalize absolute sequence abundances
how.dists <- sweep(how.dists, 1, rowSums(how.dists), "/")

# All within-colony pairwise distances
all.dists <- as.matrix(vegdist(how.dists, method = "bray")) %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(!from, names_to = "to", values_to = "dist") %>%
  left_join(relevant.meta, by = c("from" = "sample_name")) %>%
  left_join(relevant.meta, by = c("to" = "sample_name")) %>%
  filter(COLONY_ID.x == COLONY_ID.y) %>%
  group_by(COLONY_ID.x) %>%
  mutate(
    date.rank.x = dense_rank(DATE.x),
    date.rank.y = dense_rank(DATE.y)
  ) %>%
  ungroup() %>%
  distinct()

# Distance between consecutive observations
sequential <- all.dists %>%
  filter(date.rank.y == date.rank.x + 1) %>%
  select(COLONY_ID.x, DATE.x, DATE.y, dist) %>%
  distinct()

# Distance from initial community
initial <- all.dists %>%
  filter(date.rank.x == 1, date.rank.y > 1) %>%
  select(COLONY_ID.x, DATE.y, DIST_INITIAL = dist) %>%
  distinct()

# Combine and format
howcolonydists <- sequential %>%
  left_join(initial, by = c("COLONY_ID.x", "DATE.y")) %>%
  arrange(COLONY_ID.x, DATE.x) %>%
  rename(
    COLONY_ID = COLONY_ID.x,
    FROM_DATE = DATE.x,
    TO_DATE = DATE.y,
    DIST = dist
  ) %>%
  left_join(colony.meta, by = "COLONY_ID") %>%
  select(
    COLONY_ID, FROM_DATE, TO_DATE,
    DIST, DIST_INITIAL,
    LATITUDE, LONGITUDE, SPECIES
  )

howcolonydists

##------------ JUNG ET AL 2021 ------------------

jung.divs <- readxl::read_excel(
  "data/Jung_2021_Symportal_post_med_absolut abundances_meta.xlsx"
)

# Prepare metadata
jung.meta.full <- jung.divs %>%
  filter(!is.na(`Sample ID`)) %>%
  select(`Sample ID`) %>%
  separate(
    `Sample ID`,
    into = c("Mon", "Tidal", "Bleached", "Individual"),
    sep = " ",
    remove = FALSE
  ) %>%
  mutate(
    sample_name = `Sample ID`,
    COLONY_ID = paste0("Jung_AASP_", Individual),
    SPECIES = "ACROPORA_ASPERA",
    DATE = case_when(
      Mon == "Apr" ~ as.Date("2016-04-01"),
      Mon == "Nov" ~ as.Date("2016-11-01")
    ),
    LATITUDE = -16.479,
    LONGITUDE = 123.0435
  ) %>%
  distinct()

relevant.meta <- jung.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- jung.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare DIV community composition matrix
jung.dists <- jung.divs %>%
  filter(!is.na(`Sample ID`)) %>%
  select(!sample_name) %>%
  rename(sample_name = `Sample ID`) %>%
  select(sample_name, `11843_A`:`11852_C`) %>%
  filter(rowSums(across(-sample_name)) > 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Convert absolute abundances to relative abundances
jung.dists <- sweep(
  jung.dists, 1, rowSums(jung.dists), "/"
)

# All within-colony pairwise distances
all.dists <- as.matrix(vegdist(jung.dists, method = "bray")) %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(!from, names_to = "to", values_to = "dist") %>%
  left_join(relevant.meta, by = c("from" = "sample_name")) %>%
  left_join(relevant.meta, by = c("to" = "sample_name")) %>%
  filter(COLONY_ID.x == COLONY_ID.y) %>%
  group_by(COLONY_ID.x) %>%
  mutate(
    date.rank.x = dense_rank(DATE.x),
    date.rank.y = dense_rank(DATE.y)
  ) %>%
  ungroup() %>%
  distinct()

# Distance between sequential observations
sequential <- all.dists %>%
  filter(date.rank.y == date.rank.x + 1) %>%
  select(COLONY_ID.x, DATE.x, DATE.y, dist) %>%
  distinct()

# Distance from initial community
initial <- all.dists %>%
  filter(date.rank.x == 1, date.rank.y > 1) %>%
  select(COLONY_ID.x, DATE.y, DIST_INITIAL = dist) %>%
  distinct()

# Combine and format
jungcolonydists <- sequential %>%
  left_join(initial, by = c("COLONY_ID.x", "DATE.y")) %>%
  arrange(COLONY_ID.x, DATE.x) %>%
  rename(
    COLONY_ID = COLONY_ID.x,
    FROM_DATE = DATE.x,
    TO_DATE = DATE.y,
    DIST = dist
  ) %>%
  left_join(colony.meta, by = "COLONY_ID") %>%
  select(
    COLONY_ID, FROM_DATE, TO_DATE,
    DIST, DIST_INITIAL,
    LATITUDE, LONGITUDE, SPECIES
  )

jungcolonydists

#-----------ALESSI ET AL 2024 -------------

alessi.divs <- read.csv("data/Alessi_ITS2_results.csv")

# Prepare metadata
alessi.meta.full <- alessi.divs %>%
  mutate(
    sample_name = paste0("ALESSI_SAMPLE_", row_number()),
    COLONY_ID = paste("ALESSI", "AMUR", ID, sep = "_"),
    SPECIES = "ACROPORA_MURICATA",
    DATE = case_when(
      Time == "T1" ~ as.Date("2021-01-18"),
      Time == "T2" ~ as.Date("2021-05-05")
    ),
    LATITUDE = ifelse(Site == "Bouraké", -21.952, -21.965),
    LONGITUDE = ifelse(Site == "Bouraké", 165.990, 165.945)
  ) %>%
  filter(!is.na(DATE), !is.na(ID)) %>%
  distinct()

relevant.meta <- alessi.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- alessi.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare DIV community composition matrix
alessi.dists <- alessi.divs %>%
  mutate(sample_name = paste0("ALESSI_SAMPLE_", row_number())) %>%
  semi_join(relevant.meta, by = "sample_name") %>%
  select(sample_name, C1:C21bh) %>%
  filter(rowSums(across(-sample_name)) > 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Normalize to ensure each sample sums to one
alessi.dists <- sweep(
  alessi.dists, 1, rowSums(alessi.dists), "/"
)

# All within-colony pairwise distances
all.dists <- as.matrix(vegdist(alessi.dists, method = "bray")) %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(!from, names_to = "to", values_to = "dist") %>%
  left_join(relevant.meta, by = c("from" = "sample_name")) %>%
  left_join(relevant.meta, by = c("to" = "sample_name")) %>%
  filter(COLONY_ID.x == COLONY_ID.y) %>%
  group_by(COLONY_ID.x) %>%
  mutate(
    date.rank.x = dense_rank(DATE.x),
    date.rank.y = dense_rank(DATE.y)
  ) %>%
  ungroup() %>%
  distinct()

# Distance between sequential observations
sequential <- all.dists %>%
  filter(date.rank.y == date.rank.x + 1) %>%
  select(COLONY_ID.x, DATE.x, DATE.y, dist) %>%
  distinct()

# Distance from initial community
initial <- all.dists %>%
  filter(date.rank.x == 1, date.rank.y > 1) %>%
  select(COLONY_ID.x, DATE.y, DIST_INITIAL = dist) %>%
  distinct()

# Combine and format
alessicolonydists <- sequential %>%
  left_join(initial, by = c("COLONY_ID.x", "DATE.y")) %>%
  arrange(COLONY_ID.x, DATE.x) %>%
  rename(
    COLONY_ID = COLONY_ID.x,
    FROM_DATE = DATE.x,
    TO_DATE = DATE.y,
    DIST = dist
  ) %>%
  left_join(colony.meta, by = "COLONY_ID") %>%
  select(
    COLONY_ID, FROM_DATE, TO_DATE,
    DIST, DIST_INITIAL,
    LATITUDE, LONGITUDE, SPECIES
  )

alessicolonydists


#-----------Glass et al., 2026 -------------

glass.divs <- read.csv(
  "data/Glass_2026/PostMED_relative_abundance.csv"
)

# Prepare metadata
glass.meta.full <- glass.divs %>%
  mutate(
    DATE = case_when(
      Year == 2022 ~ as.Date("2022-05-01"),
      Year == 2023 ~ as.Date("2023-08-15"),
      Year == 2024 ~ as.Date("2024-04-28")
    ),
    COLONY_ID = paste("GLASS_SSID", Coral_ID, sep = "_"),
    SPECIES = "SIDERASTREA_SIDEREA",
    LATITUDE = case_when(
      Site == "Cristobal Island" ~ 9.268270,
      Site == "Hospital Point" ~ 9.333190,
      Site == "Punta Donato" ~ 9.35845,
      Site == "STRI Point" ~ 9.34868
    ),
    LONGITUDE = case_when(
      Site == "Cristobal Island" ~ -82.242280,
      Site == "Hospital Point" ~ -82.217650,
      Site == "Punta Donato" ~ -82.36804,
      Site == "STRI Point" ~ -82.2627
    ),
    sample_name = paste0("GLASS_SAMPLE_", row_number())
  ) %>%
  filter(!is.na(DATE), !is.na(Coral_ID))

relevant.meta <- glass.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- glass.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare DIV community composition matrix
glass.dists <- glass.meta.full %>%
  select(sample_name, A3:D2u) %>%
  filter(rowSums(across(-sample_name)) > 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Normalize relative abundances to sum to one
glass.dists <- sweep(
  glass.dists, 1, rowSums(glass.dists), "/"
)

# All within-colony pairwise distances
all.dists <- as.matrix(vegdist(glass.dists, method = "bray")) %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(!from, names_to = "to", values_to = "dist") %>%
  left_join(relevant.meta, by = c("from" = "sample_name")) %>%
  left_join(relevant.meta, by = c("to" = "sample_name")) %>%
  filter(COLONY_ID.x == COLONY_ID.y) %>%
  group_by(COLONY_ID.x) %>%
  mutate(
    date.rank.x = dense_rank(DATE.x),
    date.rank.y = dense_rank(DATE.y)
  ) %>%
  ungroup() %>%
  distinct()

# Distance between sequential observations
sequential <- all.dists %>%
  filter(date.rank.y == date.rank.x + 1) %>%
  select(COLONY_ID.x, DATE.x, DATE.y, dist) %>%
  distinct()

# Distance from initial community to each subsequent observation
initial <- all.dists %>%
  filter(date.rank.x == 1, date.rank.y > 1) %>%
  select(COLONY_ID.x, DATE.y, DIST_INITIAL = dist) %>%
  distinct()

# Combine and format
glasscolonydists <- sequential %>%
  left_join(initial, by = c("COLONY_ID.x", "DATE.y")) %>%
  arrange(COLONY_ID.x, DATE.x) %>%
  rename(
    COLONY_ID = COLONY_ID.x,
    FROM_DATE = DATE.x,
    TO_DATE = DATE.y,
    DIST = dist
  ) %>%
  left_join(colony.meta, by = "COLONY_ID") %>%
  select(
    COLONY_ID, FROM_DATE, TO_DATE,
    DIST, DIST_INITIAL,
    LATITUDE, LONGITUDE, SPECIES
  )

glasscolonydists

##----------- POWELL ET AL 2026 ---------------

powell.divs <- read.csv(
  "data/symportal_profiles/Powell_DIVs.csv",
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

powell.meta <- read.csv(
  "data/symportal_profiles/Powell_sample_info_curacao_2020.csv",
  fileEncoding = "latin1",
  colClasses = c(number = "character")
)

# Prepare metadata
powell.meta.full <- powell.meta %>%
  filter(
    !is.na(host_genus),
    !is.na(host_species),
    !is.na(number),
    !is.na(site)
  ) %>%
  mutate(
    sample_name = as.character(id),
    
    DATE = case_when(
      grepl("^[0-9]{6}$", collection_date) ~
        as.Date(paste0(collection_date, "01"), format = "%Y%m%d"),
      TRUE ~ as.Date(collection_date, format = "%m/%d/%y")
    ),
    
    COLONY_ID = paste0(
      "POWELL_",
      toupper(substr(host_genus, 1, 1)),
      toupper(substr(host_species, 1, 3)),
      "_", site, "_", number
    ),
    
    SPECIES = toupper(paste(host_genus, host_species, sep = "_"))
    
  ) %>%
  filter(!is.na(DATE)) %>%
  distinct()

# Use site-level coordinates to avoid inconsistent
# coordinates in individual metadata records
powell.sites <- powell.meta.full %>%
  group_by(site) %>%
  summarise(
    LATITUDE = as.numeric(
      names(sort(table(collection_latitude), decreasing = TRUE))[1]
    ),
    LONGITUDE = -abs(as.numeric(
      names(sort(table(collection_longitude), decreasing = TRUE))[1]
    )),
    .groups = "drop"
  )

powell.meta.full <- powell.meta.full %>%
  select(-collection_latitude, -collection_longitude) %>%
  left_join(powell.sites, by = "site")

# Prepare DIV matrix and normalize raw counts
powell.dists <- powell.divs %>%
  rename(sample_name = id) %>%
  semi_join(powell.meta.full, by = "sample_name") %>%
  select(sample_name, where(is.numeric)) %>%
  filter(rowSums(across(-sample_name)) > 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

powell.dists <- sweep(
  powell.dists, 1, rowSums(powell.dists), "/"
)

# Average relative abundance profiles within colony-date
# when multiple samples represent the same observation
powell.averaged <- powell.dists %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample_name") %>%
  left_join(
    powell.meta.full %>%
      select(sample_name, COLONY_ID, DATE),
    by = "sample_name"
  ) %>%
  select(-sample_name) %>%
  group_by(COLONY_ID, DATE) %>%
  summarise(
    across(where(is.numeric), mean),
    .groups = "drop"
  ) %>%
  mutate(
    sample_name = paste(
      COLONY_ID, format(DATE, "%Y%m%d"), sep = "__"
    )
  )

relevant.meta <- powell.averaged %>%
  select(sample_name, COLONY_ID, DATE)

colony.meta <- powell.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Community matrix: one row per colony and date
powell.dists <- powell.averaged %>%
  select(-COLONY_ID, -DATE) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Re-normalize averaged profiles
powell.dists <- sweep(
  powell.dists, 1, rowSums(powell.dists), "/"
)

# All pairwise Bray-Curtis distances
all.dists <- as.matrix(
  vegdist(powell.dists, method = "bray")
) %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(!from, names_to = "to", values_to = "dist") %>%
  left_join(relevant.meta, by = c("from" = "sample_name")) %>%
  left_join(relevant.meta, by = c("to" = "sample_name")) %>%
  filter(COLONY_ID.x == COLONY_ID.y) %>%
  group_by(COLONY_ID.x) %>%
  mutate(
    date.rank.x = dense_rank(DATE.x),
    date.rank.y = dense_rank(DATE.y)
  ) %>%
  ungroup() %>%
  distinct()

# Distance between sequential observations
sequential <- all.dists %>%
  filter(date.rank.y == date.rank.x + 1) %>%
  select(COLONY_ID.x, DATE.x, DATE.y, dist) %>%
  distinct()

# Distance from initial community
initial <- all.dists %>%
  filter(date.rank.x == 1, date.rank.y > 1) %>%
  select(COLONY_ID.x, DATE.y, DIST_INITIAL = dist) %>%
  distinct()

# Combine and format
powellcolonydists <- sequential %>%
  left_join(initial, by = c("COLONY_ID.x", "DATE.y")) %>%
  arrange(COLONY_ID.x, DATE.x) %>%
  rename(
    COLONY_ID = COLONY_ID.x,
    FROM_DATE = DATE.x,
    TO_DATE = DATE.y,
    DIST = dist
  ) %>%
  left_join(colony.meta, by = "COLONY_ID") %>%
  select(
    COLONY_ID, FROM_DATE, TO_DATE,
    DIST, DIST_INITIAL,
    LATITUDE, LONGITUDE, SPECIES
  )

powellcolonydists


## --------- HUFFMYER ET AL 2026 --------------


huff.divs <- read.table(
  "data/Huffmyer-Symportal.txt",
  sep = "\t", header = TRUE, check.names = FALSE
)

# Prepare metadata
huff.meta.full <- huff.divs %>%
  select(
    sample_name, collection_longitude, collection_latitude,
    collection_date, collection_depth, host_genus, host_species
  ) %>%
  filter(
    !is.na(collection_longitude),
    collection_longitude != 999
  ) %>%
  mutate(
    DATE = as.Date(
      as.character(collection_date),
      format = "%Y%m%d"
    )
  ) %>%
  separate(
    sample_name,
    into = c("spec", "ID", "tp"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    COLONY_ID = paste("HUFF", spec, ID, sep = "_"),
    SPECIES = toupper(paste(host_genus, host_species, sep = "_")),
    LATITUDE = as.numeric(collection_latitude),
    LONGITUDE = as.numeric(collection_longitude)
  ) %>%
  filter(
    !is.na(DATE),
    !is.na(COLONY_ID),
    !is.na(LATITUDE),
    !is.na(LONGITUDE),
    LATITUDE != 999,
    LONGITUDE != 999
  ) %>%
  distinct()

relevant.meta <- huff.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- huff.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare DIV community composition matrix
huff.dists <- huff.divs %>%
  semi_join(relevant.meta, by = "sample_name") %>%
  select(sample_name, A1ee:G3ab) %>%
  filter(rowSums(across(-sample_name)) > 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Normalize abundances
huff.dists <- sweep(
  huff.dists, 1, rowSums(huff.dists), "/"
)

# All within-colony pairwise distances
all.dists <- as.matrix(
  vegdist(huff.dists, method = "bray")
) %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(!from, names_to = "to", values_to = "dist") %>%
  left_join(relevant.meta, by = c("from" = "sample_name")) %>%
  left_join(relevant.meta, by = c("to" = "sample_name")) %>%
  filter(COLONY_ID.x == COLONY_ID.y) %>%
  group_by(COLONY_ID.x) %>%
  mutate(
    date.rank.x = dense_rank(DATE.x),
    date.rank.y = dense_rank(DATE.y)
  ) %>%
  ungroup() %>%
  distinct()

# Distance between consecutive observations
sequential <- all.dists %>%
  filter(date.rank.y == date.rank.x + 1) %>%
  select(COLONY_ID.x, DATE.x, DATE.y, dist) %>%
  distinct()

# Distance from initial community
initial <- all.dists %>%
  filter(date.rank.x == 1, date.rank.y > 1) %>%
  select(COLONY_ID.x, DATE.y, DIST_INITIAL = dist) %>%
  distinct()

# Combine and format
huffcolonydists <- sequential %>%
  left_join(initial, by = c("COLONY_ID.x", "DATE.y")) %>%
  arrange(COLONY_ID.x, DATE.x) %>%
  rename(
    COLONY_ID = COLONY_ID.x,
    FROM_DATE = DATE.x,
    TO_DATE = DATE.y,
    DIST = dist
  ) %>%
  left_join(colony.meta, by = "COLONY_ID") %>%
  select(
    COLONY_ID, FROM_DATE, TO_DATE,
    DIST, DIST_INITIAL,
    LATITUDE, LONGITUDE, SPECIES
  )

huffcolonydists

# -------- CUNNING ET AL 2024 -------


cun.divs <- read.table(
  "data/Cunning_Symportal_111_20200624_2020-06-24_03-40-50.238182.seqs.absolute.abund_and_meta.txt",
  sep = "\t", header = TRUE, check.names = FALSE
)

# Prepare metadata
cun.meta.full <- cun.divs %>%
  select(
    sample_name, collection_latitude, collection_longitude,
    collection_depth, collection_date,
    host_genus, host_species
  ) %>%
  separate(
    sample_name,
    into = c("SITE", "SPEC", "ID", "YEAR"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    COLONY_ID = paste("CUN", SITE, SPEC, ID, sep = "_"),
    DATE = as.Date(
      as.character(collection_date),
      format = "%Y.%m.%d"
    ),
    LATITUDE = as.numeric(collection_latitude),
    LONGITUDE = -abs(as.numeric(collection_longitude)),
    SPECIES = toupper(paste(host_genus, host_species, sep = "_"))
  ) %>%
  filter(
    !is.na(DATE),
    !is.na(COLONY_ID),
    !is.na(LATITUDE),
    !is.na(LONGITUDE)
  ) %>%
  distinct()

cun.meta.full <- cun.meta.full %>%
  mutate(
    host_genus = if_else(
      sample_name == "ET_Ofra_59_2017",
      "Orbicella", host_genus
    ),
    host_species = if_else(
      sample_name == "ET_Ofra_59_2017",
      "franksi", host_species
    ),
    SPECIES = toupper(paste(host_genus, host_species, sep = "_"))
  )

cun.meta.full <- cun.meta.full %>%
  group_by(COLONY_ID, DATE) %>%
  filter(n() == 1) %>%
  ungroup()

relevant.meta <- cun.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- cun.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare DIV community composition matrix
cun.dists <- cun.divs %>%
  semi_join(relevant.meta, by = "sample_name") %>%
  select(sample_name, A4:G3b) %>%
  filter(rowSums(across(-sample_name)) > 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Convert absolute abundances to relative abundances
cun.dists <- sweep(
  cun.dists, 1, rowSums(cun.dists), "/"
)

# All within-colony pairwise Bray-Curtis distances
all.dists <- as.matrix(
  vegdist(cun.dists, method = "bray")
) %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(!from, names_to = "to", values_to = "dist") %>%
  left_join(relevant.meta, by = c("from" = "sample_name")) %>%
  left_join(relevant.meta, by = c("to" = "sample_name")) %>%
  filter(COLONY_ID.x == COLONY_ID.y) %>%
  group_by(COLONY_ID.x) %>%
  mutate(
    date.rank.x = dense_rank(DATE.x),
    date.rank.y = dense_rank(DATE.y)
  ) %>%
  ungroup() %>%
  distinct()

# Distance between consecutive observations
sequential <- all.dists %>%
  filter(date.rank.y == date.rank.x + 1) %>%
  select(COLONY_ID.x, DATE.x, DATE.y, dist) %>%
  distinct()

# Distance from initial community
initial <- all.dists %>%
  filter(date.rank.x == 1, date.rank.y > 1) %>%
  select(
    COLONY_ID.x,
    DATE.y,
    DIST_INITIAL = dist
  ) %>%
  distinct()

# Combine and format
cuncolonydists <- sequential %>%
  left_join(initial, by = c("COLONY_ID.x", "DATE.y")) %>%
  arrange(COLONY_ID.x, DATE.x) %>%
  rename(
    COLONY_ID = COLONY_ID.x,
    FROM_DATE = DATE.x,
    TO_DATE = DATE.y,
    DIST = dist
  ) %>%
  left_join(colony.meta, by = "COLONY_ID") %>%
  select(
    COLONY_ID, FROM_DATE, TO_DATE,
    DIST, DIST_INITIAL,
    LATITUDE, LONGITUDE, SPECIES
  )

cuncolonydists

##--------------- de Souza et al., 2023 -------------

desouza.divs <- readxl::read_excel(
  "data/deSouza_Symportal.xlsx"
)

# Prepare metadata
desouza.meta.full <- desouza.divs %>%
  select(sample_name, host_genus, host_species) %>%
  separate(
    sample_name,
    into = c("Date", "ColonyID"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    DATE = case_when(
      Date == "Clo2017" ~ as.Date("2018-01-01"),
      Date == "Clo2019" ~ as.Date("2019-10-01")
    ),
    COLONY_ID = paste0("deSouz_MCAP_", ColonyID),
    SPECIES = toupper(paste(host_genus, host_species, sep = "_")),
    LATITUDE = 21.47098,
    LONGITUDE = -157.82333
  ) %>%
  filter(!is.na(DATE), !is.na(ColonyID)) %>%
  distinct()

relevant.meta <- desouza.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- desouza.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare DIV community matrix
desouza.dists <- desouza.divs %>%
  semi_join(relevant.meta, by = "sample_name") %>%
  select(sample_name, B1:F3.1) %>%
  filter(rowSums(across(-sample_name)) > 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Convert abundances to relative abundances
desouza.dists <- sweep(
  desouza.dists, 1, rowSums(desouza.dists), "/"
)

# All within-colony pairwise distances
all.dists <- as.matrix(
  vegdist(desouza.dists, method = "bray")
) %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(!from, names_to = "to", values_to = "dist") %>%
  left_join(relevant.meta, by = c("from" = "sample_name")) %>%
  left_join(relevant.meta, by = c("to" = "sample_name")) %>%
  filter(COLONY_ID.x == COLONY_ID.y) %>%
  group_by(COLONY_ID.x) %>%
  mutate(
    date.rank.x = dense_rank(DATE.x),
    date.rank.y = dense_rank(DATE.y)
  ) %>%
  ungroup() %>%
  distinct()

# Distance between sequential observations
sequential <- all.dists %>%
  filter(date.rank.y == date.rank.x + 1) %>%
  select(COLONY_ID.x, DATE.x, DATE.y, dist) %>%
  distinct()

# Distance from initial community
initial <- all.dists %>%
  filter(date.rank.x == 1, date.rank.y > 1) %>%
  select(
    COLONY_ID.x,
    DATE.y,
    DIST_INITIAL = dist
  ) %>%
  distinct()

# Combine and format
desouzacolonydists <- sequential %>%
  left_join(
    initial,
    by = c("COLONY_ID.x", "DATE.y"),
    relationship = "one-to-one"
  ) %>%
  arrange(COLONY_ID.x, DATE.x) %>%
  rename(
    COLONY_ID = COLONY_ID.x,
    FROM_DATE = DATE.x,
    TO_DATE = DATE.y,
    DIST = dist
  ) %>%
  left_join(
    colony.meta,
    by = "COLONY_ID",
    relationship = "many-to-one"
  ) %>%
  select(
    COLONY_ID, FROM_DATE, TO_DATE,
    DIST, DIST_INITIAL,
    LATITUDE, LONGITUDE, SPECIES
  )

desouzacolonydists

## --------- Starko et al., 2023 --------
library(phyloseq)
load("data/Starko_Porites_DIV.RData")

# Extract metadata
stark.meta.full <- data.frame(sample_data(Porites_DIV)) %>% dplyr::select(!sample_name) %>%
  rownames_to_column(var = "sample_name") %>%
  mutate(
    DATE = as.Date(
      as.character(collection_date),
      format = "%Y%m%d"
    ),
    COLONY_ID = paste("STARK_PLOB", coral_tag, sep = "_"),
    SPECIES = toupper(paste(host_genus, host_species, sep = "_")),
    LATITUDE = as.numeric(collection_latitude),
    LONGITUDE = as.numeric(collection_longitude)
  ) %>%
  filter(
    !is.na(DATE),
    !is.na(coral_tag),
    !is.na(LATITUDE),
    !is.na(LONGITUDE),
    LATITUDE != 999,
    LONGITUDE != 999
  ) %>%
  distinct() %>%
  filter(
    !is.na(coral_tag),
    !trimws(as.character(coral_tag)) %in% c("", "NA", "NaN", "999")
  )

stark.meta.full <- stark.meta.full %>%
  filter(!COLONY_ID %in% c("STARK_PLOB_NA", "STARK_PLOB_", "STARK_PLOB_999"))

relevant.meta <- stark.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- stark.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Remove samples with missing colony identifiers
stark.meta.full <- stark.meta.full %>%
  filter(
    !is.na(coral_tag),
    !trimws(as.character(coral_tag)) %in% c("", "NA", "NaN", "999")
  )

# Extract full DIV abundance matrix
stark.dists <- as(otu_table(Porites_DIV), "matrix")

# Ensure samples are rows
if (taxa_are_rows(Porites_DIV)) {
  stark.dists <- t(stark.dists)
}

# Retain samples with valid metadata
stark.dists <- stark.dists[
  rownames(stark.dists) %in% stark.meta.full$sample_name,
  ,
  drop = FALSE
]

# Remove samples with zero total abundance
stark.dists <- stark.dists[
  rowSums(stark.dists) > 0,
  ,
  drop = FALSE
]

# Normalize each sample to relative abundance
stark.dists <- sweep(
  stark.dists, 1, rowSums(stark.dists), "/"
)

# Average relative DIV abundances for replicate samples
stark.averaged <- stark.dists %>%
  as.data.frame() %>%
  rownames_to_column("sample_name") %>%
  left_join(
    stark.meta.full %>%
      select(sample_name, COLONY_ID, DATE),
    by = "sample_name",
    relationship = "many-to-one"
  ) %>%
  select(-sample_name) %>%
  group_by(COLONY_ID, DATE) %>%
  summarise(
    across(where(is.numeric), mean),
    .groups = "drop"
  ) %>%
  mutate(
    sample_name = paste(
      COLONY_ID,
      format(DATE, "%Y%m%d"),
      sep = "__"
    )
  )

# Update sample metadata
relevant.meta <- stark.averaged %>%
  select(sample_name, COLONY_ID, DATE)

# Rebuild community matrix
stark.dists <- stark.averaged %>%
  select(-COLONY_ID, -DATE) %>%
  as.data.frame() %>%
  column_to_rownames("sample_name")


# All within-colony pairwise Bray-Curtis distances
all.dists <- as.matrix(
  vegan::vegdist(stark.dists, method = "bray")
) %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(!from, names_to = "to", values_to = "dist") %>%
  left_join(
    relevant.meta,
    by = c("from" = "sample_name")
  ) %>%
  left_join(
    relevant.meta,
    by = c("to" = "sample_name")
  ) %>%
  filter(COLONY_ID.x == COLONY_ID.y) %>%
  group_by(COLONY_ID.x) %>%
  mutate(
    date.rank.x = dense_rank(DATE.x),
    date.rank.y = dense_rank(DATE.y)
  ) %>%
  ungroup() %>%
  distinct()

# Distance between consecutive observations
sequential <- all.dists %>%
  filter(date.rank.y == date.rank.x + 1) %>%
  select(COLONY_ID.x, DATE.x, DATE.y, dist) %>%
  distinct()

# Distance from initial community
initial <- all.dists %>%
  filter(date.rank.x == 1, date.rank.y > 1) %>%
  select(
    COLONY_ID.x,
    DATE.y,
    DIST_INITIAL = dist
  ) %>%
  distinct()

# Combine and format
starkcolonydists <- sequential %>%
  left_join(
    initial,
    by = c("COLONY_ID.x", "DATE.y"),
    relationship = "one-to-one"
  ) %>%
  arrange(COLONY_ID.x, DATE.x) %>%
  rename(
    COLONY_ID = COLONY_ID.x,
    FROM_DATE = DATE.x,
    TO_DATE = DATE.y,
    DIST = dist
  ) %>%
  left_join(
    colony.meta,
    by = "COLONY_ID",
    relationship = "many-to-one"
  ) %>%
  select(
    COLONY_ID, FROM_DATE, TO_DATE,
    DIST, DIST_INITIAL,
    LATITUDE, LONGITUDE, SPECIES
  )

starkcolonydists

##-------- Nunez Ponz 2023 -----------------
#symportal did not run like it should and almost no samples passed screening
#However, authors report results, so what's up? The symportal entry also doesn't match the region or timeframe it's supposed to
#I think I know why I used DADA2 now; it didn't assemble contigs with symportal


library(dplyr)
library(tidyr)
library(tibble)
library(vegan)


load("data/get.taxa.nunez.Rdata")
load("data/Nunez_OTUtable.ITS2.Rdata")

otus <- as.matrix(seqtab.nochim)

# Retain samples with >= 5000 reads
otus <- otus[rowSums(otus) >= 5000, , drop = FALSE]

# Retain ASVs with >= 10 reads across retained samples
otus <- otus[, colSums(otus) >= 10, drop = FALSE]


tax.assign <- data.frame(get.taxa) 
tax.assign %>%
  rownames_to_column("ASV") -> tax.assign

# Keep ASVs with genus-level assignments, as in
# your original pipeline, but do not collapse them
# into genus/species abundance categories.

valid.asvs <- tax.assign %>%
  filter(
    !is.na(Genus),
    !Genus %in% c("", "g__", "NA")
  ) %>%
  pull(ASV)

otus <- otus[
  ,
  colnames(otus) %in% valid.asvs,
  drop = FALSE
]

# Check taxonomy-to-ASV matching
stopifnot(ncol(otus) > 0)


srr_2_samn <- read.csv(
  "data/PRJNA794042_meta.csv",
  stringsAsFactors = FALSE
) %>%
  select(Run, BioSample, SampleName) %>%
  distinct()

nunez.meta <- read.table(
  "data/NunezMeta.tsv",
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE
)

otu.meta <- tibble(ID = rownames(otus)) %>%
  left_join(
    srr_2_samn,
    by = c("ID" = "Run"),
    relationship = "one-to-one"
  ) %>%
  left_join(
    nunez.meta,
    by = c("BioSample" = "V1"),
    relationship = "one-to-one"
  ) %>%
  filter(
    !is.na(V7),
    V7 != "BLANK"
  ) %>%
  mutate(
    SPECIES = gsub(" ", "_", toupper(V7)),
    DATE = as.Date(V11, format = "%Y-%m-%d"),
    COLONY_ID = paste(
      "NUNEZ",
      substr(SPECIES, 1, 3),
      gsub("_t[0-9]+$", "", SampleName),
      sep = "_"
    ),
    LATITUDE = 21.461,
    LONGITUDE = -157.823
  ) %>%
  filter(
    !is.na(DATE),
    !is.na(SampleName),
    !is.na(SPECIES)
  ) %>%
  select(
    ID, COLONY_ID, DATE,
    SPECIES, LATITUDE, LONGITUDE
  )

# Check metadata consistency within colonies
otu.meta %>%
  group_by(COLONY_ID) %>%
  summarise(
    n_species = n_distinct(SPECIES),
    n_lat = n_distinct(LATITUDE),
    n_lon = n_distinct(LONGITUDE),
    .groups = "drop"
  ) %>%
  filter(n_species > 1 | n_lat > 1 | n_lon > 1)

# Check for multiple samples per colony-date
otu.meta %>%
  count(COLONY_ID, DATE) %>%
  filter(n > 1)


# Keep samples with valid metadata
otus <- otus[
  rownames(otus) %in% otu.meta$ID,
  ,
  drop = FALSE
]

# Exclude samples with no remaining Symbiodiniaceae reads
otus <- otus[
  rowSums(otus) > 0,
  ,
  drop = FALSE
]

# Normalize to relative abundance
otus.rel <- sweep(otus, 1, rowSums(otus), "/")


# Each biological sample contributes equally,
# irrespective of sequencing depth.

otu.averaged <- as.data.frame(
  otus.rel,
  check.names = FALSE
) %>%
  rownames_to_column("ID") %>%
  left_join(
    otu.meta %>%
      select(ID, COLONY_ID, DATE),
    by = "ID",
    relationship = "one-to-one"
  ) %>%
  group_by(COLONY_ID, DATE) %>%
  summarise(
    across(all_of(colnames(otus.rel)), mean),
    .groups = "drop"
  ) %>%
  mutate(
    sample_name = paste(
      COLONY_ID,
      format(DATE, "%Y%m%d"),
      sep = "__"
    )
  )

# New sample metadata
relevant.meta <- otu.averaged %>%
  select(sample_name, COLONY_ID, DATE)

# Rebuild ASV matrix
otus.rel <- otu.averaged %>%
  select(sample_name, all_of(colnames(otus.rel))) %>%
  as.data.frame(check.names = FALSE) %>%
  column_to_rownames("sample_name") %>%
  as.matrix()

# Verify alignment and uniqueness
stopifnot(
  !anyDuplicated(relevant.meta$sample_name),
  !anyDuplicated(relevant.meta[c("COLONY_ID", "DATE")]),
  all(rownames(otus.rel) %in% relevant.meta$sample_name)
)


# Calculate distances within each colony rather
# than creating a full study-wide distance matrix.

colony.groups <- relevant.meta %>%
  arrange(COLONY_ID, DATE) %>%
  group_split(COLONY_ID)

nunez.dists <- lapply(colony.groups, function(m) {
  
  # Only longitudinally sampled colonies
  if (nrow(m) < 2) return(NULL)
  
  x <- otus.rel[
    m$sample_name,
    ,
    drop = FALSE
  ]
  
  d <- as.matrix(
    vegan::vegdist(x, method = "bray")
  )
  
  tibble(
    COLONY_ID = m$COLONY_ID[-1],
    FROM_DATE = m$DATE[-nrow(m)],
    TO_DATE = m$DATE[-1],
    DIST = d[cbind(
      seq_len(nrow(m) - 1),
      seq_len(nrow(m) - 1) + 1
    )],
    DIST_INITIAL = d[cbind(
      rep(1, nrow(m) - 1),
      seq.int(2, nrow(m))
    )]
  )
  
}) %>%
  bind_rows()


colony.meta <- otu.meta %>%
  select(
    COLONY_ID, LATITUDE,
    LONGITUDE, SPECIES
  ) %>%
  distinct()

# Ensure one metadata record per colony
stopifnot(!anyDuplicated(colony.meta$COLONY_ID))

nunezcolonydists <- nunez.dists %>%
  left_join(
    colony.meta,
    by = "COLONY_ID",
    relationship = "many-to-one"
  ) %>%
  arrange(COLONY_ID, FROM_DATE) %>%
  select(
    COLONY_ID, FROM_DATE, TO_DATE,
    DIST, DIST_INITIAL,
    LATITUDE, LONGITUDE, SPECIES
  )


nunezcolonydists %>%
  group_by(SPECIES) %>%
  summarise(
    n_colonies = n_distinct(COLONY_ID),
    n_comparisons = n(),
    mean_dist = mean(DIST),
    .groups = "drop"
  )


#### ----------- COMBINE DATASETS --------------
all.colonydists <- list(
  LEE = lee.colonydists,
  TERRANEO = terrcolonydists,
  CLAAR = claarcolonydists,
  BUZZONI = buzzcolonydists,
  VAN_NYNATTEN = vncolonydists,
  LEINBACH = lbachcolonydists,
  HOWELLS = howcolonydists,
  JUNG = jungcolonydists,
  ALESSI = alessicolonydists,
  GLASS = glasscolonydists,
  POWELL = powellcolonydists,
  HUFFMYER = huffcolonydists,
  CUNNING = cuncolonydists,
  DE_SOUZA = desouzacolonydists,
  STARKO = starkcolonydists,
  NUNEZ_PONS = nunezcolonydists
) %>%
  purrr::imap_dfr(
    ~ .x %>%
      rename_with(~ rep("COLONY_ID", length(.x)), any_of("COLONYID"))%>%
      mutate(STUDY = .y)
  ) %>%
  select(
    STUDY, COLONY_ID, SPECIES,
    FROM_DATE, TO_DATE,
    DIST, DIST_INITIAL,
    LATITUDE, LONGITUDE
  ) %>%
  arrange(STUDY, COLONY_ID, FROM_DATE)

#Now pull environmental data for the species-level. 
save(all.colonydists, file = "data/ColonyDistances.Rdata")
