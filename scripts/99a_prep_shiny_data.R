#########################################################################
# PREP DATA FOR SHINY APP (delta_PC / STUDY-SITE-COLONY version)
#########################################################################

#==============================================================================
#Development Note: I used Claude's Sonnet 5 model to develop the Shiny app
#described in this git repository. I manually generated and revised all scripts
#leading up to this point, and this script draws directly from those curated outputs.
#That said, this Shiny app data should be considered exploratory and an extension
#of the results of the main analysis; the use of AI was intended to increase 
#interpretability of results and to bridge my knowledge gap in building a Shiny app.
#===============================================================================

# Delta-model analog of 99a_prep_shiny_data.R. Gathers everything the
# Shiny app needs for the DELTA-PC model -- fitted model draws, phylo/
# trait trees, environmental PCA rotation, and species/site summary
# tables -- into data/shiny_app_data_delta.RDS.
#
# Run this AFTER fitting the delta model, i.e. after
# data/stan_dirichlet_explicit_delta.RDS exists (fit_dirichlet_partition_within
# from your delta-model script).
#
# CHANGES FROM 99a_prep_shiny_data.R (flagged inline with "DELTA CHANGE:"):
#   1. delta_PC1/delta_PC2/delta_PC3 are computed per-row (consecutive-
#      observation differences), the same way delta_t already was.
#   2. The model formula/random-effects structure now matches
#      coral_dirichlet_partition_site.stan: (1 | STUDY/SITE/COLONY_ID)
#      (3 nested levels) instead of (1 | STUDY/COLONY_ID) (2 levels), and
#      delta_PC1/2/3 replace PC1/2/3 as both fixed effects and random slopes.
#   3. spec_summary is now built from mean delta_PC1/2/3 per species
#      (still named PC1_mean/PC2_mean/PC3_mean for compatibility with
#      compute_species_switching() in helpers.R, which hardcodes those names).
#   4. species_pc_range is now built from observed delta_PC1/2/3 ranges,
#      since "outside observed range" should now mean "outside the range
#      of environmental CHANGE seen," not raw environmental state.
#
# ASSUMPTION FLAGGED BELOW: this script assumes a `SITE` column exists in
# genlevel/dat for the new STUDY:SITE nesting level. If your actual column
# is named differently (e.g. Site_ID, REEF), update the two lines marked
# "# ASSUMPTION: SITE column" accordingly.
#
# Everything NOT flagged as a DELTA CHANGE below is copied unchanged from
# 99a_prep_shiny_data.R -- tree construction, trait PCA, genus handling,
# and site/env lookups for display purposes don't depend on raw-vs-delta.
#########################################################################

library(tidyverse)
library(ape)
library(phytools)
library(vegan)
library(posterior)
library(brms)

source("scripts/helpers.R")

#### ----------- 1. ENVIRONMENTAL PCA (rotation + standardization refs) ----------- ####
# UNCHANGED from 99a -- env.pc itself doesn't depend on raw-vs-delta.

if (file.exists("data/EnvPCA_model.Rdata")) {

  load("data/EnvPCA_model.Rdata")  # expects: env.pc, env_raw_means, env_raw_sds

} else {

  message("data/EnvPCA_model.Rdata not found -- rebuilding env PCA from data/envdata/*.txt ",
          "(fast, since these are cached; see note at bottom of this script to avoid this step next time)")

  sst.files <- list.files("data/envdata", full.names = TRUE, pattern = "sst")
  sst.data  <- do.call(rbind, lapply(sst.files, read.table, header = FALSE, nrows = 1))
  colnames(sst.data) <- c("lat", "lon", "date",
                          "min_SST", "max_SST", "sd_SST", "mean_SST",
                          "max_DHW", "min_DHW", "mean_DHW",
                          "mean_ANOM", "sd_ANOM", "min_ANOM", "max_ANOM", "cum_ANOM")

  precip.files <- list.files("data/envdata/", full.names = TRUE, pattern = "precip")
  precip.data  <- do.call(rbind, lapply(precip.files, read.table, header = FALSE, nrows = 1))
  colnames(precip.data) <- c("lat", "lon", "date",
                             "min_precip", "max_precip", "sd_precip", "mean_precip")

  color.files <- list.files("data/envdata/", full.names = TRUE, pattern = "color")
  color.data  <- do.call(rbind, lapply(color.files, read.table, header = FALSE, nrows = 1))
  colnames(color.data) <- c("lat", "lon", "date",
                            "min_chla", "max_chla", "sd_chla", "mean_chla",
                            "min_turb", "max_turb", "sd_turb", "mean_turb")

  join.data <- left_join(left_join(sst.data, precip.data), color.data)

  means.only <- join.data %>%
    select(mean_ANOM, min_ANOM, max_ANOM, cum_ANOM, mean_DHW, sd_ANOM,
           mean_SST, min_SST, max_SST, sd_SST,
           mean_precip, min_precip, max_precip, sd_precip,
           mean_chla, min_chla, max_chla, sd_chla,
           mean_turb, min_turb, max_turb, sd_turb) %>%
    mutate(across(everything(), ~ if_else(is.infinite(.), NA, .))) %>%
    mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

  env_raw_means <- sapply(means.only, mean, na.rm = TRUE)
  env_raw_sds   <- sapply(means.only, sd,   na.rm = TRUE)

  mean.stand <- decostand(means.only, method = "standardize")
  env.pc     <- prcomp(mean.stand, center = TRUE, scale = TRUE)

  save(env.pc, env_raw_means, env_raw_sds, join.data, file = "data/EnvPCA_model.Rdata")
}

means_only_for_bounds <- join.data %>%
  select(all_of(names(env_raw_means))) %>%
  mutate(across(everything(), ~ if_else(is.infinite(.), NA, .))) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

env_obs_min <- sapply(means_only_for_bounds, min, na.rm = TRUE)
env_obs_max <- sapply(means_only_for_bounds, max, na.rm = TRUE)
env_range_pad <- 0.2 * (env_obs_max - env_obs_min)
env_input_min <- env_obs_min - env_range_pad
env_input_max <- env_obs_max + env_range_pad
non_negative_domain <- env_obs_min >= 0
env_input_min[non_negative_domain] <- pmax(0, env_input_min[non_negative_domain])
#### ----------- 1b. REDUCE TO A NON-REDUNDANT SET OF INPUT VARIABLES ----------- ####
# UNCHANGED from 99a.

CORR_THRESHOLD <- 0.7
n_pcs_for_ranking <- 3

loadings <- env.pc$rotation[, seq_len(n_pcs_for_ranking), drop = FALSE]
var_importance <- apply(abs(loadings), 1, max)
var_importance <- sort(var_importance, decreasing = TRUE)

cor_mat <- cor(means_only_for_bounds, use = "pairwise.complete.obs")

selected_env_vars <- character(0)
redundant_map <- list()

for (v in names(var_importance)) {
  if (length(selected_env_vars) == 0) {
    selected_env_vars <- v
    next
  }
  cors_to_selected <- abs(cor_mat[v, selected_env_vars])
  names(cors_to_selected) <- selected_env_vars
  if (max(cors_to_selected) > CORR_THRESHOLD) {
    proxy <- names(which.max(cors_to_selected))
    fit <- lm(means_only_for_bounds[[v]] ~ means_only_for_bounds[[proxy]])
    redundant_map[[v]] <- list(
      proxy       = proxy,
      intercept   = unname(coef(fit)[1]),
      slope       = unname(coef(fit)[2]),
      r_squared   = summary(fit)$r.squared,
      correlation = cor_mat[v, proxy]
    )
  } else {
    selected_env_vars <- c(selected_env_vars, v)
  }
}

redundant_env_vars <- tibble(
  variable    = names(redundant_map),
  proxy       = map_chr(redundant_map, "proxy"),
  intercept   = map_dbl(redundant_map, "intercept"),
  slope       = map_dbl(redundant_map, "slope"),
  r_squared   = map_dbl(redundant_map, "r_squared"),
  correlation = map_dbl(redundant_map, "correlation")
)

weakest_fit_msg <- if (nrow(redundant_env_vars) > 0) {
  redundant_env_vars %>% slice_min(r_squared, n = 1) %>%
    {sprintf("%s ~ %s (R\u00b2 = %.2f)", .$variable, .$proxy, .$r_squared)}
} else {
  "none (no redundant variables found -- all kept as inputs)"
}

message(length(selected_env_vars), " of ", length(var_importance), " environmental variables kept as ",
        "user-facing inputs (max pairwise |r| among them: ",
        round(max(abs(cor_mat[selected_env_vars, selected_env_vars])[
          lower.tri(cor_mat[selected_env_vars, selected_env_vars])]), 2), "); the rest are auto-set ",
        "from their most-correlated selected variable. Weakest of these auto-set relationships ",
        "(lowest R\u00b2 -- worth a look if it's surprisingly low): ", weakest_fit_msg)

#### ----------- 2. RECONSTRUCT MODEL INPUTS (mirrors delta-model fitting script) ----------- ####

load("data/LongZooxInputData_wEnv.Rdata")  # genlevel, all_datapoints

dat <- genlevel |>
  filter(GENUS != "NA") %>%
  mutate(
    Bleached = ifelse(is.na(BLEACHING_OBS) | BLEACHING_OBS == "NA", "N", BLEACHING_OBS),
    Bleached = factor(Bleached, levels = c("N", "Y"))
  ) %>%
  mutate(SITE = paste(LONGITUDE, LATITUDE, sep ="_")) %>%
  distinct(STUDY, HOST_SPEC_REVISED, SITE, SPECIES, COLONY_ID, DATE, .keep_all = TRUE) |>
  arrange(STUDY, HOST_SPEC_REVISED, SPECIES, SITE,COLONY_ID, DATE) |>
  group_by(STUDY, HOST_SPEC_REVISED, SPECIES, SITE, COLONY_ID) |>
  mutate(
    shifted     = case_when(
      GenShift %in% c("MAJOR", "MINOR", "SUBTYPE_SWITCH") ~ 1,
      GenShift %in% c("NONE")                              ~ 0,
      is.na(GenShift)                                      ~ NA_real_
    ),
    delta_t     = as.numeric(DATE - lag(DATE)),
    log_delta_t = log(delta_t),
    last_genus  = factor(lag(GENUS)),
    # DELTA CHANGE: consecutive-observation differences in the environmental
    # PCs, computed with the SAME group_by/arrange/lag() as delta_t above --
    # this must stay in the same mutate() block so the row alignment (which
    # "previous observation" each delta is measured against) is identical
    # for delta_t and delta_PC1/2/3.
    delta_PC1   = PC1 - lag(PC1),
    delta_PC2   = PC2 - lag(PC2),
    delta_PC3   = PC3 - lag(PC3)
  ) |>
  filter(!is.na(delta_t), !is.na(shifted)) %>%
  filter(is.finite(max_ANOM), is.finite(max_DHW), is.finite(max_SST)) %>%
  separate(SPECIES, into = c("HOST_GEN", "HOST_SPEC"), sep = "_", remove = FALSE) %>%
  ungroup()

dat$last_genus <- relevel(dat$last_genus, ref = "C")

dat$HOST_SPEC <- ifelse(dat$HOST_SPEC == "INTERMEDIA(NOBILIS)", "INTERMEDIA", dat$HOST_SPEC)
dat$HOST_GEN  <- case_when(
  dat$HOST_SPEC == "HYSTRIX"            ~ "SERIATOPORA",
  dat$HOST_GEN  == "GONIASTREA (FAVIA)" ~ "GONIASTREA",
  dat$HOST_SPEC == "PALIFERA"           ~ "ISOPORA",
  .default = dat$HOST_GEN
)

dat <- dat %>%
  mutate(clean.spec = paste(str_to_title(HOST_GEN), tolower(HOST_SPEC), sep = "_")) %>%
  mutate(clean.spec = case_when(
    clean.spec == "Platygyra_daedalea" ~ "Platygyra_daedalea",
    clean.spec == "Porites_sp."        ~ "Porites_lobata",
    clean.spec == "Pachyseris_rugosa"  ~ "Pachyseris_speciosa",
    .default = clean.spec
  )) %>%
  filter(clean.spec != "Acropora_sp.")

# DELTA CHANGE: filter out rows where the delta PCs couldn't be computed
# (first observation per colony has no lag() to difference against --
# same reason delta_t/shifted rows get dropped above; this just makes the
# equivalent drop explicit for the delta columns too, in case any slipped
# through with an NA delta_t edge case)
dat <- dat %>% filter(!is.na(delta_PC1), !is.na(delta_PC2), !is.na(delta_PC3))

# ---- Phylogenetic tree (same manual grafts as 99a) ----
trees <- read.nexus("data/mtree.trees")
tree  <- maxCladeCred(trees)

poc_dam_node <- getMRCA(tree, c("Pocillopora_damicornis", "Pocillopora_effusa"))
bl <- tree$edge.length[which(tree$edge[, 2] == poc_dam_node)]
tree <- bind.tip(tree, "Pocillopora_acuta", where = poc_dam_node, edge.length = bl, position = 0.5 * bl)

asp_node <- which(tree$tip.label == "Acropora_aspera")
bl <- tree$edge.length[which(tree$edge[, 2] == asp_node)]
tree <- bind.tip(tree, "Acropora_downingi", where = asp_node, edge.length = bl, position = 0.5 * bl)

etp_node <- getMRCA(tree, c("Pocillopora_inflata", "Pocillopora_elegans"))
bl <- tree$edge.length[which(tree$edge[, 2] == etp_node)]
tree <- bind.tip(tree, "Pocillopora_capitata", where = etp_node, edge.length = bl, position = 0.5 * bl)

tree_final <- keep.tip(tree, unique(dat$clean.spec))
tree_final <- force.ultrametric(tree_final, method = "extend")
A          <- ape::vcv.phylo(tree_final, corr = TRUE)

# ---- Trait covariance + trait dendrogram (unchanged) ----
load("data/McMillianTraitPCs.Rdata")  # mcmillian.sc

mcmillian.sc <- mcmillian.sc %>% mutate(species = gsub(" ", "_", species))
mcmillian.sc$species[mcmillian.sc$species == "Montastrea_cavernosa"] <- "Montastraea_cavernosa"
mcmillian.sc$species[mcmillian.sc$species == "Pocillopora_eydouxi"]  <- "Pocillopora_grandis"
rownames(mcmillian.sc)  <- mcmillian.sc$species
colnames(mcmillian.sc)  <- gsub(" ", "", colnames(mcmillian.sc))

phy.dat <- dat %>% filter(clean.spec %in% rownames(A))
phy.dat <- left_join(phy.dat, mcmillian.sc, by = c("clean.spec" = "species"))
colnames(phy.dat) <- gsub(" ", "", colnames(phy.dat))
phy.dat$clean.spec       <- as.factor(phy.dat$clean.spec)
phy.dat$clean.spec_resid <- phy.dat$clean.spec
phy.dat$clean.spec_trait <- phy.dat$clean.spec
phy.dat <- phy.dat |>
  mutate(shift_ord = factor(
    case_when(
      GenShift == "NONE"           ~ 0,
      GenShift == "SUBTYPE_SWITCH" ~ 1,
      GenShift == "MINOR"          ~ 2,
      GenShift == "MAJOR"          ~ 3
    ),
    levels = c(0, 1, 2, 3), ordered = TRUE
  ))

trait.sub   <- mcmillian.sc[as.character(unique(phy.dat$clean.spec)), ]
trait.dists <- as.matrix(dist(trait.sub, method = "euclidean"))
T_mat       <- 1 - (trait.dists / max(trait.dists))
diag(T_mat) <- 1
T_mat       <- as.matrix(Matrix::nearPD(T_mat)$mat)

trait_dist_mat <- as.dist(1 - T_mat)
trait_hclust   <- hclust(trait_dist_mat, method = "average")
trait_tree     <- as.phylo(trait_hclust)

species_levels <- levels(phy.dat$clean.spec)
n_spec         <- length(species_levels)

#### ----------- 3. LOAD FITTED DELTA MODEL & EXTRACT DRAWS ----------- ####

fit_dirichlet_partition_delta <- readRDS("data/stan_dirichlet_explicit_delta.RDS")
draws_df <- as_draws_df(fit_dirichlet_partition_delta$draws())

# DELTA CHANGE: formula updated to delta_PC1/2/3 and the 3-level
# STUDY/SITE/COLONY_ID nesting, matching coral_dirichlet_partition_site.stan
# and the stan_dat_delta object from your model-fitting script.
#
# ASSUMPTION: SITE column -- this assumes phy.dat has a `SITE` column.
# If your actual site-identifier column is named differently, change
# STUDY/SITE/COLONY_ID below to match.
stan_dat_delta <- brms::standata(
  shift_ord ~ log_delta_t + delta_PC1 + delta_PC2 + delta_PC3 + last_genus +
    (1 | STUDY/SITE/COLONY_ID) +
    (1 + delta_PC1 + delta_PC2 + delta_PC3 | gr(clean.spec,       cov = A)) +
    (1 + delta_PC1 + delta_PC2 + delta_PC3 | gr(clean.spec_trait, cov = T_mat)) +
    (1 + delta_PC1 + delta_PC2 + delta_PC3 | clean.spec_resid),
  data   = phy.dat,
  data2  = list(A = A, T_mat = T_mat),
  family = cumulative(link = "logit")
)

X_col_means <- colMeans(stan_dat_delta$X)
X_cols      <- colnames(stan_dat_delta$X)

# For cumulative()/ordinal families, brms's stan_dat$X does NOT include an
# Intercept column (same as the original model) -- confirm that still
# holds for the delta model before trusting the b_names indexing below.
if ("Intercept" %in% X_cols) {
  stop("Unexpected: 'Intercept' found in stan_dat_delta$X for a cumulative() family model. ",
       "This function assumed ordinal models never include it in X -- if your brms/cmdstanr ",
       "version does something different, adjust the b_names/X_cols indexing below accordingly.")
}

b_names <- paste0("b[", seq_along(X_cols), "]")
missing_b <- setdiff(b_names, names(draws_df))
if (length(missing_b) > 0) {
  stop("Expected columns ", paste(missing_b, collapse = ", "), " not found in the fitted delta model's ",
       "draws. This means the reconstructed formula/data above doesn't match what was actually fit -- ",
       "check for any formula changes since this prep script was written.")
}
b_mat   <- as.matrix(draws_df[, b_names])
colnames(b_mat) <- X_cols

int_mat <- as.matrix(draws_df[, c("Intercept[1]", "Intercept[2]", "Intercept[3]")])

message("Delta model design matrix reconstructed and verified. Columns: ",
        paste(colnames(b_mat), collapse = ", "))

extract_re_local <- function(block, coef_slot) {
  cols <- paste0("r_", block, "_", coef_slot, "[", 1:n_spec, "]")
  missing <- setdiff(cols, names(draws_df))
  if (length(missing) > 0) stop("Missing columns: ", paste(missing[1:min(3, length(missing))], collapse = ", "), "...")
  as.matrix(draws_df[, cols])
}

# NOTE: field names kept as r1_pc1/r2_pc1/r3_pc1 etc. for compatibility
# with existing downstream scripts (helpers.R, the reviewer-response
# analyses) -- these now represent sensitivity to delta_PC1/2/3, not raw
# PC1/2/3. See the "DELTA CHANGE" comments in downstream scripts.
r1_int <- extract_re_local(1, 1); r1_pc1 <- extract_re_local(1, 2)
r1_pc2 <- extract_re_local(1, 3); r1_pc3 <- extract_re_local(1, 4)
r2_int <- extract_re_local(2, 1); r2_pc1 <- extract_re_local(2, 2)
r2_pc2 <- extract_re_local(2, 3); r2_pc3 <- extract_re_local(2, 4)
r3_int <- extract_re_local(3, 1); r3_pc1 <- extract_re_local(3, 2)
r3_pc2 <- extract_re_local(3, 3); r3_pc3 <- extract_re_local(3, 4)

category_values <- 0:3  # NONE, SUBTYPE_SWITCH, MINOR, MAJOR

mean_log_delta_t <- X_col_means["log_delta_t"]
genus_col_means  <- X_col_means[grepl("^last_genus", names(X_col_means))]

genus_reference_level <- "C"
genus_levels_available <- c(genus_reference_level, gsub("^last_genus", "", names(genus_col_means)))

#### ----------- 4a. PER-SPECIES DEFAULT (MOST COMMON) SYMBIONT GENUS ----------- ####
# UNCHANGED from 99a -- genus handling doesn't depend on raw-vs-delta.

species_default_genus <- dat %>%
  filter(clean.spec %in% species_levels, last_genus %in% genus_levels_available) %>%
  count(clean.spec, last_genus, name = "n_obs") %>%
  group_by(clean.spec) %>%
  slice_max(n_obs, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(clean.spec, default_genus = last_genus, n_obs) %>%
  mutate(default_genus = as.character(default_genus))

missing_default <- setdiff(species_levels, species_default_genus$clean.spec)
if (length(missing_default) > 0) {
  species_default_genus <- bind_rows(
    species_default_genus,
    tibble(clean.spec = missing_default, default_genus = genus_reference_level, n_obs = 0)
  )
  message(length(missing_default), " species had no single-letter genus observations -- ",
          "defaulted to reference level '", genus_reference_level, "': ",
          paste(missing_default, collapse = ", "))
}

#### ----------- 4b. PER-SPECIES SUMMARY (delta-baseline version) ----------- ####
# DELTA CHANGE: spec_summary now holds each species' MEAN delta_PC1/2/3
# (typical magnitude/direction of environmental change it has experienced),
# rather than mean raw PC1/2/3. Column names kept as PC1_mean/PC2_mean/
# PC3_mean since compute_species_switching() in helpers.R hardcodes those
# names as the "hold constant" baseline for whichever PC isn't being
# varied. If you want a strict zero (\"no change\") baseline instead of
# each species' own mean observed change, override these to 0 downstream
# (as done in the 4f tree-plotting script) rather than changing them here,
# since the Shiny app itself may still want each species' typical/observed
# delta as a sensible default.

spec_summary <- phy.dat |>
  group_by(clean.spec) |>
  summarise(
    PC1_mean = mean(delta_PC1, na.rm = TRUE),
    PC2_mean = mean(delta_PC2, na.rm = TRUE),
    PC3_mean = mean(delta_PC3, na.rm = TRUE),
    n_obs    = n(),
    .groups  = "drop"
  )

#### ----------- 5. SITE-LEVEL ENV PCA LOOKUP (named reef + month/year) ----------- ####
# UNCHANGED from 99a -- this section shows RAW/absolute environmental PC
# conditions at a site/year for display purposes (e.g. "what were
# conditions like at this reef in this month"), which is independent of
# whether the SWITCHING model itself uses raw or delta PCs.

site_year_lookup <- join.data %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         year = year(date)) %>%
  bind_cols(as.data.frame(predict(env.pc, newdata = decostand(
    join.data %>%
      select(mean_ANOM, min_ANOM, max_ANOM, cum_ANOM, mean_DHW, sd_ANOM,
             mean_SST, min_SST, max_SST, sd_SST,
             mean_precip, min_precip, max_precip, sd_precip,
             mean_chla, min_chla, max_chla, sd_chla,
             mean_turb, min_turb, max_turb, sd_turb) %>%
      mutate(across(everything(), ~ if_else(is.infinite(.), NA, .))) %>%
      mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))),
    method = "standardize"
  )))) %>%
  select(lat, lon, year, date, PC1, PC2, PC3) %>%
  mutate(month = month(date, label = TRUE, abbr = FALSE))

#### ---- 5a. Cluster nearby coordinates into the same "reef system" ---- ####
# UNCHANGED from 99a.

CLUSTER_KM <- 25

unique_sites <- site_year_lookup %>% distinct(lat, lon)

haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371; to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat / 2)^2 + cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  2 * R * asin(sqrt(pmin(1, a)))
}

n_sites <- nrow(unique_sites)
if (n_sites > 400) {
  message("NOTE: ", n_sites, " unique site coordinates -- pairwise distance clustering ",
          "may take a little while (O(n^2)).")
}
dist_mat <- matrix(0, n_sites, n_sites)
for (i in seq_len(n_sites)) {
  dist_mat[i, ] <- haversine_km(unique_sites$lat[i], unique_sites$lon[i], unique_sites$lat, unique_sites$lon)
}
unique_sites$cluster_id <- cutree(hclust(as.dist(dist_mat), method = "complete"), h = CLUSTER_KM)

cluster_centroids <- unique_sites %>%
  group_by(cluster_id) %>%
  summarise(lat = mean(lat), lon = mean(lon), n_subsites = n(), .groups = "drop")

#### ---- 5b. Broad region: point-in-polygon join against AIMS reef regions ---- ####
# UNCHANGED from 99a.

aims_region_path <- "data/World_AIMS_Reef-regions-main/finaldata/World_AIMS_Reef-regions_V1.shp"

if (file.exists(aims_region_path)) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for the region join. Install it with install.packages('sf').")
  }
  aims_regions <- sf::st_read(aims_region_path, quiet = TRUE)
  pts <- sf::st_as_sf(cluster_centroids, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

  joined <- suppressWarnings(sf::st_join(pts, aims_regions[, "Region"]))

  missing <- is.na(joined$Region)
  if (any(missing)) {
    nearest_idx <- sf::st_nearest_feature(pts[missing, ], aims_regions)
    joined$Region[missing] <- aims_regions$Region[nearest_idx]
    fallback_dist_km <- as.numeric(sf::st_distance(pts[missing, ], aims_regions[nearest_idx, ], by_element = TRUE)) / 1000
    message(sum(missing), " site cluster(s) fell outside all AIMS region polygons -- ",
            "assigned to nearest region instead (max fallback distance: ",
            round(max(fallback_dist_km), 1), " km). Worth spot-checking these.")
  }

  cluster_centroids$basin_region <- joined$Region
  message("Joined ", nrow(cluster_centroids), " site clusters to AIMS reef regions.")
} else {
  message(aims_region_path, " not found -- basin_region will be NA for all sites. ",
          "Clone https://github.com/eatlas/World_AIMS_Reef-regions and copy finaldata/* to ",
          "data/World_AIMS_Reef-regions/ to enable this.")
  cluster_centroids$basin_region <- NA_character_
}

cluster_centroids$reef_name <- NA_character_

#### ---- 5c. Editable override file (auto-generated once, never overwritten) ---- ####
# UNCHANGED from 99a.

site_name_file <- "data/site_name_lookup.csv"

if (!file.exists(site_name_file)) {
  site_name_table <- cluster_centroids %>%
    select(cluster_id, lat, lon, n_subsites, basin_region, reef_name)
  write_csv(site_name_table, site_name_file)
  message("Created ", site_name_file, " with AIMS-derived basin_region for each site cluster ",
          "(reef_name left blank pending the WCMC polygon join). Please review before deploying -- ",
          "you can hand-edit either column here, and edits won't be overwritten on rerun.")
} else {
  site_name_table <- read_csv(site_name_file, show_col_types = FALSE)
}

unique_sites <- unique_sites %>%
  left_join(site_name_table %>% select(cluster_id, basin_region, reef_name), by = "cluster_id")

site_year_lookup <- site_year_lookup %>%
  left_join(unique_sites, by = c("lat", "lon")) %>%
  mutate(reef_name = ifelse(is.na(reef_name), basin_region, reef_name))

site_year_lookup <- site_year_lookup %>%
  filter(!is.na(basin_region)) %>%
  group_by(cluster_id, reef_name, basin_region, year, month) %>%
  summarise(
    PC1 = median(PC1, na.rm = TRUE),
    PC2 = median(PC2, na.rm = TRUE),
    PC3 = median(PC3, na.rm = TRUE),
    lat = median(lat, na.rm = TRUE),
    lon = median(lon, na.rm = TRUE),
    n_subsites = n(),
    .groups = "drop"
  ) %>%
  mutate(site_label = ifelse(
    reef_name != basin_region,
    paste0(reef_name, ", ", basin_region, " (", month, " ", year, ")"),
    paste0(basin_region, " (", month, " ", year, ")")
  ))

#### ----------- 5d. PER-SPECIES ENV MEANS & OBSERVED delta-CONDITION RANGE ----------- ####
# species_raw_means (absolute env conditions, for display/auto-fill) is
# UNCHANGED. species_pc_range is a DELTA CHANGE: rebuilt from observed
# delta_PC1/2/3 ranges rather than raw PC1/2/3 ranges, since flagging
# "outside observed range" should now mean "outside the range of
# environmental CHANGE the model actually saw."

species_raw_means <- phy.dat %>%
  select(clean.spec, all_of(names(env_raw_means))) %>%
  mutate(across(all_of(names(env_raw_means)), ~ ifelse(is.infinite(.), NA, .))) %>%
  group_by(clean.spec) %>%
  summarise(across(everything(), ~ mean(., na.rm = TRUE)), .groups = "drop")

nan_mask <- is.na(as.matrix(species_raw_means %>% select(-clean.spec)))
if (any(nan_mask)) {
  affected <- which(nan_mask, arr.ind = TRUE)
  affected_desc <- sprintf("%s (%s)",
                           species_raw_means$clean.spec[affected[, 1]],
                           colnames(species_raw_means %>% select(-clean.spec))[affected[, 2]])
  message(sum(nan_mask), " species/variable combination(s) had zero non-missing observations ",
          "(likely an incomplete environmental join for that species -- worth checking data/envdata/ ",
          "coverage for it): ", paste(unique(affected_desc), collapse = ", "),
          ". Falling back to the global dataset mean for those specific combinations.")
  for (i in seq_len(nrow(affected))) {
    var_name <- colnames(species_raw_means %>% select(-clean.spec))[affected[i, 2]]
    species_raw_means[[var_name]][affected[i, 1]] <- env_raw_means[var_name]
  }
}

# DELTA CHANGE: species_site_key now carries delta_PC1/2/3 (not PC1/2/3)
# for the range-widening logic below.
# ASSUMPTION: SITE column -- LATITUDE/LONGITUDE join keys unchanged, but if
# your site identifier differs from lat/lon matching, adjust accordingly.
species_site_key <- phy.dat %>%
  select(clean.spec, LATITUDE, LONGITUDE, delta_PC1, delta_PC2, delta_PC3) %>%
  left_join(unique_sites %>% select(lat, lon, cluster_id),
           by = c("LATITUDE" = "lat", "LONGITUDE" = "lon"))

unmatched_pct <- round(100 * mean(is.na(species_site_key$cluster_id)), 1)
if (unmatched_pct > 0) {
  message(unmatched_pct, "% of phy.dat rows didn't match a site cluster by exact lat/lon ",
          "(rounding/precision mismatch between dat and join.data) -- those rows fall back to ",
          "species' own range only, without neighborhood widening.")
}

species_own_range <- species_site_key %>%
  group_by(clean.spec) %>%
  summarise(
    n_obs      = n(),
    PC1_own_min = min(delta_PC1, na.rm = TRUE), PC1_own_max = max(delta_PC1, na.rm = TRUE),
    PC2_own_min = min(delta_PC2, na.rm = TRUE), PC2_own_max = max(delta_PC2, na.rm = TRUE),
    PC3_own_min = min(delta_PC3, na.rm = TRUE), PC3_own_max = max(delta_PC3, na.rm = TRUE),
    .groups = "drop"
  )

species_clusters <- species_site_key %>%
  filter(!is.na(cluster_id)) %>%
  distinct(clean.spec, cluster_id)

neighborhood_range <- species_clusters %>%
  left_join(species_site_key %>% filter(!is.na(cluster_id)) %>% select(cluster_id, delta_PC1, delta_PC2, delta_PC3),
           by = "cluster_id", relationship = "many-to-many") %>%
  group_by(clean.spec) %>%
  summarise(
    PC1_nbhd_min = min(delta_PC1, na.rm = TRUE), PC1_nbhd_max = max(delta_PC1, na.rm = TRUE),
    PC2_nbhd_min = min(delta_PC2, na.rm = TRUE), PC2_nbhd_max = max(delta_PC2, na.rm = TRUE),
    PC3_nbhd_min = min(delta_PC3, na.rm = TRUE), PC3_nbhd_max = max(delta_PC3, na.rm = TRUE),
    .groups = "drop"
  )

# NOTE: field/column names kept as PC1_min/PC1_max etc. (not
# delta_PC1_min) for compatibility with any existing app.R code that reads
# species_pc_range -- but the VALUES are now delta ranges. Update any
# labels in app.R accordingly so the UI doesn't say "PC1" when it means
# "\u0394PC1".
species_pc_range <- species_own_range %>%
  left_join(neighborhood_range, by = "clean.spec") %>%
  mutate(
    PC1_min = pmin(PC1_own_min, PC1_nbhd_min, na.rm = TRUE), PC1_max = pmax(PC1_own_max, PC1_nbhd_max, na.rm = TRUE),
    PC2_min = pmin(PC2_own_min, PC2_nbhd_min, na.rm = TRUE), PC2_max = pmax(PC2_own_max, PC2_nbhd_max, na.rm = TRUE),
    PC3_min = pmin(PC3_own_min, PC3_nbhd_min, na.rm = TRUE), PC3_max = pmax(PC3_own_max, PC3_nbhd_max, na.rm = TRUE)
  ) %>%
  select(clean.spec, n_obs, PC1_min, PC1_max, PC2_min, PC2_max, PC3_min, PC3_max)

#### ----------- 6. TRAIT PCA FROM RAW McWILLIAM ET AL. 2018 DATA ----------- ####
# UNCHANGED from 99a -- trait PCA doesn't depend on raw-vs-delta.

mcwilliam_raw_path <- "data/McWilliam_2018/data/traitbiogeography.csv"

if (file.exists(mcwilliam_raw_path)) {

  trait_raw  <- read.csv(mcwilliam_raw_path, row.names = 1)
  trait_vars <- c("cat_growthrate", "cat_skeletaldensity", "cat_corallitesize",
                  "cat_colonydiameter", "cat_colonyheight", "cat_SA_vol", "cat_spacesize")

  traits_num <- na.omit(trait_raw[, trait_vars])
  traits_sp  <- na.omit(trait_raw[, c(trait_vars, "species")])

  trait_pca    <- prcomp(traits_num, center = TRUE, scale. = TRUE)
  trait_varexp <- round((trait_pca$sdev^2 / sum(trait_pca$sdev^2)) * 100, 1)

  trait_pca_lookup <- tibble(
    clean.spec = gsub(" ", "_", traits_sp$species),
    PC1 = trait_pca$x[, 1], PC2 = -trait_pca$x[, 2],
    PC3 = trait_pca$x[, 3], PC4 = trait_pca$x[, 4]
  ) %>%
    mutate(clean.spec = case_when(
      clean.spec == "Montastrea_cavernosa" ~ "Montastraea_cavernosa",
      clean.spec == "Pocillopora_eydouxi"  ~ "Pocillopora_grandis",
      .default = clean.spec
    ))

  trait_labels <- c(
    cat_growthrate = "Growth rate", cat_skeletaldensity = "Skeletal density",
    cat_corallitesize = "Corallite size", cat_colonydiameter = "Colony diameter",
    cat_colonyheight = "Colony height", cat_SA_vol = "SA:vol ratio", cat_spacesize = "Space size"
  )

  trait_vecs_df <- tibble(
    varname = rownames(trait_pca$rotation),
    PC1 = trait_pca$rotation[, 1],
    PC2 = -trait_pca$rotation[, 2]
  ) %>%
    mutate(label = trait_labels[varname])

} else {

  message(mcwilliam_raw_path, " not found -- falling back to mcmillian.sc scores only ",
          "(no biplot arrows available). Place McWilliam et al. 2018's raw data folder at ",
          "data/McWilliam_2018/ (so the trait csv is at data/McWilliam_2018/data/traitbiogeography.csv) ",
          "to enable the full trait biplot with arrows.")

  trait_pc_cols <- grep("^trait_?PC[0-9]+$", colnames(mcmillian.sc), value = TRUE)
  trait_pca_lookup <- mcmillian.sc %>%
    as.data.frame() %>%
    rownames_to_column("clean.spec") %>%
    select(clean.spec, all_of(trait_pc_cols)) %>%
    rename_with(~ gsub("trait_?", "", .x), all_of(trait_pc_cols))
  trait_vecs_df <- tibble(varname = character(), PC1 = numeric(), PC2 = numeric(), label = character())
  trait_varexp  <- c(NA_real_, NA_real_)
}

#### ----------- 7. BUNDLE + SAVE ----------- ####

shiny_app_data_delta <- list(
  env.pc            = env.pc,
  env_raw_means     = env_raw_means,
  env_raw_sds       = env_raw_sds,
  env_input_min     = env_input_min,
  env_input_max     = env_input_max,
  selected_env_vars = selected_env_vars,
  redundant_env_vars = redundant_env_vars,
  species_raw_means = species_raw_means,
  species_pc_range  = species_pc_range,   # now delta_PC1/2/3 ranges -- see 5d note
  site_year_lookup  = site_year_lookup,
  trait_pca_lookup  = trait_pca_lookup,
  trait_vecs_df     = trait_vecs_df,
  trait_varexp      = trait_varexp,
  tree_final        = tree_final,
  trait_tree        = trait_tree,
  species_levels    = species_levels,
  n_spec            = n_spec,
  spec_summary      = spec_summary,       # now mean delta_PC1/2/3 -- see 4b note
  species_default_genus  = species_default_genus,
  genus_levels_available = genus_levels_available,
  genus_reference_level  = genus_reference_level,
  genus_col_means        = genus_col_means,
  category_values   = category_values,
  mean_log_delta_t  = mean_log_delta_t,
  b_mat             = b_mat,
  int_mat           = int_mat,
  r1_int = r1_int, r1_pc1 = r1_pc1, r1_pc2 = r1_pc2, r1_pc3 = r1_pc3,
  r2_int = r2_int, r2_pc1 = r2_pc1, r2_pc2 = r2_pc2, r2_pc3 = r2_pc3,
  r3_int = r3_int, r3_pc1 = r3_pc1, r3_pc2 = r3_pc2, r3_pc3 = r3_pc3
)

saveRDS(shiny_app_data_delta, "data/shiny_app_data_delta.RDS")
message("Saved data/shiny_app_data_delta.RDS (",
        format(file.size("data/shiny_app_data_delta.RDS") / 1e6, digits = 3), " MB)")
