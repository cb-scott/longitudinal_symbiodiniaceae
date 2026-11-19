############################################################
##### BAYESIAN TRAIT PC INTERACTION SENSITIVITY MODEL #####
############################################################
# Tests whether host TraitPC1-TraitPC3 predict species-level sensitivity to
# environmental change in the full genus-level ordinal dataset. Trait PCs enter
# as explicit fixed effects and interactions with delta_PC1-delta_PC3. The
# trait-covariance random-effect block used by script 03a is intentionally
# omitted to avoid testing trait similarity with slopes already structured by
# the same trait matrix. Phylogenetic and idiosyncratic species-level random
# intercepts and environmental slopes are retained.


#-------- Load required packages ------------
library(tidyverse)
library(lubridate)
library(brms)
library(cmdstanr)
library(ape)
library(phangorn)
library(phytools)
library(posterior)
library(Matrix)
library(ggdist)
library(writexl)

results_dir <- "results/03h_Bayesian_GenusOrdinal_TraitPCInteractions"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

source("scripts/helpers.R")


##########################################
##### PREPARE GENUS LEVEL MODEL DATA #####
##########################################
# Rebuild the complete genus-level transition dataset used by script 03a and
# attach one standardized set of TraitPC1-TraitPC3 scores to each host species.

#-------- Build sequential genus-level outcomes ------------
load("data/LongZooxInputData_wEnv.Rdata")
# Contains: genlevel (symbiont genus-level summary)
#           all_datapoints (all observations with environmental merge)

dat <- genlevel |>
  filter(GENUS != "NA") |>
  mutate(
    Bleached = ifelse(is.na(BLEACHING_OBS) | BLEACHING_OBS == "NA", "N", BLEACHING_OBS),
    Bleached = factor(Bleached, levels = c("N", "Y"))
  ) |>
  distinct(STUDY, HOST_SPEC_REVISED, SPECIES, COLONY_ID, DATE, .keep_all = TRUE) |>
  arrange(STUDY, HOST_SPEC_REVISED, SPECIES, COLONY_ID, DATE) |>
  group_by(STUDY, HOST_SPEC_REVISED, SPECIES, COLONY_ID) |>
  mutate(
    shifted = case_when(
      GenShift %in% c("MAJOR", "MINOR", "SUBTYPE_SWITCH") ~ 1,
      GenShift == "NONE" ~ 0,
      is.na(GenShift) ~ NA_real_
    ),
    delta_t = as.numeric(DATE - lag(DATE)),
    log_delta_t = log(delta_t),
    delta_SST = max_SST - lag(min_SST, 1),
    delta_PC1 = PC1 - lag(PC1, 1),
    delta_PC2 = PC2 - lag(PC2, 1),
    delta_PC3 = PC3 - lag(PC3, 1),
    dist_PC = sqrt((PC1 - lag(PC1, 1))^2 + (PC2 - lag(PC2, 1))^2),
    last_genus = factor(lag(GENUS))
  ) |>
  filter(!is.na(delta_t), !is.na(shifted)) |>
  filter(is.finite(max_ANOM), is.finite(max_DHW), is.finite(max_SST)) |>
  separate(SPECIES, into = c("HOST_GEN", "HOST_SPEC"), sep = "_", remove = FALSE) |>
  ungroup() |>
  mutate(
    ocean_basin = assign_ocean_basin(LONGITUDE, LATITUDE),
    Year = year(DATE),
    Years_Since_2000 = Year - 2000
  )

dat$last_genus <- relevel(dat$last_genus, ref = "C")


#-------- Standardize host species names ------------
dat$HOST_SPEC <- ifelse(dat$HOST_SPEC == "INTERMEDIA(NOBILIS)", "INTERMEDIA", dat$HOST_SPEC)
dat$HOST_GEN <- case_when(
  dat$HOST_SPEC == "HYSTRIX" ~ "SERIATOPORA",
  dat$HOST_GEN == "GONIASTREA (FAVIA)" ~ "GONIASTREA",
  dat$HOST_SPEC == "PALIFERA" ~ "ISOPORA",
  .default = dat$HOST_GEN
)

dat <- dat |>
  mutate(clean.spec = paste(str_to_title(HOST_GEN), tolower(HOST_SPEC), sep = "_")) |>
  mutate(clean.spec = case_when(
    clean.spec == "Platygyra_daedalea" ~ "Platygyra_daedalea",
    clean.spec == "Porites_sp." ~ "Porites_lobata",
    clean.spec == "Pachyseris_rugosa" ~ "Pachyseris_speciosa",
    clean.spec == "Dipastraea_pallida" ~ "Dipsastraea_pallida",
    clean.spec == "Porites_harrisoni" ~ "Porites_lobata",
    clean.spec == "Cyphastrea_microphtalma" ~ "Cyphastrea_microphthalma",
    clean.spec == "Pocillopora_spp" ~ "Pocillopora_capitata",
    .default = clean.spec
  )) |>
  filter(clean.spec != "Acropora_sp.")


#-------- Build and graft the host phylogeny ------------
# Source: Gault et al. (2021), https://doi.org/10.1126/sciadv.abh4243
trees <- read.nexus("data/mtree.trees")
tree <- maxCladeCred(trees)

poc_dam_node <- getMRCA(tree, c("Pocillopora_damicornis", "Pocillopora_effusa"))
bl <- tree$edge.length[which(tree$edge[, 2] == poc_dam_node)]
tree <- bind.tip(
  tree, "Pocillopora_acuta",
  where = poc_dam_node, edge.length = bl, position = 0.5 * bl
)

asp_node <- which(tree$tip.label == "Acropora_aspera")
bl <- tree$edge.length[which(tree$edge[, 2] == asp_node)]
tree <- bind.tip(
  tree, "Acropora_downingi",
  where = asp_node, edge.length = bl, position = 0.5 * bl
)

etp_node <- getMRCA(tree, c("Pocillopora_inflata", "Pocillopora_elegans"))
bl <- tree$edge.length[which(tree$edge[, 2] == etp_node)]
tree <- bind.tip(
  tree, "Pocillopora_capitata",
  where = etp_node, edge.length = bl, position = 0.5 * bl
)


#-------- Prepare equal-weighted species trait scores ------------
# Standardize trait PCs across the modeled host species, rather than across
# observations, so heavily sampled species do not determine the scale.
load("data/McMillianTraitPCs.Rdata")

mcmillian.sc <- mcmillian.sc |>
  mutate(species = gsub(" ", "_", species))
mcmillian.sc$species[mcmillian.sc$species == "Montastrea_cavernosa"] <- "Montastraea_cavernosa"
mcmillian.sc$species[mcmillian.sc$species == "Pocillopora_eydouxi"] <- "Pocillopora_grandis"
mcmillian.sc$species[mcmillian.sc$species == "Dipsastrea_pallida"] <- "Dipsastraea_pallida"
colnames(mcmillian.sc) <- gsub(" ", "", colnames(mcmillian.sc))

trait_scores <- mcmillian.sc |>
  select(species, trait_PC1, trait_PC2, trait_PC3) |>
  distinct(species, .keep_all = TRUE) |>
  filter(species %in% unique(dat$clean.spec)) |>
  drop_na() |>
  mutate(across(starts_with("trait_PC"), ~ as.numeric(scale(.x))))

missing_traits <- setdiff(unique(dat$clean.spec), trait_scores$species)
cat("Modeled species missing TraitPC scores:", length(missing_traits), "\n")
print(missing_traits)

# Retain species represented in both the trait dataset and host phylogeny.
phy.dat <- dat |>
  filter(clean.spec %in% trait_scores$species, clean.spec %in% tree$tip.label) |>
  left_join(trait_scores, by = c("clean.spec" = "species")) |>
  mutate(
    shift_ord = factor(
      case_when(
        GenShift == "NONE" ~ 0,
        GenShift == "SUBTYPE_SWITCH" ~ 0,
        GenShift == "MINOR" ~ 1,
        GenShift == "MAJOR" ~ 2
      ),
      levels = c(0, 1, 2),
      ordered = TRUE
    ),
    SITE = paste(LONGITUDE, LATITUDE, sep = "_")
  )

stopifnot(!anyNA(phy.dat[c("trait_PC1", "trait_PC2", "trait_PC3")]))

# Prune and scale the phylogenetic covariance to the retained species.
tree_final <- keep.tip(tree, unique(phy.dat$clean.spec))
tree_final <- force.ultrametric(tree_final, method = "extend")
A <- ape::vcv.phylo(tree_final, corr = TRUE)

phy.dat$clean.spec <- factor(phy.dat$clean.spec, levels = rownames(A))
phy.dat$clean.spec_resid <- phy.dat$clean.spec

analysis_sample_size <- phy.dat |>
  summarise(
    n_studies = n_distinct(STUDY),
    n_species = n_distinct(clean.spec),
    n_colonies = n_distinct(interaction(STUDY, SPECIES, COLONY_ID, drop = TRUE)),
    n_intervals = n()
  )
write_csv(analysis_sample_size, file.path(results_dir, "Analysis_Sample_Size.csv"))
write_csv(trait_scores, file.path(results_dir, "Standardized_TraitPC_Scores.csv"))


##########################################
##### FIT TRAIT PC INTERACTION MODEL #####
##########################################
# Estimate direct TraitPC effects on baseline turnover and on sensitivity to
# each environmental PC while retaining unexplained phylogenetic and residual
# species-level variation.

#-------- Construct and verify Stan data ------------
model_formula <-
  shift_ord ~ log_delta_t + delta_PC1 + delta_PC2 + delta_PC3 + last_genus +
  trait_PC1 + trait_PC2 + trait_PC3 +
  (delta_PC1 + delta_PC2 + delta_PC3):(trait_PC1 + trait_PC2 + trait_PC3) +
  (1 | STUDY/SITE/COLONY_ID) +
  (1 + delta_PC1 + delta_PC2 + delta_PC3 | gr(clean.spec, cov = A)) +
  (1 + delta_PC1 + delta_PC2 + delta_PC3 | clean.spec_resid)

stan_dat_trait <- standata(
  model_formula,
  data = phy.dat,
  data2 = list(A = A),
  family = cumulative(link = "logit")
)

interaction_terms <- as.vector(outer(
  paste0("delta_PC", 1:3),
  paste0("trait_PC", 1:3),
  paste,
  sep = ":"
))

if (!all(interaction_terms %in% colnames(stan_dat_trait$X))) {
  stop(
    "TraitPC interaction columns are missing from the model matrix: ",
    paste(setdiff(interaction_terms, colnames(stan_dat_trait$X)), collapse = ", ")
  )
}

stopifnot(
  stan_dat_trait$M_1 == 4,
  stan_dat_trait$M_2 == 4,
  stan_dat_trait$M_3 == 1,
  stan_dat_trait$M_4 == 1,
  stan_dat_trait$M_5 == 1
)

stan_dat_trait$prior_only <- 0L


#-------- Load the saved model or sample it when absent ------------
model_file <- "data/stan_dirichlet_traitpc_interactions_ordinal.RDS"

if (file.exists(model_file)) {
  fit_traitpc_interactions <- readRDS(model_file)
} else {
  mod_trait <- cmdstan_model("scripts/coral_traitpc_interactions_ordinal.stan")

  fit_traitpc_interactions <- mod_trait$sample(
    data = stan_dat_trait,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    seed = 42,
    adapt_delta = 0.95,
    refresh = 50
  )

  fit_traitpc_interactions$save_object(file = model_file)
}


##########################################
##### SUMMARIZE TRAIT PC ASSOCIATIONS #####
##########################################
# Summarize fixed TraitPC effects and interaction coefficients. Interaction
# coefficients directly test whether each trait axis predicts species-level
# sensitivity to each environmental axis.

#-------- Extract labeled posterior draws ------------
draws_df <- posterior::as_draws_df(fit_traitpc_interactions$draws())
predictor_names <- colnames(stan_dat_trait$X)
b_names <- paste0("b[", seq_len(stan_dat_trait$Kc), "]")
names(predictor_names) <- b_names

stopifnot(all(b_names %in% names(draws_df)))

b_long <- draws_df |>
  select(all_of(b_names)) |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  mutate(term = predictor_names[variable])


#-------- Save all fixed effects and ordinal thresholds ------------
fixed_summary <- b_long |>
  group_by(term) |>
  summarise(
    median = median(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    p_positive = mean(value > 0),
    p_negative = mean(value < 0),
    .groups = "drop"
  )

threshold_summary <- draws_df |>
  select(starts_with("Intercept[")) |>
  pivot_longer(everything(), names_to = "term", values_to = "value") |>
  group_by(term) |>
  summarise(
    median = median(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    p_positive = mean(value > 0),
    p_negative = mean(value < 0),
    .groups = "drop"
  )

writexl::write_xlsx(
  bind_rows(fixed_summary, threshold_summary),
  file.path(results_dir, "Fixed_Effects_TraitPC_Interaction_Model.xlsx")
)


#-------- Summarize TraitPC by environmental-PC interactions ------------
interaction_summary <- b_long |>
  filter(term %in% interaction_terms) |>
  separate(term, into = c("environment", "trait"), sep = ":", remove = FALSE) |>
  group_by(term, environment, trait) |>
  summarise(
    median = median(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    p_positive = mean(value > 0),
    p_negative = mean(value < 0),
    directional_support = pmax(p_positive, p_negative),
    .groups = "drop"
  ) |>
  mutate(
    environment = recode(
      environment,
      delta_PC1 = "Delta PC1",
      delta_PC2 = "Delta PC2",
      delta_PC3 = "Delta PC3"
    ),
    trait = recode(
      trait,
      trait_PC1 = "Trait PC1",
      trait_PC2 = "Trait PC2",
      trait_PC3 = "Trait PC3"
    )
  )

writexl::write_xlsx(
  interaction_summary,
  file.path(results_dir, "TraitPC_by_EnvironmentalPC_Interactions.xlsx")
)

interaction_plot_data <- b_long |>
  filter(term %in% interaction_terms) |>
  separate(term, into = c("environment", "trait"), sep = ":") |>
  mutate(
    environment = factor(environment, levels = paste0("delta_PC", 1:3), labels = paste("Delta PC", 1:3)),
    trait = factor(trait, levels = paste0("trait_PC", 1:3), labels = paste("Trait PC", 1:3))
  )

interaction_plot <- ggplot(interaction_plot_data, aes(x = value, y = trait)) +
  stat_halfeye(
    .width = c(0.5, 0.95),
    point_interval = median_qi,
    fill = "#276C73",
    alpha = 0.8,
    point_size = 1
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~environment, nrow = 1, scales = "free_x") +
  labs(
    x = "Trait PC interaction effect on cumulative log-odds",
    y = NULL
  ) +
  theme_bw(base_size = 10)

ggsave(
  file.path(results_dir, "TraitPC_by_EnvironmentalPC_Interactions.svg"),
  plot = interaction_plot,
  width = 6.5,
  height = 3
)


#-------- Summarize remaining species-level variance ------------
# After explicit TraitPC effects are included, h2 partitions unexplained
# species-level variance between phylogenetic and idiosyncratic components.
h2_draws <- draws_df |>
  select(starts_with("h2[")) |>
  pivot_longer(everything(), names_to = "parameter", values_to = "proportion") |>
  mutate(
    coefficient_index = as.integer(str_extract(parameter, "(?<=\\[)\\d(?=,)")),
    source_index = as.integer(str_extract(parameter, "(?<=,)\\d(?=\\])")),
    coefficient = factor(
      coefficient_index,
      levels = 1:4,
      labels = c("Baseline turnover", "Delta PC1 sensitivity", "Delta PC2 sensitivity", "Delta PC3 sensitivity")
    ),
    source = factor(source_index, levels = 1:2, labels = c("Phylogenetic", "Idiosyncratic"))
  )

h2_summary <- h2_draws |>
  group_by(coefficient, source) |>
  summarise(
    median = median(proportion),
    lower = quantile(proportion, 0.025),
    upper = quantile(proportion, 0.975),
    .groups = "drop"
  )

h2_probabilities <- h2_draws |>
  group_by(coefficient, source) |>
  summarise(
    p_larger_half = mean(proportion > 0.5),
    .groups = "drop"
  )

writexl::write_xlsx(
  list(
    variance_partition = h2_summary,
    posterior_probabilities = h2_probabilities
  ),
  file.path(results_dir, "Residual_Species_Variance_Partition.xlsx")
)


##########################################
##### CALCULATE SPECIES LEVEL SLOPES #####
##########################################
# Combine population environmental effects, TraitPC interactions, and the two
# remaining species random-effect components to recover each species' total
# response to Delta PC1-Delta PC3.

#-------- Map Stan species indices to host species ------------
model_vars <- c(
  "shift_ord", "log_delta_t", "delta_PC1", "delta_PC2", "delta_PC3",
  "last_genus", "trait_PC1", "trait_PC2", "trait_PC3",
  "STUDY", "SITE", "COLONY_ID", "clean.spec", "clean.spec_resid"
)

model_dat <- phy.dat |>
  select(all_of(model_vars)) |>
  drop_na()

if (
  nrow(model_dat) != stan_dat_trait$N ||
  !identical(as.integer(model_dat$shift_ord), as.integer(stan_dat_trait$Y))
) {
  stop("Model observations do not match standata; check filtering and row order.")
}

species_map <- model_dat |>
  mutate(j_phylo = stan_dat_trait$J_1, j_resid = stan_dat_trait$J_2) |>
  transmute(
    species = as.character(clean.spec),
    trait_PC1, trait_PC2, trait_PC3,
    j_phylo, j_resid
  ) |>
  distinct() |>
  arrange(species)

stopifnot(
  !anyDuplicated(species_map$species),
  !anyDuplicated(species_map$j_phylo),
  !anyDuplicated(species_map$j_resid)
)


#-------- Extract fixed and random slope draws ------------
b_draws <- as.matrix(draws_df[, b_names, drop = FALSE])
colnames(b_draws) <- colnames(stan_dat_trait$X)

extract_random_slope <- function(group, slot, indices) {
  n_levels <- stan_dat_trait[[paste0("N_", group)]]
  columns <- paste0("r_", group, "[", seq_len(n_levels), ",", slot, "]")
  if (!all(columns %in% names(draws_df))) {
    stop("Missing random-effect draws: ", paste(head(setdiff(columns, names(draws_df))), collapse = ", "))
  }
  as.matrix(draws_df[, columns, drop = FALSE])[, indices, drop = FALSE]
}

find_interaction_column <- function(environment, trait) {
  candidates <- c(paste(environment, trait, sep = ":"), paste(trait, environment, sep = ":"))
  found <- intersect(candidates, colnames(b_draws))
  if (length(found) != 1) {
    stop("Could not uniquely identify interaction for ", environment, " and ", trait)
  }
  found
}


#-------- Summarize total slopes for every species ------------
species_slope_summaries <- map_dfr(seq_len(3), function(environment_index) {
  environment <- paste0("delta_PC", environment_index)
  slot <- environment_index + 1

  fixed_slopes <- matrix(
    b_draws[, environment],
    nrow = nrow(b_draws),
    ncol = nrow(species_map)
  )

  for (trait_index in seq_len(3)) {
    trait <- paste0("trait_PC", trait_index)
    interaction <- find_interaction_column(environment, trait)
    fixed_slopes <- fixed_slopes + outer(
      b_draws[, interaction],
      species_map[[trait]],
      "*"
    )
  }

  total_slopes <- fixed_slopes +
    extract_random_slope(1, slot, species_map$j_phylo) +
    extract_random_slope(2, slot, species_map$j_resid)

  tibble(
    species = species_map$species,
    trait_PC1 = species_map$trait_PC1,
    trait_PC2 = species_map$trait_PC2,
    trait_PC3 = species_map$trait_PC3,
    environmental_PC = paste0("Delta PC", environment_index),
    median = apply(total_slopes, 2, median),
    lower = apply(total_slopes, 2, quantile, probs = 0.025),
    upper = apply(total_slopes, 2, quantile, probs = 0.975),
    p_positive = colMeans(total_slopes > 0),
    p_negative = colMeans(total_slopes < 0)
  )
})

writexl::write_xlsx(
  species_slope_summaries,
  file.path(results_dir, "Species_Environmental_Slopes_with_TraitPCs.xlsx")
)


#-------- Plot species slopes against each TraitPC axis ------------
slope_plot_data <- species_slope_summaries |>
  pivot_longer(
    cols = starts_with("trait_PC"),
    names_to = "trait_PC",
    values_to = "trait_score"
  ) |>
  mutate(
    trait_PC = factor(trait_PC, levels = paste0("trait_PC", 1:3), labels = paste("Trait PC", 1:3)),
    environmental_PC = factor(environmental_PC, levels = paste0("Delta PC", 1:3))
  )

slope_trait_plot <- ggplot(slope_plot_data, aes(x = trait_score, y = median)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, alpha = 0.25) +
  geom_point(color = "#276C73", size = 1.4, alpha = 0.8) +
  facet_grid(environmental_PC ~ trait_PC, scales = "free_x") +
  labs(
    x = "Standardized host trait score",
    y = "Species-specific environmental slope\n(cumulative logit scale)"
  ) +
  theme_bw(base_size = 9)

ggsave(
  file.path(results_dir, "Species_Slopes_by_TraitPC.svg"),
  plot = slope_trait_plot,
  width = 7,
  height = 6.5
)


#-------- Save sampler diagnostics ------------
sampler_diagnostics <- fit_traitpc_interactions$diagnostic_summary()

# CmdStan returns separate tables for divergences, tree depth, and E-BFMI.
if (is.data.frame(sampler_diagnostics)) {
  sampler_diagnostics <- list(summary = sampler_diagnostics)
}

writexl::write_xlsx(
  lapply(sampler_diagnostics, as.data.frame),
  file.path(results_dir, "Sampler_Diagnostics.xlsx")
)
