######################################
# BAYESIAN GENUS-LEVEL ORDINAL MODEL #
######################################
# Builds phylogenetic and trait covariance matrices for the coral host species,
# then fits the main ordinal model testing symbiont switching against
# environmental PCs and prior symbiont state. Species-level variance is
# partitioned into phylogenetic, trait-associated, and residual components.
# The script produces convergence diagnostics, fixed/random-effect summaries,
# variance-partition ternary plots, and phylogeny/trait-tree response maps.
# The fitted model is provided as an .RDS object so sampling can be skipped.


#-------- Load required packages ------------
library(tidyverse)
library(brms)
library(cmdstanr)
library(ggtree)
library(patchwork)
library(ape)
library(phytools)
library(posterior)
library(Matrix)
library(ggdist)
library(Ternary)
library(viridisLite)
library(writexl)

# Create the script-specific results directory before writing outputs.
dir.create("results/03a_Bayesian_GenusOrdinal", recursive = TRUE, showWarnings = FALSE)

source("scripts/helpers.R")

##################################
# PREPARE GENUS-LEVEL MODEL DATA #
##################################
# Build sequential genus-level outcomes and align the environmental, host,
# site, and colony identifiers required by the hierarchical model.

# Load the genus-level interval dataset created by script 02a.
load("data/LongZooxInputData_wEnv.Rdata")
# Contains: genlevel (symbiont genus-level summary)
#           all_datapoints (all observations with environmental merge)

dat <- genlevel |>
  filter(GENUS != "NA") %>%
  mutate(
    Bleached = ifelse(is.na(BLEACHING_OBS) | BLEACHING_OBS == "NA", "N", BLEACHING_OBS),
    Bleached = factor(Bleached, levels = c("N", "Y"))
  ) %>%
  distinct(STUDY, HOST_SPEC_REVISED, SPECIES, COLONY_ID, DATE, .keep_all = TRUE) |>
  arrange(STUDY, HOST_SPEC_REVISED, SPECIES, COLONY_ID, DATE) |>
  group_by(STUDY, HOST_SPEC_REVISED, SPECIES, COLONY_ID) |>
  mutate(
    shifted     = case_when(
      GenShift %in% c("MAJOR", "MINOR", "SUBTYPE_SWITCH") ~ 1, #actually no longer count these subtype switches as per reviewer suggestions
      GenShift %in% c("NONE")                              ~ 0,
      is.na(GenShift)                                      ~ NA_real_
    ),
    delta_t     = as.numeric(DATE - lag(DATE)),
    log_delta_t = log(delta_t),
    delta_SST   = max_SST - lag(min_SST, 1),
    delta_PC1   = PC1 - lag(PC1, 1),
    delta_PC2 = PC2 - lag(PC2, 1),
    delta_PC3 = PC3 - lag(PC3, 1), 
    dist_PC     = sqrt((PC1 - lag(PC1, 1))^2 + (PC2 - lag(PC2, 1))^2),
    last_genus  = factor(lag(GENUS))
  ) |>
  filter(!is.na(delta_t), !is.na(shifted)) %>%
  filter(is.finite(max_ANOM)) %>%
  filter(is.finite(max_DHW))  %>%
  filter(is.finite(max_SST))  %>%
  separate(SPECIES, into = c("HOST_GEN", "HOST_SPEC"), sep = "_", remove = FALSE) %>%
  ungroup() %>%
  mutate(
    ocean_basin    = assign_ocean_basin(LONGITUDE, LATITUDE),
    Year           = year(DATE),
    Years_Since_2000 = Year - 2000
  )

dat$last_genus <- relevel(dat$last_genus, ref = "C")


#-------- Standardize host species names ------------
dat$HOST_SPEC <- ifelse(dat$HOST_SPEC == "INTERMEDIA(NOBILIS)", "INTERMEDIA", dat$HOST_SPEC)
dat$HOST_GEN  <- case_when(
  dat$HOST_SPEC == "HYSTRIX"                     ~ "SERIATOPORA",
  dat$HOST_GEN  == "GONIASTREA (FAVIA)"          ~ "GONIASTREA",
  dat$HOST_SPEC == "PALIFERA"                    ~ "ISOPORA",
  .default = dat$HOST_GEN
)

dat <- dat %>%
  mutate(clean.spec = paste(str_to_title(HOST_GEN), tolower(HOST_SPEC), sep = "_")) %>%
  mutate(clean.spec = case_when(
    clean.spec == "Platygyra_daedalea"  ~ "Platygyra_daedalea",
    clean.spec == "Porites_sp."         ~ "Porites_lobata",
    clean.spec == "Pachyseris_rugosa"   ~ "Pachyseris_speciosa",
    clean.spec == "Dipastraea_pallida" ~ "Dipsastraea_pallida",
    clean.spec == "Porites_harrisoni" ~ "Porites_lobata",
    clean.spec == "Cyphastrea_microphtalma" ~ "Cyphastrea_microphthalma",
    clean.spec == "Pocillopora_spp" ~ "Pocillopora_capitata",
    .default = clean.spec
  )) %>%
  filter(clean.spec != "Acropora_sp.")


#-------- Build the host phylogenetic covariance matrix ------------
# Source: Gault, J. A., Bentlage, B., Huang, D., & Kerr, A. M. (2021).
# Lineage-specific variation in the evolutionary stability of coral
# photosymbiosis. Science Advances. https://doi.org/10.1126/sciadv.abh4243
# Data:   https://datadryad.org/dataset/doi:10.5061/dryad.tdz08kq0p

trees     <- read.nexus("data/mtree.trees")
tree      <- maxCladeCred(trees)

# Inspect coverage
missing <- setdiff(unique(dat$clean.spec), tree$tip.label)
cat("Species missing from tree:", length(missing), "\n")
print(missing)

#-------- Graft missing species onto the host phylogeny ------------

# 1. Pocillopora acuta — sister to P. damicornis / P. effusa clade
poc_dam_node <- getMRCA(tree, c("Pocillopora_damicornis", "Pocillopora_effusa"))
bl <- tree$edge.length[which(tree$edge[, 2] == poc_dam_node)]
tree <- bind.tip(tree, "Pocillopora_acuta",
                 where = poc_dam_node, edge.length = bl, position = 0.5 * bl)

# 2. Acropora downingi — sister to A. aspera
asp_node <- which(tree$tip.label == "Acropora_aspera")
bl <- tree$edge.length[which(tree$edge[, 2] == asp_node)]
tree <- bind.tip(tree, "Acropora_downingi",
                 where = asp_node, edge.length = bl, position = 0.5 * bl)

# 3. Pocillopora capitata — sister to ETP clade (inflata + elegans)
etp_node <- getMRCA(tree, c("Pocillopora_inflata", "Pocillopora_elegans"))
bl <- tree$edge.length[which(tree$edge[, 2] == etp_node)]
tree <- bind.tip(tree, "Pocillopora_capitata",
                 where = etp_node, edge.length = bl, position = 0.5 * bl)

stopifnot(length(setdiff(unique(dat$clean.spec), tree$tip.label)) == 0)

#-------- Prune the phylogeny and calculate host covariance ------------
tree_final <- keep.tip(tree, unique(dat$clean.spec))
tree_final <- force.ultrametric(tree_final, method = "extend")
A          <- ape::vcv.phylo(tree_final, corr = TRUE)

plot(tree_final)
#-------- Build the host trait-similarity covariance matrix ------------
# McWilliam, Hoogenboom & Baird (2018) trait PCs
# Source: https://github.com/mikemcwilliam/biogeographical_disparity

load("data/McMillianTraitPCs.Rdata")

mcmillian.sc <- mcmillian.sc %>%
  mutate(species = gsub(" ", "_", species))
mcmillian.sc$species[mcmillian.sc$species == "Montastrea_cavernosa"]  <- "Montastraea_cavernosa"
mcmillian.sc$species[mcmillian.sc$species == "Pocillopora_eydouxi"]   <- "Pocillopora_grandis"
mcmillian.sc$species[mcmillian.sc$species == "Dipsastrea_pallida"] <- "Dipsastraea_pallida"
rownames(mcmillian.sc) <- mcmillian.sc$species
colnames(mcmillian.sc) <- gsub(" ", "", colnames(mcmillian.sc))

#-------- Assemble the model-ready analysis dataset ------------
phy.dat <- dat %>% filter(clean.spec %in% rownames(A))
phy.dat <- left_join(phy.dat, mcmillian.sc, by = c("clean.spec" = "species"))
colnames(phy.dat) <- gsub(" ", "", colnames(phy.dat))

phy.dat$clean.spec       <- as.factor(phy.dat$clean.spec)
phy.dat$clean.spec_resid <- phy.dat$clean.spec

phy.dat <- phy.dat |>
  mutate(
    shift_ord = factor(
      case_when(
        GenShift == "NONE"          ~ 0,
        GenShift == "SUBTYPE_SWITCH" ~ 0, ###CHANGED HERE TO BETTER REFLECT REVIEWRE COMMENTS - ONLY LOOK AT GENUS LEVEL SHIFTS FIRST
        GenShift == "MINOR"          ~ 1,
        GenShift == "MAJOR"          ~ 2
      ),
      levels  = c(0,1, 2),
      ordered = TRUE
    )
  )

category_values <- as.numeric(levels(factor(phy.dat$shift_ord)))

#-------- Calculate trait distances and covariance ------------
trait.sub   <- mcmillian.sc[as.character(unique(phy.dat$clean.spec)), ] %>% drop_na()
trait.dists <- as.matrix(dist(trait.sub, method = "euclidean"))
T_mat       <- 1 - (trait.dists / max(trait.dists))
diag(T_mat) <- 1
T_mat       <- as.matrix(nearPD(T_mat)$mat)   # ensure positive-definite

phy.dat$clean.spec_trait <- phy.dat$clean.spec

#-------- Build the trait dendrogram used in figures ------------
trait_dist_mat <- as.dist(1 - T_mat)
trait_hclust   <- hclust(trait_dist_mat, method = "average")   # UPGMA
trait_tree     <- as.phylo(trait_hclust)

plot(trait_tree)

phy.dat <- phy.dat %>% mutate(SITE = paste(LONGITUDE, LATITUDE, sep = "_"))

#####################################
#FIT THE GENUS-LEVEL ORDINAL MODEL #
####################################
# Translate the brms formula into Stan data, compile the custom
# variance-partition model, sample the posterior, and save the fitted object.

#-------- Prepare Stan data and verify model-matrix names ------------
stan_dat_delta <- standata(
  shift_ord ~ log_delta_t + delta_PC1 + delta_PC2 + delta_PC3 + last_genus +
    (1 | STUDY/SITE/COLONY_ID) +
    (1 + delta_PC1 + delta_PC2 + delta_PC3 | gr(clean.spec, cov = A)) +
    (1 + delta_PC1 + delta_PC2 + delta_PC3 | gr(clean.spec_trait, cov = T_mat)) +
    (1 + delta_PC1 + delta_PC2 + delta_PC3 | clean.spec_resid),
  data  = phy.dat,
  data2 = list(A = A, T_mat = T_mat),
  family = cumulative(link = "logit")
)

names(stan_dat_delta)[grep("^Z_[123]", names(stan_dat_delta))]

#-------- Compile and sample the custom ordinal Stan model ------------
mod_within <- cmdstan_model("scripts/coral_dirichlet_partition_site.stan", force_recompile = T)

stan_dat_delta$prior_only <- 0L

fit_dirichlet_partition_within <- mod_within$sample(
  data            = stan_dat_delta,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  seed            = 42,
  adapt_delta     = 0.95,
  refresh         = 50
)

fit_dirichlet_partition_within$save_object(
  file = "data/stan_dirichlet_explicit_delta.RDS"
)

fit_dirichlet_partition_within = readRDS("data/stan_dirichlet_explicit_delta.RDS")

#########################################
# SUMMARIZE THE ORDINAL MODEL POSTERIOR #
#########################################

# Convert posterior draws to labeled summaries, inspect fixed and random
# effects, and quantify the phylogenetic, trait, and residual variance partition.

#-------- Extract posterior draws and reconstruct parameter names ------------
draws    <- fit_dirichlet_partition_within$draws()
draws_df <- as_draws_df(draws)

# Reconstruct centered design matrix (for PPC)
Xc      <- sweep(stan_dat_delta$X, 2, colMeans(stan_dat_delta$X), "-")
y_obs   <- stan_dat_delta$Y
b_names <- paste0("b[", 1:stan_dat_delta$Kc, "]")
b_draws <- as.matrix(draws_df[, b_names])
intercept_draws <- as.matrix(draws_df[, c("Intercept[1]", "Intercept[2]")])
predictor_labels <- colnames(stan_dat_delta$X)
names(predictor_labels) <- b_names


#-------- Summarize fixed effects ------------
fes = draws_df |>
  select(all_of(b_names)) |>
  pivot_longer(everything(), names_to = "param", values_to = "value") |>
  mutate(predictor = predictor_labels[param],
         predictor = factor(predictor, levels = rev(predictor_labels))) |>
  mutate(category = case_when(
    grepl("PC", predictor) == T ~ "Environmental",
    grepl("last_genus", predictor) == T~"Symbiont",
    .default = "other"
  )) %>% mutate(predictor = as.character(predictor)) %>%
  mutate(predictor = ifelse(category == "Symbiont",
                            gsub("last_genus", "", predictor), predictor)) %>%
  mutate(predictor = gsub("delta_", "\u0394", predictor)) %>%
  mutate(is_mixed = ifelse(grepl("_", predictor) == T, "Y", "N")) %>%
  arrange(is_mixed, predictor)

fes$predictor = factor(fes$predictor, levels = rev(unique(fes$predictor)))

fe.p <- fes |> filter(category != "other") |>
  filter(!(predictor %in% c("A_C_D", "B_C_D", "D_G", "C_G", "A_C_G"))) |>
  ggplot(aes(x = value, y = predictor)) +
  stat_halfeye(.width = c(0.5, 0.95), point_interval = median_qi,
               aes(fill = interaction(category, is_mixed)), alpha = 0.8, point_size =1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(x = "Log-odds of genus turnover (cumulative)", y = NULL) +
  labs(fill = "Codominance") + 
  #    title = "Fixed effect posteriors",
  #    subtitle = "Median with 50% and 95% credible intervals") +
  theme_bw(base_size = 10) + 
  xlim(c(-2, 3)) + 
  theme(legend.position = "none")+
  scale_fill_manual(values = c("#c3d781", "#81C3D7", "#d7c081"))+
  facet_wrap(~category, nrow = 2, scales = "free_y")
fe.p

ggsave(fe.p, file = "results/03a_Bayesian_GenusOrdinal/fes_zoox_enviro.svg", width =3.25, height = 4)

#-------- Summarize fixed and random effects ------------
fixed_summary <- b_draws |>
  as_draws_df() |>
  summarize_draws(median, ~quantile2(.x, probs = c(0.025, 0.975))) |>
  mutate(term = predictor_labels[variable]) |>
  select(term, median, q2.5, q97.5)

# thresholds (cumulative logit intercepts)
threshold_summary <- intercept_draws |>
  as_draws_df() |>
  summarize_draws(median, ~quantile2(.x, probs = c(0.025, 0.975))) |>
  mutate(term = paste0("Threshold ", variable)) |>
  select(term, median, q2.5, q97.5)

fixed_effects_table <- bind_rows(fixed_summary, threshold_summary)
writexl::write_xlsx(fixed_effects_table, "results/03a_Bayesian_GenusOrdinal/FE_DiriOrdinal.xlsx")


# vector-valued groups (Intercept, PC1, PC2, PC3 each)
sd_labels_vec <- list(
  sd_1 = c("HOST.SPECIES_Phylo_Intercept", "HOST.SPECIES_Phylo_PC1",
           "HOST.SPECIES_Phylo_PC2",       "HOST.SPECIES_Phylo_PC3"),
  sd_2 = c("HOST.SPECIES_Idio_Intercept", "HOST.SPECIES_Idio_PC1",
           "HOST.SPECIES_Idio_PC2",       "HOST.SPECIES_Idio_PC3"),
  sd_3 = c("HOST.SPECIES_Trait_Intercept", "HOST.SPECIES_Trait_PC1",
           "HOST.SPECIES_Trait_PC2",       "HOST.SPECIES_Trait_PC3")
)

vec_summary <- purrr::imap_dfr(sd_labels_vec, function(labels, sd_var) {
  draws <- fit_dirichlet_partition_within$draws(variables = sd_var, format = "draws_df") |>
    subset_draws(variable = sd_var)
  
  summarize_draws(draws, median, ~quantile2(.x, probs = c(0.025, 0.975))) |>
    mutate(term = labels[match(variable, paste0(sd_var, "[", seq_along(labels), "]"))]) |>
    select(term, median, q2.5, q97.5)
})

# scalar groups (intercept-only)
scalar_labels <- c(
  sd_4 = "STUDY_Intercept",
  sd_5 = "STUDY:SITE_Intercept",
  sd_6 = "STUDY:SITE:COLONY_ID_Intercept"
)

scalar_summary <- purrr::imap_dfr(scalar_labels, function(label, sd_var) {
  draws <- fit_dirichlet_partition_within$draws(variables = sd_var, format = "draws_df")
  
  summarize_draws(draws, median, ~quantile2(.x, probs = c(0.025, 0.975))) |>
    mutate(term = label) |>
    select(term, median, q2.5, q97.5)
})

random_summary <- bind_rows(vec_summary, scalar_summary)
writexl::write_xlsx(random_summary, "results/03a_Bayesian_GenusOrdinal/RE_DiriOrdinal.xlsx")

#-------- Inspect random-effect standard deviations ------------

all_vars <- c(scalar_labels, sd_labels_vec)

raw_draws_long <- purrr::imap_dfr(all_vars, function(labels, sd_var) {
  draws <- fit_dirichlet_partition_within$draws(variables = sd_var, format = "draws_df") |>
    as_draws_df()
  
  # pivot the variable column(s) into long format, keep draw/chain ids
  long <- draws |>
    pivot_longer(
      cols = starts_with(sd_var),
      names_to = "variable",
      values_to = "value"
    )
  
  if (length(labels) == 1) {
    # scalar group: sd_4, sd_5 -- one label, no bracket indexing needed
    long <- long |> mutate(term = labels)
  } else {
    
    long <- long |>
      mutate(term = labels[match(variable, paste0(sd_var, "[", seq_along(labels), "]"))])
  }
  
  long |> select(.chain, .iteration, .draw, term, value)
})

raw_draws_long %>%
  ggplot(aes(x = value, y = term)) +
  stat_halfeye(.width = c(0.5, 0.95), point_interval = median_qi,
               normalize = "panels", alpha = 0.85) +
  xlab("sd(cumulative log odds)") + ylab("Random Effect") +
  theme_bw() +
  theme(text=element_text(size = 8))

ggsave("results/03a_Bayesian_GenusOrdinal/Posterior_REs.svg", width = 4, height = 6)


#-------- Calculate and plot the variance partition ------------

h2_draws <- draws_df |>
  select(starts_with("h2")) |>
  pivot_longer(everything(), names_to = "param", values_to = "value") |>
  mutate(
    coef_idx   = as.integer(str_extract(param, "(?<=\\[)\\d(?=,)")),
    source_idx = as.integer(str_extract(param, "(?<=,)\\d(?=\\])")),
    coefficient = factor(coef_idx, levels = 1:4,
                         labels = c("Baseline rate", "PC1 sensitivity",
                                    "PC2 sensitivity", "PC3 sensitivity")),
    source = factor(source_idx, levels = 1:3,
                    labels = c("Phylogenetic", "Idiosyncratic", "Trait"))
  )

h2_draws |>
  group_by(coefficient, source) |>
  mutate(.draw = row_number()) |>
  ungroup() |>
  pivot_wider(
    id_cols    = c(coefficient, .draw),
    names_from = source,
    values_from = value
  ) |>
  group_by(coefficient) |>
  summarise(
    p_phylo_above_uniform = mean(Phylogenetic > 1/3),
    p_trait_above_uniform = mean(Trait > 1/3),
    p_resid_above_uniform = mean(Idiosyncratic > 1/3),
    p_trait_leads         = mean(Trait > Phylogenetic & Trait > Idiosyncratic),
    p_phylo_leads         = mean(Phylogenetic > Trait & Phylogenetic > Idiosyncratic),
    p_resid_leads         = mean(Idiosyncratic > Trait & Idiosyncratic > Phylogenetic)
  )

h2_draws |>
  group_by(coefficient, source) |>
  mutate(.draw = row_number()) |>
  ungroup() |>
  pivot_wider(
    id_cols     = c(coefficient, .draw),
    names_from  = source,
    values_from = value
  ) |>
  group_by(coefficient) |>
  summarise(
    # probability of dominating
    p_trait_leads         = mean(Trait > Phylogenetic & Trait > Idiosyncratic),
    p_phylo_leads         = mean(Phylogenetic > Trait & Phylogenetic > Idiosyncratic),
    p_resid_leads         = mean(Idiosyncratic > Trait & Idiosyncratic > Phylogenetic),
    # mean proportion when dominating
    trait_median_when_leads = median(Trait[Trait > Phylogenetic & Trait > Idiosyncratic]),
    phylo_median_when_leads = median(Phylogenetic[Phylogenetic > Trait & Phylogenetic > Idiosyncratic]),
    resid_median_when_leads = median(Idiosyncratic[Idiosyncratic > Trait & Idiosyncratic > Phylogenetic])
  )

draws_df |>
  select(starts_with("total_var")) |>
  pivot_longer(everything(), names_to = "param", values_to = "value") |>
  mutate(coefficient = factor(
    str_extract(param, "\\d+"),
    levels = 1:4,
    labels = c("Baseline rate", "PC1 sensitivity", 
               "PC2 sensitivity", "PC3 sensitivity")
  )) |>
  group_by(coefficient) |>
  summarise(
    median = median(value),
    lower  = quantile(value, 0.025),
    upper  = quantile(value, 0.975)
  )

#look at univariate density histograms
ggplot(h2_draws, aes(x = value, y = source, fill = source)) +
  stat_halfeye(.width = c(0.5, 0.95), point_interval = median_qi,
               normalize = "panels", alpha = 0.85) +
  facet_wrap(~coefficient, nrow = 1) +
  scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_fill_manual(values = c("Phylogenetic"  = "#2166ac",
                               "Trait"         = "#d6604d",
                               "Idiosyncratic" = "#4dac26")) +
  labs(x = "Fraction of species-level variance", y = NULL,
       title    = "Variance decomposition of species-level switching",
       subtitle = "Direct Dirichlet(1,1,1) partition prior | median with 50% and 95% intervals") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

# Ternary plots
dens.draw_full <- h2_draws |>
  group_by(source, coefficient) |>
  mutate(.draw = row_number()) |>
  ungroup() |>
  select(value, coefficient, source, .draw) |>
  pivot_wider(names_from = source, values_from = value)


partition_probabilities <- h2_wide |>
  group_by(coefficient) |>
  summarise(
    p_phylo_above_third = mean(Phylogenetic > 1/3),
    p_trait_above_third = mean(Trait > 1/3),
    p_resid_above_third = mean(Idiosyncratic > 1/3),
    p_phylo_leads = mean(Phylogenetic > Trait & Phylogenetic > Idiosyncratic),
    p_trait_leads = mean(Trait > Phylogenetic & Trait > Idiosyncratic),
    p_resid_leads = mean(Idiosyncratic > Trait & Idiosyncratic > Phylogenetic),
    .groups = "drop"
  )
print(partition_probabilities)

writexl::write_xlsx(partition_probabilities,
                    "results/03a_Bayesian_GenusOrdinal/Variance_Partition_Probabilities_DiriOrdinal.xlsx")



h2_summary <- h2_draws |>
  group_by(coefficient, source) |>
  summarise(median = median(value), mean = mean(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975), .groups = "drop")
print(h2_summary)
writexl::write_xlsx(h2_summary, "results/03a_Bayesian_GenusOrdinal/Variance_Partition_DiriOrdinal.xlsx")

trait_when_dominant <- h2_wide %>%
  group_by(coefficient) %>%
  summarise(
    p_trait_dominant = mean(
      Trait > Phylogenetic & Trait > Idiosyncratic
    ),
    median_trait_when_dominant = median(
      Trait[Trait > Phylogenetic & Trait > Idiosyncratic]
    ),
    .groups = "drop"
  )

print(trait_when_dominant)


phylo_when_dominant <- h2_wide %>%
  group_by(coefficient) %>%
  summarise(
    p_phylo_dominant = mean(
      Phylogenetic > Trait & Phylogenetic > Idiosyncratic
    ),
    median_phylo_when_dominant = median(
      Phylogenetic[Phylogenetic > Trait & Phylogenetic > Idiosyncratic]
    ),
    .groups = "drop"
  )

phylo_when_dominant
coef_labels <- c("Baseline rate", "\u0394PC1", "\u0394PC2", "\u0394PC3")
names(coef_labels) = c("Baseline rate", "PC1 sensitivity", "PC2 sensitivity", "PC3 sensitivity")

svg("results/03a_Bayesian_GenusOrdinal/Variance_Partitioned_delta.svg", width = 6.5, height = 3)
#dev.off()
par(mfrow = c(1, 4), mar = c(1, 1, 2, 1), cex = 8/12)
for (cf in names(coef_labels)) {
  sub <- dens.draw_full |>
    filter(coefficient == cf) |>
    select(Trait, Phylogenetic, Idiosyncratic) |>
    as.matrix()
  TernaryPlot(
    alab = "Trait \u2192", blab = "Phylo. \u2192", clab = "\u2190 Idiosyn.",
    main = coef_labels[cf],
    grid.minor.lines = 0, lab.offset = 0.11,
    axis.labels = FALSE, axis.tick = FALSE,
    axis.cex = .8, lab.cex = 1
  )
  dens <- TernaryDensity(sub, resolution = 20)
  ColourTernary(dens, spectrum = mako(256))
  TernaryLines(list(c(0, 1, 0), c(1, 0, 0)), col = "gray70", lwd = 0.5)
}
dev.off()

h2_summary <- h2_draws |>
  group_by(coefficient, source) |>
  summarise(median = median(value), mean = mean(value),
            lower = quantile(value, 0.025), upper = quantile(value, 0.975),
            .groups = "drop")
print(h2_summary)



###########################################
# PLOT SPECIES-SPECIFIC ORDINAL RESPONSES #
###########################################

# Recover species-specific environmental slopes and display total responses
# and deviations alongside host phylogeny and trait similarity.

#-------- Extract posterior draws for species-level responses ------------

fit_ord <- fit_dirichlet_partition_within
sd_ord <- stan_dat_delta

draws_ord <- posterior::as_draws_df(
  fit_ord$draws()
)

# Identify fixed-effect coefficients from the model matrix
predictor_names <- colnames(sd_ord$X)

b_names <- paste0(
  "b[", seq_len(sd_ord$Kc), "]"
)

stopifnot(
  length(predictor_names) == length(b_names),
  all(b_names %in% names(draws_ord))
)

b_draws_ord <- as.matrix(
  draws_ord[, b_names, drop = FALSE]
)

colnames(b_draws_ord) <- predictor_names

# Confirm environmental predictors are present
pc_names <- paste0("delta_PC", 1:3)

stopifnot(
  all(pc_names %in% colnames(b_draws_ord))
)

#-------- Map model random-effect indices to species ------------

# Reconstruct observations retained by brms::standata
# so the posterior random effects match the correct species.

model_vars <- c(
  "shift_ord",
  "log_delta_t",
  "delta_PC1",
  "delta_PC2",
  "delta_PC3",
  "last_genus",
  "STUDY",
  "SITE",
  "COLONY_ID",
  "clean.spec",
  "clean.spec_trait",
  "clean.spec_resid"
)

model_dat <- phy.dat %>%
  select(all_of(model_vars)) %>%
  drop_na()

# Verify that row order matches the Stan data
if (
  nrow(model_dat) != sd_ord$N ||
  !identical(
    as.integer(model_dat$shift_ord),
    as.integer(sd_ord$Y)
  )
) {
  stop(
    "Model observations do not match standata. ",
    "Check filtering and row ordering before mapping species."
  )
}

# J_1, J_2, J_3 correspond to the species-level
# phylogenetic, trait, and residual random-effect groups.

species_map_ord <- model_dat %>%
  mutate(
    j1 = sd_ord$J_1,
    j2 = sd_ord$J_2,
    j3 = sd_ord$J_3
  ) %>%
  transmute(
    species = as.character(clean.spec),
    j1, j2, j3
  ) %>%
  distinct() %>%
  arrange(species)

# Each species must map uniquely to each RE group
stopifnot(
  !anyDuplicated(species_map_ord$species),
  !anyDuplicated(species_map_ord$j1),
  !anyDuplicated(species_map_ord$j2),
  !anyDuplicated(species_map_ord$j3)
)

#-------- Extract species-level random effects ------------

extract_re_ord <- function(block, coef_slot) {
  
  n <- sd_ord[[paste0("N_", block)]]
  
  cols <- paste0(
    "r_", block, "[",
    seq_len(n), ",", coef_slot, "]"
  )
  
  if (!all(cols %in% names(draws_ord))) {
    stop(
      "Missing posterior random-effect columns: ",
      paste(head(setdiff(cols, names(draws_ord))), collapse = ", ")
    )
  }
  
  as.matrix(
    draws_ord[, cols, drop = FALSE]
  )
}

# Coefficient slots:
# 1 = intercept
# 2 = delta_PC1
# 3 = delta_PC2
# 4 = delta_PC3

# Sum the phylogenetic, trait, and residual REs
re_by_slot_ord <- lapply(1:4, function(k) {
  
  extract_re_ord(1, k)[
    , species_map_ord$j1, drop = FALSE
  ] +
    
    extract_re_ord(2, k)[
      , species_map_ord$j2, drop = FALSE
    ] +
    
    extract_re_ord(3, k)[
      , species_map_ord$j3, drop = FALSE
    ]
  
})

#-------- Define the species-level posterior summary helper ------------

summarize_species_ord <- function(mat, coefficient) {
  
  stopifnot(
    ncol(mat) == nrow(species_map_ord)
  )
  
  p_positive <- colMeans(mat > 0)
  p_negative <- colMeans(mat < 0)
  
  tibble(
    species = species_map_ord$species,
    
    median = apply(mat, 2, median),
    lower = apply(mat, 2, quantile, probs = 0.025),
    upper = apply(mat, 2, quantile, probs = 0.975),
    
    p_positive = p_positive,
    
    # Signed posterior directional support:
    # -1 = all draws negative
    #  0 = equal support for both directions
    # +1 = all draws positive
    
    signed_probability = p_positive - p_negative,
    
    coefficient = coefficient
  )
}

#-------- Calculate total species-specific environmental slopes ------------

# Total slope =
# population-level fixed effect +
# phylogenetic RE + trait RE + residual RE

total_effects_ord <- lapply(1:3, function(k) {
  
  pc <- paste0("PC", k)
  
  fixed <- b_draws_ord[, paste0("delta_", pc)]
  
  total_draws <- sweep(
    re_by_slot_ord[[k + 1]],
    1,
    fixed,
    "+"
  )
  
  summarize_species_ord(
    total_draws,
    coefficient = pc
  )
  
}) %>%
  bind_rows()

#-------- Calculate species deviations from population slopes ------------

# Deviations include the sum of the three species-level
# random effects but exclude the population-level fixed slope.

deviation_effects_ord <- lapply(1:3, function(k) {
  
  pc <- paste0("PC", k)
  
  summarize_species_ord(
    re_by_slot_ord[[k + 1]],
    coefficient = pc
  )
  
}) %>%
  bind_rows()

#-------- Prepare the host phylogenetic tree ------------

species_to_plot <- species_map_ord$species

# Prune original coral phylogeny to modeled species
phylo_tree_ord <- ape::keep.tip(
  tree_final,
  species_to_plot
)

stopifnot(
  setequal(
    phylo_tree_ord$tip.label,
    species_to_plot
  )
)

# Ladderize for presentation
phylo_tree_ord <- ape::ladderize(
  phylo_tree_ord
)

# Obtain species ordering from the actual tree
tree_plot_base <- ggtree(
  phylo_tree_ord,
  layout = "rectangular"
)

species_order <- tree_plot_base$data %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label)

n_species <- length(species_order)

#-------- Prepare effect summaries for plotting ------------

prepare_effects_ord <- function(dat) {
  
  dat %>%
    mutate(
      y_position = match(
        species,
        species_order
      ),
      
      coefficient = factor(
        coefficient,
        levels = c("PC1", "PC2", "PC3")
      )
    )
}

total_effects_ord <- prepare_effects_ord(
  total_effects_ord
)

deviation_effects_ord <- prepare_effects_ord(
  deviation_effects_ord
)

stopifnot(
  !anyNA(total_effects_ord$y_position),
  !anyNA(deviation_effects_ord$y_position)
)

#-------- Define the environmental-response plotting function ------------

plot_species_ordinal_phylo <- function(dat, deviation = FALSE) {
  
  panel_names <- c("PC1", "PC2", "PC3")
  
  #-------- Environmental PC panels ------------
  
  effect_panels <- lapply(seq_along(panel_names), function(i) {
    
    pc <- panel_names[i]
    
    subdat <- dat %>%
      filter(coefficient == pc)
    
    p <- ggplot(
      subdat,
      aes(x = median, y = y_position)
    ) +
      
      # Zero-effect reference
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "gray60",
        linewidth = 0.35
      ) +
      
      # 95% posterior credible intervals
      geom_segment(
        aes(
          x = lower,
          xend = upper,
          yend = y_position,
          color = signed_probability
        ),
        linewidth = 0.45
      ) +
      
      # Posterior medians
      geom_point(
        aes(color = signed_probability),
        size = 1.5
      ) +
      
      # Signed posterior directional support
      scale_color_gradient2(
        low = "#2166ac",
        mid = "gray85",
        high = "#b2182b",
        midpoint = 0,
        limits = c(-1, 1),
        breaks = c(-1, -0.5, 0, 0.5, 1),
        labels = c("-1", "-0.5", "0", "0.5", "1"),
        name = "Posterior \nsupport",
        guide = guide_colorbar(
          direction = "vertical",
          barheight = grid::unit(5, "cm"),
          barwidth = grid::unit(0.5, "cm")
        )
      ) +
      
      # Match species positions across panels
      scale_y_continuous(
        limits = c(0.5, n_species + 0.5),
        breaks = seq_len(n_species),
        labels = gsub("_", " ", species_order),
        expand = c(0, 0)
      ) +
      
      labs(
        title = paste0("\u0394", pc),
        x = NULL,
        y = NULL
      ) +
      
      theme_bw(base_size = 10) +
      
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 11
        ),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(5, 6, 5, 6)
      )
    
    # Species labels only on first panel
    if (i == 1) {
      
      p <- p +
        theme(
          axis.text.y = element_text(
            size = 7,
            face = "italic",
            color = "black"
          )
        )
      
    } else {
      
      p <- p +
        theme(
          axis.text.y = element_blank()
        )
    }
    
    p
  })
  
#-------- Draw the phylogenetic tree panel ------------
  
  phylo_plot <- ggtree(
    phylo_tree_ord,
    layout = "rectangular",
    branch.length = "branch.length",
    size = 0.35
  ) +
    
    # Match vertical species coordinates
    scale_y_continuous(
      limits = c(0.5, n_species + 0.5),
      expand = c(0, 0)
    ) +
    
    # Root left, tips right
    # No scale_x_reverse()
    
    theme_void() +
    
    theme(
      plot.margin = margin(5, 0, 5, 6)
    )
  
#-------- Combine tree and response panels ------------
  
  combined_plot <- wrap_plots(
    c(list(phylo_plot), effect_panels),
    nrow = 1,
    widths = c(1.2, 3.2, 2.5, 2.5),
    guides = "collect"
  ) &
    theme(
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
  
  #-------- Shared x-axis title ------------
  
  x_label <- if (deviation) {
    "Deviation from population-level slope (cumulative logit scale)"
  } else {
    "Species-specific environmental response (cumulative logit scale)"
  }
  
  final_plot <- combined_plot /
    
    wrap_elements(
      full = grid::textGrob(
        x_label,
        gp = grid::gpar(fontsize = 11)
      )
    ) +
    
    plot_layout(
      heights = c(1, 0.06)
    )
  
  return(final_plot)
}

#-------- Generate phylogeny-ordered response figures ------------

# Main: total species-specific responses
fig_total_ord_phylo <- plot_species_ordinal_phylo(
  total_effects_ord,
  deviation = FALSE
)

# Supplementary: deviations from population-level slopes
fig_deviation_ord_phylo <- plot_species_ordinal_phylo(
  deviation_effects_ord,
  deviation = TRUE
)

print(fig_total_ord_phylo)
print(fig_deviation_ord_phylo)

#-------- Save phylogeny-ordered response figures ------------


ggsave(
  "results/03a_Bayesian_GenusOrdinal/Species_Total_Slopes_Phylogeny_ORDINAL.svg",
  plot = fig_total_ord_phylo,
  width = 6.5,
  height = 8
)

ggsave(
  "results/03a_Bayesian_GenusOrdinal/Species_Deviation_Slopes_Phylogeny_ORDINAL.svg",
  plot = fig_deviation_ord_phylo,
  width = 12,
  height = 10
)


#-------- Cluster modeled species by trait similarity ------------

species_to_plot <- species_map_ord$species

trait_mat <- T_mat[
  species_to_plot,
  species_to_plot,
  drop = FALSE
]

trait_hc <- hclust(
  as.dist(1 - trait_mat),
  method = "average"
)

species_order_trait <- trait_hc$labels[trait_hc$order]
n_species_trait <- length(species_order_trait)

dend_data_trait <- ggdendro::dendro_data(
  as.dendrogram(trait_hc),
  type = "rectangle"
)

# Confirm species ordering
stopifnot(
  setequal(species_order_trait, species_to_plot)
)



#-------- Prepare trait-ordered effect summaries ------------

# Use existing posterior summaries from the ordinal model.
# Recalculate positions according to trait similarity.

prepare_effects_trait <- function(dat) {
  
  dat %>%
    mutate(
      y_position = match(
        as.character(species),
        species_order_trait
      ),
      
      coefficient = factor(
        as.character(coefficient),
        levels = c("PC1", "PC2", "PC3")
      )
    )
}

total_effects_ord_trait <- prepare_effects_trait(
  total_effects_ord
)

deviation_effects_ord_trait <- prepare_effects_trait(
  deviation_effects_ord
)

stopifnot(
  !anyNA(total_effects_ord_trait$y_position),
  !anyNA(deviation_effects_ord_trait$y_position)
)

#-------- Define the trait-ordered response plotting function ------------

plot_species_ordinal_trait <- function(dat, deviation = FALSE) {
  
  panel_names <- c("PC1", "PC2", "PC3")
  
#-------- Environmental PC panels ------------
  
  effect_panels <- lapply(seq_along(panel_names), function(i) {
    
    pc <- panel_names[i]
    
    subdat <- dat %>%
      filter(coefficient == pc)
    
    p <- ggplot(
      subdat,
      aes(x = median, y = y_position)
    ) +
      
      # Zero-effect reference
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "gray60",
        linewidth = 0.35
      ) +
      
      # 95% posterior credible intervals
      geom_segment(
        aes(
          x = lower,
          xend = upper,
          yend = y_position,
          color = signed_probability
        ),
        linewidth = 0.45
      ) +
      
      # Posterior medians
      geom_point(
        aes(color = signed_probability),
        size = 1.5
      ) +
      
      # Signed posterior directional support
      scale_color_gradient2(
        low = "#2166ac",
        mid = "gray85",
        high = "#b2182b",
        midpoint = 0,
        limits = c(-1, 1),
        breaks = c(-1, -0.5, 0, 0.5, 1),
        labels = c("-1", "-0.5", "0", "0.5", "1"),
        name = "Posterior \nsupport",
        guide = guide_colorbar(
          direction = "vertical",
          barheight = grid::unit(5, "cm"),
          barwidth = grid::unit(0.5, "cm")
        )
      ) +
      
      # Species positions ordered by trait similarity
      scale_y_continuous(
        limits = c(0.5, n_species_trait + 0.5),
        breaks = seq_len(n_species_trait),
        labels = gsub("_", " ", species_order_trait),
        expand = c(0, 0)
      ) +
      
      labs(
        title = paste0("\u0394", pc),
        x = NULL,
        y = NULL
      ) +
      
      theme_bw(base_size = 10) +
      
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 11
        ),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(5, 6, 5, 6)
      )
    
    # Species labels only on first panel
    if (i == 1) {
      
      p <- p +
        theme(
          axis.text.y = element_text(
            size = 7,
            face = "italic",
            color = "black"
          )
        )
      
    } else {
      
      p <- p +
        theme(
          axis.text.y = element_blank()
        )
    }
    
    p
  })
  
#-------- Draw the trait dendrogram ------------
  
  dend_segments <- ggdendro::segment(dend_data_trait)
  
  dend_plot <- ggplot(dend_segments) +
    
    geom_segment(
      aes(
        x = -y,
        xend = -yend,
        y = x,
        yend = xend
      ),
      linewidth = 0.35,
      color = "gray35"
    ) +
    
    scale_y_continuous(
      limits = c(0.5, n_species_trait + 0.5),
      expand = c(0, 0)
    ) +
    
    labs(
      x = NULL,
      y = NULL
    ) +
    
    theme_void() +
    
    theme(
      plot.margin = margin(5, 0, 5, 6)
    )
  
#-------- Combine trait dendrogram and response panels ------------
  
  combined_plot <- wrap_plots(
    c(list(dend_plot), effect_panels),
    nrow = 1,
    widths = c(1.2, 3.2, 2.5, 2.5),
    guides = "collect"
  ) &
    theme(
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
  
#-------- Shared x-axis title ------------
  
  x_label <- if (deviation) {
    "Deviation from population-level slope (cumulative logit scale)"
  } else {
    "Species-specific environmental response (cumulative logit scale)"
  }
  
  final_plot <- combined_plot /
    
    wrap_elements(
      full = grid::textGrob(
        x_label,
        gp = grid::gpar(fontsize = 11)
      )
    ) +
    
    plot_layout(
      heights = c(1, 0.06)
    )
  
  return(final_plot)
}

#-------- Generate trait-ordered response figures ------------

# Main figure: total species-specific slopes
fig_total_ord_trait <- plot_species_ordinal_trait(
  total_effects_ord_trait,
  deviation = FALSE
)

# Supplementary figure: deviations from population-level slopes
fig_deviation_ord_trait <- plot_species_ordinal_trait(
  deviation_effects_ord_trait,
  deviation = TRUE
)

print(fig_total_ord_trait)
print(fig_deviation_ord_trait)

#-------- Save trait-ordered response figures ------------

ggsave(
  "results/03a_Bayesian_GenusOrdinal/Species_Total_Slopes_TraitTree_ORDINAL.svg",
  plot = fig_total_ord_trait,
  width = 6.5,
  height = 8
)

ggsave(
  "results/03a_Bayesian_GenusOrdinal/Species_Deviation_Slopes_TraitTree_ORDINAL.svg",
  plot = fig_deviation_ord_trait,
  width = 12,
  height = 10
)

