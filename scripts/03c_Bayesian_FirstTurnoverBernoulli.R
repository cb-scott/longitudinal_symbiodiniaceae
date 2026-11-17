###########################################
# BAYESIAN FIRST-TURNOVER BERNOULLI MODEL #
###########################################
# Builds the genus-level analysis data, restricts each colony to observations
# through its first modeled turnover, and fits the Bernoulli sensitivity model
# used for comparison with annual transition probabilities from the Markov
# analysis. Species-level intercepts and environmental slopes are partitioned
# among phylogenetic, trait-associated, and residual components. The existing
# fitted model is loaded by default; sampling occurs only when that file is absent.


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

# Create the script-specific results directory before writing outputs.
dir.create("results/03c_Bayesian_FirstTurnoverBernoulli", recursive = TRUE, showWarnings = FALSE)

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


######################################################
# FIT THE FIRST-TURNOVER BERNOULLI SENSITIVITY MODEL #
######################################################
# Restrict each colony to observations through its first modeled turnover,
# fit the Bernoulli variance-partition model, and save it for comparison with
# the Markov analysis.
#-------- Prepare the first-turnover Bernoulli model data ------------
# Reuse the hierarchical structure above with a binary first-turnover outcome.

first.shift = phy.dat %>% group_by(COLONY_ID) %>%
  mutate(shifted = ifelse(shift_ord !=0, 1, 0)) %>%
  arrange(DATE) %>%
  filter(lag(cumsum(shifted), default = 0) < 1)

stan_dat.ber.first <- standata(
  shifted ~ log_delta_t + delta_PC1 + delta_PC2 + delta_PC3 + last_genus +
    (1 | STUDY/SITE/COLONY_ID) +
    (1 + delta_PC1 + delta_PC2 + delta_PC3 | gr(clean.spec,       cov = A)) +
    (1 + delta_PC1 + delta_PC2 + delta_PC3 | gr(clean.spec_trait, cov = T_mat)) +
    (1 + delta_PC1 + delta_PC2 + delta_PC3 | clean.spec_resid),
  data   = first.shift,
  data2  = list(A = A, T_mat = T_mat),
  family = bernoulli()
)
stan_dat.ber.first$prior_only <- 0L

#-------- Load the saved model or sample it when absent ------------
model_file <- "data/stan_dirichlet_explicit_bernoulli_first.RDS"

if (file.exists(model_file)) {
  fit_dirichlet_partition_ber_first <- readRDS(model_file)
} else {
  mod.ber <- cmdstan_model("scripts/coral_dirichlet_partition_bernoulli_site.stan")

  fit_dirichlet_partition_ber_first <- mod.ber$sample(
    data            = stan_dat.ber.first,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 1000,
    iter_sampling   = 1000,
    seed            = 42,
    adapt_delta     = 0.95,
    refresh         = 50
  )

  fit_dirichlet_partition_ber_first$save_object(file = model_file)
}




######################################################
# SUMMARIZE AND PLOT THE BERNOULLI SENSITIVITY MODEL #
######################################################


fit <- fit_dirichlet_partition_ber_first
sd_data <- stan_dat.ber.first
draws_df <- as_draws_df(fit$draws())

# Fixed-effect coefficient names and predictor labels
b_names <- paste0("b[", seq_len(sd_data$Kc), "]")

x_names <- colnames(sd_data$X)[-1]

if (is.null(x_names)) {
  stop("stan_data$X needs predictor column names.")
}

stopifnot(
  length(x_names) == length(b_names),
  all(c(b_names, "Intercept") %in% names(draws_df))
)

names(x_names) <- b_names

# Posterior matrices
b_draws <- as.matrix(draws_df[, b_names, drop = FALSE])
colnames(b_draws) <- unname(x_names)

intercept_draws <- draws_df$Intercept

model_means <- colMeans(sd_data$X[, -1, drop = FALSE])
names(model_means) <- colnames(sd_data$X)[-1]

# Posterior summary function
summarize_params <- function(vars) {
  draws_df %>%
    select(all_of(vars)) %>%
    pivot_longer(everything(),
                 names_to = "variable",
                 values_to = "value") %>%
    group_by(variable) %>%
    summarise(
      median = median(value),
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975),
      .groups = "drop"
    )
}

pretty_predictor <- function(x) {
  x <- ifelse(x == "log_delta_t", "Log elapsed time", x)
  x <- sub("^delta_", "\u0394", x)
  x
}


#-------- Summarize and plot fixed effects ------------

fe_draws <- draws_df %>%
  select(all_of(b_names)) %>%
  pivot_longer(
    everything(),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  mutate(
    predictor = pretty_predictor(unname(x_names[parameter])),
    predictor = factor(
      predictor,
      levels = rev(unique(pretty_predictor(unname(x_names))))
    )
  )

fe_summary <- summarize_params(b_names) %>%
  mutate(term = pretty_predictor(unname(x_names[variable])))

print(fe_summary)

fe.p <- ggplot(fe_draws, aes(x = value, y = predictor)) +
  stat_halfeye(
    .width = c(0.5, 0.95),
    point_interval = median_qi,
    fill = "#81C3D7",
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray50"
  ) +
  labs(
    x = "Effect on log-odds of first turnover",
    y = NULL
  ) +
  theme_bw(base_size = 10)

print(fe.p)

ggsave(
  "results/03c_Bayesian_FirstTurnoverBernoulli/Fixed_Effects.tiff",
  fe.p,
  width = 4.5,
  height = 4,
  dpi = 600,
  device = ragg::agg_tiff
)

write_xlsx(
  fe_summary,
  "results/03c_Bayesian_FirstTurnoverBernoulli/Fixed_Effects.xlsx"
)



#-------- Summarize random-effect standard deviations ------------

coef_names <- c("Intercept", "PC1", "PC2", "PC3")

sd_vars <- c(
  unlist(lapply(1:3, function(i) {
    paste0("sd_", i, "[", 1:4, "]")
  })),
  paste0("sd_", 4:6, "[1]")
)

stopifnot(all(sd_vars %in% names(draws_df)))

sd_lookup <- tibble(
  variable = sd_vars,
  source = c(
    rep("Phylogenetic", 4),
    rep("Idiosyncratic", 4),
    rep("Trait", 4),
    "Study", "Site", "Colony"
  ),
  coefficient = c(
    rep(coef_names, 3),
    rep("Intercept", 3)
  )
) %>%
  mutate(
    term = paste(source, coefficient, sep = " - ")
  )

re_summary <- summarize_params(sd_vars) %>%
  left_join(sd_lookup, by = "variable")

print(re_summary)

re_draws <- draws_df %>%
  select(all_of(sd_vars)) %>%
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "value"
  ) %>%
  left_join(sd_lookup, by = "variable")

re.p <- ggplot(
  re_draws,
  aes(x = value, y = term)
) +
  stat_halfeye(
    .width = c(0.5, 0.95),
    point_interval = median_qi,
    normalize = "panels",
    alpha = 0.85
  ) +
  labs(
    x = "Random-effect SD (log-odds scale)",
    y = NULL
  ) +
  theme_bw(base_size = 9)

print(re.p)

ggsave(
  "results/03c_Bayesian_FirstTurnoverBernoulli/Random_Effects.tiff",
  re.p,
  width = 5,
  height = 6,
  dpi = 600,
  device = ragg::agg_tiff
)

write_xlsx(
  re_summary,
  "results/03c_Bayesian_FirstTurnoverBernoulli/Random_Effects.xlsx"
)


#-------- Calculate the Bernoulli variance partition ------------

coef_labels <- c(
  "Baseline turnover",
  "PC1 sensitivity",
  "PC2 sensitivity",
  "PC3 sensitivity"
)

source_labels <- c(
  "Phylogenetic",
  "Idiosyncratic",
  "Trait"
)

h2_vars <- grep(
  "^h2\\[[0-9]+,[0-9]+\\]$",
  names(draws_df),
  value = TRUE
)

stopifnot(length(h2_vars) == 12L)

h2_draws <- draws_df %>%
  select(.draw, all_of(h2_vars)) %>%
  pivot_longer(
    all_of(h2_vars),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  mutate(
    coef_idx = as.integer(
      str_match(parameter, "^h2\\[([0-9]+),([0-9]+)\\]$")[, 2]
    ),
    source_idx = as.integer(
      str_match(parameter, "^h2\\[([0-9]+),([0-9]+)\\]$")[, 3]
    ),
    coefficient = factor(
      coef_idx,
      levels = 1:4,
      labels = coef_labels
    ),
    source = factor(
      source_idx,
      levels = 1:3,
      labels = source_labels
    )
  )

h2_summary <- h2_draws %>%
  group_by(coefficient, source) %>%
  summarise(
    median = median(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    .groups = "drop"
  )

print(h2_summary)

# Preserve joint posterior draws for comparisons
h2_wide <- h2_draws %>%
  select(.draw, coefficient, source, value) %>%
  pivot_wider(
    id_cols = c(.draw, coefficient),
    names_from = source,
    values_from = value
  )

partition_probabilities <- h2_wide %>%
  group_by(coefficient) %>%
  summarise(
    p_phylo_above_third = mean(Phylogenetic > 1/3),
    p_trait_above_third = mean(Trait > 1/3),
    p_resid_above_third = mean(Idiosyncratic > 1/3),
    
    p_phylo_leads = mean(
      Phylogenetic > Trait &
        Phylogenetic > Idiosyncratic
    ),
    p_trait_leads = mean(
      Trait > Phylogenetic &
        Trait > Idiosyncratic
    ),
    p_resid_leads = mean(
      Idiosyncratic > Trait &
        Idiosyncratic > Phylogenetic
    ),
    .groups = "drop"
  )

print(partition_probabilities)

# Posterior distributions
h2.p <- ggplot(
  h2_draws,
  aes(x = value, y = source, fill = source)
) +
  stat_halfeye(
    .width = c(0.5, 0.95),
    point_interval = median_qi,
    normalize = "panels",
    alpha = 0.85
  ) +
  facet_wrap(~coefficient, nrow = 1) +
  scale_x_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  scale_fill_manual(
    values = c(
      "Phylogenetic" = "#2166ac",
      "Trait" = "#d6604d",
      "Idiosyncratic" = "#4dac26"
    )
  ) +
  labs(
    x = "Fraction of species-level latent variance",
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

print(h2.p)

ggsave(
  "results/03c_Bayesian_FirstTurnoverBernoulli/Variance_Partition.tiff",
  h2.p,
  width = 9,
  height = 3,
  dpi = 600,
  device = ragg::agg_tiff
)



#-------- Plot the Bernoulli variance partition in ternary space ------------

if (
  requireNamespace("Ternary", quietly = TRUE) &&
  requireNamespace("viridisLite", quietly = TRUE)
) {
  
  png(
    "results/03c_Bayesian_FirstTurnoverBernoulli/Variance_Partition_Ternary.png",
    width = 8,
    height = 2.5,
    units = "in",
    res = 500
  )
  
  oldpar <- par(
    mfrow = c(1, 4),
    mar = c(1, 1, 2, 1),
    cex = 8/12
  )
  
  for (cf in levels(h2_wide$coefficient)) {
    
    sub <- h2_wide %>%
      filter(coefficient == cf) %>%
      select(Trait, Phylogenetic, Idiosyncratic) %>%
      as.matrix()
    
    Ternary::TernaryPlot(
      alab = "Trait \u2192",
      blab = "Phylo. \u2192",
      clab = "\u2190 Idiosyn.",
      main = c(
        "Baseline", "\u0394PC1", "\u0394PC2", "\u0394PC3"
      )[match(cf, coef_labels)],
      grid.minor.lines = 0,
      axis.labels = FALSE,
      axis.tick = FALSE
    )
    
    dens <- Ternary::TernaryDensity(
      sub,
      resolution = 20
    )
    
    Ternary::ColourTernary(
      dens,
      spectrum = viridisLite::mako(256)
    )
  }
  
  par(oldpar)
  dev.off()
}


write_xlsx(
  list(
    Partition = h2_summary,
    Probabilities = partition_probabilities
  ),
  "results/03c_Bayesian_FirstTurnoverBernoulli/Variance_Partition.xlsx"
)

# Total species-level variances
total_var_vars <- grep(
  "^total_var\\[[0-9]+\\]$",
  names(draws_df),
  value = TRUE
)

if (length(total_var_vars) == 4) {
  total_var_summary <- summarize_params(total_var_vars) %>%
    mutate(
      coefficient = factor(
        as.integer(str_extract(variable, "[0-9]+")),
        levels = 1:4,
        labels = coef_labels
      )
    )
  print(total_var_summary)
}



