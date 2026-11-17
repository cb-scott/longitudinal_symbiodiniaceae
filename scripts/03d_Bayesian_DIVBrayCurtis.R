####################################################
# BAYESIAN MODEL OF DIV-LEVEL BRAY-CURTIS DISTANCE #
####################################################
# Fits the main DIV-level beta model, partitions species-level variance among
# phylogenetic, trait, and residual components, and summarizes species-specific
# responses to environmental change.

#-------- Load required packages ------------
library(tidyverse)
library(brms)
library(ggtree)
library(patchwork)
library(ape)
library(phytools)
library(posterior)
library(Matrix)
library(ggdist)

# Create the script-specific results directory before writing outputs.
dir.create("results/03d_Bayesian_DIVBrayCurtis", recursive = TRUE, showWarnings = FALSE)

source("scripts/helpers.R")

######################################
# PREPARE THE BRAY-CURTIS MODEL DATA #
######################################
# Join sequential DIV distances to environmental PCs and construct the host,
# site, and colony identifiers required by the hierarchical model.
# Join environmental PCs to sequential Bray-Curtis distances. Lee observations
# use the standardized coordinates assigned during environmental extraction.
load("data/ColonyDistances.Rdata")
load("data/EnvData_PCs.Rdata")

tomerge.env = env.data.final %>%
  select(lat, lon, DATE, PC1:PC3) %>%
  mutate(lat = round(lat, 8),
         lon = round(lon, 8),
         DATE = as.Date(DATE)) %>% distinct()

all.colonydists = all.colonydists %>%
  mutate(LATITUDE = ifelse(STUDY == "LEE", 5.610409515086883, LATITUDE),
         LONGITUDE = ifelse(STUDY == "LEE", 103.06064931859011, LONGITUDE),
         LATITUDE = round(LATITUDE, 8),
         LONGITUDE = round(LONGITUDE, 8),
         FROM_DATE = as.Date(FROM_DATE),
         TO_DATE = as.Date(TO_DATE))

c.dists = all.colonydists %>% 
  left_join(tomerge.env, by = c("FROM_DATE" = "DATE", "LATITUDE" = "lat", "LONGITUDE" = "lon")) %>%
  left_join(tomerge.env, by = c("TO_DATE" = "DATE", "LATITUDE" = "lat", "LONGITUDE" = "lon")) %>%
  mutate(delta_PC1 = PC1.y - PC1.x,
         delta_PC2 = PC2.y - PC2.x,
         delta_PC3 = PC3.y - PC3.x,
         delta_t = as.numeric(TO_DATE - FROM_DATE),
         log_delta_t = log(delta_t)) %>%
  separate(SPECIES, into = c("HOST_GEN", "HOST_SPEC"), sep = "_", remove = FALSE)

#-------- Standardize host species names ------------
c.dists$HOST_SPEC <- ifelse(c.dists$HOST_SPEC == "INTERMEDIA(NOBILIS)", "INTERMEDIA", c.dists$HOST_SPEC)
c.dists$HOST_GEN  <- case_when(
  c.dists$HOST_SPEC == "HYSTRIX"                     ~ "SERIATOPORA",
  c.dists$HOST_GEN  == "GONIASTREA (FAVIA)"          ~ "GONIASTREA",
  c.dists$HOST_SPEC == "PALIFERA"                    ~ "ISOPORA",
  .default = c.dists$HOST_GEN
)

c.dists <- c.dists %>%
  mutate(clean.spec = paste(str_to_title(HOST_GEN), tolower(HOST_SPEC), sep = "_")) %>%
  mutate(clean.spec = case_when(
    clean.spec == "Platygyra_daedalea"  ~ "Platygyra_daedalea",
    clean.spec == "Porites_sp."         ~ "Porites_lobata",
    clean.spec == "Pachyseris_rugosa"   ~ "Pachyseris_speciosa",
    clean.spec == "Dipastraea_pallida" ~ "Dipsastraea_pallida",
    clean.spec == "Porites_harrisoni" ~ "Porites_lobata",
    clean.spec == "Cyphastrea_microphtalma" ~ "Cyphastrea_microphthalma",
    clean.spec == "Diploria_labrynthiformis" ~ "Diploria_labyrinthiformis",
    clean.spec == "Pocillopora_spp" ~ "Pocillopora_capitata",
    .default = clean.spec
  )) %>%
  filter(clean.spec != "Acropora_sp.") %>%
  filter(clean.spec != "Dipsastraea_pallida") #no trait data for dipsastraea

dat = c.dists
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

# The DIV dataset contains fewer host species than the genus-level dataset, so
# only the missing taxa listed above require manual grafting.

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

phy.dat <- phy.dat %>% mutate(SITE = paste(LONGITUDE, LATITUDE, sep = "_")) %>%
  filter(clean.spec != "Porites_porites")


max_below_one <- max(phy.dat$DIST[phy.dat$DIST < 1], na.rm = TRUE)

phy.dat <- phy.dat %>%
  mutate(
    DIST_beta = if_else(
      DIST == 1,
      (1 + max_below_one) / 2,
      DIST
    )
  ) %>% filter(is.na(DIST_beta) == F)

phy.dat$last_genus <- relevel(
  factor(phy.dat$last_genus),
  ref = "C"
)

######################################
# FIT THE BRAY-CURTIS DISTANCE MODEL #
######################################
# Transform boundary distances for the beta likelihood, prepare Stan data,
# sample the custom distance model, and save the fitted object.
# Use Bray-Curtis distance as the response for the primary DIV-level model.
stan_dat_distances <- standata(
  DIST_beta ~ log_delta_t + delta_PC1 + delta_PC2 + delta_PC3 + last_genus + 
    (1 | STUDY/SITE/COLONY_ID) +
    (1 + delta_PC1 + delta_PC2 + delta_PC3 | gr(clean.spec, cov = A)) +
    (1 + delta_PC1 + delta_PC2 + delta_PC3 | gr(clean.spec_trait, cov = T_mat)) +
    (1 + delta_PC1 + delta_PC2 + delta_PC3 | clean.spec_resid),
  data = phy.dat,
  data2 = list(A = A, T_mat = T_mat),
  family = Beta(link = "logit")
)


mod_dists <- cmdstanr::cmdstan_model(
  "scripts/coral_dirichlet_partition_distance_beta.stan"#,
  #force_recompile = TRUE
)

stan_dat_distances$prior_only <- 0L

fit_dirichlet_partition_distances <- mod_dists$sample(
  data            = stan_dat_distances,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  seed            = 42,
  adapt_delta     = 0.95,
  refresh         = 50
)

fit_dirichlet_partition_distances$save_object(
  file = "data/stan_dirichlet_distances_bray.RDS"
)

fit_dirichlet_partition_distances <- readRDS("data/stan_dirichlet_distances_bray.RDS")

# Summarize the fitted beta-distance model and generate posterior diagnostics.
# Posterior visualization for the continuous beta DIST model
# Requires: fit_dirichlet_partition_distances, stan_dat_distances, phy.dat,
#           tidyverse, posterior, ggdist, writexl, Ternary, viridisLite
# Phylogeny figures use the tree_final and trait_tree objects prepared above.
# Optional trait association: trait.sub


#############################################
# SUMMARIZE THE BRAY-CURTIS MODEL POSTERIOR #
#############################################
# Label posterior parameters, summarize effects, partition species variance, and calculate expected community distances.

#-------- Extract posterior draws and center the predictor matrix ------------

fit <- fit_dirichlet_partition_distances  # updated fit including last_genus
sd_data <- stan_dat_distances
draws_df <- posterior::as_draws_df(fit$draws())

# Kc changes when categorical genus predictors are added; K (= intercept + Kc) follows.
stopifnot(sd_data$Kc == ncol(sd_data$X) - 1L,
          sd_data$K == ncol(sd_data$X),
          all(c("log_delta_t", "delta_PC1", "delta_PC2", "delta_PC3") %in%
                colnames(sd_data$X)),
          any(startsWith(colnames(sd_data$X), "last_genus")))
predictor_labels <- colnames(sd_data$X)[-1]
if (is.null(predictor_labels)) {
  stop("stan_dat_distances$X needs column names to label fixed effects")
}
b_names <- paste0("b[", seq_len(sd_data$Kc), "]")
stopifnot(all(c(b_names, "Intercept") %in% names(draws_df)))
names(predictor_labels) <- b_names

Xc <- sweep(sd_data$X[, -1, drop = FALSE], 2,
            colMeans(sd_data$X[, -1, drop = FALSE]), "-")
y_obs <- sd_data$Y
b_draws <- as.matrix(draws_df[, b_names, drop = FALSE])
colnames(b_draws) <- unname(predictor_labels)
intercept_draws <- draws_df$Intercept
# b_Intercept, if needed, is the uncentered intercept in generated quantities.

#-------- Summarize and plot fixed effects ------------
pretty_predictor <- function(x) {
  x <- sub("^delta_", "\u0394", x)
  x <- ifelse(x == "log_delta_t", "Log elapsed time", x)
  x
}

# Extract fixed-effect posteriors
fes <- draws_df |>
  select(all_of(b_names)) |>
  pivot_longer(
    everything(),
    names_to = "param",
    values_to = "value"
  ) |>
  mutate(
    predictor = unname(predictor_labels[param]),
    category = case_when(
      grepl("delta_PC", predictor) ~ "Environmental",
      grepl("last_genus", predictor) ~ "Symbiont",
      TRUE ~ "other"
    ),
    predictor = ifelse(
      category == "Symbiont",
      gsub("^last_genus", "", predictor),
      predictor
    ),
    predictor = gsub("^delta_", "\u0394", predictor),
    is_mixed = ifelse(
      category == "Symbiont" & grepl("_", predictor),
      "Y", "N"
    )
  ) |>
  arrange(is_mixed, predictor)

# Match ordering of ordinal model
fes$predictor <- factor(
  fes$predictor,
  levels = rev(unique(fes$predictor))
)

# Plot fixed effects
fe.p <- fes |>
  filter(category != "other") |>
  filter(!(predictor %in% c("A_C_D", "B_C_D", "D_G"))) |>
  ggplot(aes(x = value, y = predictor)) +
  stat_halfeye(
    aes(fill = interaction(category, is_mixed)),
    .width = c(0.5, 0.95),
    point_interval = median_qi,
    alpha = 0.8, point_size = 1
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray50"
  ) +
  labs(
    x = "Logit of DIV distance", y = NULL,
    fill = "Codominance"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c("#c3d781", "#81C3D7", "#d7c081")
  ) +
  facet_wrap(
    ~category,
    nrow = 2,
    scales = "free_y"
  ) +
  xlim(c(-1.5, 2))

fe.p

ggsave(
  "results/03d_Bayesian_DIVBrayCurtis/FEs_Distance_Beta_BRAY_WITH_GENUS.svg",
  plot = fe.p,
  width = 3.25,
  height = 4
)


#-------- Summarize fixed and random effects ------------
summarize_params <- function(vars) {
  vars <- intersect(vars, names(draws_df))
  if (!length(vars)) stop("No matching draw columns")
  posterior::summarise_draws(
    posterior::as_draws_df(draws_df[, vars, drop = FALSE]),
    median, ~posterior::quantile2(.x, probs = c(0.025, 0.975))
  ) |>
    select(variable, median, q2.5, q97.5)
}

fixed_summary <- summarize_params(b_names) |>
  mutate(term = unname(predictor_labels[variable])) |>
  select(term, median, q2.5, q97.5)

intercept_summary <- summarize_params("Intercept") |>
  mutate(term = "Beta mean intercept (centered logit scale)") |>
  select(term, median, q2.5, q97.5)

precision_summary <- summarize_params("phi") |>
  mutate(term = "Beta precision (phi)") |>
  select(term, median, q2.5, q97.5)

fixed_effects_table <- bind_rows(fixed_summary, intercept_summary,
                                 precision_summary)
writexl::write_xlsx(fixed_effects_table, "results/03d_Bayesian_DIVBrayCurtis/FE_DiriBeta_BRAY_WITH_GENUS.xlsx")

coef_names <- c("Intercept", "PC1", "PC2", "PC3")
sd_labels_vec <- list(
  sd_1 = paste0("HOST.SPECIES_Phylo_", coef_names),
  sd_2 = paste0("HOST.SPECIES_Idio_", coef_names),
  sd_3 = paste0("HOST.SPECIES_Trait_", coef_names)
)
scalar_labels <- c(sd_4 = "STUDY_Intercept",
                   sd_5 = "STUDY:SITE_Intercept",
                   sd_6 = "STUDY:SITE:COLONY_ID_Intercept")

sd_vars <- c(unlist(lapply(1:3, function(i)
  paste0("sd_", i, "[", 1:4, "]"))), paste0("sd_", 4:6, "[1]"))
stopifnot(all(sd_vars %in% names(draws_df)))

sd_lookup <- tibble(
  variable = sd_vars,
  term = c(unlist(sd_labels_vec, use.names = FALSE), unname(scalar_labels))
)

random_summary <- summarize_params(sd_vars) |>
  left_join(sd_lookup, by = "variable") |>
  select(term, median, q2.5, q97.5)
writexl::write_xlsx(random_summary, "results/03d_Bayesian_DIVBrayCurtis/RE_DiriBeta_BRAY_WITH_GENUS.xlsx")

raw_draws_long <- draws_df |>
  select(.chain, .iteration, .draw, all_of(sd_vars)) |>
  pivot_longer(all_of(sd_vars), names_to = "variable", values_to = "value") |>
  left_join(sd_lookup, by = "variable")

re.p <- ggplot(raw_draws_long, aes(x = value, y = term)) +
  stat_halfeye(.width = c(0.5, 0.95), point_interval = median_qi,
               normalize = "panels", alpha = 0.85) +
  labs(x = "SD (logit expected community distance)", y = "Random effect") +
  theme_bw(base_size = 8)
print(re.p)

ggsave("results/03d_Bayesian_DIVBrayCurtis/Posterior_REs_Distance_BRAY_WITH_GENUS.svg", plot = re.p,
       width = 4, height = 6)

#-------- Calculate and plot the variance partition ------------
coef_labels <- c("Baseline dissimilarity", "PC1 sensitivity",
                 "PC2 sensitivity", "PC3 sensitivity")
source_labels <- c("Phylogenetic", "Idiosyncratic", "Trait")
h2_vars <- grep("^h2\\[[0-9]+,[0-9]+\\]$", names(draws_df), value = TRUE)
stopifnot(length(h2_vars) == 12L)

h2_draws <- draws_df |>
  select(.draw, all_of(h2_vars)) |>
  pivot_longer(all_of(h2_vars), names_to = "param", values_to = "value") |>
  mutate(
    coef_idx = as.integer(str_match(param, "^h2\\[([0-9]+),([0-9]+)\\]$")[, 2]),
    source_idx = as.integer(str_match(param, "^h2\\[([0-9]+),([0-9]+)\\]$")[, 3]),
    coefficient = factor(coef_idx, levels = 1:4, labels = coef_labels),
    source = factor(source_idx, levels = 1:3, labels = source_labels)
  )

# Use genuine .draw identifiers to preserve the joint posterior.
h2_wide <- h2_draws |>
  select(.draw, coefficient, source, value) |>
  pivot_wider(id_cols = c(.draw, coefficient),
              names_from = source, values_from = value)

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
                    "results/03d_Bayesian_DIVBrayCurtis/Variance_Partition_Probabilities_Beta_BRAY_WITH_GENUS.xlsx")


h2_summary <- h2_draws |>
  group_by(coefficient, source) |>
  summarise(median = median(value), mean = mean(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975), .groups = "drop")
print(h2_summary)
writexl::write_xlsx(h2_summary, "results/03d_Bayesian_DIVBrayCurtis/Variance_Partition_Beta_BRAY_WITH_GENUS.xlsx")

total_var_vars <- grep("^total_var\\[[0-9]+\\]$", names(draws_df), value = TRUE)
total_var_summary <- summarize_params(total_var_vars) |>
  mutate(coefficient = factor(as.integer(str_extract(variable, "[0-9]+")),
                              levels = 1:4, labels = coef_labels))
print(total_var_summary)

h2.p <- ggplot(h2_draws, aes(x = value, y = source, fill = source)) +
  stat_halfeye(.width = c(0.5, 0.95), point_interval = median_qi,
               normalize = "panels", alpha = 0.85) +
  facet_wrap(~coefficient, nrow = 1) +
  scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_fill_manual(values = c("Phylogenetic" = "#2166ac",
                               "Trait" = "#d6604d",
                               "Idiosyncratic" = "#4dac26")) +
  labs(x = "Fraction of species-level latent variance", y = NULL,
       title = "Variance decomposition of community dissimilarity") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))
print(h2.p)

ggsave("results/03d_Bayesian_DIVBrayCurtis/Variance_Partition_Beta_BRAY_WITH_GENUS.svg", plot = h2.p,
       width = 9, height = 3)



#-------- Compare trait-associated and phylogenetic variance ------------

trait_vs_phylo <- h2_wide %>%
  group_by(coefficient) %>%
  summarise(
    p_trait_gt_phylo = mean(Trait > Phylogenetic),
    
    median_difference = median(Trait - Phylogenetic),
    lower = quantile(Trait - Phylogenetic, 0.025),
    upper = quantile(Trait - Phylogenetic, 0.975),
    
    .groups = "drop"
  ) %>%
  mutate(
    # Symmetric Dirichlet prior implies prior odds = 1
    BF_trait_gt_phylo =
      p_trait_gt_phylo / (1 - p_trait_gt_phylo)
  )

print(trait_vs_phylo)

#-------- Plot posterior differences among variance components ------------

trait_phylo_draws <- h2_wide %>%
  mutate(
    difference = Trait - Phylogenetic
  )

trait_phylo.p <- ggplot(
  trait_phylo_draws,
  aes(x = difference, y = coefficient)
) +
  ggdist::stat_halfeye(
    .width = c(0.5, 0.95),
    point_interval = ggdist::median_qi,
    fill = "#81C3D7",
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray50"
  ) +
  labs(
    x = "Trait-associated − phylogenetic variance fraction",
    y = NULL
  ) +
  theme_bw(base_size = 10)

trait_phylo.p

ggsave(
  "results/03d_Bayesian_DIVBrayCurtis/Trait_vs_Phylo_Variance_DIV.svg",
  plot = trait_phylo.p,
  width = 5,
  height = 3
)




# Ternary plots, joint posterior draws, same palette as original analysis.
# Requires package Ternary (with TernaryPlot, TernaryDensity, ColourTernary).
if (requireNamespace("Ternary", quietly = TRUE) &&
    requireNamespace("viridisLite", quietly = TRUE)) {
  svg("results/03d_Bayesian_DIVBrayCurtis/Variance_Partitioned_Distance_BRAY_WITH_GENUS.svg",
      width = 6.5, height = 3)
  oldpar <- par(mfrow = c(1, 4), mar = c(1, 1, 2, 1), cex = 8/12)
  for (cf in levels(h2_wide$coefficient)) {
    sub <- h2_wide |>
      filter(coefficient == cf) |>
      select(Trait, Phylogenetic, Idiosyncratic) |>
      as.matrix()
    Ternary::TernaryPlot(
      alab = "Trait \u2192", blab = "Phylo. \u2192", clab = "\u2190 Idiosyn.",
      main = c("Baseline", "\u0394PC1", "\u0394PC2", "\u0394PC3")[match(cf, coef_labels)],
      grid.minor.lines = 0, lab.offset = 0.11,
      axis.labels = FALSE, axis.tick = FALSE, axis.cex = .8, lab.cex = 1
    )
    dens <- Ternary::TernaryDensity(sub, resolution = 20)
    Ternary::ColourTernary(dens, spectrum = viridisLite::mako(256))
    Ternary::TernaryLines(list(c(0, 1, 0), c(1, 0, 0)), col = "gray70", lwd = 0.5)
  }
  par(oldpar)
  dev.off()
}

#-------- Calculate species-level expected Bray-Curtis distance ------------
# This section replaces ordinal expected-switching functions. It aligns
# the actual species indices through standata J_1 / J_2 / J_3.
extract_re <- function(block, coef_slot) {
  n <- sd_data[[paste0("N_", block)]]
  cols <- paste0("r_", block, "[", seq_len(n), ",", coef_slot, "]")
  missing <- setdiff(cols, names(draws_df))
  if (length(missing)) stop("Missing columns: ", paste(head(missing, 3), collapse = ", "))
  as.matrix(draws_df[, cols, drop = FALSE])
}

# Reconstruct the exact observation rows retained by brms::standata.
# STOP on a mismatch rather than silently mislabeling species.
model_vars <- c("DIST_beta", "log_delta_t", "delta_PC1", "delta_PC2", "delta_PC3",
                "last_genus", "STUDY", "SITE", "COLONY_ID", "clean.spec",
                "clean.spec_resid", "clean.spec_trait")
model_dat <- phy.dat |>
  select(all_of(model_vars)) |>
  drop_na()
if (nrow(model_dat) != sd_data$N ||
    !isTRUE(all.equal(as.numeric(model_dat$DIST_beta), as.numeric(sd_data$Y),
                      tolerance = 1e-10, check.attributes = FALSE))) {
  stop("Model-row order differs from phy.dat complete cases. Reconstruct the retained rows before mapping J indices to species.")
}

species_map <- model_dat |>
  mutate(j1 = sd_data$J_1, j2 = sd_data$J_2, j3 = sd_data$J_3) |>
  transmute(species = as.character(clean.spec), j1, j2, j3) |>
  distinct() |>
  arrange(species)

if (anyDuplicated(species_map$species) ||
    anyDuplicated(species_map$j1) ||
    anyDuplicated(species_map$j2) ||
    anyDuplicated(species_map$j3)) {
  stop("Species/group indices are not one-to-one; inspect the three clean.spec variables")
}
print(species_map)

# Predict conditional mean with study, site and colony REs set to zero.
# This is an average-level prediction, not a posterior predictive draw from Beta.
# X is already centered by the Stan model; use the same means.
model_means <- colMeans(sd_data$X[, -1, drop = FALSE])
names(model_means) <- colnames(sd_data$X)[-1]

re_by_slot <- lapply(1:4, function(k) {
  extract_re(1, k)[, species_map$j1, drop = FALSE] +
    extract_re(2, k)[, species_map$j2, drop = FALSE] +
    extract_re(3, k)[, species_map$j3, drop = FALSE]
})

# A one-hot contrast for last_genus is included in X, with C as reference.
# New predictions hold the initial genus state at C unless specified otherwise.
genus_levels <- levels(droplevels(model_dat$last_genus))
stopifnot("C" %in% genus_levels)

predict_species_distance <- function(pc = "delta_PC1", value = 0,
                                     log_time = median(model_dat$log_delta_t, na.rm = TRUE),
                                     genus = "C") {
  if (!pc %in% c("delta_PC1", "delta_PC2", "delta_PC3")) stop("Unknown PC")
  if (!genus %in% genus_levels) stop("Unknown last_genus: ", genus)
  
  new_x <- setNames(rep(0, length(model_means)), names(model_means))
  new_x["log_delta_t"] <- log_time
  new_x[pc] <- value
  # Match treatment contrasts generated by brms::standata.
  # All last_genus columns are zero for reference genus C.
  genus_cols <- grep("^last_genus", names(model_means), value = TRUE)
  if (genus != "C") {
    active_col <- paste0("last_genus", genus)
    if (!active_col %in% genus_cols) stop("Missing dummy column: ", active_col)
    new_x[active_col] <- 1
  }
  
  # Same fixed-effect centering as Stan: Xc = X[,-1] - column means.
  fixed <- as.vector(b_draws %*% (new_x - model_means))
  nd <- nrow(b_draws)
  ns <- nrow(species_map)
  eta <- matrix(intercept_draws + fixed, nrow = nd, ncol = ns)
  eta <- eta + re_by_slot[[1]] +
    value * re_by_slot[[match(pc, c("delta_PC1", "delta_PC2", "delta_PC3")) + 1]]
  mu <- plogis(eta)
  tibble(species = species_map$species,
         median = apply(mu, 2, median),
         lower = apply(mu, 2, quantile, probs = 0.025),
         upper = apply(mu, 2, quantile, probs = 0.975),
         varied_pc = pc, pc_value = value, last_genus = genus)
}

# Equal and opposite environmental changes: 90th percentile of absolute values.
q90_abs_pc1 <- unname(quantile(abs(model_dat$delta_PC1), 0.90, na.rm = TRUE))
q90_abs_pc3 <- unname(quantile(abs(model_dat$delta_PC3), 0.90, na.rm = TRUE))

pred_pc1_summary <- bind_rows(
  predict_species_distance("delta_PC1", -q90_abs_pc1) |> mutate(level = "Low"),
  predict_species_distance("delta_PC1", 0) |> mutate(level = "Baseline"),
  predict_species_distance("delta_PC1", q90_abs_pc1) |> mutate(level = "High")
)
pred_pc3_summary <- bind_rows(
  predict_species_distance("delta_PC3", -q90_abs_pc3) |> mutate(level = "Low"),
  predict_species_distance("delta_PC3", 0) |> mutate(level = "Baseline"),
  predict_species_distance("delta_PC3", q90_abs_pc3) |> mutate(level = "High")
)

writexl::write_xlsx(list(PC1 = pred_pc1_summary, PC3 = pred_pc3_summary),
                    "results/03d_Bayesian_DIVBrayCurtis/Species_Distance_Predictions_BRAY_WITH_GENUS.xlsx")

# Basic visualization that does not require phylogenetic plotting helpers.
species.p <- bind_rows(pred_pc1_summary, pred_pc3_summary) |>
  ggplot(aes(x = pc_value, y = median, group = species)) +
  geom_line(alpha = 0.15, linewidth = 0.35) +
  geom_point(alpha = 0.3, size = 0.5) +
  facet_wrap(~varied_pc, scales = "free_x") +
  labs(x = "Environmental change (PC units)",
       y = "Expected community distance (conditional beta mean)") +
  theme_bw(base_size = 10)
print(species.p)

ggsave("results/03d_Bayesian_DIVBrayCurtis/Species_Distance_Responses_BRAY_WITH_GENUS.tiff", plot = species.p,
       width = 6, height = 3.5, dpi = 600, device = ragg::agg_tiff)

###############################################
# PLOT SPECIES-SPECIFIC BRAY-CURTIS RESPONSES #
###############################################
# Recover total species slopes and deviations, then display them in trait-similarity and phylogenetic order.

#-------- Cluster modeled species by trait similarity ------------

species_to_plot <- species_map$species

trait_mat <- T_mat[
  species_to_plot,
  species_to_plot,
  drop = FALSE
]

trait_hc <- hclust(
  as.dist(1 - trait_mat),
  method = "average"
)

species_order <- trait_hc$labels[trait_hc$order]
n_species <- length(species_order)

dend_data <- ggdendro::dendro_data(
  as.dendrogram(trait_hc),
  type = "rectangle"
)

#-------- Define the species-level posterior summary helper ------------

summarize_species_draws <- function(mat, coefficient) {
  
  stopifnot(ncol(mat) == nrow(species_map))
  
  p_positive <- colMeans(mat > 0)
  p_negative <- colMeans(mat < 0)
  
  tibble(
    species = species_map$species,
    
    median = apply(mat, 2, median),
    lower = apply(mat, 2, quantile, probs = 0.025),
    upper = apply(mat, 2, quantile, probs = 0.975),
    
    p_positive = p_positive,
    
    # Signed posterior directional support:
    # -1 = all posterior draws negative
    #  0 = equal support for positive and negative
    # +1 = all posterior draws positive
    signed_probability = p_positive - p_negative,
    
    coefficient = coefficient
  )
}

#-------- Calculate total species-specific environmental slopes ------------

# Total slope = population-level fixed effect +
# phylogenetic + trait + idiosyncratic random effects

total_effects <- lapply(1:3, function(k) {
  
  pc <- paste0("PC", k)
  
  fixed <- b_draws[, paste0("delta_", pc)]
  
  total_draws <- sweep(
    re_by_slot[[k + 1]],
    1,
    fixed,
    "+"
  )
  
  summarize_species_draws(
    total_draws,
    coefficient = pc
  )
  
}) %>%
  bind_rows()

#-------- Calculate species deviations from population slopes ------------

# Deviation = phylogenetic + trait + idiosyncratic random effects
# No population-level fixed effect included

deviation_effects <- lapply(1:3, function(k) {
  
  pc <- paste0("PC", k)
  
  summarize_species_draws(
    re_by_slot[[k + 1]],
    coefficient = pc
  )
  
}) %>%
  bind_rows()

#-------- Prepare trait-ordered effect summaries ------------

prepare_effects <- function(dat) {
  
  dat %>%
    mutate(
      y_position = match(species, species_order),
      
      coefficient = factor(
        coefficient,
        levels = c("PC1", "PC2", "PC3")
      )
    )
}

total_effects <- prepare_effects(total_effects)
deviation_effects <- prepare_effects(deviation_effects)

#-------- Define the trait-ordered response plotting function ------------

plot_species_with_trait_tree <- function(dat, deviation = FALSE) {
  
  panel_names <- c("PC1", "PC2", "PC3")
  
  #-------- Individual environmental PC panels ------------
  
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
      
      # Continuous signed posterior directional support
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
      
      # Identical species positions across panels
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
    
    # Species labels appear only on the first effect panel
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
  
  dend_segments <- ggdendro::segment(dend_data)
  
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
      limits = c(0.5, n_species + 0.5),
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
  
  #-------- One shared x-axis title ------------
  
  x_label <- if (deviation) {
    "Deviation from population-level slope (logit scale)"
  } else {
    "Species-specific environmental response (logit scale)"
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

# Main figure: total species-specific environmental responses
fig_total <- plot_species_with_trait_tree(
  total_effects,
  deviation = FALSE
)

# Supplementary figure: deviations from population-level slopes
fig_deviation <- plot_species_with_trait_tree(
  deviation_effects,
  deviation = TRUE
)

print(fig_total)
print(fig_deviation)

#-------- Save trait-ordered response figures ------------


ggsave(
  "results/03d_Bayesian_DIVBrayCurtis/Species_Total_Slopes_TraitTree_BRAY_WITH_GENUS.svg",
  plot = fig_total,
  width = 6.5,
  height = 6
)

ggsave(
  "results/03d_Bayesian_DIVBrayCurtis/Species_Deviation_Slopes_TraitTree_BRAY_WITH_GENUS.svg",
  plot = fig_deviation,
  width = 12,
  height = 10
)



#-------- Plot species responses ordered by coral phylogeny ------------


#-------- Prepare the coral phylogeny ------------

# Use the original phylogenetic tree rather than clustering covariance matrix A

species_to_plot <- as.character(species_map$species)

# Keep only modeled species
phylo_tree <- ape::keep.tip(
  tree_final,
  species_to_plot
)

# Check coverage
stopifnot(
  setequal(phylo_tree$tip.label, species_to_plot)
)

# Ladderize for clearer presentation
phylo_tree <- ape::ladderize(phylo_tree)

# Extract the vertical ordering of species from the plotted phylogeny
tree_plot_base <- ggtree(
  phylo_tree,
  layout = "rectangular"
)

tree_plot_data <- tree_plot_base$data

tip_data <- tree_plot_data %>%
  filter(isTip) %>%
  arrange(y)

species_order <- tip_data$label
n_species <- length(species_order)

#-------- Reapply the species-level posterior summary helper ------------

summarize_species_draws <- function(mat, coefficient) {
  
  stopifnot(ncol(mat) == nrow(species_map))
  
  p_positive <- colMeans(mat > 0)
  p_negative <- colMeans(mat < 0)
  
  tibble(
    species = species_map$species,
    
    median = apply(mat, 2, median),
    lower = apply(mat, 2, quantile, probs = 0.025),
    upper = apply(mat, 2, quantile, probs = 0.975),
    
    p_positive = p_positive,
    
    # Signed posterior directional support:
    # -1 = all posterior draws negative
    #  0 = equal posterior support for both directions
    # +1 = all posterior draws positive
    signed_probability = p_positive - p_negative,
    
    coefficient = coefficient
  )
}

#-------- Calculate phylogeny-ordered total species slopes ------------

total_effects <- lapply(1:3, function(k) {
  
  pc <- paste0("PC", k)
  
  fixed <- b_draws[, paste0("delta_", pc)]
  
  total_draws <- sweep(
    re_by_slot[[k + 1]],
    1,
    fixed,
    "+"
  )
  
  summarize_species_draws(
    total_draws,
    coefficient = pc
  )
  
}) %>%
  bind_rows()

#-------- Calculate phylogeny-ordered species deviations ------------

deviation_effects <- lapply(1:3, function(k) {
  
  pc <- paste0("PC", k)
  
  summarize_species_draws(
    re_by_slot[[k + 1]],
    coefficient = pc
  )
  
}) %>%
  bind_rows()

#-------- Prepare phylogeny-ordered effect summaries ------------

prepare_effects <- function(dat) {
  
  dat %>%
    mutate(
      y_position = match(as.character(species), species_order),
      
      coefficient = factor(
        coefficient,
        levels = c("PC1", "PC2", "PC3")
      )
    )
}

total_effects <- prepare_effects(total_effects)
deviation_effects <- prepare_effects(deviation_effects)

stopifnot(
  !anyNA(total_effects$y_position),
  !anyNA(deviation_effects$y_position)
)

#-------- Define the phylogeny-ordered plotting function ------------

plot_species_with_phylogeny <- function(dat, deviation = FALSE) {
  
  panel_names <- c("PC1", "PC2", "PC3")
  
  #-------- Individual environmental PC panels ------------
  
  effect_panels <- lapply(seq_along(panel_names), function(i) {
    
    pc <- panel_names[i]
    
    subdat <- dat %>%
      filter(coefficient == pc)
    
    p <- ggplot(
      subdat,
      aes(x = median, y = y_position)
    ) +
      
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "gray60",
        linewidth = 0.35
      ) +
      
      geom_segment(
        aes(
          x = lower,
          xend = upper,
          yend = y_position,
          color = signed_probability
        ),
        linewidth = 0.45
      ) +
      
      geom_point(
        aes(color = signed_probability),
        size = 1.5
      ) +
      
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
    phylo_tree,
    layout = "rectangular",
    branch.length = "branch.length",
    size = 0.35
  ) +
    
    # Match response panels vertically
    scale_y_continuous(
      limits = c(0.5, n_species + 0.5),
      expand = c(0, 0)
    ) +
    
    # Reverse x-axis so tips face the species labels
    #scale_x_reverse() +
    
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
    "Deviation from population-level slope (logit scale)"
  } else {
    "Species-specific environmental response (logit scale)"
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

fig_total_phylo <- plot_species_with_phylogeny(
  total_effects,
  deviation = FALSE
)

fig_deviation_phylo <- plot_species_with_phylogeny(
  deviation_effects,
  deviation = TRUE
)

print(fig_total_phylo)
print(fig_deviation_phylo)

#-------- Save phylogeny-ordered response figures ------------

ggsave(
  "results/03d_Bayesian_DIVBrayCurtis/Species_Total_Slopes_Phylogeny_BRAY_WITH_GENUS.svg",
  plot = fig_total_phylo,
  width = 6.5,
  height = 6
)

ggsave(
  "results/03d_Bayesian_DIVBrayCurtis/Species_Deviation_Slopes_Phylogeny_BRAY_WITH_GENUS.svg",
  plot = fig_deviation_phylo,
  width = 12,
  height = 10
)

