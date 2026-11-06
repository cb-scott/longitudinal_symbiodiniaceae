source("scripts/helpers.R")
library(tidyverse)
library(brms)
library(cmdstanr)
library(ggtree)
library(patchwork)
#library(rotl)
library(ape)
library(phytools)
library(phangorn)
library(posterior)
library(Matrix)
library(posterior)
library(bayesplot)
library(tidybayes)
library(ggdist)
library(Ternary)
library(viridisLite)


# ==============================================================================
# 1. DATA PROCESSING
# ==============================================================================
#need to merge this with the environmetnal data 
#lee and adig didn't have environmental dta
load("data/ColonyDistances.Rdata")
load("data/EnvData_PCs.Rdata")
tomerge.env = env.data.final %>% select(lat, lon, DATE, PC1:PC3)
c.dists = all.colonydists %>% 
  left_join(tomerge.env, by = c("FROM_DATE" = "DATE", "LATITUDE" = "lat", "LONGITUDE" = "lon")) %>%
  left_join(tomerge.env, by = c("TO_DATE" = "DATE", "LATITUDE" = "lat", "LONGITUDE" = "lon")) %>%
  mutate(delta_PC1 = PC1.y - PC1.x,
         delta_PC2 = PC2.y - PC2.x,
         delta_PC3 = PC3.y - PC3.x,
         delta_t = as.numeric(TO_DATE - FROM_DATE),
         log_delta_t = log(delta_t)) %>%
  separate(SPECIES, into = c("HOST_GEN", "HOST_SPEC"), sep = "_", remove = FALSE)

# ---- Species name standardisation ----
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
# ==============================================================================
# 2. PHYLOGENETIC TREE & COVARIANCE MATRIX
# ==============================================================================
# Source:Gault, Jordan A.; Bentlage, Bastian; Huang, Danwei; Kerr, Alexander M. (2021). Lineage-specific variation in the evolutionary stability of coral photosymbiosis. Science Advances. https://doi.org/10.1126/sciadv.abh4243
# Data:   https://datadryad.org/dataset/doi:10.5061/dryad.tdz08kq0p

trees     <- read.nexus("data/mtree.trees")
tree      <- maxCladeCred(trees)

# Inspect coverage
missing <- setdiff(unique(dat$clean.spec), tree$tip.label)
cat("Species missing from tree:", length(missing), "\n")
print(missing)

#less species are contained in the ITS data, so we don't need to graft as many nodes here. 

# ---- Manual tip grafts for missing species ----

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

# ---- Prune to study species and build covariance matrix ----
tree_final <- keep.tip(tree, unique(dat$clean.spec))
tree_final <- force.ultrametric(tree_final, method = "extend")
A          <- ape::vcv.phylo(tree_final, corr = TRUE)

plot(tree_final)
# ==============================================================================
# 3. TRAIT COVARIANCE MATRIX
# ==============================================================================
# McWilliam, Hoogenboom & Baird (2018) trait PCs
# Source: https://github.com/mikemcwilliam/biogeographical_disparity

load("data/McMillianTraitPCs.Rdata")

mcmillian.sc <- mcmillian.sc %>%
  mutate(species = gsub(" ", "_", species))
mcmillian.sc$species[mcmillian.sc$species == "Montastrea_cavernosa"]  <- "Montastraea_cavernosa"
mcmillian.sc$species[mcmillian.sc$species == "Pocillopora_eydouxi"]   <- "Pocillopora_grandis"
rownames(mcmillian.sc) <- mcmillian.sc$species
colnames(mcmillian.sc) <- gsub(" ", "", colnames(mcmillian.sc))

# ---- Assemble analysis dataset ----
phy.dat <- dat %>% filter(clean.spec %in% rownames(A))
phy.dat <- left_join(phy.dat, mcmillian.sc, by = c("clean.spec" = "species"))
colnames(phy.dat) <- gsub(" ", "", colnames(phy.dat))

phy.dat$clean.spec       <- as.factor(phy.dat$clean.spec)
phy.dat$clean.spec_resid <- phy.dat$clean.spec

# ---- Trait distance matrix ----
trait.sub   <- mcmillian.sc[as.character(unique(phy.dat$clean.spec)), ] %>% drop_na()
trait.dists <- as.matrix(dist(trait.sub, method = "euclidean"))
T_mat       <- 1 - (trait.dists / max(trait.dists))
diag(T_mat) <- 1
T_mat       <- as.matrix(nearPD(T_mat)$mat)   # ensure positive-definite

phy.dat$clean.spec_trait <- phy.dat$clean.spec

# ---- Trait dendrogram (for tree visualisation) ----
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
  )
###########===================================###########
## BUILD NEW STAN MODEL 
###=====================================================

stan_dat_distances <- standata(
  DIST_beta ~ log_delta_t + delta_PC1 + delta_PC2 + delta_PC3 +
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
  file = "data/stan_dirichlet_distances.RDS"
)

##
# Posterior visualization for the continuous beta DIST model
# Requires: fit_dirichlet_partition_distances, stan_dat_distances, phy.dat,
#           tidyverse, posterior, ggdist, writexl, Ternary, viridisLite
# Optional phylogeny figures: tree_final, trait_tree, plot_switching_tree()
# Optional trait association: trait.sub

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdist)
library(posterior)
library(stringr)
library(purrr)
library(writexl)

#dir.create("results", showWarnings = FALSE, recursive = TRUE)

# ---- 4c. Extract draws; centered predictor matrix ----
fit <- fit_dirichlet_partition_distances
sd_data <- stan_dat_distances
draws_df <- posterior::as_draws_df(fit$draws())

stopifnot(sd_data$Kc == 4L, sd_data$K == 5L)
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

# ---- 4d. Fixed effects ----
pretty_predictor <- function(x) {
  x <- sub("^delta_", "\u0394", x)
  x <- ifelse(x == "log_delta_t", "Log elapsed time", x)
  x
}

fes <- draws_df |>
  select(all_of(b_names)) |>
  pivot_longer(everything(), names_to = "param", values_to = "value") |>
  mutate(predictor = pretty_predictor(unname(predictor_labels[param])),
         predictor = factor(predictor,
                            levels = rev(pretty_predictor(unname(predictor_labels)))))

fe.p <- ggplot(fes, aes(x = value, y = predictor)) +
  stat_halfeye(.width = c(0.5, 0.95), point_interval = median_qi,
               fill = "#81C3D7", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(x = "Effect on logit(expected community distance)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")
print(fe.p)

ggsave("results/FEs_Distance_Beta.tiff", plot = fe.p,
       width = 4, height = 3.5, units = "in", dpi = 600, device = ragg::agg_tiff)

# ---- 4e. Fixed and random-effect summaries ----
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
writexl::write_xlsx(fixed_effects_table, "results/FE_DiriBeta.xlsx")

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
writexl::write_xlsx(random_summary, "results/RE_DiriBeta.xlsx")

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

ggsave("results/Posterior_REs_Distance.tiff", plot = re.p,
       width = 4, height = 6, units = "in", dpi = 600, device = ragg::agg_tiff)

# ---- 4f. Dirichlet variance partition ----
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
                    "results/Variance_Partition_Probabilities_Beta.xlsx")

h2_summary <- h2_draws |>
  group_by(coefficient, source) |>
  summarise(median = median(value), mean = mean(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975), .groups = "drop")
print(h2_summary)
writexl::write_xlsx(h2_summary, "results/Variance_Partition_Beta.xlsx")

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
ggsave("results/Variance_Partition_Beta.tiff", plot = h2.p,
       width = 9, height = 3, dpi = 600, device = ragg::agg_tiff)

# Ternary plots, joint posterior draws, same palette as original analysis.
# Requires package Ternary (with TernaryPlot, TernaryDensity, ColourTernary).
if (requireNamespace("Ternary", quietly = TRUE) &&
    requireNamespace("viridisLite", quietly = TRUE)) {
  png("results/Variance_Partitioned_Distance.png",
      width = 6.5, height = 3, units = "in", res = 500)
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

# ---- 4g. Species-level expected beta distance (conditional mu) ----
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
                "STUDY", "SITE", "COLONY_ID", "clean.spec",
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

predict_species_distance <- function(pc = "delta_PC1", value = 0,
                                     log_time = median(phy.dat$log_delta_t, na.rm = TRUE)) {
  new_x <- setNames(rep(0, length(model_means)), names(model_means))
  new_x["log_delta_t"] <- log_time
  if (!pc %in% c("delta_PC1", "delta_PC2", "delta_PC3")) stop("Unknown PC")
  new_x[pc] <- value
  
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
         varied_pc = pc, pc_value = value)
}

# Equal and opposite environmental changes: 90th percentile of absolute values.
q90_abs_pc1 <- unname(quantile(abs(phy.dat$delta_PC1), 0.90, na.rm = TRUE))
q90_abs_pc3 <- unname(quantile(abs(phy.dat$delta_PC3), 0.90, na.rm = TRUE))

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
                    "results/Species_Distance_Predictions.xlsx")

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
ggsave("results/Species_Distance_Responses.tiff", plot = species.p,
       width = 6, height = 3.5, dpi = 600, device = ragg::agg_tiff)



# ==============================================================
# 4h. Species-level predictions on phylogenetic and trait trees
# No custom plotting or prediction functions
# ==============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtree)
library(viridis)

# ---- 1. Extract species-specific intercept predictions ----
# Includes the population intercept and all three species-level
# random intercepts (phylogenetic, residual, trait).
# Does not include PC slopes or study/site/colony random effects.

eta_intercept <- matrix(
  intercept_draws,
  nrow = length(intercept_draws),
  ncol = nrow(species_map)
) + re_by_slot[[1]]

mu_intercept <- plogis(eta_intercept)

pred_intercept <- tibble(
  species = species_map$species,
  median = apply(mu_intercept, 2, median),
  lower = apply(mu_intercept, 2, quantile, probs = 0.025),
  upper = apply(mu_intercept, 2, quantile, probs = 0.975),
  level = "Intercept"
)

# ---- 2. Combine predictions for each tree ----
# Existing Low / Baseline / High predictions use median
# elapsed time, with the other environmental PCs set to zero.
# Intercept is the model's species-specific intercept component.

pred_phylo <- bind_rows(
  pred_intercept,
  pred_pc3_summary %>%
    select(species, median, lower, upper, level)
) %>%
  mutate(
    condition = factor(
      level,
      levels = c("Intercept", "Low", "Baseline", "High"),
      labels = c("Intercept", "Low PC3", "Zero change", "High PC3")
    )
  )

pred_trait <- bind_rows(
  pred_intercept,
  pred_pc1_summary %>%
    select(species, median, lower, upper, level)
) %>%
  mutate(
    condition = factor(
      level,
      levels = c("Intercept", "Low", "Baseline", "High"),
      labels = c("Intercept", "Low PC1", "Zero change", "High PC1")
    )
  )

# Shared color scale across all trees
fill_limits <- range(
  c(pred_phylo$median, pred_trait$median),
  na.rm = TRUE
)

# ==============================================================
# 3. Phylogenetic tree: Intercept and PC3 predictions
# ==============================================================

phylo_base <- ggtree(tree_final, layout = "rectangular")

phylo_tips <- phylo_base$data %>%
  filter(isTip) %>%
  select(label, y)

phylo_max_x <- max(phylo_base$data$x, na.rm = TRUE)
phylo_span <- diff(range(phylo_base$data$x, na.rm = TRUE))

phylo_plot_data <- pred_phylo %>%
  rename(label = species) %>%
  inner_join(phylo_tips, by = "label") %>%
  mutate(
    x_plot = phylo_max_x +
      case_when(
        condition == "Intercept"   ~ 0.08 * phylo_span,
        condition == "Low PC3"     ~ 0.15 * phylo_span,
        condition == "Zero change" ~ 0.22 * phylo_span,
        condition == "High PC3"    ~ 0.29 * phylo_span
      )
  )

# Add labels to the right of the prediction columns
phylo_labels <- phylo_tips %>%
  mutate(
    x_label = phylo_max_x + 0.37 * phylo_span,
    short_name = gsub("^(\\w)\\w+_", "\\1. ", label)
  )

p_phylo <- phylo_base +
  geom_point(
    data = phylo_plot_data,
    aes(x = x_plot, y = y, fill = median),
    shape = 22,
    size = 2.5,
    color = "gray30",
    stroke = 0.25,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = phylo_labels,
    aes(x = x_label, y = y, label = short_name),
    hjust = 0,
    size = 2.3,
    fontface = "italic",
    inherit.aes = FALSE
  ) +
  scale_fill_viridis_c(
    limits = fill_limits,
    option = "viridis",
    name = "E[DIST]"
  ) +
  expand_limits(x = phylo_max_x + 0.75 * phylo_span) +
  coord_cartesian(clip = "off") +
  labs(title = "Phylogenetic similarity") +
  theme_tree2() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 10)
  )

print(p_phylo)

ggsave(
  "results/Expected_Distance_Phylo.tiff",
  plot = p_phylo,
  width = 5,
  height = 6,
  dpi = 600,
  device = ragg::agg_tiff
)

# ==============================================================
# 4. Trait tree: Intercept and PC1 predictions
# ==============================================================

trait_base <- ggtree(trait_tree, layout = "rectangular")

trait_tips <- trait_base$data %>%
  filter(isTip) %>%
  select(label, y)

trait_max_x <- max(trait_base$data$x, na.rm = TRUE)
trait_span <- diff(range(trait_base$data$x, na.rm = TRUE))

trait_plot_data <- pred_trait %>%
  rename(label = species) %>%
  inner_join(trait_tips, by = "label") %>%
  mutate(
    x_plot = trait_max_x +
      case_when(
        condition == "Intercept"   ~ 0.08 * trait_span,
        condition == "Low PC1"     ~ 0.15 * trait_span,
        condition == "Zero change" ~ 0.22 * trait_span,
        condition == "High PC1"    ~ 0.29 * trait_span
      )
  )

trait_labels <- trait_tips %>%
  mutate(
    x_label = trait_max_x + 0.37 * trait_span,
    short_name = gsub("^(\\w)\\w+_", "\\1. ", label)
  )

p_trait <- trait_base +
  geom_point(
    data = trait_plot_data,
    aes(x = x_plot, y = y, fill = median),
    shape = 22,
    size = 2.5,
    color = "gray30",
    stroke = 0.25,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = trait_labels,
    aes(x = x_label, y = y, label = short_name),
    hjust = 0,
    size = 2.3,
    fontface = "italic",
    inherit.aes = FALSE
  ) +
  scale_fill_viridis_c(
    limits = fill_limits,
    option = "viridis",
    name = "E[DIST]"
  ) +
  expand_limits(x = trait_max_x + 0.75 * trait_span) +
  coord_cartesian(clip = "off") +
  labs(title = "Trait similarity") +
  theme_tree2() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 10)
  )

print(p_trait)

ggsave(
  "results/Expected_Distance_Trait.tiff",
  plot = p_trait,
  width = 5,
  height = 6,
  dpi = 600,
  device = ragg::agg_tiff
)

# ==============================================================
# 4i. Trait associations with expected community distance
# ==============================================================

high_pc1 <- pred_pc1_summary %>%
  filter(level == "High") %>%
  rename(
    clean.spec = species,
    expected_dist = median
  ) %>%
  left_join(
    trait.sub,
    by = c("clean.spec" = "species")
  )

# Predicted distance vs trait PC1
p_trait_pc1 <- ggplot(
  high_pc1,
  aes(x = trait_PC1, y = expected_dist)
) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Trait PC1",
    y = "Expected distance at high environmental PC1"
  ) +
  theme_bw()

print(p_trait_pc1)

# Baseline/intercept prediction vs trait PC1
baseline_trait <- pred_intercept %>%
  rename(
    clean.spec = species,
    expected_dist = median
  ) %>%
  left_join(
    trait.sub,
    by = c("clean.spec" = "species")
  )

p_baseline_trait <- ggplot(
  baseline_trait,
  aes(x = trait_PC1, y = expected_dist)
) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Trait PC1",
    y = "Species intercept expected distance"
  ) +
  theme_bw()

print(p_baseline_trait)

ggsave(
  "results/Baseline_Distance_vs_TraitPC1.tiff",
  plot = p_baseline_trait,
  width = 4,
  height = 3,
  dpi = 600,
  device = ragg::agg_tiff
)

# ==============================================================
# 4j. Align predicted distances with the two covariance matrices
# ==============================================================

species <- Reduce(
  intersect,
  list(
    rownames(A),
    rownames(T_mat),
    baseline_trait$clean.spec
  )
)

A_ord <- A[species, species, drop = FALSE]
T_mat_ord <- T_mat[species, species, drop = FALSE]

pred.vals <- baseline_trait$expected_dist[
  match(species, baseline_trait$clean.spec)
]

names(pred.vals) <- species

vals <- as.matrix(dist(pred.vals))

# A_ord, T_mat_ord, and vals have identical species ordering.

