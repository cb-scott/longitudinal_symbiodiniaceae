################################################################################
# BRMS/STAN MODELING OF SYMBIONT SWITCHING WITH DIRICHLET VARIANCE PARTITION
################################################################################
# Builds phylogenetic and trait covariance matrices for the coral host species,
# then fits three Stan models (ordinal severity, binary, and first-switch-only)
# testing symbiont switching against environmental PCs and prior symbiont
# state, partitioning species-level variance into phylogenetic/trait/residual
# components. Produces convergence diagnostics, fixed/random effect summaries,
# variance-partition ternary plots, and phylo/trait-tree switching-probability
# maps. Model fits (~20-45 min runtime each) are provided as .RDS objects so  
# sampling can be optionally skipped.
################################################################################

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

source("scripts/helpers.R")

#REDUCED MODEL TO COMPARE DATES

# ==============================================================================
# 1. DATA PROCESSING
# ==============================================================================

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

#reduced model -- change me on rerun
dat <- dat %>% filter(DATE >= as.Date("2012-08-30") & DATE <= as.Date("2024-04-28"))

# ---- Species name standardisation ----
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

# 3. Pocillopora capitata — sister to ETP clade (inflata + elegans)
etp_node <- getMRCA(tree, c("Pocillopora_inflata", "Pocillopora_elegans"))
bl <- tree$edge.length[which(tree$edge[, 2] == etp_node)]
tree <- bind.tip(tree, "Pocillopora_capitata",
                 where = etp_node, edge.length = bl, position = 0.5 * bl)

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
mcmillian.sc$species[mcmillian.sc$species == "Dipsastrea_pallida"] <- "Dipsastraea_pallida"
rownames(mcmillian.sc) <- mcmillian.sc$species
colnames(mcmillian.sc) <- gsub(" ", "", colnames(mcmillian.sc))

# ---- Assemble analysis dataset ----
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

phy.dat <- phy.dat %>% mutate(SITE = paste(LONGITUDE, LATITUDE, sep = "_"))

# ==============================================================================
# 4. (new)  DIRILECHT VAR PARTITION WITH DELTA PCS
# ==============================================================================
#ah I've probably already dropped the baseline observation, which is a problem for this centering
#-----4a. Set up standata and check variable naming conventions------
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

#--- 4b. read in the stan model that has a RE for site in addition to study and colony -----
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
  file = "data/stan_dirichlet_explicit_delta_REDUCED.RDS"
)

#fit_dirichlet_partition_within = readRDS("data/stan_dirichlet_explicit_delta.RDS")

# ---- 4c. Extract draws and reconstruct variables names  ----
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


# ----- 4d. Extract fixed effects from draws ------- 
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
  filter(!(predictor %in% c("A_C_D", "B_C_D", "D_G"))) |>
  ggplot(aes(x = value, y = predictor)) +
  stat_halfeye(.width = c(0.5, 0.95), point_interval = median_qi,
               aes(fill = interaction(category, is_mixed)), alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(x = "Log-odds of switching (cumulative)", y = NULL) +
  labs(fill = "Codominance") + 
  #    title = "Fixed effect posteriors",
  #    subtitle = "Median with 50% and 95% credible intervals") +
  theme_bw(base_size = 10) + 
  theme(legend.position = "none")+
  scale_fill_manual(values = c("#c3d781", "#81C3D7", "#d7c081"))+
  facet_wrap(~category, nrow = 2, scales = "free_y")
fe.p

ggsave(fe.p, file = "results/fes_zoox_enviro_REDUCED.svg", width =4, height = 4)

#------ 4e. Fixed and Random Effects Summary ----- 
library(posterior)
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
writexl::write_xlsx(fixed_effects_table, "results/FE_DiriOrdinal_REDUCED.xlsx")


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
#library(writexl)
writexl::write_xlsx(random_summary, "results/RE_DiriOrdinal_REDUCED.xlsx")

###histogram the SD, but not the % that was estimated

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

ggsave("results/Posterior_REs_REDUCED.svg", width = 4, height = 6)


# ---- 4f. Variance partition (h2) from parameter estimation----

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


#quantify posterior joints
#This code isn't working right

# h2_draws |>
#   pivot_wider(names_from = source, values_from = value) |>
#   group_by(coefficient) |>
#   summarise(
#     p_phylo_dominant = mean(Phylogenetic > 1/3),
#     p_trait_dominant = mean(Trait > 1/3),
#     p_resid_dominant = mean(Idiosyncratic > 1/3),
#     p_trait_leads    = mean(Trait > Phylogenetic & Trait > Idiosyncratic),
#     p_phylo_leads    = mean(Phylogenetic > Trait & Phylogenetic > Idiosyncratic)
#   )

coef_labels <- c("Baseline rate", "\u0394PC1", "\u0394PC2", "\u0394PC3")
names(coef_labels) = c("Baseline rate", "PC1 sensitivity", "PC2 sensitivity", "PC3 sensitivity")

svg("results/Variance_Partitioned_delta_REDUCED.svg", width = 6.5, height = 3)
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

