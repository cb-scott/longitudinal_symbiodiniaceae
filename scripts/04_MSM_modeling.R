load("data/LongZooxInputData_wEnv.Rdata")
source("scripts/helpers.R")
library(msm)
library(tidyverse)
library(igraph)
set.seed(123)


dat <- genlevel |>
  filter(GENUS != "NA") %>%
  group_by(COLONY_ID) |>
  distinct(STUDY, HOST_SPEC_REVISED, SPECIES, COLONY_ID, DATE, .keep_all = TRUE) |>
  arrange(DATE) |>
  mutate(
    baseline_genus = first(GENUS),
    state = case_when(
     # is.na(GENUS)           ~ 4L,  # lost/dead
      GENUS == baseline_genus ~ 1L,  # original
      GENUS != baseline_genus ~ 2L   # switched
    )
  ) %>%
  ungroup() %>% dplyr::mutate(ocean_basin = assign_ocean_basin(LONGITUDE, LATITUDE)) %>%
  mutate(Year = year(DATE),
         Years_Since_2000 = Year - 2000) %>%
  arrange(COLONY_ID, DATE)

dat <- dat |>
  group_by(COLONY_ID) |>
  filter(n() > 1) |>
  mutate(numeric_years = (as.numeric(as.Date(DATE))-min(as.numeric(as.Date(DATE))))/365) |>
  ungroup()

Q_init <- rbind(
  c(0,   0.1),  # Original → Switched, Absorbed
  c(0.5, 0)  # Switched → Original, Absorbed
)

fit <- msm(state ~ numeric_years, subject = COLONY_ID,
           data = dat,
           qmatrix = Q_init)
           #covariates = ~ PC1 + PC2 + PC3)  # robust SEs for study clustering
fit




fit.pcs <- msm(state ~ numeric_years, subject = COLONY_ID,
           data = dat,
           qmatrix = Q_init,
          covariates = ~ PC1 + PC2 + PC3)  # robust SEs for study clustering
fit.pcs #this is for switching and reversion



msm.out <- data.frame(
  transition   = rownames(msm.form.qoutput(fit.pcs)),
  base.Estimate = msm.form.qoutput(fit.pcs)$base.Estimate,
  base.L        = msm.form.qoutput(fit.pcs)$base.L,
  base.U        = msm.form.qoutput(fit.pcs)$base.U,
  PC1.Estimate  = msm.form.qoutput(fit.pcs)$PC1[, "Estimate"],
  PC1.L         = msm.form.qoutput(fit.pcs)$PC1[, "L"],
  PC1.U         = msm.form.qoutput(fit.pcs)$PC1[, "U"],
  PC2.Estimate  = msm.form.qoutput(fit.pcs)$PC2[, "Estimate"],
  PC2.L         = msm.form.qoutput(fit.pcs)$PC2[, "L"],
  PC2.U         = msm.form.qoutput(fit.pcs)$PC2[, "U"],
  PC3.Estimate  = msm.form.qoutput(fit.pcs)$PC3[, "Estimate"],
  PC3.L         = msm.form.qoutput(fit.pcs)$PC3[, "L"],
  PC3.U         = msm.form.qoutput(fit.pcs)$PC3[, "U"]
)

hr_long <- msm.out |>
  pivot_longer(
    cols = -transition,
    names_to = c("covariate", ".value"),
    names_sep = "\\."
  ) |>
  rename(L = L, U = U, HR = Estimate) |>
  filter(!is.na(HR), HR != 0)

hrs = hr_long %>% filter(transition != "State 1 - State 1") %>%
  filter(transition != "State 2 - State 2") %>%
  filter(grepl("^State 2", transition) == F)

hrs = hrs %>% filter(covariate != "base") %>% mutate(sig = ifelse(L < 1 & U > 1, "ns", "significant"))
hrs$covariate = factor(hrs$covariate, levels = rev(sort(hrs$covariate)))

hrs %>% ggplot(aes(x = HR, y=covariate, col = sig)) +
  geom_errorbarh(aes(xmin = L, xmax = U), width = .1) + 
  geom_vline(xintercept = 1) + 
  geom_point(size = 3) + 
  scale_color_manual(values = c( "gray60", "black")) + 
  theme_classic() + 
  theme(legend.position = "none",
        text = element_text(size = 8),
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5)) +
  ylab("") + 
  xlab("HR Initial Switch")
ggsave("results/Switch_fxn_PCs.png", width = 1.5, height = 2, units = "in", dpi = 500)

#### Just plot the baseline fits


# Extract Q matrix estimates
Q_est <- qmatrix.msm(fit.pcs)$estimates
Q_L <- qmatrix.msm(fit.pcs)$L
Q_H <- qmatrix.msm(fit.pcs)$U
#Q_est <- pmatrix.msm(fit, t = 1)

# State names
state_names <- c("Base", "Switch")
rownames(Q_est) <- state_names
colnames(Q_est) <- state_names
rownames(Q_L) <- state_names
colnames(Q_L) <- state_names
rownames(Q_H) <- state_names
colnames(Q_H) <- state_names
# Build edge list
edges <- data.frame(from = character(), to = character(), rate = numeric())

for (i in 1:nrow(Q_est)) {
  for (j in 1:ncol(Q_est)) {
    if (i != j && Q_est[i, j] > 0) {
      edges <- rbind(edges, data.frame(
        from = state_names[i],
        to   = state_names[j],
        rate = Q_est[i, j]
      ))
    }
  }
}

# Build graph
g <- graph_from_data_frame(edges, directed = TRUE)

# Curve bidirectional edges, leave unidirectional straight
mutual_edges <- which_mutual(g)
edge_curves  <- ifelse(mutual_edges, 0.3, 0)

# Fix layout
layout_coords <- layout_in_circle(g)

# Plot
png("results/Baseline_Switch.png", width = 2, height = 2, units = "in", res = 500)
par(mar = c(.5, .5, .5, .5), bg = "white", family = "Arial", ps = 12)
plot(g,
     layout             = layout_coords,
     edge.width         = edges$rate * 5,
     edge.arrow.size    = 0.5,
     edge.arrow.width   = 2,
     edge.curved        = edge_curves,
     edge.label         = round(edges$rate, 3),
     edge.label.cex     = 8/12,
     edge.label.color   = "black",
     edge.label.family  = "Arial",
     edge.color         = "#81C3D7",
     vertex.shape       = "circle",
     vertex.size        = 60,
     vertex.color       = "gray95",
     vertex.frame.color = "black",
     vertex.label.color = "black",
     vertex.label.cex   = 8/12,
     vertex.label.family = "Arial")
dev.off()


####### PREDICT SPECIES SPECIFIC RATES FROM THE MODEL ##############

# Get species-mean covariate values
species_covs <- dat |>
  group_by(SPECIES) |>
  summarise(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    PC3 = mean(PC3, na.rm = TRUE),
    n_colonies = n_distinct(COLONY_ID),
    .groups = "drop"
  ) |>
  filter(n_colonies >= 3)  # drop species with very few colonies

# Predict transition intensities at species-mean covariates
species_rates <- species_covs |>
  rowwise() |>
  mutate(
    qmat        = list(qmatrix.msm(fit.pcs,
                                   covariates = list(PC1 = PC1, PC2 = PC2, PC3 = PC3))),
    rate_switch = qmat$estimates[1, 2],
    rate_revert = qmat$estimates[2, 1]
  ) |>
  mutate(
    mean_time_to_switch = 1 / rate_switch,
    mean_time_to_revert = 1 / rate_revert
  ) |>
  select(-qmat) |>
  ungroup()

print(species_rates)

species_rates |>
  tidyr::pivot_longer(cols = c(rate_switch, rate_revert),
                      names_to = "transition", values_to = "rate") |>
  mutate(transition = recode(transition,
                             rate_switch = "Switching (1→2)",
                             rate_revert = "Reversion (2→1)")) |>
  ggplot(aes(x = rate, y = reorder(SPECIES, rate), color = transition)) +
  geom_point(aes(size = n_colonies)) +
  scale_color_manual(values = c("Switching (1→2)" = "#E57C3A",
                                "Reversion (2→1)" = "#4A90BF")) +
  labs(x = "Transition intensity (per year)", y = NULL,
       size = "N colonies", color = "Transition") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))

#######################################
##### MODEL IT WITH THESE OTEHR STATE TRANSITIONS BUT GET RID OF RARE ONES

dat <- genlevel |>
  filter(GENUS != "NA") %>%
  group_by(COLONY_ID) |>
  distinct(STUDY, HOST_SPEC_REVISED, SPECIES, COLONY_ID, DATE, .keep_all = TRUE) |>
  arrange(DATE) |>
  filter(GENUS %in% c("A", "A_C", "A_D", "B", "B_C", "C", "C_D", "D")) %>%
  ungroup() %>% dplyr::mutate(ocean_basin = assign_ocean_basin(LONGITUDE, LATITUDE)) %>%
  mutate(Year = year(DATE),
         Years_Since_2000 = Year - 2000) %>%
  arrange(COLONY_ID, DATE)

dat <- dat |>
  group_by(COLONY_ID) |>
  filter(n() > 1) |>
  mutate(numeric_years = (as.numeric(as.Date(DATE))-min(as.numeric(as.Date(DATE))))/365) |>
  ungroup()

#genus -> states 
stateids = data.frame(GENUS = sort(unique(dat$GENUS)), state = 1:length(unique(dat$GENUS)))
dat = left_join(dat, stateids)

mixed_states = c("A_C", "A_D", "B_C", "C_D")
dat <- dat |>
  mutate(state_censored = ifelse(GENUS %in% mixed_states, NA, state))

stateids_pure <- data.frame(
  GENUS = c("A", "B", "C", "D"),
  state = 1:4
)
stateids_censored <- data.frame(
  GENUS=mixed_states,
  state = 900:903
)
new_stateids = rbind(stateids_pure, stateids_censored)

dat <- dat |>
  select(-state) |>
  left_join(new_stateids, by = "GENUS") #|>
 # mutate(state_new = ifelse(GENUS %in% mixed_states, NA, state))

observed <- statetable.msm(state, COLONY_ID, data = dat[dat$state < 100,])
Q_init <- matrix(0, nrow = 4, ncol = 4)
Q_init[observed > 0] <- 0.005
diag(Q_init) <- 0
Q_init
#any transition is possible 
# Check observed transitions first
statetable.msm(state, COLONY_ID, data = dat)


fit <- msm(state ~ numeric_years, subject = COLONY_ID,
           data = dat,
           qmatrix = Q_init,
           control = list(fnscale = 8000, maxit = 3000),
           censor = c(900, 901, 902, 903),
           censor.states = list(c(1,3), c(1,4), c(2,3), c(3,4)))
#covariates = ~ PC1 + PC2 + PC3)  # robust SEs for study clustering
fit


# fit.pcs <- msm(state ~ numeric_years, subject = COLONY_ID,
#                data = dat,
#                qmatrix = Q_init,
#                control = list(fnscale = 1000, maxit = 5000),
#                censor = c(900, 901, 902, 903),
#                censor.states = list(c(1,3), c(1,4), c(2,3), c(3,4)),
#                covariates = ~ PC1 + PC2 + PC3)  # robust SEs for study clustering
# fit.pcs #this super doesn't want to converge, perhaps just ignore it!
# save(fit.pcs, file = "data/fit_msm_wPCs.Rdata")
# 

library(ggplot2)
library(dplyr)
#load("data/fit_msm_wPCs.Rdata")
# Extract hazard ratios from model
hr_table <- hazard.msm(fit)
#need to make a forest plot from this.



library(igraph)

# Extract Q matrix estimates
Q_est <- qmatrix.msm(fit)$estimates
Q_L <- qmatrix.msm(fit)$L
Q_H <- qmatrix.msm(fit)$U
#Q_est <- pmatrix.msm(fit, t = 1)

# State names
state_names <- c("A", "B", "C", "D")
rownames(Q_est) <- state_names
colnames(Q_est) <- state_names
rownames(Q_L) <- state_names
colnames(Q_L) <- state_names
rownames(Q_H) <- state_names
colnames(Q_H) <- state_names
# Build edge list
edges <- data.frame(from = character(), to = character(), rate = numeric())

for (i in 1:nrow(Q_est)) {
  for (j in 1:ncol(Q_est)) {
    if (i != j && Q_est[i, j] > 0) {
      edges <- rbind(edges, data.frame(
        from = state_names[i],
        to   = state_names[j],
        rate = Q_est[i, j]
      ))
    }
  }
}

# Build graph
g <- graph_from_data_frame(edges, directed = TRUE)

# Curve bidirectional edges, leave unidirectional straight
mutual_edges <- which_mutual(g)
edge_curves  <- ifelse(mutual_edges, 0.3, 0)

# Fix layout
layout_coords <- layout_with_fr(g)

# Plot
png("results/ZooxNetwork.png", width = 3.5, height = 3.5, units = "in", res = 500)
par(mar = c(.5, .5, .5, .5), bg = "white", family = "Arial", ps = 12)
plot(g,
     layout             = layout_coords,
     edge.width         = edges$rate * 50,
     edge.arrow.size    = 0.5,
     edge.arrow.width  = 2,
     edge.curved        = edge_curves,
     edge.label         = round(edges$rate, 3),
     edge.label.cex     = 8/12,
     edge.label.color   = "black",
     edge.label.family  = "Arial",
     edge.color         = "#81C3D7",
     vertex.size        = 35,
     vertex.color       = "gray95",
     vertex.frame.color = "black",
     vertex.label.color = "black",
     vertex.label.cex   = 1.0,
     vertex.label.family = "Arial")
dev.off()

#forest plot
df.est <- as.data.frame(Q_est) |>
  tibble::rownames_to_column("from") |>
  pivot_longer(-from, names_to = "to", values_to = "rate")
df.L <-  as.data.frame(Q_L) |>
  tibble::rownames_to_column("from") |>
  pivot_longer(-from, names_to = "to", values_to = "L")
df.U <- as.data.frame(Q_H) |>
  tibble::rownames_to_column("from") |>
  pivot_longer(-from, names_to = "to", values_to = "U")
df.all = left_join(left_join(df.est, df.L), df.U) %>%
  filter(from != to) %>%
  mutate(label = paste(from, to, sep = " \u2192 ")) %>%
  mutate(pair = paste(pmin(from, to), pmax(from, to), sep = "_")) |>
  arrange(pair, from)

df.overlap <- df.all |>
  group_by(pair) |>
  summarise(
    overlap = max(L) < min(U),
    .groups = "drop"
  ) |>
  right_join(df.all, by = "pair")

df.overlap$label = factor(df.overlap$label, levels = df.overlap$label)

df.overlap %>%
  ggplot(aes(x = rate, y = label)) + 
  geom_point(aes(color = overlap)) + 
  geom_errorbarh(aes(xmin = L, xmax = U, col = overlap), width = .2) + 
  theme_classic() +
  scale_color_manual(values = c("black", "gray60")) +
  ylab("") + 
  xlab("Annual Switching Rate") + 
  facet_wrap(~pair, ncol = 1, scales = "free_y") +
  theme(
    legend.position = "none",
    text=element_text(size = 8),
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = .5)
  )

ggsave("results/Ann_Switch_Rate.png", width = 3, height = 3.5, units = "in", dpi = 500)


#################################################
### SPECIES LEVEL SWITCHING AND REVERSION ########
#################################################
# Check how many species would pass each filter
species_summary <- dat |>
  group_by(SPECIES) |>
  summarise(
    n_colonies   = n_distinct(COLONY_ID),
    n_obs        = n(),
    n_switches   = sum(state == 2 & lag(state, default = 1) == 1, na.rm = TRUE),
    n_reversions = sum(state == 1 & lag(state, default = 1) == 2, na.rm = TRUE),
    .groups = "drop"
  ) 

print(species_summary)


# ============================================================
# Species-level multi-state Markov models
# Fits switching/reversion model with PC covariates
# for species with >= 10 colonies
# ============================================================

library(msm)
library(tidyverse)

# ---- Identify qualifying species ----
qualifying_species <- dat |>
  group_by(SPECIES) |>
  summarise(
    n_colonies = n_distinct(COLONY_ID),
    n_obs      = n(),
    n_switches = sum(diff(state) != 0, na.rm = TRUE),
    .groups    = "drop"
  ) |>
  filter(n_colonies >= 3)

cat("Species qualifying (>= 3 colonies):", nrow(qualifying_species), "\n")
print(qualifying_species)

# ---- Split data by species ----
species_dat   <- dat |>
  filter(SPECIES %in% qualifying_species$SPECIES) |>
  group_split(SPECIES)

species_names <- dat |>
  filter(SPECIES %in% qualifying_species$SPECIES) |>
  pull(SPECIES) |>
  unique() |>
  sort()

# ---- Fit per species ----
species_fits <- map2(species_dat, species_names, function(d, sp) {
  tryCatch(
    msm(state ~ numeric_years,
        subject    = COLONY_ID,
        data       = d,
        qmatrix    = Q_init,
        covariates = ~ PC1 + PC2 + PC3),
    error   = function(e) { message("Failed: ", sp, " — ", e$message);   NULL },
    warning = function(w) { message("Warning: ", sp, " — ", w$message); NULL }
  )
})
names(species_fits) <- species_names

# Remove failed fits
failed        <- names(species_fits)[sapply(species_fits, is.null)]
species_fits  <- Filter(Negate(is.null), species_fits)
cat("Successful fits:", length(species_fits), "of", length(species_names), "\n")
if (length(failed) > 0) cat("Failed species:", paste(failed, collapse = ", "), "\n")

# ---- Extract hazard ratios for each species ----
extract_hrs <- function(fit, sp) {
  tryCatch({
    out <- msm.form.qoutput(fit)
    df  <- data.frame(
      SPECIES       = sp,
      transition    = rownames(out),
      base.Estimate = out$base.Estimate,
      base.L        = out$base.L,
      base.U        = out$base.U,
      PC1.Estimate  = out$PC1[, "Estimate"],
      PC1.L         = out$PC1[, "L"],
      PC1.U         = out$PC1[, "U"],
      PC2.Estimate  = out$PC2[, "Estimate"],
      PC2.L         = out$PC2[, "L"],
      PC2.U         = out$PC2[, "U"],
      PC3.Estimate  = out$PC3[, "Estimate"],
      PC3.L         = out$PC3[, "L"],
      PC3.U         = out$PC3[, "U"]
    )
    df
  }, error = function(e) { message("Extract failed: ", sp); NULL })
}

all_hrs_raw <- map2_dfr(species_fits, names(species_fits), extract_hrs)

# ---- Reshape to long format ----
all_hrs <- all_hrs_raw |>
  pivot_longer(
    cols      = -c(SPECIES, transition),
    names_to  = c("covariate", ".value"),
    names_sep = "\\."
  ) |>
  rename(HR = Estimate) |>
  filter(!is.na(HR), HR != 0) |>
  # Keep only switching transition (State 1 -> State 2) and PC covariates
  filter(transition == "State 1 - State 2") |>
  filter(covariate != "base") |>
  mutate(
    sig        = ifelse(L < 1 & U > 1, "ns", "significant"),
    covariate  = factor(covariate, levels = c("PC1", "PC2", "PC3")),
    SPECIES    = factor(SPECIES)
  ) |>
  left_join(qualifying_species, by = "SPECIES")

# ---- Forest plot: HR per species per PC ----
all_hrs |>
  ggplot(aes(x = HR, y = reorder(SPECIES, HR), col = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = L, xmax = U), height = 0.2) +
  geom_point(aes(size = n_colonies)) +
  scale_color_manual(values = c("ns" = "gray60", "significant" = "black")) +
  scale_size_continuous(range = c(1, 4)) +
  facet_wrap(~ covariate, nrow = 1, scales = "free_x") +
  labs(
    x     = "Hazard ratio for switching",
    y     = NULL,
    size  = "N colonies",
    color = NULL,
    title = "Species-level covariate effects on switching rate"
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.y    = element_text(face = "italic"),
    legend.position = "bottom",
    panel.border   = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.text     = element_text(face = "bold")
  )
ggsave("results/Species_HR_Switch.png", width = 7, height = 6, dpi = 500)

# ---- Also extract baseline switching and reversion rates ----
baseline_rates <- map_dfr(names(species_fits), function(sp) {
  fit  <- species_fits[[sp]]
  qmat <- qmatrix.msm(fit)
  tibble(
    SPECIES         = sp,
    rate_switch     = qmat$estimates[1, 2],
    rate_switch_L   = qmat$L[1, 2],
    rate_switch_U   = qmat$U[1, 2],
    rate_revert     = qmat$estimates[2, 1],
    rate_revert_L   = qmat$L[2, 1],
    rate_revert_U   = qmat$U[2, 1]
  )
}) |>
  left_join(qualifying_species, by = "SPECIES") |>
  mutate(
    mean_time_to_switch = 1 / rate_switch,
    mean_time_to_revert = 1 / rate_revert
  )

print(baseline_rates)

# ---- Plot baseline switching and reversion rates ----
baseline_rates |>
  pivot_longer(
    cols      = c(rate_switch, rate_revert),
    names_to  = "transition",
    values_to = "rate"
  ) |>
  mutate(
    L = ifelse(transition == "rate_switch", rate_switch_L, rate_revert_L),
    U = ifelse(transition == "rate_switch", rate_switch_U, rate_revert_U),
    transition = recode(transition,
                        rate_switch = "Switching (1\u21922)",
                        rate_revert = "Reversion (2\u21921)")
  ) |>
  ggplot(aes(x = rate, y = reorder(SPECIES, rate), color = transition)) +
  geom_errorbarh(aes(xmin = L, xmax = U), height = 0.2) +
  geom_point(aes(size = n_colonies)) +
  scale_color_manual(values = c(
    "Switching (1\u21922)" = "#E57C3A",
    "Reversion (2\u21921)" = "#4A90BF"
  )) +
  scale_size_continuous(range = c(1, 4)) +
  labs(
    x     = "Transition intensity (per year)",
    y     = NULL,
    size  = "N colonies",
    color = "Transition",
    title = "Species-level baseline switching and reversion rates"
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.y     = element_text(face = "italic"),
    legend.position = "bottom",
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
ggsave("results/Species_Baseline_Rates.png", width = 4, height = 6, dpi = 500)

# ---- Save results ----
save(species_fits,    file = "data/msm_species_fits.Rdata")
write_csv(all_hrs,        "results/species_msm_hrs.csv")
write_csv(baseline_rates, "results/species_msm_baseline_rates.csv")