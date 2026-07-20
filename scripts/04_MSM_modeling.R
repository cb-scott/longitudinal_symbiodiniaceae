################################################################################
# MULTISTATE MARKOV (MSM) MODELS OF SYMBIONT SWITCHING RATES
################################################################################
# Fits msm models to genus-level symbiont data: first a two-state (baseline
# vs. switched) model of annual switching/reversion probability as a function
# of environmental PCs, then a four-state (A/B/C/D) pairwise transition rate
# model with mixed-genotype states censored. Produces hazard ratio plots,
# transition-rate network diagrams, and pairwise switching-rate forest plots.
################################################################################



load("data/LongZooxInputData_wEnv.Rdata")
source("scripts/helpers.R")
library(msm)
library(tidyverse)
library(igraph)
set.seed(123)


#===============================================================
# MODEL ANNUAL SWITCHING AND REVERSION PROBABILITIES 
#===============================================================

# ------ Clean data determine if a switch event happened at any level -------
dat <- genlevel |>
  filter(GENUS != "NA") %>%
  group_by(COLONY_ID) |>
  distinct(STUDY, HOST_SPEC_REVISED, SPECIES, COLONY_ID, DATE, .keep_all = TRUE) |>
  arrange(DATE) |>
  mutate(
    baseline_genus = first(GENUS),
    state = case_when(
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
  group_by(STUDY, HOST_SPEC_REVISED, SPECIES, COLONY_ID) |>
  mutate(
    delta_PC1   = PC1 - lag(PC1, 1),
    delta_PC2 = PC2 - lag(PC2, 1),
    delta_PC3 = PC3 - lag(PC3, 1)
  ) %>%
  ungroup()

#------ Supply initial transition probabilities ------- 

Q_init <- rbind(
  c(0,   0.1),  # Original → Switched
  c(0.5, 0)  # Switched → Original
)

#----- Fit and check an initial model without PC covariates -----

fit <- msm(state ~ numeric_years, subject = COLONY_ID,
           data = dat,
           qmatrix = Q_init)
fit 

#----- Fit a final model with PC covariates ------
fit.pcs <- msm(state ~ numeric_years, subject = COLONY_ID,
           data = dat,
           qmatrix = Q_init,
          covariates = ~ delta_PC1 + delta_PC2 + delta_PC3)
fit.pcs 


#------ Extract transition rates and CIs ------ 
msm.out <- data.frame(
  transition   = rownames(msm.form.qoutput(fit.pcs)),
  base.Estimate = msm.form.qoutput(fit.pcs)$base.Estimate,
  base.L        = msm.form.qoutput(fit.pcs)$base.L,
  base.U        = msm.form.qoutput(fit.pcs)$base.U,
  PC1.Estimate  = msm.form.qoutput(fit.pcs)$delta_PC1[, "Estimate"],
  PC1.L         = msm.form.qoutput(fit.pcs)$delta_PC1[, "L"],
  PC1.U         = msm.form.qoutput(fit.pcs)$delta_PC1[, "U"],
  PC2.Estimate  = msm.form.qoutput(fit.pcs)$delta_PC2[, "Estimate"],
  PC2.L         = msm.form.qoutput(fit.pcs)$delta_PC2[, "L"],
  PC2.U         = msm.form.qoutput(fit.pcs)$delta_PC2[, "U"],
  PC3.Estimate  = msm.form.qoutput(fit.pcs)$delta_PC3[, "Estimate"],
  PC3.L         = msm.form.qoutput(fit.pcs)$delta_PC3[, "L"],
  PC3.U         = msm.form.qoutput(fit.pcs)$delta_PC3[, "U"]
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

# ------ Plot hazard ratios --------
hrs %>% ggplot(aes(x = HR, y=covariate, col = sig)) +
  geom_errorbarh(aes(xmin = L, xmax = U), width = .1) + 
  geom_vline(xintercept = 1) + 
  geom_point(size = 3) + 
  scale_color_manual(values = c( "gray60", "black")) + 
  scale_x_continuous(breaks = c(1, 1.1, 1.2)) +
  theme_classic() + 
  theme(legend.position = "none",
        text = element_text(size = 8),
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5)) +
  ylab("") + 
  xlab("HR Initial Switch")
ggsave("results/Switch_fxn_PCs.png", width = 1.5, height = 2, units = "in", dpi = 500)

# ------ Plot switch and revert graph ----- 

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


#===============================================================
# MODEL ANNUAL SYMBIONT GENUS PAIRWISE TRANSITION RATES 
#===============================================================

# ------ Clean data to get rid of rare transition states -------- 

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

stateids = data.frame(GENUS = sort(unique(dat$GENUS)), state = 1:length(unique(dat$GENUS)))
dat = left_join(dat, stateids)

#---- censor mixed states ----- 
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
  left_join(new_stateids, by = "GENUS") 

#----- initialize transtions which were observed in the data at 0.005 ------ 
observed <- statetable.msm(state, COLONY_ID, data = dat[dat$state < 100,])
Q_init <- matrix(0, nrow = 4, ncol = 4)
Q_init[observed > 0] <- 0.005
diag(Q_init) <- 0
Q_init
#any transition is possible 

#---- fit msm model with censored states ----- 
fit <- msm(state ~ numeric_years, subject = COLONY_ID,
           data = dat,
           qmatrix = Q_init,
           control = list(fnscale = 8000, maxit = 3000),
           censor = c(900, 901, 902, 903),
           censor.states = list(c(1,3), c(1,4), c(2,3), c(3,4)))
fit

# Extract hazard ratios from model
hr_table <- hazard.msm(fit)

# ------ Plot pairwise switching graph ----- 

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

#------ Plot pairwise switching forest plot------- 
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
