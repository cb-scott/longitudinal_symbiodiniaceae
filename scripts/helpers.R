# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# ------------------------------------------------------------------------------
# assign_ocean_basin()
# Map longitude/latitude coordinates to named ocean basins.
# ------------------------------------------------------------------------------
assign_ocean_basin <- function(lon, lat) {
  lon <- ifelse(lon > 180, lon - 360, lon)
  dplyr::case_when(
    lon >= 48  & lon <= 57  & lat >= 22 & lat <= 30 ~ "Persian Gulf",
    lon >= -85 & lon <= -55 & lat >= 10 & lat <= 35 ~ "Gulf of Mexico / Caribbean / Western Atlantic",
    lon >= -155 & lon <= -70 & lat >= -10 & lat <= 30 ~ "Eastern Tropical Pacific",
    (lon >= 30 & lon <= 180) | (lon >= -180 & lon < -130) ~ "Indo-Pacific",
    TRUE ~ "Other / Unassigned"
  )
}

# ------------------------------------------------------------------------------
# plot_species_effects_phylo()
# Caterpillar plot of species-level random effects from a brms model.
# group: grouping term name (e.g., "clean.spec", "clean.spec_trait")
# param: coefficient name (e.g., "PC1", "Intercept")
# ------------------------------------------------------------------------------
plot_species_effects_phylo <- function(model, param, group = "clean.spec_resid") {
  est_col <- paste0("Estimate.", param)
  q2_col  <- paste0("Q2.5.",    param)
  q97_col <- paste0("Q97.5.",   param)
  
  species_output <- coef(model)[[group]] |>
    as.data.frame() |>
    rownames_to_column("SPECIES") |>
    mutate(GENUS = word(SPECIES, 1, sep = "_"))
  
  order_plot <- species_output |>
    group_by(SPECIES) |>
    mutate(mean_effect = mean(.data[[est_col]])) |>
    arrange(desc(mean_effect))
  
  species_output$SPECIES <- factor(species_output$SPECIES,
                                   levels = unique(order_plot$SPECIES))
  
  species_output <- species_output |>
    mutate(signif = case_when(
      .data[[q2_col]]  > 0 ~ "more flexible",
      .data[[q97_col]] < 0 ~ "less flexible",
      .default = "ns"
    ))
  
  species_output |>
    ggplot(aes(y = SPECIES, x = .data[[est_col]], col = signif)) +
    geom_errorbarh(aes(xmin = .data[[q2_col]], xmax = .data[[q97_col]]), height = 0.3) +
    geom_point(size = 2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    labs(x = paste0(param, " (log-odds)"), y = NULL,
         title = paste0("Species-level ", param, " effects on symbiont switching"),
         col = NULL) +
    scale_color_manual(values = c("more flexible" = "red3",
                                  "less flexible" = "blue3",
                                  "ns"            = "gray60")) +
    scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
    theme_bw() +
    theme(axis.text.y = element_text(face = "italic"))
}

# ------------------------------------------------------------------------------
# plot_species_effects()
# Simpler version using HOST_SPEC_REVISED grouping (legacy; kept for
# compatibility with older model objects).
# ------------------------------------------------------------------------------
plot_species_effects <- function(model, param) {
  est_col <- paste0("Estimate.", param)
  q2_col  <- paste0("Q2.5.",    param)
  q97_col <- paste0("Q97.5.",   param)
  
  ordinal_output <- coef(model)$`HOST_SPEC_REVISED` |>
    as.data.frame() |>
    rownames_to_column("SPECIES")
  
  order_plot <- ordinal_output |>
    group_by(SPECIES) |>
    mutate(mean_effect = mean(.data[[est_col]])) |>
    arrange(desc(mean_effect))
  
  ordinal_output$SPECIES <- factor(ordinal_output$SPECIES,
                                   levels = unique(order_plot$SPECIES))
  
  ordinal_output <- ordinal_output |>
    mutate(signif = case_when(
      .data[[q2_col]]  > 0 ~ "more flexible",
      .data[[q97_col]] < 0 ~ "less flexible",
      .default = "ns"
    ))
  
  ordinal_output |>
    ggplot(aes(y = SPECIES, x = .data[[est_col]], col = signif)) +
    geom_errorbarh(aes(xmin = .data[[q2_col]], xmax = .data[[q97_col]])) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = paste0(param, " Slope (log-odds)"), y = NULL,
         title = paste0("Species-specific ", param, " effects on symbiont switching")) +
    scale_color_manual(values = c("more flexible" = "red3",
                                  "less flexible" = "blue3",
                                  "ns"            = "gray")) +
    theme_bw()
}

# ------------------------------------------------------------------------------
# plot_species_effects_grid()
# Wrap multiple plot_species_effects() calls into a patchwork grid.
# ------------------------------------------------------------------------------
plot_species_effects_grid <- function(model, params) {
  plots <- lapply(params, function(p) plot_species_effects(model, p))
  wrap_plots(plots, ncol = length(params)) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

# ------------------------------------------------------------------------------
# graft_missing_species()
# Graft taxa absent from `tree` onto it using `ref_tree` for distance-informed
# placement. Falls back to genus crown if the species is not in ref_tree.
# ------------------------------------------------------------------------------
graft_missing_species <- function(tree, missing, ref_tree) {
  ref_dist <- cophenetic(ref_tree)
  shared   <- intersect(tree$tip.label, ref_tree$tip.label)
  
  for (sp in missing) {
    genus              <- sub("_.*$", "", sp)
    congeners_in_tree  <- tree$tip.label[grepl(paste0("^", genus, "_"), tree$tip.label)]
    
    if (!sp %in% ref_tree$tip.label) {
      message("Not in reference tree: ", sp, " — grafting to genus crown")
      graft_node <- if (length(congeners_in_tree) > 1) {
        getMRCA(tree, congeners_in_tree)
      } else if (length(congeners_in_tree) == 1) {
        which(tree$tip.label == congeners_in_tree[1])
      } else {
        length(tree$tip.label) + 1
      }
      bl   <- tree$edge.length[which(tree$edge[, 2] == graft_node)]
      tree <- bind.tip(tree, sp, where = graft_node,
                       edge.length = bl, position = 0.5 * bl)
      next
    }
    
    sp_dists       <- ref_dist[sp, shared]
    nearest        <- names(sort(sp_dists))[1:min(2, length(shared))]
    nearest_genera <- sub("_.*$", "", nearest)
    
    if (length(unique(nearest_genera)) > 1) {
      nearest <- nearest[1]
      message("Grafted ", sp, " near ", nearest, " (single nearest — cross-genus avoided)")
    } else {
      message("Grafted ", sp, " near ", paste(nearest, collapse = " & "))
    }
    if (!genus %in% sub("_.*$", "", nearest))
      message("  Note: cross-genus placement for ", sp)
    
    graft_node <- if (length(nearest) == 1) {
      which(tree$tip.label == nearest[1])
    } else {
      getMRCA(tree, nearest)
    }
    
    bl_ref    <- ref_tree$edge.length[which(ref_tree$edge[, 2] ==
                                              which(ref_tree$tip.label == sp))]
    bl_target <- tree$edge.length[which(tree$edge[, 2] == graft_node)]
    bl        <- ifelse(length(bl_ref) > 0, bl_ref[1], bl_target)
    tree      <- bind.tip(tree, sp, where = graft_node,
                          edge.length = bl, position = 0.5 * bl_target)
  }
  return(tree)
}

# ------------------------------------------------------------------------------
# dirichlet_soft_stanvars()
# Build a brms stanvar block that adds a Dirichlet(alpha_conc, alpha_conc,
# alpha_conc) soft penalty on the implied variance-partition simplex for
# coefficient coef_idx.
# sd ordering: sd_1 = phylogenetic, sd_2 = residual, sd_3 = trait.
# ------------------------------------------------------------------------------
dirichlet_soft_stanvars <- function(coef_idx, alpha_conc = 1) {
  stanvar(scode = glue::glue("
    {{
      // sd_1 = phylogenetic (clean.spec), sd_2 = residual (clean.spec_resid), sd_3 = trait (clean.spec_trait)
      real total_h2_{coef_idx} = square(sd_1[{coef_idx}]) + square(sd_2[{coef_idx}]) + square(sd_3[{coef_idx}]);
      vector[3] pi_implied_{coef_idx};
      pi_implied_{coef_idx}[1] = square(sd_1[{coef_idx}]) / total_h2_{coef_idx};  // phylogenetic
      pi_implied_{coef_idx}[2] = square(sd_2[{coef_idx}]) / total_h2_{coef_idx};  // residual
      pi_implied_{coef_idx}[3] = square(sd_3[{coef_idx}]) / total_h2_{coef_idx};  // trait
      target += dirichlet_lpdf(pi_implied_{coef_idx} | rep_vector({alpha_conc}, 3));
    }}
  "), block = "model")
}

# ------------------------------------------------------------------------------
# compute_var_partition_slopes()
# Sample-based variance contribution of a grouping term (with random slopes)
# to the linear predictor variance across observations, per posterior draw.
# ------------------------------------------------------------------------------
compute_var_partition_slopes <- function(re_array, group_var, data,
                                         slope_vars = c("PC1", "PC2", "PC3")) {
  n_draws   <- dim(re_array)[1]
  group_idx <- match(data[[group_var]], dimnames(re_array)[[2]])
  X_slopes  <- as.matrix(data[, slope_vars])
  purrr::map_dbl(1:n_draws, function(d) {
    b0   <- re_array[d, group_idx, "Intercept"]
    bslp <- matrix(re_array[d, group_idx, slope_vars],
                   nrow = length(group_idx), ncol = length(slope_vars))
    var(b0 + rowSums(bslp * X_slopes))
  })
}

# ------------------------------------------------------------------------------
# compute_var_partition_intercept()
# Sample-based variance contribution of an intercept-only grouping term
# (e.g., STUDY, COLONY_ID) per posterior draw.
# ------------------------------------------------------------------------------
compute_var_partition_intercept <- function(re_array, group_var, data) {
  n_draws   <- dim(re_array)[1]
  group_idx <- match(data[[group_var]], dimnames(re_array)[[2]])
  purrr::map_dbl(1:n_draws, function(d) var(re_array[d, group_idx, "Intercept"]))
}

# ------------------------------------------------------------------------------
# compute_var_partition_single_coef()
# Variance contribution of a single coefficient (intercept or slope) from one
# grouping term, per posterior draw.
# ------------------------------------------------------------------------------
compute_var_partition_single_coef <- function(re_array, group_var, data, coef_name) {
  n_draws   <- dim(re_array)[1]
  group_idx <- match(data[[group_var]], dimnames(re_array)[[2]])
  if (coef_name == "Intercept") {
    purrr::map_dbl(1:n_draws, function(d) var(re_array[d, group_idx, "Intercept"]))
  } else {
    x_vals <- data[[coef_name]]
    purrr::map_dbl(1:n_draws, function(d) var(re_array[d, group_idx, coef_name] * x_vals))
  }
}

# ------------------------------------------------------------------------------
# compute_species_switching()
# Compute expected switching intensity (E[shift_ord]) for each species under
# min/mean/max values of a chosen PC, holding all others at species mean.
# Requires objects in the calling environment:
#   b_mat, int_mat, r1_int, r1_pc1, r1_pc2, r1_pc3,
#   r2_int, r2_pc1, r2_pc2, r2_pc3, r3_int, r3_pc1, r3_pc2, r3_pc3,
#   stan_dat, species_levels, category_values
# ------------------------------------------------------------------------------
compute_species_switching <- function(spec_df, pc_col, pc_level_name,
                                      vary_pc = c("PC1", "PC2", "PC3"),
                                      ndraws = 500) {
  vary_pc           <- match.arg(vary_pc)
  draw_idx          <- sample(nrow(b_mat), ndraws)
  b_log_delta_t_idx <- which(colnames(stan_dat_delta$X) == "log_delta_t")
  
  purrr::map_dfr(seq_along(species_levels), function(sp_idx) {
    sp     <- species_levels[sp_idx]
    sp_row <- spec_df[spec_df$clean.spec == sp, ]
    if (nrow(sp_row) == 0) return(NULL)
    
    pc1 <- sp_row$PC1_mean
    pc2 <- sp_row$PC2_mean
    pc3 <- sp_row$PC3_mean
    
    if (vary_pc == "PC1") pc1 <- sp_row[[pc_col]]
    if (vary_pc == "PC2") pc2 <- sp_row[[pc_col]]
    if (vary_pc == "PC3") pc3 <- sp_row[[pc_col]]
    
    mu <- r1_int[draw_idx, sp_idx] + r1_pc1[draw_idx, sp_idx] * pc1 +
      r1_pc2[draw_idx, sp_idx] * pc2 + r1_pc3[draw_idx, sp_idx] * pc3 +
      r2_int[draw_idx, sp_idx] + r2_pc1[draw_idx, sp_idx] * pc1 +
      r2_pc2[draw_idx, sp_idx] * pc2 + r2_pc3[draw_idx, sp_idx] * pc3 +
      r3_int[draw_idx, sp_idx] + r3_pc1[draw_idx, sp_idx] * pc1 +
      r3_pc2[draw_idx, sp_idx] * pc2 + r3_pc3[draw_idx, sp_idx] * pc3 +
      b_mat[draw_idx, b_log_delta_t_idx] *
      (log(360) - mean(stan_dat_delta$X[, b_log_delta_t_idx]))
    
    thres    <- int_mat[draw_idx, ]
    e_switch <- purrr::map_dbl(seq_along(draw_idx), function(i) {
      cum_p <- plogis(as.numeric(thres[i, ]) - mu[i])
      cat_p <- c(cum_p, 1) - c(0, cum_p)
      sum(cat_p * category_values)
    })
    
    tibble(clean.spec = sp, pc_level = pc_level_name,
           e_switch = e_switch, .draw = draw_idx)
  })
}

# ------------------------------------------------------------------------------
# summarise_switching()
# Summarise compute_species_switching() output to posterior median + 95% CrI.
# ------------------------------------------------------------------------------
summarise_switching <- function(pred_df) {
  pred_df |>
    group_by(clean.spec, pc_level) |>
    summarise(
      prob     = mean(e_switch),
      prob_low = quantile(e_switch, 0.025),
      prob_hi  = quantile(e_switch, 0.975),
      .groups  = "drop"
    )
}

# ------------------------------------------------------------------------------
# plot_switching_tree()
# Overlay predicted switching intensity onto a phylogenetic or trait tree.
# pred_pc3: output of summarise_switching() with columns clean.spec, pc_level,
#           prob (pc_level used for x-axis offsets via `levs`).
# ------------------------------------------------------------------------------
plot_switching_tree <- function(tree, pred_pc3, fill_limits, title = "",
                                levs = c(10, 20, 30), textoffset = 40) {
  tree_layout <- ggtree(tree)$data |>
    filter(isTip) |>
    select(label, x, y)
  
  pred_positioned <- pred_pc3 |>
    rename(label = clean.spec) |>
    left_join(tree_layout, by = "label") |>
    mutate(x_offset = x + case_when(
      pc_level == "min"  ~ levs[1],
      pc_level == "mean" ~ levs[2],
      pc_level == "max"  ~ levs[3]
    )) 
  
  ggtree(tree, layout = "rectangular") +
    geom_point(data = pred_positioned,
               aes(x = x_offset, y = y, fill = prob),
               size = 2.7, color = "gray30", stroke = 0.3,
               shape = 22, inherit.aes = FALSE) +
    scale_fill_viridis_c(trans = "sqrt", option = "viridis",
                         limits = fill_limits) +
    geom_tiplab(aes(label = gsub("(\\w)\\w+_(\\w+)", "\\1. \\2", label)),
                fontface = "italic", size = 2.5, offset = textoffset) +
    ggtitle(title) +
    theme(legend.position = "right") +
    hexpand(0.3)
}

# ------------------------------------------------------------------------------
# add_epred_draws_manual()
# Manual posterior prediction for the cmdstanr Dirichlet model, replicating
# add_epred_draws() from tidybayes but using the raw Stan parameter arrays.
# Requires stan_dat to be in the calling environment.
# ------------------------------------------------------------------------------
add_epred_draws_manual <- function(newdata, b_draws, intercept_draws,
                                   re_draws, means_X, species_levels,
                                   category_values, ndraws = 500) {
  X_new  <- as.matrix(newdata[, colnames(stan_dat_delta$X)])
  Xc_new <- sweep(X_new, 2, means_X, "-")
  draw_idx <- sample(nrow(b_draws), ndraws)
  n_obs    <- nrow(newdata)
  
  purrr::map_dfr(seq_along(draw_idx), function(i) {
    d  <- draw_idx[i]
    mu <- as.numeric(Xc_new %*% b_draws[d, ])
    
    for (n in seq_len(n_obs)) {
      sp_name <- as.character(newdata$clean.spec[n])
      sp_idx  <- which(species_levels == sp_name)
      if (length(sp_idx) == 0) next
      
      pc1_val <- newdata$PC1[n]; pc2_val <- newdata$PC2[n]; pc3_val <- newdata$PC3[n]
      
      mu[n] <- mu[n] +
        as.numeric(re_draws[d, paste0("r_1_1[", sp_idx, "]")]) +
        as.numeric(re_draws[d, paste0("r_1_2[", sp_idx, "]")]) * pc1_val +
        as.numeric(re_draws[d, paste0("r_1_3[", sp_idx, "]")]) * pc2_val +
        as.numeric(re_draws[d, paste0("r_1_4[", sp_idx, "]")]) * pc3_val +
        as.numeric(re_draws[d, paste0("r_2_1[", sp_idx, "]")]) +
        as.numeric(re_draws[d, paste0("r_2_2[", sp_idx, "]")]) * pc1_val +
        as.numeric(re_draws[d, paste0("r_2_3[", sp_idx, "]")]) * pc2_val +
        as.numeric(re_draws[d, paste0("r_2_4[", sp_idx, "]")]) * pc3_val +
        as.numeric(re_draws[d, paste0("r_3_1[", sp_idx, "]")]) +
        as.numeric(re_draws[d, paste0("r_3_2[", sp_idx, "]")]) * pc1_val +
        as.numeric(re_draws[d, paste0("r_3_3[", sp_idx, "]")]) * pc2_val +
        as.numeric(re_draws[d, paste0("r_3_4[", sp_idx, "]")]) * pc3_val
    }
    
    thres        <- as.numeric(intercept_draws[d, ])
    cum_prob_mat <- plogis(outer(thres, mu, "-"))
    cat_prob_mat <- rbind(cum_prob_mat, 1) - rbind(0, cum_prob_mat)
    e_val        <- colSums(cat_prob_mat * category_values)
    
    tibble(.draw = d, obs_idx = seq_len(n_obs),
           clean.spec = as.character(newdata$clean.spec),
           .epred = e_val,
           .epred_cat1 = cat_prob_mat[1, ], .epred_cat2 = cat_prob_mat[2, ],
           .epred_cat3 = cat_prob_mat[3, ], .epred_cat4 = cat_prob_mat[4, ]) |>
      bind_cols(newdata |> select(any_of(c("PC1", "PC2", "PC3",
                                           "log_delta_t", "last_genus"))))
  })
}

# ------------------------------------------------------------------------------
# rordered_logistic_vec()
# Vectorised draw from an ordered logistic distribution.
# Used for posterior predictive checks with the cmdstanr Dirichlet model.
# ------------------------------------------------------------------------------
rordered_logistic_vec <- function(mu, thres) {
  n_thres      <- length(thres)
  cum_prob_mat <- plogis(outer(mu, thres, function(m, t) t - m))
  cat_prob_mat <- cbind(cum_prob_mat, 1) - cbind(0, cum_prob_mat)
  apply(cat_prob_mat, 1, function(p) sample(n_thres + 1, 1, prob = p))
}


# ------------------------------------------------------------------------------
# extract_re()
#  Extract n_spec-column matrix for one RE block + coefficient slot
# ------------------------------------------------------------------------------
extract_re <- function(block, coef_slot) {
  cols <- paste0("r_", block, "_", coef_slot, "[", 1:n_spec, "]")
  missing <- setdiff(cols, names(draws_df))
  if (length(missing) > 0) stop("Missing columns: ", paste(missing[1:3], collapse = ", "), "...")
  as.matrix(draws_df[, cols])
}


# ------------------------------------------------------------------------------
# compute_species_switching_genus()
# Same as compute_species_switching(), but adds each species' modal last_genus
# fixed effect into mu, rather than silently using the reference level for
# every species.
# ------------------------------------------------------------------------------
compute_species_switching_genus <- function(spec_df, pc_col, pc_level_name,
                                            vary_pc = c("PC1", "PC2", "PC3"),
                                            ndraws = 500) {
  vary_pc           <- match.arg(vary_pc)
  draw_idx          <- sample(nrow(b_mat), ndraws)
  b_log_delta_t_idx <- which(colnames(stan_dat_delta$X) == "log_delta_t")
  
  purrr::map_dfr(seq_along(species_levels), function(sp_idx) {
    sp     <- species_levels[sp_idx]
    sp_row <- spec_df[spec_df$clean.spec == sp, ]
    if (nrow(sp_row) == 0) return(NULL)
    
    pc1 <- sp_row$PC1_mean
    pc2 <- sp_row$PC2_mean
    pc3 <- sp_row$PC3_mean
    
    if (vary_pc == "PC1") pc1 <- sp_row[[pc_col]]
    if (vary_pc == "PC2") pc2 <- sp_row[[pc_col]]
    if (vary_pc == "PC3") pc3 <- sp_row[[pc_col]]
    
    # Look up the b_mat column for this species' modal genus, if it isn't
    # the reference level (reference level contributes 0, column won't exist)
    genus_col <- paste0("last_genus", sp_row$modal_genus)
    genus_effect <- if (genus_col %in% colnames(b_mat)) {
      b_mat[draw_idx, genus_col]
    } else {
      0
    }
    
    mu <- r1_int[draw_idx, sp_idx] + r1_pc1[draw_idx, sp_idx] * pc1 +
      r1_pc2[draw_idx, sp_idx] * pc2 + r1_pc3[draw_idx, sp_idx] * pc3 +
      r2_int[draw_idx, sp_idx] + r2_pc1[draw_idx, sp_idx] * pc1 +
      r2_pc2[draw_idx, sp_idx] * pc2 + r2_pc3[draw_idx, sp_idx] * pc3 +
      r3_int[draw_idx, sp_idx] + r3_pc1[draw_idx, sp_idx] * pc1 +
      r3_pc2[draw_idx, sp_idx] * pc2 + r3_pc3[draw_idx, sp_idx] * pc3 +
      b_mat[draw_idx, b_log_delta_t_idx] *
      (log(360) - mean(stan_dat_delta$X[, b_log_delta_t_idx])) +
      genus_effect
    
    thres    <- int_mat[draw_idx, ]
    e_switch <- purrr::map_dbl(seq_along(draw_idx), function(i) {
      cum_p <- plogis(as.numeric(thres[i, ]) - mu[i])
      cat_p <- c(cum_p, 1) - c(0, cum_p)
      sum(cat_p * category_values)
    })
    
    tibble(clean.spec = sp, pc_level = pc_level_name,
           e_switch = e_switch, .draw = draw_idx)
  })
}


compute_pc1_switching_genus <- function(spec_df, pc_col, pc_level_name,
                                        ndraws = 500) {
  draw_idx          <- sample(nrow(b_mat), ndraws)
  b_log_delta_t_idx <- which(colnames(stan_dat_delta$X) == "log_delta_t")
  
  purrr::map_dfr(seq_along(species_levels), function(sp_idx) {
    sp     <- species_levels[sp_idx]
    sp_row <- spec_df[spec_df$clean.spec == sp, ]
    if (nrow(sp_row) == 0) return(NULL)
    
    pc1 <- sp_row[[pc_col]]
    pc2 <- sp_row$PC2_mean
    pc3 <- sp_row$PC3_mean
    
    genus_col <- paste0("last_genus", sp_row$modal_genus)
    genus_effect <- if (genus_col %in% colnames(b_mat)) {
      b_mat[draw_idx, genus_col]
    } else {
      0
    }
    
    mu <- r1_int[draw_idx, sp_idx] + r1_pc1[draw_idx, sp_idx] * pc1 +
      r1_pc2[draw_idx, sp_idx] * pc2 + r1_pc3[draw_idx, sp_idx] * pc3 +
      r2_int[draw_idx, sp_idx] + r2_pc1[draw_idx, sp_idx] * pc1 +
      r2_pc2[draw_idx, sp_idx] * pc2 + r2_pc3[draw_idx, sp_idx] * pc3 +
      r3_int[draw_idx, sp_idx] + r3_pc1[draw_idx, sp_idx] * pc1 +
      r3_pc2[draw_idx, sp_idx] * pc2 + r3_pc3[draw_idx, sp_idx] * pc3 +
      b_mat[draw_idx, b_log_delta_t_idx] *
      (log(360) - mean(stan_dat_delta$X[, b_log_delta_t_idx])) +
      genus_effect
    
    thres    <- int_mat[draw_idx, ]
    e_switch <- purrr::map_dbl(seq_along(draw_idx), function(i) {
      cum_p <- plogis(as.numeric(thres[i, ]) - mu[i])
      cat_p <- c(cum_p, 1) - c(0, cum_p)
      sum(cat_p * category_values)
    })
    
    tibble(clean.spec = sp, pc_level = pc_level_name,
           e_switch = e_switch, .draw = draw_idx)
  })
}
