#########################################################################
# SHINY APP: CORAL SYMBIONT SWITCHING EXPLORER
#########################################################################


#==============================================================================
#Development Note: I used Claude's Sonnet 5 model to develop the Shiny app
#described in this git repository. I manually generated and revised all scripts
#leading up to this point, and this script draws directly from those curated outputs.
#That said, this Shiny app data should be considered exploratory and an extension
#of the results of the main analysis; the use of AI was intended to increase 
#interpretability of results and to bridge my knowledge gap in building a Shiny app.
#===============================================================================

# Interactive front-end for the coral symbiont switching meta-analysis
# DELTA model. Lets users:
#   1. Enter how much each of several NON-REDUNDANT environmental
#      variables changed (delta), which get converted directly to
#      delta_PC1/2/3 -- the actual inputs the fitted delta model takes --
#      via a single linear projection (see compute_delta_pc() below), and
#      get predicted switching probability for one or many species.
#   2. Prune the phylogenetic and trait trees to species of interest,
#      optionally colored by predicted switching probability.
#   3. Explore where a given site+year (environmental PCA) or species
#      (trait PCA) falls in PCA space, and see the user's entered
#      scenario as an arrow from the dataset average (delta = 0) outward.
#
# DELTA-MODEL ARCHITECTURE NOTE: earlier drafts of this app computed
# delta_PC1/2/3 by projecting a full STARTING environmental vector and a
# full ENDING vector (each filled in with fallback means for anything not
# specified) and then subtracting the two projections. Since PCA
# projection is linear, this is mathematically equivalent to projecting
# the vector of DIFFERENCES directly -- and doing it that way means:
#   - No fallback means are needed at all for delta_PC (a variable held
#     constant across "start" and "end" always contributes exactly zero
#     to a delta regardless of which constant it's held at -- so a
#     species-mean vs. dataset-mean fallback toggle can never do
#     anything different for this calculation. Removed accordingly.)
#   - Auto-filling a correlated variable only ever needs its linear
#     SLOPE against whichever selected variable it's proxying (any
#     intercept cancels exactly in a difference).
#   - The user-facing inputs are exactly the NON-REDUNDANT variables
#     (app_data$selected_env_vars, built by the prep script's greedy
#     0.7-correlation reduction) -- anything correlated with one of
#     those above the same 0.7 threshold (app_data$redundant_env_vars)
#     is auto-set from it; everything else is held at delta = 0 (i.e.
#     genuinely unspecified/unchanged), which is now explicit rather
#     than an emergent side effect of a fallback-mean cancellation.
#
# Requires data/shiny_app_data_delta.RDS, produced by
# scripts/99b_prep_shiny_data_delta.R (run that script once after the
# delta model has been fit). No changes to the prep script are needed
# for this version -- selected_env_vars/redundant_env_vars already exist
# from the original 0.7-threshold reduction.
#########################################################################

library(shiny)
library(tidyverse)
library(vegan)
library(ape)
library(ggtree)
library(bslib)

source("scripts/helpers.R")  # for plot_switching_tree()

app_data <- readRDS("data/shiny_app_data_delta.RDS")

env_vars       <- names(app_data$env_raw_means)  # ALL raw variables (needed for the full delta vector)
env_vars_shown <- app_data$selected_env_vars      # non-redundant set (|r| < 0.7 with each other) -- user-facing

# Friendly labels (with units) for the raw environmental inputs.
# NOTE on precip units: assumed mm/day based on the IMERG_monthly_global_precip
# product used in 01_EnvironmentalData.R -- double check against your rerddap
# query if you want to confirm the exact unit that product reports in.
env_var_labels <- c(
  mean_SST    = "Mean SST (\u00b0C)",         min_SST    = "Min SST (\u00b0C)",
  max_SST     = "Max SST (\u00b0C)",          sd_SST     = "SST SD (\u00b0C)",
  mean_DHW    = "Mean Degree Heating Weeks (\u00b0C-weeks)",
  mean_ANOM   = "Mean SST Anomaly (\u00b0C)", sd_ANOM    = "SST Anomaly SD (\u00b0C)",
  min_ANOM    = "Min SST Anomaly (\u00b0C)",  max_ANOM   = "Max SST Anomaly (\u00b0C)",
  cum_ANOM    = "Cumulative SST Anomaly (\u00b0C)",
  mean_precip = "Mean Precipitation (mm/day)", min_precip = "Min Precipitation (mm/day)",
  max_precip  = "Max Precipitation (mm/day)",  sd_precip  = "Precipitation SD (mm/day)",
  mean_chla   = "Mean Chlorophyll-a (mg/m\u00b3)", min_chla = "Min Chlorophyll-a (mg/m\u00b3)",
  max_chla    = "Max Chlorophyll-a (mg/m\u00b3)",  sd_chla  = "Chlorophyll-a SD (mg/m\u00b3)",
  mean_turb   = "Mean Turbidity, Kd490 (m\u207b\u00b9)", min_turb = "Min Turbidity (m\u207b\u00b9)",
  max_turb    = "Max Turbidity (m\u207b\u00b9)",         sd_turb  = "Turbidity SD (m\u207b\u00b9)"
)
env_var_labels <- env_var_labels[env_vars]  # keep only vars actually in the model

#### ----------- CORE PREDICTION FUNCTIONS ----------- ####

# Projects a vector of RAW-VARIABLE CHANGES (deltas) directly to
# delta_PC1/2/3, using the linearity of PCA: delta_PC = rotation^T *
# (delta_x / env_raw_sds). No fallback means are needed -- a variable
# that isn't specified and isn't auto-filled from a correlated selected
# variable simply gets delta_x[v] = 0 (held constant / unchanged), which
# contributes nothing to delta_PC regardless of what its absolute value
# would have been.
compute_delta_pc <- function(driving_deltas, app_data) {
  # driving_deltas: named list/vector over env_vars_shown, e.g.
  # list(mean_SST = 1.5, sd_precip = -5, max_ANOM = 0.8, ...). Any entry
  # left NA/blank is treated as "no change" (0).
  full_vars <- names(app_data$env_raw_means)
  delta_x <- setNames(rep(0, length(full_vars)), full_vars)
  
  provided <- intersect(names(driving_deltas), full_vars)
  provided <- provided[!is.na(driving_deltas[provided]) & driving_deltas[provided] != ""]
  if (length(provided) > 0) {
    delta_x[provided] <- as.numeric(unlist(driving_deltas[provided]))
  }
  
  # Auto-fill correlated variables: delta_v = slope * delta_proxy. Only
  # the slope is needed -- the intercept in "v ~ proxy" cancels exactly
  # in a difference. Anything neither provided nor auto-filled stays at
  # delta_x[v] = 0 (explicitly held constant).
  redundant <- app_data$redundant_env_vars
  used_proxy_fill <- character(0)
  if (!is.null(redundant) && nrow(redundant) > 0) {
    for (i in seq_len(nrow(redundant))) {
      v <- redundant$variable[i]
      if (v %in% provided) next  # never overwrite a variable the user explicitly set
      proxy <- redundant$proxy[i]
      delta_x[v] <- redundant$slope[i] * delta_x[proxy]
      used_proxy_fill <- c(used_proxy_fill, v)
    }
  }
  
  held_constant <- setdiff(full_vars, c(provided, used_proxy_fill))
  
  delta_z  <- delta_x / app_data$env_raw_sds
  delta_pc <- as.numeric(delta_z %*% app_data$env.pc$rotation)
  names(delta_pc) <- colnames(app_data$env.pc$rotation)
  
  list(delta_pc = delta_pc, used_proxy_fill = used_proxy_fill, held_constant = held_constant)
}

# NOTE ON SCOPE: this adds the species-level random effects (phylogenetic +
# residual + trait) to a log_delta_t fixed-effect term, following the same
# approach as compute_species_switching() in helpers.R. It does NOT add the
# population-level (fixed-effect) coefficients for PC1/PC2/PC3 themselves --
# only the species-specific deviations. This matches your existing
# E_Switch_Phylo/Trait_global.tiff figures for consistency, but flagging
# again here in case that omission in the original helper was unintentional.
#
# DELTA-MODEL NOTE: pc1/pc2/pc3 below are delta_PC1/delta_PC2/delta_PC3
# (from compute_delta_pc() above) rather than raw PC1/2/3 -- the function
# body itself is unchanged from earlier versions, since the underlying
# math (species' random slope times whatever PC value it's handed) is
# identical either way; only the caller's interpretation has changed.

predict_switching <- function(pc1, pc2, pc3, species, app_data, delta_days = 365, genus_choice = "default") {
  
  sp_idx <- match(species, app_data$species_levels)
  valid  <- !is.na(sp_idx)
  
  if (!any(valid)) return(NULL)
  
  species  <- species[valid]
  sp_idx   <- sp_idx[valid]
  n_draws  <- nrow(app_data$b_mat)
  
  # pc1/pc2/pc3 will normally be a single shared scenario here (one delta
  # scenario applies to every species) -- recycled to per-species length
  # so the loop below always indexes explicitly.
  if (length(pc1) == 1) pc1 <- rep(pc1, length(species))
  if (length(pc2) == 1) pc2 <- rep(pc2, length(species))
  if (length(pc3) == 1) pc3 <- rep(pc3, length(species))
  
  log_dt_offset <- log(delta_days) - app_data$mean_log_delta_t
  b_log_dt      <- if ("log_delta_t" %in% colnames(app_data$b_mat)) app_data$b_mat[, "log_delta_t"] else 0
  
  # Which genus to condition on for each species: either that species'
  # most-observed genus (genus_choice == "default"), or one fixed genus
  # (A/B/C/D) applied to every species.
  genus_per_species <- if (identical(genus_choice, "default")) {
    lookup <- app_data$species_default_genus
    lookup$default_genus[match(species, lookup$clean.spec)]
  } else {
    rep(genus_choice, length(species))
  }
  
  # brms centers the ENTIRE non-intercept design matrix (Xc <- sweep(X, 2,
  # colMeans(X), "-")), including all last_genus dummy columns together --
  # not just "the selected one". So the genus contribution to mu has to be
  # the dot product across ALL dummy columns (raw - colMean), which is
  # nonzero even for the reference level (its dummy is 0 everywhere, but
  # colMean isn't, so centered value is -colMean). Handling only the
  # "selected" column and zeroing the rest would silently bias predictions.
  genus_dummy_names <- names(app_data$genus_col_means)
  
  genus_contribution <- function(genus) {
    raw_dummy <- setNames(rep(0, length(genus_dummy_names)), genus_dummy_names)
    target_col <- paste0("last_genus", genus)
    if (target_col %in% genus_dummy_names) {
      raw_dummy[target_col] <- 1
    } else if (!identical(as.character(genus), app_data$genus_reference_level)) {
      warning("Genus '", genus, "' not among fitted last_genus levels (",
              paste(app_data$genus_levels_available, collapse = ", "), ") -- treating as reference level.")
    }
    centered <- raw_dummy - app_data$genus_col_means[genus_dummy_names]
    as.numeric(app_data$b_mat[, genus_dummy_names, drop = FALSE] %*% centered)
  }
  
  # Flag a species' scenario as outside what the model actually saw for it
  # if PC1/2/3 falls outside that species' observed range (its own
  # observations unioned with its 25km reef-neighborhood's observations,
  # from the prep script -- see species_pc_range there for why).
  flag_out_of_range <- function(sp, p1, p2, p3) {
    rng <- app_data$species_pc_range[app_data$species_pc_range$clean.spec == sp, ]
    if (nrow(rng) == 0) return("No range data for this species")
    outside <- c(
      if (p1 < rng$PC1_min || p1 > rng$PC1_max) "PC1",
      if (p2 < rng$PC2_min || p2 > rng$PC2_max) "PC2",
      if (p3 < rng$PC3_min || p3 > rng$PC3_max) "PC3"
    )
    if (length(outside) == 0) "" else paste0("Outside observed range (", paste(outside, collapse = ", "), ")")
  }
  
  purrr::map_dfr(seq_along(species), function(i) {
    sp <- species[i]; si <- sp_idx[i]
    b_genus <- genus_contribution(genus_per_species[i])
    p1 <- pc1[i]; p2 <- pc2[i]; p3 <- pc3[i]
    
    mu <- app_data$r1_int[, si] + app_data$r1_pc1[, si] * p1 +
      app_data$r1_pc2[, si] * p2 + app_data$r1_pc3[, si] * p3 +
      app_data$r2_int[, si] + app_data$r2_pc1[, si] * p1 +
      app_data$r2_pc2[, si] * p2 + app_data$r2_pc3[, si] * p3 +
      app_data$r3_int[, si] + app_data$r3_pc1[, si] * p1 +
      app_data$r3_pc2[, si] * p2 + app_data$r3_pc3[, si] * p3 +
      b_log_dt * log_dt_offset + b_genus
    
    # cumulative logit: P(Y <= k) = plogis(threshold_k - mu), per draw
    thres  <- app_data$int_mat
    cum_p  <- plogis(thres - matrix(mu, nrow = n_draws, ncol = ncol(thres)))
    cat_p  <- cbind(cum_p, 1) - cbind(0, cum_p)
    e_switch <- as.numeric(cat_p %*% app_data$category_values)
    
    tibble(
      clean.spec = sp,
      genus_used = genus_per_species[i],
      flag       = flag_out_of_range(sp, p1, p2, p3),
      draw       = seq_len(n_draws),
      e_switch   = e_switch,
      p_none          = cat_p[, 1],
      p_subtype_switch = cat_p[, 2],
      p_minor         = cat_p[, 3],
      p_major         = cat_p[, 4]
    )
  })
}

summarise_predictions <- function(pred_df) {
  pred_df |>
    group_by(clean.spec, genus_used, flag) |>
    summarise(
      e_switch_median = median(e_switch),
      e_switch_low    = quantile(e_switch, 0.025),
      e_switch_high   = quantile(e_switch, 0.975),
      p_none          = mean(p_none),
      p_subtype_switch = mean(p_subtype_switch),
      p_minor         = mean(p_minor),
      p_major         = mean(p_major),
      .groups = "drop"
    ) |>
    arrange(desc(e_switch_median))
}

#### ----------- UI ----------- ####

ui <- page_navbar(
  title = "Coral Symbiont Switching Explorer",
  theme = bs_theme(bootswatch = "flatly"),
  
  nav_panel(
    "Switching Calculator",
    layout_sidebar(
      sidebar = sidebar(
        width = 350,
        h5("Environmental change scenario"),
        p(class = "text-muted small",
          "Enter how much each variable below changed (end minus start) -- the model responds ",
          "directly to this change, not to any absolute starting or ending value. Leave a field ",
          "blank to assume no change (\u0394 = 0) for that variable. Only variables NOT strongly ",
          "correlated with each other (|r| < 0.7 in the historical dataset) are shown here -- ",
          "closely related variables (e.g. other SST metrics) shift automatically based on these ",
          "(see 'Show variables that are auto-set' below). For scale, each field's typical ",
          "dataset-wide standard deviation is shown as a rough reference for what a 'substantial' ",
          "change might look like."),
        purrr::imap(env_var_labels[env_vars_shown], function(lab, var) {
          sd_ref <- app_data$env_raw_sds[var]
          div(
            numericInput(paste0("delta_", var), paste0("\u0394 ", lab), value = NA),
            p(class = "text-muted", style = "font-size: 0.75em; margin-top: -12px;",
              sprintf("(dataset-wide SD for reference: %.2f)", sd_ref))
          )
        }),
        hr(),
        checkboxInput("show_autoset", "Show variables that are auto-set (and why)", value = FALSE),
        conditionalPanel("input.show_autoset == true", tableOutput("autoset_table")),
        hr(),
        numericInput("delta_days", "Prediction window (days)", value = 365, min = 1),
        hr(),
        h5("Symbiont genus"),
        radioButtons("genus_mode", NULL,
                     choices = c("Default (each species' most common genus)" = "default",
                                 "Fixed genus for all species" = "fixed"),
                     selected = "default"),
        conditionalPanel(
          "input.genus_mode == 'fixed'",
          selectInput("genus_fixed", "Genus", choices = app_data$genus_levels_available)
        ),
        hr(),
        selectizeInput("species_calc", "Species (leave blank = all)",
                       choices = app_data$species_levels, multiple = TRUE),
        actionButton("calc_btn", "Calculate switching probability", class = "btn-primary")
      ),
      card(
        card_header("Predicted switching probability"),
        textOutput("defaults_used_msg"),
        div(style = "max-height: 600px; overflow-y: auto;",
            plotOutput("switch_plot")),
        div(style = "max-height: 350px; overflow-y: auto; margin-top: 10px;",
            tableOutput("switch_table"))
      )
    )
  ),
  
  nav_panel(
    "Trees",
    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        selectizeInput("species_tree", "Species to include (leave blank = all)",
                       choices = app_data$species_levels, multiple = TRUE),
        checkboxInput("color_by_switch", "Color tips by predicted switching (from Calculator tab)", value = TRUE),
        p(class = "text-muted small",
          "Run a calculation on the Switching Calculator tab first if you want tip coloring.")
      ),
      layout_columns(
        card(card_header("Phylogenetic tree"), plotOutput("phylo_tree_plot", height = "600px")),
        card(card_header("Trait dendrogram"), plotOutput("trait_tree_plot", height = "600px"))
      )
    )
  ),
  
  nav_panel(
    "PCA Explorer",
    layout_columns(
      col_widths = c(8, 4),
      card(
        card_header("Environmental PCA (site \u00d7 year)"),
        selectizeInput("region_filter", "Filter by region",
                       choices = c("All regions", sort(unique(app_data$site_year_lookup$basin_region))),
                       selected = "All regions"),
        selectizeInput("site_year_select", "Highlight a site + year", choices = NULL, selected = "(none)"),
        checkboxInput("show_custom_point", "Show my \u0394 scenario (Calculator tab)", value = FALSE),
        plotOutput("env_pca_plot", height = "750px")
      ),
      card(
        card_header("Trait PCA (species)"),
        selectizeInput("trait_species_select", "Highlight a species",
                       choices = c("(none)", app_data$species_levels), selected = "(none)"),
        plotOutput("trait_pca_plot", height = "500px")
      )
    )
  )
)

#### ----------- SERVER ----------- ####

server <- function(input, output, session) {
  
  last_prediction <- reactiveVal(NULL)
  last_delta_pc   <- reactiveVal(NULL)
  
  observeEvent(input$calc_btn, {
    species <- if (length(input$species_calc) == 0) app_data$species_levels else input$species_calc
    genus_choice <- if (identical(input$genus_mode, "fixed")) input$genus_fixed else "default"
    
    driving_deltas <- setNames(
      lapply(env_vars_shown, function(v) input[[paste0("delta_", v)]]),
      env_vars_shown
    )
    
    d <- compute_delta_pc(driving_deltas, app_data)
    last_delta_pc(d$delta_pc)
    
    pred <- predict_switching(
      pc1 = d$delta_pc["PC1"], pc2 = d$delta_pc["PC2"], pc3 = d$delta_pc["PC3"],
      species = species, app_data = app_data, delta_days = input$delta_days,
      genus_choice = genus_choice
    )
    last_prediction(summarise_predictions(pred))
    
    output$defaults_used_msg <- renderText({
      n_specified <- sum(!vapply(driving_deltas, function(x) is.null(x) || is.na(x), logical(1)))
      n_proxy     <- length(d$used_proxy_fill)
      n_constant  <- length(d$held_constant)
      if (n_specified == 0) {
        base_msg <- "No changes entered -- this is a \u0394PC \u2248 0 scenario (no environmental change)."
      } else {
        base_msg <- paste0(n_specified, " of ", length(env_vars_shown), " variables specified.")
      }
      paste0(base_msg, " (", n_proxy, " correlated variable", if (n_proxy != 1) "s" else "",
             " auto-set from these -- see 'Show variables that are auto-set' above; ",
             n_constant, " variables held constant/unspecified. Check the Flag column in the ",
             "results for any species evaluated outside the range of environmental CHANGE it's ",
             "actually been observed at.)")
    })
  })
  
  output$autoset_table <- renderTable({
    req(app_data$redundant_env_vars)
    app_data$redundant_env_vars %>%
      transmute(
        "Auto-set variable" = variable,
        "Based on" = proxy,
        "Correlation" = round(correlation, 2),
        "R\u00b2" = round(r_squared, 2)
      )
  }, digits = 2)
  
  output$switch_table <- renderTable({
    req(last_prediction())
    last_prediction() |>
      mutate(flag = ifelse(flag == "", "OK", flag)) |>
      rename(Species = clean.spec, Genus = genus_used, Flag = flag,
             `E[Switch] (median)` = e_switch_median,
             `2.5%` = e_switch_low, `97.5%` = e_switch_high,
             `P(none)` = p_none, `P(subtype switch)` = p_subtype_switch,
             `P(minor)` = p_minor, `P(major)` = p_major)
  }, digits = 2)
  
  output$switch_plot <- renderPlot({
    req(last_prediction())
    df <- last_prediction() |>
      mutate(clean.spec = factor(clean.spec, levels = rev(clean.spec)),
             flagged = flag != "")
    ggplot(df, aes(x = e_switch_median, y = clean.spec)) +
      geom_errorbarh(aes(xmin = e_switch_low, xmax = e_switch_high, color = flagged), height = 0.3, alpha = 0.7) +
      geom_point(aes(color = flagged), size = 3) +
      scale_color_manual(values = c(`FALSE` = "#2166ac", `TRUE` = "#c0392b"),
                         labels = c(`FALSE` = "Within observed conditions", `TRUE` = "Outside observed conditions"),
                         name = NULL) +
      labs(x = "Predicted E[switching severity] (0 = none, 3 = major)", y = NULL) +
      theme_bw(base_size = 12) +
      theme(axis.text.y = element_text(face = "italic"), legend.position = "top")
  }, height = function() {
    n <- if (is.null(last_prediction())) 1 else nrow(last_prediction())
    # ~22px per species row, floor of 300px (a handful of species) so the
    # plot never goes stubby, no hard ceiling so it can grow to fit
    # everything -- the outer div below makes that scrollable rather than
    # pushing the whole page down.
    max(300, 22 * n)
  })
  
  #### ---- Trees tab ---- ####
  
  pruned_trees <- reactive({
    sp <- if (length(input$species_tree) == 0) app_data$species_levels else input$species_tree
    sp <- intersect(sp, app_data$species_levels)
    req(length(sp) >= 2)  # need at least 2 tips for a meaningful tree
    list(
      phylo = keep.tip(app_data$tree_final, intersect(sp, app_data$tree_final$tip.label)),
      trait = keep.tip(app_data$trait_tree, intersect(sp, app_data$trait_tree$tip.label))
    )
  })
  
  tip_colors <- reactive({
    if (!isTRUE(input$color_by_switch)) return(NULL)
    pred <- last_prediction()
    if (is.null(pred)) return(NULL)
    pred |> select(clean.spec, prob = e_switch_median)
  })
  
  output$phylo_tree_plot <- renderPlot({
    trees <- pruned_trees()
    cols  <- tip_colors()
    if (!is.null(cols)) {
      pred_df <- cols |> rename(clean.spec = clean.spec) |> mutate(pc_level = "mean")
      tree_depth   <- max(node.depth.edgelength(trees$phylo))
      box_offset   <- tree_depth * 0.05   # colored square: just past the tip
      label_offset <- tree_depth * 0.15   # species label: further out, past the box
      plot_switching_tree(trees$phylo, pred_df, fill_limits = range(cols$prob, na.rm = TRUE),
                          title = "Phylogeny", levs = rep(box_offset, 3), textoffset = label_offset)
    } else {
      plot(trees$phylo, main = "Phylogeny (pruned)", cex = 0.8)
    }
  })
  
  output$trait_tree_plot <- renderPlot({
    trees <- pruned_trees()
    cols  <- tip_colors()
    if (!is.null(cols)) {
      pred_df <- cols |> rename(clean.spec = clean.spec) |> mutate(pc_level = "mean")
      tree_depth   <- max(node.depth.edgelength(trees$trait))
      box_offset   <- tree_depth * 0.05
      label_offset <- tree_depth * 0.15
      plot_switching_tree(trees$trait, pred_df, fill_limits = range(cols$prob, na.rm = TRUE),
                          title = "Trait similarity", levs = rep(box_offset, 3), textoffset = label_offset)
    } else {
      plot(trees$trait, main = "Trait dendrogram (pruned)", cex = 0.8)
    }
  })
  
  #### ---- PCA explorer tab ---- ####
  
  filtered_site_year <- reactive({
    d <- app_data$site_year_lookup
    if (!is.null(input$region_filter) && input$region_filter != "All regions") {
      d <- d |> filter(basin_region == input$region_filter)
    }
    d
  })
  
  observeEvent(input$region_filter, {
    choices <- c("(none)", sort(filtered_site_year()$site_label))
    updateSelectizeInput(session, "site_year_select", choices = choices, selected = "(none)")
  }, ignoreNULL = FALSE)
  
  output$env_pca_plot <- renderPlot({
    all_d      <- app_data$site_year_lookup
    d          <- filtered_site_year()
    is_filtered <- !is.null(input$region_filter) && input$region_filter != "All regions"
    
    p <- ggplot() +
      theme_bw(base_size = 12) +
      labs(x = "PC1", y = "PC2",
           title = if (is_filtered) paste0("Environmental PCA -- ", input$region_filter) else "Environmental PCA -- all sampling site \u00d7 years")
    
    if (is_filtered) {
      # gray out everything outside the selected region for context, then
      # draw the selected region's points on top in color
      p <- p + geom_point(data = all_d, aes(PC1, PC2), color = "gray90", alpha = 0.5) +
        geom_point(data = d, aes(PC1, PC2), color = "#2166ac", alpha = 0.7)
    } else {
      p <- p + geom_point(data = all_d, aes(PC1, PC2), color = "gray70", alpha = 0.6)
    }
    
    if (!is.null(input$site_year_select) && input$site_year_select != "(none)") {
      hl <- d |> filter(site_label == input$site_year_select)
      p <- p + geom_point(data = hl, aes(PC1, PC2), color = "#c0392b", size = 4) +
        ggrepel::geom_label_repel(data = hl, aes(PC1, PC2, label = site_label), color = "#c0392b")
    }
    
    # DELTA CHANGE: there's no absolute "start" point anymore -- the
    # scenario is a pure difference, so it's drawn as an arrow FROM the
    # dataset average (PC = 0, since PCA scores are mean-centered by
    # construction) OUT BY the entered delta_PC1/delta_PC2. This is
    # arguably more honest than the old start/end-point version, since
    # the delta model has no real notion of "your specific starting
    # point" -- only of the size of the shift.
    if (isTRUE(input$show_custom_point) && !is.null(last_delta_pc())) {
      dp <- last_delta_pc()
      origin <- tibble(PC1 = 0, PC2 = 0)
      endpt  <- tibble(PC1 = dp["PC1"], PC2 = dp["PC2"])
      p <- p +
        geom_segment(data = bind_cols(origin, endpt %>% rename(xend = PC1, yend = PC2)),
                     aes(x = PC1, y = PC2, xend = xend, yend = yend),
                     arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
                     color = "#2e5f74", linewidth = 0.8, inherit.aes = FALSE) +
        geom_point(data = origin, aes(PC1, PC2), color = "#2e5f74", size = 4, shape = 1, stroke = 1.2) +
        geom_point(data = endpt,  aes(PC1, PC2), color = "#2e5f74", size = 4, shape = 17) +
        annotate("text", x = 0, y = 0, label = "no change", vjust = -1, color = "#2e5f74") +
        annotate("text", x = dp["PC1"], y = dp["PC2"], label = "your \u0394 scenario", vjust = -1, color = "#2e5f74")
    }
    p
  })
  
  output$trait_pca_plot <- renderPlot({
    d <- app_data$trait_pca_lookup
    if (!all(c("PC1", "PC2") %in% colnames(d))) {
      plot.new(); text(0.5, 0.5, "Trait PCA columns not found -- check trait_pca_lookup in shiny_app_data.RDS")
      return(invisible(NULL))
    }
    
    vars <- app_data$trait_varexp
    xlab <- if (!is.na(vars[1])) paste0("PC1 (", vars[1], "%)") else "PC1"
    ylab <- if (!is.na(vars[2])) paste0("PC2 (", vars[2], "%)") else "PC2"
    
    p <- ggplot(d, aes(PC1, PC2)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray85") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray85") +
      geom_point(color = "gray70", alpha = 0.6) +
      theme_bw(base_size = 12) +
      labs(title = "Trait PCA -- all species", x = xlab, y = ylab)
    
    vecs <- app_data$trait_vecs_df
    if (nrow(vecs) > 0) {
      vec_scale <- max(abs(c(d$PC1, d$PC2)), na.rm = TRUE) / max(abs(c(vecs$PC1, vecs$PC2)), na.rm = TRUE) * 0.8
      vecs <- vecs %>% mutate(PC1 = PC1 * vec_scale, PC2 = PC2 * vec_scale)
      p <- p +
        geom_segment(data = vecs, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                     arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                     color = "#c0392b", linewidth = 0.5, inherit.aes = FALSE) +
        geom_text(data = vecs, aes(x = PC1 * 1.1, y = PC2 * 1.1, label = label),
                  color = "#c0392b", size = 3, fontface = "bold", inherit.aes = FALSE)
    }
    
    if (!is.null(input$trait_species_select) && input$trait_species_select != "(none)") {
      hl <- d |> filter(clean.spec == input$trait_species_select)
      p <- p + geom_point(data = hl, color = "#2e5f74", size = 4) +
        ggrepel::geom_label_repel(data = hl, aes(label = clean.spec), color = "#2e5f74", fontface = "italic")
    }
    p
  })
}

shinyApp(ui, server)