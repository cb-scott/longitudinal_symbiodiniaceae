
########################################################
# DIV-LEVEL BRAY-CURTIS AND WEIGHTED UNIFRAC DISTANCES #
########################################################


#-------- Load required packages ------------
library(tidyverse)
library(readxl)
library(phyloseq)

#-------- Load phylogeny ------------

div.tree <- ape::read.tree("data/SymPortal_DIVs_tree.nwk")
div.tree <- phangorn::midpoint(div.tree)

stopifnot(
  !anyDuplicated(div.tree$tip.label),
  !is.null(div.tree$edge.length),
  all(is.finite(div.tree$edge.length)),
  all(div.tree$edge.length >= 0)
)

# Collect coverage and calibration data from each study
coverage.by.study <- new.env(parent = emptyenv())
calibration.by.study <- new.env(parent = emptyenv())


##################################################
# DEFINE COMMUNITY-DISTANCE AND TAXONOMY HELPERS #
##################################################
# These functions standardize DIV taxonomy, calculate community distances,
# and provide a calibrated fallback when direct UniFrac cannot be computed.

#-------- Extract genus labels from DIV identifiers ------------

get_div_genus <- function(ids) {
  
  genus <- rep(NA_character_, length(ids))
  
  # Standard SymPortal: A4, C3, D1, etc.
  standard <- grepl("^[A-GX][0-9]", ids)
  
  genus[standard] <- substr(ids[standard], 1, 1)
  
  # Numeric SymPortal: 2855659_A, 11843_C, etc.
  suffix <- grepl("_[A-GX]$", ids)
  
  genus[suffix] <- sub(
    ".*_([A-GX])$", "\\1", ids[suffix]
  )
  
  genus
}


#-------- Classify sample-level dominant and co-dominant genera ------------

get_last_genus <- function(comm, taxonomy = NULL,
                           threshold = 0.10) {
  
  comm <- as.matrix(comm)
  storage.mode(comm) <- "numeric"
  
  stopifnot(
    !is.null(rownames(comm)),
    all(is.finite(comm)),
    all(comm >= 0),
    all(rowSums(comm) > 0)
  )
  
  if (is.null(taxonomy)) {
    genus <- get_div_genus(colnames(comm))
  } else {
    genus <- unname(taxonomy[colnames(comm)])
  }
  
  if (any(is.na(genus) & colSums(comm) > 0)) {
    warning("Observed DIVs without genus assignments")
  }
  
  comm <- sweep(comm, 1, rowSums(comm), "/")
  
  genera <- sort(unique(na.omit(genus)))
  
  if (!length(genera)) {
    return(tibble(
      sample_name = rownames(comm),
      last_genus = NA_character_
    ))
  }
  
  genus.comm <- vapply(genera, function(g) {
    rowSums(comm[, !is.na(genus) & genus == g,
                 drop = FALSE])
  }, numeric(nrow(comm)))
  
  genus.comm <- matrix(
    genus.comm,
    nrow = nrow(comm),
    dimnames = list(rownames(comm), genera)
  )
  
  categories <- apply(genus.comm, 1, function(x) {
    
    present <- names(x)[x >= threshold]
    
    if (!length(present)) return(NA_character_)
    
    paste(sort(present), collapse = "_")
  })
  
  tibble(
    sample_name = rownames(comm),
    last_genus = unname(categories)
  )
}


#-------- Calculate direct weighted UniFrac when phylogenetic coverage is available ------------
# Flag samples if:
#  - any unmatched DIV >= 10% abundance
#  - total unmatched abundance > 50%
#
# Remaining matched DIVs are renormalized.

get_unifrac <- function(comm, tree = div.tree,
                        min_abund = 0.10,
                        max_missing = 0.50) {
  
  comm <- as.matrix(comm)
  storage.mode(comm) <- "numeric"
  
  stopifnot(
    !anyDuplicated(rownames(comm)),
    all(is.finite(comm)),
    all(comm >= 0),
    all(rowSums(comm) > 0)
  )
  
  comm <- sweep(comm, 1, rowSums(comm), "/")
  
  shared <- intersect(colnames(comm), tree$tip.label)
  missing <- setdiff(colnames(comm), tree$tip.label)
  
  missing_abund <- if (length(missing)) {
    rowSums(comm[, missing, drop = FALSE])
  } else {
    rep(0, nrow(comm))
  }
  
  max_missing_div <- if (length(missing)) {
    apply(comm[, missing, drop = FALSE], 1, max)
  } else {
    rep(0, nrow(comm))
  }
  
  flagged <- missing_abund > max_missing |
    max_missing_div >= min_abund
  
  coverage <- tibble(
    sample_name = rownames(comm),
    missing_abundance = as.numeric(missing_abund),
    max_missing_div = as.numeric(max_missing_div),
    flagged = as.logical(flagged)
  )
  
  result <- matrix(
    NA_real_,
    nrow(comm), nrow(comm),
    dimnames = list(rownames(comm), rownames(comm))
  )
  
  if (length(shared) >= 2) {
    
    valid <- !flagged &
      rowSums(comm[, shared, drop = FALSE]) > 0
    
    if (sum(valid) >= 2) {
      
      x <- comm[valid, shared, drop = FALSE]
      x <- x[, colSums(x) > 0, drop = FALSE]
      
      if (ncol(x) >= 2) {
        
        x <- sweep(x, 1, rowSums(x), "/")
        
        tree.sub <- ape::keep.tip(tree, colnames(x))
        x <- x[, tree.sub$tip.label, drop = FALSE]
        
        ps <- phyloseq::phyloseq(
          phyloseq::otu_table(
            x, taxa_are_rows = FALSE
          ),
          phyloseq::phy_tree(tree.sub)
        )
        
        d <- as.matrix(
          phyloseq::UniFrac(
            ps,
            weighted = TRUE,
            normalized = TRUE,
            parallel = FALSE
          )
        )
        
        result[rownames(d), colnames(d)] <- d
      }
    }
  }
  
  attr(result, "coverage") <- coverage
  
  message(
    sum(flagged), "/", nrow(comm),
    " samples flagged (tree coverage)"
  )
  
  result
}


#-------- Construct a genus-level proxy for incomplete DIV coverage ------------
# Use all known genera regardless of abundance.
# Sum the abundance of all DIVs within each genus.
# Distribute that genus's abundance equally across its tree tips.
#
# This captures genus-level turnover but not within-genus turnover.

get_genus_proxy_unifrac <- function(
    comm, tree = div.tree,
    min_genus_coverage = 0.90) {
  
  comm <- as.matrix(comm)
  storage.mode(comm) <- "numeric"
  
  comm <- sweep(comm, 1, rowSums(comm), "/")
  
  comm.genus <- get_div_genus(colnames(comm))
  tree.genus <- get_div_genus(tree$tip.label)
  
  genera <- intersect(
    unique(na.omit(comm.genus)),
    unique(na.omit(tree.genus))
  )
  
  result <- matrix(
    NA_real_,
    nrow(comm), nrow(comm),
    dimnames = list(rownames(comm), rownames(comm))
  )
  
  if (!length(genera)) return(result)
  
  genus.comm <- vapply(genera, function(g) {
    rowSums(
      comm[, !is.na(comm.genus) &
             comm.genus == g, drop = FALSE]
    )
  }, numeric(nrow(comm)))
  
  genus.comm <- matrix(
    genus.comm,
    nrow = nrow(comm),
    dimnames = list(rownames(comm), genera)
  )
  
  coverage <- rowSums(genus.comm)
  valid <- coverage >= min_genus_coverage
  
  attr(result, "genus_coverage") <- coverage
  
  if (sum(valid) < 2) return(result)
  
  tips <- tree$tip.label[
    !is.na(tree.genus) &
      tree.genus %in% genera
  ]
  
  if (length(tips) < 2) return(result)
  
  tree.sub <- ape::keep.tip(tree, tips)
  
  proxy <- matrix(
    0,
    nrow = sum(valid),
    ncol = length(tree.sub$tip.label),
    dimnames = list(
      rownames(comm)[valid],
      tree.sub$tip.label
    )
  )
  
  tip.genus <- get_div_genus(tree.sub$tip.label)
  
  for (g in genera) {
    
    idx <- which(tip.genus == g)
    
    if (length(idx)) {
      proxy[, idx] <- genus.comm[valid, g] / length(idx)
    }
  }
  
  proxy <- sweep(proxy, 1, rowSums(proxy), "/")
  
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(
      proxy, taxa_are_rows = FALSE
    ),
    phyloseq::phy_tree(tree.sub)
  )
  
  d <- as.matrix(
    phyloseq::UniFrac(
      ps,
      weighted = TRUE,
      normalized = TRUE,
      parallel = FALSE
    )
  )
  
  result[rownames(d), colnames(d)] <- d
  
  attr(result, "genus_coverage") <- coverage
  
  result
}


#-------- Partition Bray-Curtis turnover within and among genera ------------
# BC_within = BC_total - BC_genus
#
# This measures the DIV-level turnover not explained by changes
# in genus relative abundance.

get_within_genus_bray <- function(comm) {
  
  comm <- as.matrix(comm)
  storage.mode(comm) <- "numeric"
  
  comm <- sweep(comm, 1, rowSums(comm), "/")
  
  genus <- get_div_genus(colnames(comm))
  
  # Cannot accurately partition BC if observed taxa have
  # unknown genus assignments.
  if (any(is.na(genus) & colSums(comm) > 0)) {
    warning("Unknown genus: within-genus BC unavailable")
    
    return(matrix(
      NA_real_, nrow(comm), nrow(comm),
      dimnames = list(rownames(comm), rownames(comm))
    ))
  }
  
  genera <- sort(unique(na.omit(genus)))
  
  genus.comm <- vapply(genera, function(g) {
    rowSums(comm[, genus == g, drop = FALSE])
  }, numeric(nrow(comm)))
  
  genus.comm <- matrix(
    genus.comm,
    nrow = nrow(comm),
    dimnames = list(rownames(comm), genera)
  )
  
  total.bc <- as.matrix(
    vegan::vegdist(comm, method = "bray")
  )
  
  if (length(genera) == 1) {
    genus.bc <- matrix(
      0, nrow(comm), nrow(comm),
      dimnames = dimnames(total.bc)
    )
  } else {
    genus.bc <- as.matrix(
      vegan::vegdist(genus.comm, method = "bray")
    )
  }
  
  within.bc <- total.bc - genus.bc
  
  # Prevent negative values caused by floating-point precision
  within.bc[within.bc < 0] <- 0
  
  within.bc
}


#-------- Select direct or calibrated hybrid UniFrac estimates ------------
# alpha is provisional until calibrated across all studies.

assign_unifrac <- function(dist, uf_exact, uf_genus,
                           bc_within, alpha = 0.5,
                           bray_threshold = 0.10) {
  
  # Hybrid estimate
  uf_hybrid <- uf_genus +
    (1 - uf_genus) * alpha * bc_within
  
  uf_hybrid <- pmin(1, pmax(0, uf_hybrid))
  
  tibble(
    unifrac = case_when(
      !is.na(uf_exact) ~ uf_exact,
      !is.na(uf_hybrid) ~ uf_hybrid,
      dist < bray_threshold ~ 0,
      TRUE ~ NA_real_
    ),
    
    unifrac_method = case_when(
      !is.na(uf_exact) ~ "sequence",
      !is.na(uf_hybrid) ~ "hybrid_proxy",
      dist < bray_threshold ~ "bray_zero",
      TRUE ~ "unavailable"
    )
  )
}

#-------- Calculate sequential and baseline colony distances ------------
# Preserves original study workflow:
# all.dists -> sequential -> initial -> colony output
#
# Also retains the raw components needed for global calibration.

calculate_colony_distances <- function(
    comm, relevant.meta, colony.meta,
    study, colony_col = "COLONY_ID",
    date_col = "DATE", alpha = 0.5) {
  
  comm <- as.matrix(comm)
  storage.mode(comm) <- "numeric"
  
  meta <- relevant.meta %>%
    transmute(
      sample_name = as.character(sample_name),
      COLONY_ID = as.character(.data[[colony_col]]),
      DATE = as.Date(.data[[date_col]])
    ) %>%
    distinct()
  
  if (anyDuplicated(meta$sample_name)) {
    stop(study, ": duplicate sample metadata")
  }
  
  # Exclude accession rows and any samples absent from metadata
  keep <- !is.na(rownames(comm)) &
    rownames(comm) %in% meta$sample_name
  
  comm <- comm[keep, , drop = FALSE]
  
  if (anyDuplicated(rownames(comm))) {
    stop(study, ": duplicate sample IDs")
  }
  
  if (any(rowSums(comm) <= 0)) {
    stop(study, ": zero-abundance samples")
  }
  
  if (anyDuplicated(colony.meta$COLONY_ID)) {
    stop(study, ": duplicate colony metadata")
  }
  
  # Calculate component matrices
  genus <- get_last_genus(comm)
  
  uf <- get_unifrac(comm)
  
  uf.proxy <- get_genus_proxy_unifrac(comm)
  
  within.bc <- get_within_genus_bray(comm)
  
  assign(
    study,
    attr(uf, "coverage"),
    envir = coverage.by.study
  )
  
  # All pairwise comparisons
  all.dists <- as.matrix(
    vegan::vegdist(comm, method = "bray")
  ) %>%
    as.data.frame() %>%
    rownames_to_column("from") %>%
    pivot_longer(
      !from, names_to = "to", values_to = "dist"
    ) %>%
    mutate(
      uf_exact = uf[
        cbind(
          match(from, rownames(uf)),
          match(to, colnames(uf))
        )
      ],
      
      uf_genus = uf.proxy[
        cbind(
          match(from, rownames(uf.proxy)),
          match(to, colnames(uf.proxy))
        )
      ],
      
      bc_within = within.bc[
        cbind(
          match(from, rownames(within.bc)),
          match(to, colnames(within.bc))
        )
      ]
    ) %>%
    bind_cols(
      assign_unifrac(
        dist = as.vector(
          as.matrix(vegan::vegdist(comm, method = "bray"))
        )[match(
          paste(.$from, .$to),
          as.vector(outer(
            rownames(comm), rownames(comm),
            paste, sep = " "
          ))
        )],
        uf_exact = .$uf_exact,
        uf_genus = .$uf_genus,
        bc_within = .$bc_within,
        alpha = alpha
      )
    ) %>%
    left_join(genus, by = c("from" = "sample_name")) %>%
    left_join(meta, by = c("from" = "sample_name")) %>%
    left_join(meta, by = c("to" = "sample_name")) %>%
    filter(
      !is.na(COLONY_ID.x),
      COLONY_ID.x == COLONY_ID.y
    ) %>%
    group_by(COLONY_ID.x) %>%
    mutate(
      date.rank.x = dense_rank(DATE.x),
      date.rank.y = dense_rank(DATE.y)
    ) %>%
    ungroup() %>%
    distinct()
  
  # Collect pairs for calibration.
  # Use each pair once and avoid same-date comparisons.
  calibration <- all.dists %>%
    filter(
      DATE.x < DATE.y,
      dist >= 0.05,
      !is.na(uf_exact),
      !is.na(uf_genus),
      !is.na(bc_within)
    ) %>%
    transmute(
      STUDY = study,
      from, to,
      dist, uf_exact, uf_genus, bc_within
    )
  
  assign(
    study, calibration,
    envir = calibration.by.study
  )
  
  # Sequential community comparisons
  sequential <- all.dists %>%
    filter(date.rank.y == date.rank.x + 1) %>%
    select(
      COLONY_ID.x, DATE.x, DATE.y,
      dist, unifrac, unifrac_method,
      last_genus,
      UF_EXACT = uf_exact,
      UF_GENUS = uf_genus,
      BC_WITHIN = bc_within
    ) %>%
    distinct()
  
  # Comparisons relative to initial community
  initial <- all.dists %>%
    filter(date.rank.x == 1, date.rank.y > 1) %>%
    select(
      COLONY_ID.x, DATE.y,
      DIST_INITIAL = dist,
      UNIFRAC_INITIAL = unifrac,
      UNIFRAC_INITIAL_METHOD = unifrac_method,
      initial_genus = last_genus,
      UF_EXACT_INITIAL = uf_exact,
      UF_GENUS_INITIAL = uf_genus,
      BC_WITHIN_INITIAL = bc_within
    ) %>%
    distinct()
  
  sequential %>%
    left_join(
      initial, by = c("COLONY_ID.x", "DATE.y")
    ) %>%
    arrange(COLONY_ID.x, DATE.x) %>%
    rename(
      COLONY_ID = COLONY_ID.x,
      FROM_DATE = DATE.x,
      TO_DATE = DATE.y,
      DIST = dist,
      UNIFRAC = unifrac,
      UNIFRAC_METHOD = unifrac_method
    ) %>%
    left_join(colony.meta, by = "COLONY_ID") %>%
    select(
      COLONY_ID, FROM_DATE, TO_DATE,
      DIST, DIST_INITIAL,
      UNIFRAC, UNIFRAC_INITIAL,
      UNIFRAC_METHOD, UNIFRAC_INITIAL_METHOD,
      last_genus, initial_genus,
      UF_EXACT, UF_GENUS, BC_WITHIN,
      UF_EXACT_INITIAL, UF_GENUS_INITIAL,
      BC_WITHIN_INITIAL,
      LATITUDE, LONGITUDE, SPECIES
    )
}


#-------- Estimate the cross-study hybrid-distance calibration ------------

calibrate_unifrac_alpha <- function(
    calibration_env = calibration.by.study) {
  
  studies <- ls(calibration_env)
  
  calibration <- purrr::map_dfr(
    studies,
    ~ get(.x, envir = calibration_env)
  )
  
  if (!nrow(calibration)) {
    warning("No calibration pairs. Using alpha = 0.5.")
    return(0.5)
  }
  
  # Rearrange:
  # UF_exact - UF_genus =
  # alpha * (1 - UF_genus) * BC_within
  
  x <- with(
    calibration,
    (1 - uf_genus) * bc_within
  )
  
  y <- with(
    calibration,
    uf_exact - uf_genus
  )
  
  valid <- is.finite(x) & is.finite(y) & x > 0
  
  if (!any(valid)) {
    warning("No informative calibration pairs.")
    return(0.5)
  }
  
  # Least-squares estimate through origin
  alpha <- sum(x[valid] * y[valid]) /
    sum(x[valid]^2)
  
  # Constrain to 0-1
  alpha <- max(0, min(1, alpha))
  
  message(
    "Calibrated alpha = ", round(alpha, 4),
    " using ", sum(valid), " sample pairs"
  )
  
  alpha
}


#-------- Apply the calibrated hybrid distance to incomplete observations ------------
# Apply after all studies have been processed.
# Nunez Pons is unchanged if these components are absent.

recalibrate_colony_distances <- function(df, alpha) {
  
  required <- c(
    "UF_EXACT", "UF_GENUS", "BC_WITHIN",
    "UF_EXACT_INITIAL", "UF_GENUS_INITIAL",
    "BC_WITHIN_INITIAL"
  )
  
  if (!all(required %in% names(df))) return(df)
  
  sequential <- assign_unifrac(
    df$DIST,
    df$UF_EXACT,
    df$UF_GENUS,
    df$BC_WITHIN,
    alpha
  )
  
  initial <- assign_unifrac(
    df$DIST_INITIAL,
    df$UF_EXACT_INITIAL,
    df$UF_GENUS_INITIAL,
    df$BC_WITHIN_INITIAL,
    alpha
  )
  
  df %>%
    mutate(
      UNIFRAC = sequential$unifrac,
      UNIFRAC_METHOD = sequential$unifrac_method,
      UNIFRAC_INITIAL = initial$unifrac,
      UNIFRAC_INITIAL_METHOD =
        initial$unifrac_method
    )
}



#####################################
# HARMONIZE STUDY-SPECIFIC DIV DATA #
#####################################
# Each subsection reshapes one study's SymPortal output and metadata into a
# common sample-by-DIV table before calling the shared distance functions.

#-------- Lee et al. 2026 ------------
# Lee observations are represented in the longitudinal database; this block
# reconstructs the corresponding DIV-level intervals.

lee.divs <- read.table("data/symportal_profiles/Lee2026_DIVs_corr.txt", sep = "\t", header = T, check.names = FALSE)
lee.meta <- read_xlsx("data/symportal_profiles/Lee_Symportal_Metadata.xlsx", skip=1)
head(lee.meta)
relevant.meta = lee.meta %>% select(sample_name, colonyID, collection_date) %>%
  mutate(collection_date = as.Date(collection_date, format = "%Y%m%d"))

# coordinates for study
#5.610409515086883, 103.06064931859011

#first do colony level calls to include in database
lee.dists <- lee.divs %>%
  filter(!is.na(sample_name), sample_name != "") %>%
  semi_join(relevant.meta, by = "sample_name") %>%
  select(sample_name, A6b:G3x) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")


# Standardize Lee sample metadata to the same schema used below.
relevant.meta <- relevant.meta %>%
  rename(COLONY_ID = colonyID, DATE = collection_date) %>%
  mutate(COLONY_ID = paste("LEE", "ADIG", COLONY_ID, sep = "_"))
colony.meta <- relevant.meta %>%
  distinct(COLONY_ID) %>%
  mutate(LATITUDE = 5.610409515086883,
         LONGITUDE = 103.06064931859011,
         SPECIES = "ACROPORA_DIGITIFERA")

# Bray-Curtis + weighted UniFrac + genus classification
lee.colonydists <- calculate_colony_distances(
  comm = lee.dists, relevant.meta = relevant.meta,
  colony.meta = colony.meta, study = "LEE",
  colony_col = "COLONY_ID", date_col = "DATE"
)

lee.colonydists

#-------- Terraneo et al., 2026 ------------
terr.divs <- read.table("data/symportal_profiles/Terraneo_DIVs.txt", sep = "\t", header = T, check.names = FALSE)
terr.meta <- read_xlsx("data/symportal_profiles/Terraneo_SymPortal_Metadata.xlsx", skip=1)
terr.meta2 <- read_excel("data/Tarraneo_Meta.xlsx")
terr.meta2 = terr.meta2[-1,]
terr.meta2 = janitor::clean_names(terr.meta2)

terr.meta.full = terr.meta %>% select(colonyID, collection_date, sample_name) %>% 
  left_join(terr.meta2, by= c("colonyID" = "tagged_colony")) %>% 
  select(colonyID, collection_date, sample_name, genus, x4, sampling_locality) %>% distinct() %>%
  mutate(SPECIES = toupper(paste(genus, x4, sep = "_"))) %>%
  mutate(DATE = as.Date(paste(collection_date, "01", sep = ""), format = "%Y%m%d")) %>%
  mutate(COLONY_ID = paste("TERR", toupper(substr(genus, 1, 1)),
                           toupper(substr(x4, 1, 3)), colonyID, sep = "_")) 
terr.meta.full
lonlat.key = data.frame(sampling_locality = sort(unique(terr.meta.full$sampling_locality)), 
                        LATITUDE = c(-22.477683, -22.2295633, -22.285833, -22.285833, -22.312217),
                        LONGITUDE = c(166.463033, 166.435, 166.183167, 166.183167,166.633833))
terr.meta.full = terr.meta.full %>% left_join(lonlat.key)
relevant.meta = terr.meta.full %>% select(sample_name, COLONY_ID, DATE) %>% distinct()
colony.meta = terr.meta.full %>% select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>% distinct()


terr.dists <- terr.divs %>%
  dplyr::select(sample_name, A1:`28138_D`) %>% data.frame() %>% column_to_rownames(var = "sample_name")

# Bray-Curtis + weighted UniFrac + genus classification
terrcolonydists <- calculate_colony_distances(
  comm = terr.dists, relevant.meta = relevant.meta,
  colony.meta = colony.meta, study = "TERRANEO",
  colony_col = "COLONY_ID", date_col = "DATE"
)

terrcolonydists

#-------- Buzzoni et al., 2026 ------------
# Buzzoni observations are represented in the longitudinal database; this
# block reconstructs the corresponding DIV-level intervals.


buzz.divs <- read.csv("data/symportal_profiles/Buzzoni2026_DIVs.csv", check.names = FALSE)
buzz.meta <- read.csv("data/symportal_profiles/Buzzoni2026_metadata.csv")

# Only retain samples with both metadata and DIVs
shared.samps <- intersect(buzz.divs$sample_name, buzz.meta$sample_name)

# Prepare metadata
buzz.meta.full <- buzz.meta %>%
  filter(sample_name %in% shared.samps) %>%
  mutate(
    SPECIES = "PLATYGYRA_RYUKYUENSIS",
    DATE = as.Date(as.character(collection_date), format = "%Y%m%d"),
    COLONY_ID = paste("BUZZ", "PRYUK", coral_tag, sep = "_"),
    LATITUDE = as.numeric(collection_latitude),
    LONGITUDE = as.numeric(collection_longitude)
  ) %>%
  distinct()

relevant.meta <- buzz.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- buzz.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare community composition matrix
buzz.dists <- buzz.divs %>%
  filter(sample_name %in% shared.samps) %>%
  select(sample_name, A5b:`1269165_G`) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Convert counts to relative abundances
buzz.dists <- sweep(buzz.dists, 1, rowSums(buzz.dists), "/")

# Bray-Curtis + weighted UniFrac + genus classification
buzzcolonydists <- calculate_colony_distances(
  comm = buzz.dists, relevant.meta = relevant.meta,
  colony.meta = colony.meta, study = "BUZZONI",
  colony_col = "COLONY_ID", date_col = "DATE"
)

buzzcolonydists

#-------- Claar et al., 2020 ------------
#some overlap with Buzzoni; keep only Favites, worried about double counting platygyra

claar.divs <- read.table(
  "data/symportal_profiles/Claar_DIVs.txt",
  sep = "\t", header = TRUE, check.names = FALSE
)

claar.meta <- read_xlsx(
  "data/symportal_profiles/Claar_Symportal_BioPRJ_Meta.xlsx",
  skip = 1
)

# Prepare metadata, excluding Platygyra
claar.meta.full <- claar.meta %>%
  filter(tolower(host_genus) != "platygyra") %>%
  mutate(
    SPECIES = toupper(paste(host_genus, host_species, sep = "_")),
    DATE = as.Date(
      paste(collection_date, "01", sep = ""),
      format = "%Y%m%d"
    ),
    COLONY_ID = paste(
      "CLAAR",
      toupper(substr(host_genus, 1, 1)),
      toupper(substr(host_species, 1, 3)),
      colonyID, sep = "_"
    ),
    LATITUDE = as.numeric(collection_latitude),
    LONGITUDE = as.numeric(collection_longitude)
  ) %>%
  distinct()

relevant.meta <- claar.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- claar.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare community composition matrix
claar.dists <- claar.divs %>%
  filter(!is.na(sample_name)) %>%
  semi_join(relevant.meta, by = "sample_name") %>%
  select(sample_name, A5b:D17c) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Bray-Curtis + weighted UniFrac + genus classification
claarcolonydists <- calculate_colony_distances(
  comm = claar.dists, relevant.meta = relevant.meta,
  colony.meta = colony.meta, study = "CLAAR",
  colony_col = "COLONY_ID", date_col = "DATE"
)

claarcolonydists

#-------- Van Nynatten et al., 2025 ------------

vn.divs <- read.table(
  "data/symportal_profiles/van_nynatten_DIVs_corr.txt",
  sep = "\t", header = TRUE, check.names = FALSE
)

vn.seqid <- read.csv("data/van_nynatten_sample_metadata.csv")
vn.sites <- read.csv("data/van_nynatten_survey_metadata.csv")

vn.exps <- read.csv("data/van_nynatten_expeditions.csv") %>%
  select(Expedition.code, Mean.date) %>%
  rename(
    expedition = Expedition.code,
    Mean.date = Mean.date
  ) %>%
  mutate(
    expedition = ifelse(expedition == "2023", "2023b", expedition),
    DATE = as.Date(Mean.date, format = "%m/%d/%Y")
  )

# Only retain samples with both metadata and DIVs
shared.samps <- intersect(
  vn.divs$sample_name,
  vn.seqid$library_ID
)

# Prepare metadata
vn.meta.full <- vn.seqid %>%
  filter(library_ID %in% shared.samps) %>%
  left_join(vn.exps, by = "expedition") %>%
  left_join(vn.sites, by = "site_name") %>%
  mutate(
    sample_name = library_ID,
    SPECIES = toupper(paste(host_genus, host_species, sep = "_")),
    COLONY_ID = paste0(
      "VN_",
      substr(toupper(host_genus), 1, 1),
      substr(toupper(host_species), 1, 3),
      "_", coral_tag
    ),
    LATITUDE = lat,
    LONGITUDE = lon
  ) %>%
  distinct()

relevant.meta <- vn.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- vn.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare community composition matrix
vn.dists <- vn.divs %>%
  filter(sample_name %in% shared.samps) %>%
  select(sample_name, A5z:G3aq) %>%
  filter(rowSums(across(-sample_name)) > 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Convert counts to relative abundances
vn.dists <- sweep(vn.dists, 1, rowSums(vn.dists), "/")

# Bray-Curtis + weighted UniFrac + genus classification
vncolonydists <- calculate_colony_distances(
  comm = vn.dists, relevant.meta = relevant.meta,
  colony.meta = colony.meta, study = "VAN_NYNATTEN",
  colony_col = "COLONY_ID", date_col = "DATE"
)

vncolonydists

#-------- Leinbach et al., 2023 ------------


lbach.meta <- read.csv("data/Leinbach_profs.csv")

lbach.divs <- read.table(
  "data/symportal_profiles/Leinbach_DIVs.txt",
  sep = "\t", header = TRUE, check.names = FALSE
)

# Prepare metadata
lbach.meta.full <- lbach.meta %>%
  select(Sample, Individual, month, reef_zone) %>%
  distinct() %>%
  mutate(
    DATE = as.Date(
      paste(month, "01", "2019"),
      format = "%B %d %Y"
    ),
    COLONY_ID = paste("LBACH_AHYC", Individual, sep = "_"),
    SPECIES = "ACROPORA_HYACINTHUS",
    LATITUDE = case_when(
      reef_zone == "crest" ~ -17.4741,
      reef_zone == "forereef" ~ -17.4731
    ),
    LONGITUDE = case_when(
      reef_zone == "crest" ~ -149.8173,
      reef_zone == "forereef" ~ -149.8176
    )
  )

relevant.meta <- lbach.meta.full %>%
  transmute(
    sample_name = as.character(Sample),
    COLONY_ID,
    DATE
  ) %>%
  distinct()

colony.meta <- lbach.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare DIV community composition matrix
lbach.dists <- lbach.divs %>%
  semi_join(relevant.meta, by = "sample_name") %>%
  select(sample_name, A1:`1767177_G`) %>%
  filter(rowSums(across(-sample_name)) > 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Bray-Curtis + weighted UniFrac + genus classification
lbachcolonydists <- calculate_colony_distances(
  comm = lbach.dists, relevant.meta = relevant.meta,
  colony.meta = colony.meta, study = "LEINBACH",
  colony_col = "COLONY_ID", date_col = "DATE"
)

lbachcolonydists

#-------- Howells et al. 2020 ------------

how.divs <- read_xlsx(
  "data/symportal_profiles/Howells_DIVs_Meta.xlsx",
  skip = 1
)

# Prepare metadata

# Published sampling coordinates: Howells et al. (2020)
how.sites <- tibble(
  SITE = c("Delma", "Saadiyat", "Ras Ghanada", "Dibba", "Al Aqah"),
  LATITUDE = c(24.5208, 24.5986, 24.8482, 25.6031, 25.4929),
  LONGITUDE = c(52.2781, 54.4201, 54.6903, 56.3485, 56.3635)
)

# Prepare metadata
how.meta.full <- how.divs %>%
  transmute(
    SITE = trimws(Site),
    SPECIES = toupper(gsub(" ", "_", trimws(Species))),
    ORIGINAL_ID = trimws(`Colony ID`),
    DATE = as.Date(
      paste0("01-", sub("^DNA\\s+", "", Timepoint)),
      format = "%d-%b-%Y"
    )
  ) %>%
  filter(
    !is.na(SITE), SITE != "",
    !is.na(ORIGINAL_ID), ORIGINAL_ID != "",
    !is.na(DATE)
  ) %>%
  mutate(
    SITE = sub("^.*,\\s*", "", SITE),
    COLONY_ID = paste("HOW", ORIGINAL_ID, sep = "_"),
    sample_name = paste(COLONY_ID, format(DATE, "%Y%m"), sep = "_")
  ) %>%
  left_join(how.sites, by = "SITE") %>%
  distinct()

relevant.meta <- how.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- how.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare community composition matrix
# Columns 5 onward contain ITS2 DIV abundances
how.dists <- how.divs %>%
  mutate(
    ORIGINAL_ID = trimws(`Colony ID`),
    COLONY_ID = paste("HOW", ORIGINAL_ID, sep = "_"),
    DATE = as.Date(
      paste0("01-", sub("^DNA\\s+", "", Timepoint)),
      format = "%d-%b-%Y"
    ),
    sample_name = paste(COLONY_ID, format(DATE, "%Y%m"), sep = "_")
  ) %>%
  semi_join(relevant.meta, by = "sample_name") %>%
  select(sample_name, all_of(names(how.divs)[5:ncol(how.divs)])) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Normalize absolute sequence abundances
how.dists <- sweep(how.dists, 1, rowSums(how.dists), "/")

# Bray-Curtis + weighted UniFrac + genus classification
howcolonydists <- calculate_colony_distances(
  comm = how.dists, relevant.meta = relevant.meta,
  colony.meta = colony.meta, study = "HOWELLS",
  colony_col = "COLONY_ID", date_col = "DATE"
)

howcolonydists

#-------- Jung et al. 2021 ------------

jung.divs <- readxl::read_excel(
  "data/Jung_2021_Symportal_post_med_absolut abundances_meta.xlsx"
)

# Prepare metadata
jung.meta.full <- jung.divs %>%
  filter(!is.na(`Sample ID`)) %>%
  select(`Sample ID`) %>%
  separate(
    `Sample ID`,
    into = c("Mon", "Tidal", "Bleached", "Individual"),
    sep = " ",
    remove = FALSE
  ) %>%
  mutate(
    sample_name = `Sample ID`,
    COLONY_ID = paste0("Jung_AASP_", Individual),
    SPECIES = "ACROPORA_ASPERA",
    DATE = case_when(
      Mon == "Apr" ~ as.Date("2016-04-01"),
      Mon == "Nov" ~ as.Date("2016-11-01")
    ),
    LATITUDE = -16.479,
    LONGITUDE = 123.0435
  ) %>%
  distinct()

relevant.meta <- jung.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- jung.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare DIV community composition matrix
jung.dists <- jung.divs %>%
  filter(!is.na(`Sample ID`)) %>%
  select(!sample_name) %>%
  rename(sample_name = `Sample ID`) %>%
  select(sample_name, `11843_A`:`11852_C`) %>%
  filter(rowSums(across(-sample_name)) > 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Convert absolute abundances to relative abundances
jung.dists <- sweep(
  jung.dists, 1, rowSums(jung.dists), "/"
)

# Bray-Curtis + weighted UniFrac + genus classification
jungcolonydists <- calculate_colony_distances(
  comm = jung.dists, relevant.meta = relevant.meta,
  colony.meta = colony.meta, study = "JUNG",
  colony_col = "COLONY_ID", date_col = "DATE"
)

jungcolonydists

#-------- Alessi et al. 2024 ------------

alessi.divs <- read.csv("data/Alessi_ITS2_results.csv", check.names = FALSE)

# Prepare metadata
alessi.meta.full <- alessi.divs %>%
  mutate(
    sample_name = paste0("ALESSI_SAMPLE_", row_number()),
    COLONY_ID = paste("ALESSI", "AMUR", ID, sep = "_"),
    SPECIES = "ACROPORA_MURICATA",
    DATE = case_when(
      Time == "T1" ~ as.Date("2021-01-18"),
      Time == "T2" ~ as.Date("2021-05-05")
    ),
    LATITUDE = ifelse(Site == "Bouraké", -21.952, -21.965),
    LONGITUDE = ifelse(Site == "Bouraké", 165.990, 165.945)
  ) %>%
  filter(!is.na(DATE), !is.na(ID)) %>%
  distinct()

relevant.meta <- alessi.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- alessi.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare DIV community composition matrix
alessi.dists <- alessi.divs %>%
  mutate(sample_name = paste0("ALESSI_SAMPLE_", row_number())) %>%
  semi_join(relevant.meta, by = "sample_name") %>%
  select(sample_name, C1:C21bh) %>%
  filter(rowSums(across(-sample_name)) > 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Normalize to ensure each sample sums to one
alessi.dists <- sweep(
  alessi.dists, 1, rowSums(alessi.dists), "/"
)

# Bray-Curtis + weighted UniFrac + genus classification
alessicolonydists <- calculate_colony_distances(
  comm = alessi.dists, relevant.meta = relevant.meta,
  colony.meta = colony.meta, study = "ALESSI",
  colony_col = "COLONY_ID", date_col = "DATE"
)

alessicolonydists

#-------- Glass et al., 2026 ------------

glass.divs <- read.csv(
  "data/Glass_2026/PostMED_relative_abundance.csv", check.names = FALSE
)

# Prepare metadata
glass.meta.full <- glass.divs %>%
  mutate(
    DATE = case_when(
      Year == 2022 ~ as.Date("2022-05-01"),
      Year == 2023 ~ as.Date("2023-08-15"),
      Year == 2024 ~ as.Date("2024-04-28")
    ),
    COLONY_ID = paste("GLASS_SSID", Coral_ID, sep = "_"),
    SPECIES = "SIDERASTREA_SIDEREA",
    LATITUDE = case_when(
      Site == "Cristobal Island" ~ 9.268270,
      Site == "Hospital Point" ~ 9.333190,
      Site == "Punta Donato" ~ 9.35845,
      Site == "STRI Point" ~ 9.34868
    ),
    LONGITUDE = case_when(
      Site == "Cristobal Island" ~ -82.242280,
      Site == "Hospital Point" ~ -82.217650,
      Site == "Punta Donato" ~ -82.36804,
      Site == "STRI Point" ~ -82.2627
    ),
    sample_name = paste0("GLASS_SAMPLE_", row_number())
  ) %>%
  filter(!is.na(DATE), !is.na(Coral_ID))

relevant.meta <- glass.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- glass.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare DIV community composition matrix
glass.dists <- glass.meta.full %>%
  select(sample_name, A3:D2u) %>%
  filter(rowSums(across(-sample_name)) > 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Normalize relative abundances to sum to one
glass.dists <- sweep(
  glass.dists, 1, rowSums(glass.dists), "/"
)

# Bray-Curtis + weighted UniFrac + genus classification
glasscolonydists <- calculate_colony_distances(
  comm = glass.dists, relevant.meta = relevant.meta,
  colony.meta = colony.meta, study = "GLASS",
  colony_col = "COLONY_ID", date_col = "DATE"
)

glasscolonydists

#-------- Powell et al. 2026 ------------

powell.divs <- read.csv(
  "data/symportal_profiles/Powell_DIVs.csv",
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

powell.meta <- read.csv(
  "data/symportal_profiles/Powell_sample_info_curacao_2020.csv",
  fileEncoding = "latin1",
  colClasses = c(number = "character")
)

# Prepare metadata
powell.meta.full <- powell.meta %>%
  filter(
    !is.na(host_genus),
    !is.na(host_species),
    !is.na(number),
    !is.na(site)
  ) %>%
  mutate(
    sample_name = as.character(id),
    
    DATE = case_when(
      grepl("^[0-9]{6}$", collection_date) ~
        as.Date(paste0(collection_date, "01"), format = "%Y%m%d"),
      TRUE ~ as.Date(collection_date, format = "%m/%d/%y")
    ),
    
    COLONY_ID = paste0(
      "POWELL_",
      toupper(substr(host_genus, 1, 1)),
      toupper(substr(host_species, 1, 3)),
      "_", site, "_", number
    ),
    
    SPECIES = toupper(
      paste(host_genus, host_species, sep = "_")
    )
  ) %>%
  filter(!is.na(DATE)) %>%
  distinct()

# Use site-level coordinates to avoid inconsistent
# coordinates in individual metadata records
powell.sites <- powell.meta.full %>%
  group_by(site) %>%
  summarise(
    LATITUDE = as.numeric(
      names(
        sort(table(collection_latitude), decreasing = TRUE)
      )[1]
    ),
    LONGITUDE = -abs(as.numeric(
      names(
        sort(table(collection_longitude), decreasing = TRUE)
      )[1]
    )),
    .groups = "drop"
  )

powell.meta.full <- powell.meta.full %>%
  select(-collection_latitude, -collection_longitude) %>%
  left_join(powell.sites, by = "site")

# Prepare DIV matrix and normalize raw counts
powell.dists <- powell.divs %>%
  rename(sample_name = id) %>%
  semi_join(powell.meta.full, by = "sample_name") %>%
  select(sample_name, where(is.numeric)) %>%
  filter(rowSums(across(-sample_name)) > 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

powell.dists <- sweep(
  powell.dists, 1, rowSums(powell.dists), "/"
)

# Average relative abundance profiles within colony-date
# when multiple samples represent the same observation
powell.averaged <- powell.dists %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample_name") %>%
  left_join(
    powell.meta.full %>%
      select(sample_name, COLONY_ID, DATE),
    by = "sample_name"
  ) %>%
  select(-sample_name) %>%
  group_by(COLONY_ID, DATE) %>%
  summarise(
    across(where(is.numeric), mean),
    .groups = "drop"
  ) %>%
  mutate(
    sample_name = paste(
      COLONY_ID, format(DATE, "%Y%m%d"), sep = "__"
    )
  )

relevant.meta <- powell.averaged %>%
  select(sample_name, COLONY_ID, DATE)

colony.meta <- powell.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Community matrix: one row per colony and date
powell.dists <- powell.averaged %>%
  select(-COLONY_ID, -DATE) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")


#-------- Correct Powell DIV identifiers ------------
# Powell uses genus prefixes:
#   A_A4       -> A4
#   A_A4a      -> A4a
#   A_2855659  -> 2855659_A
#   C_123456   -> 123456_C
#
# For numeric IDs, move genus to the end.
# For standard named DIVs, simply remove the prefix.

old_names <- colnames(powell.dists)

genus <- sub("^([A-GX])_.*$", "\\1", old_names)
div_id <- sub("^[A-GX]_", "", old_names)

colnames(powell.dists) <- ifelse(
  grepl("^[0-9]+$", div_id),
  paste0(div_id, "_", genus),
  div_id
)

# Check for duplicated names after conversion
stopifnot(!anyDuplicated(colnames(powell.dists)))

# Check conversion
head(
  data.frame(
    original = old_names,
    corrected = colnames(powell.dists)
  ),
  20
)

#-------- Normalize averaged profiles ------------

powell.dists <- sweep(
  powell.dists, 1, rowSums(powell.dists), "/"
)

#-------- Check phylogenetic coverage ------------

powell.taxa <- colnames(powell.dists)[
  colSums(powell.dists) > 0
]

cat(
  "Observed DIVs:", length(powell.taxa), "\n",
  "Matching tree tips:",
  length(intersect(powell.taxa, div.tree$tip.label)), "\n",
  "Missing DIVs:",
  length(setdiff(powell.taxa, div.tree$tip.label)), "\n"
)

#-------- Calculate Bray-Curtis, weighted UniFrac, and genus classification ------------

powellcolonydists <- calculate_colony_distances(
  comm = powell.dists,
  relevant.meta = relevant.meta,
  colony.meta = colony.meta,
  study = "POWELL",
  colony_col = "COLONY_ID",
  date_col = "DATE"
)

powellcolonydists

#-------- Huffmyer et al. 2025 ------------


huff.divs <- read.table(
  "data/Huffmyer-Symportal.txt",
  sep = "\t", header = TRUE, check.names = FALSE
)

# Prepare metadata
huff.meta.full <- huff.divs %>%
  select(
    sample_name, collection_longitude, collection_latitude,
    collection_date, collection_depth, host_genus, host_species
  ) %>%
  filter(
    !is.na(collection_longitude),
    collection_longitude != 999
  ) %>%
  mutate(
    DATE = as.Date(
      as.character(collection_date),
      format = "%Y%m%d"
    )
  ) %>%
  separate(
    sample_name,
    into = c("spec", "ID", "tp"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    COLONY_ID = paste("HUFF", spec, ID, sep = "_"),
    SPECIES = toupper(paste(host_genus, host_species, sep = "_")),
    LATITUDE = as.numeric(collection_latitude),
    LONGITUDE = as.numeric(collection_longitude)
  ) %>%
  filter(
    !is.na(DATE),
    !is.na(COLONY_ID),
    !is.na(LATITUDE),
    !is.na(LONGITUDE),
    LATITUDE != 999,
    LONGITUDE != 999
  ) %>%
  distinct()

relevant.meta <- huff.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- huff.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare DIV community composition matrix
huff.dists <- huff.divs %>%
  semi_join(relevant.meta, by = "sample_name") %>%
  select(sample_name, A1ee:G3ab) %>%
  filter(rowSums(across(-sample_name)) > 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Normalize abundances
huff.dists <- sweep(
  huff.dists, 1, rowSums(huff.dists), "/"
)

# Bray-Curtis + weighted UniFrac + genus classification
huffcolonydists <- calculate_colony_distances(
  comm = huff.dists, relevant.meta = relevant.meta,
  colony.meta = colony.meta, study = "HUFFMYER",
  colony_col = "COLONY_ID", date_col = "DATE"
)

huffcolonydists


#-------- Cunning et al. 2024 ------------


cun.divs <- read.table(
  "data/Cunning_Symportal_111_20200624_2020-06-24_03-40-50.238182.seqs.absolute.abund_and_meta.txt",
  sep = "\t", header = TRUE, check.names = FALSE
)

# Prepare metadata
cun.meta.full <- cun.divs %>%
  select(
    sample_name, collection_latitude, collection_longitude,
    collection_depth, collection_date,
    host_genus, host_species
  ) %>%
  separate(
    sample_name,
    into = c("SITE", "SPEC", "ID", "YEAR"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    COLONY_ID = paste("CUN", SITE, SPEC, ID, sep = "_"),
    DATE = as.Date(
      as.character(collection_date),
      format = "%Y.%m.%d"
    ),
    LATITUDE = as.numeric(collection_latitude),
    LONGITUDE = -abs(as.numeric(collection_longitude)),
    SPECIES = toupper(paste(host_genus, host_species, sep = "_"))
  ) %>%
  filter(
    !is.na(DATE),
    !is.na(COLONY_ID),
    !is.na(LATITUDE),
    !is.na(LONGITUDE)
  ) %>%
  distinct()

cun.meta.full <- cun.meta.full %>%
  mutate(
    host_genus = if_else(
      sample_name == "ET_Ofra_59_2017",
      "Orbicella", host_genus
    ),
    host_species = if_else(
      sample_name == "ET_Ofra_59_2017",
      "franksi", host_species
    ),
    SPECIES = toupper(paste(host_genus, host_species, sep = "_"))
  )

cun.meta.full <- cun.meta.full %>%
  group_by(COLONY_ID, DATE) %>%
  filter(n() == 1) %>%
  ungroup()

relevant.meta <- cun.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- cun.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare DIV community composition matrix
cun.dists <- cun.divs %>%
  semi_join(relevant.meta, by = "sample_name") %>%
  select(sample_name, A4:G3b) %>%
  filter(rowSums(across(-sample_name)) > 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Convert absolute abundances to relative abundances
cun.dists <- sweep(
  cun.dists, 1, rowSums(cun.dists), "/"
)

# Bray-Curtis + weighted UniFrac + genus classification
cuncolonydists <- calculate_colony_distances(
  comm = cun.dists, relevant.meta = relevant.meta,
  colony.meta = colony.meta, study = "CUNNING",
  colony_col = "COLONY_ID", date_col = "DATE"
)

cuncolonydists

#-------- de Souza et al., 2023 ------------

desouza.divs <- readxl::read_excel(
  "data/deSouza_Symportal.xlsx"
)

# Prepare metadata
desouza.meta.full <- desouza.divs %>%
  select(sample_name, host_genus, host_species) %>%
  separate(
    sample_name,
    into = c("Date", "ColonyID"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    DATE = case_when(
      Date == "Clo2017" ~ as.Date("2018-01-01"),
      Date == "Clo2019" ~ as.Date("2019-10-01")
    ),
    COLONY_ID = paste0("deSouz_MCAP_", ColonyID),
    SPECIES = toupper(paste(host_genus, host_species, sep = "_")),
    LATITUDE = 21.47098,
    LONGITUDE = -157.82333
  ) %>%
  filter(!is.na(DATE), !is.na(ColonyID)) %>%
  distinct()

relevant.meta <- desouza.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- desouza.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Prepare DIV community matrix
desouza.dists <- desouza.divs %>%
  semi_join(relevant.meta, by = "sample_name") %>%
  select(sample_name, B1:F3.1) %>%
  filter(rowSums(across(-sample_name)) > 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "sample_name")

# Convert abundances to relative abundances
desouza.dists <- sweep(
  desouza.dists, 1, rowSums(desouza.dists), "/"
)

# Bray-Curtis + weighted UniFrac + genus classification
desouzacolonydists <- calculate_colony_distances(
  comm = desouza.dists, relevant.meta = relevant.meta,
  colony.meta = colony.meta, study = "DE_SOUZA",
  colony_col = "COLONY_ID", date_col = "DATE"
)

desouzacolonydists

#-------- Starko et al., 2023 ------------
load("data/Starko_Porites_DIV.RData")

# Extract metadata
stark.meta.full <- data.frame(sample_data(Porites_DIV)) %>% dplyr::select(!sample_name) %>%
  rownames_to_column(var = "sample_name") %>%
  mutate(
    DATE = as.Date(
      as.character(collection_date),
      format = "%Y%m%d"
    ),
    COLONY_ID = paste("STARK_PLOB", coral_tag, sep = "_"),
    SPECIES = toupper(paste(host_genus, host_species, sep = "_")),
    LATITUDE = as.numeric(collection_latitude),
    LONGITUDE = as.numeric(collection_longitude)
  ) %>%
  filter(
    !is.na(DATE),
    !is.na(coral_tag),
    !is.na(LATITUDE),
    !is.na(LONGITUDE),
    LATITUDE != 999,
    LONGITUDE != 999
  ) %>%
  distinct() %>%
  filter(
    !is.na(coral_tag),
    !trimws(as.character(coral_tag)) %in% c("", "NA", "NaN", "999")
  )

stark.meta.full <- stark.meta.full %>%
  filter(!COLONY_ID %in% c("STARK_PLOB_NA", "STARK_PLOB_", "STARK_PLOB_999"))

relevant.meta <- stark.meta.full %>%
  select(sample_name, COLONY_ID, DATE) %>%
  distinct()

colony.meta <- stark.meta.full %>%
  select(COLONY_ID, LATITUDE, LONGITUDE, SPECIES) %>%
  distinct()

# Remove samples with missing colony identifiers
stark.meta.full <- stark.meta.full %>%
  filter(
    !is.na(coral_tag),
    !trimws(as.character(coral_tag)) %in% c("", "NA", "NaN", "999")
  )

# Extract full DIV abundance matrix
stark.dists <- as(otu_table(Porites_DIV), "matrix")

# Ensure samples are rows
if (taxa_are_rows(Porites_DIV)) {
  stark.dists <- t(stark.dists)
}

# Retain samples with valid metadata
stark.dists <- stark.dists[
  rownames(stark.dists) %in% stark.meta.full$sample_name,
  ,
  drop = FALSE
]

# Remove samples with zero total abundance
stark.dists <- stark.dists[
  rowSums(stark.dists) > 0,
  ,
  drop = FALSE
]

# Normalize each sample to relative abundance
stark.dists <- sweep(
  stark.dists, 1, rowSums(stark.dists), "/"
)

# Average relative DIV abundances for replicate samples
stark.averaged <- stark.dists %>%
  as.data.frame() %>%
  rownames_to_column("sample_name") %>%
  left_join(
    stark.meta.full %>%
      select(sample_name, COLONY_ID, DATE),
    by = "sample_name",
    relationship = "many-to-one"
  ) %>%
  select(-sample_name) %>%
  group_by(COLONY_ID, DATE) %>%
  summarise(
    across(where(is.numeric), mean),
    .groups = "drop"
  ) %>%
  mutate(
    sample_name = paste(
      COLONY_ID,
      format(DATE, "%Y%m%d"),
      sep = "__"
    )
  )

# Update sample metadata
relevant.meta <- stark.averaged %>%
  select(sample_name, COLONY_ID, DATE)

# Rebuild community matrix
stark.dists <- stark.averaged %>%
  select(-COLONY_ID, -DATE) %>%
  as.data.frame() %>%
  column_to_rownames("sample_name")

# Bray-Curtis + weighted UniFrac + genus classification
starkcolonydists <- calculate_colony_distances(
  comm = stark.dists, relevant.meta = relevant.meta,
  colony.meta = colony.meta, study = "STARKO",
  colony_col = "COLONY_ID", date_col = "DATE"
)

starkcolonydists

#-------- Nunez-Pons et al. 2023 ------------
# SymPortal screening retained too few samples, failed to assemble the expected
# contigs, and returned metadata inconsistent with the reported region and
# sampling period. Use the DADA2-derived objects for this study instead.




load("data/get.taxa.nunez.Rdata")
load("data/Nunez_OTUtable.ITS2.Rdata")

otus <- as.matrix(seqtab.nochim)

# Retain samples with >= 5000 reads
otus <- otus[rowSums(otus) >= 5000, , drop = FALSE]

# Retain ASVs with >= 10 reads across retained samples
otus <- otus[, colSums(otus) >= 10, drop = FALSE]


tax.assign <- data.frame(get.taxa) 
tax.assign %>%
  rownames_to_column("ASV") -> tax.assign

# Keep ASVs with genus-level assignments, as in
# your original pipeline, but do not collapse them
# into genus/species abundance categories.

valid.asvs <- tax.assign %>%
  filter(
    !is.na(Genus),
    !Genus %in% c("", "g__", "NA")
  ) %>%
  pull(ASV)

otus <- otus[
  ,
  colnames(otus) %in% valid.asvs,
  drop = FALSE
]

# Check taxonomy-to-ASV matching
stopifnot(ncol(otus) > 0)


srr_2_samn <- read.csv(
  "data/PRJNA794042_meta.csv",
  stringsAsFactors = FALSE
) %>%
  select(Run, BioSample, SampleName) %>%
  distinct()

nunez.meta <- read.table(
  "data/NunezMeta.tsv",
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE
)

otu.meta <- tibble(ID = rownames(otus)) %>%
  left_join(
    srr_2_samn,
    by = c("ID" = "Run"),
    relationship = "one-to-one"
  ) %>%
  left_join(
    nunez.meta,
    by = c("BioSample" = "V1"),
    relationship = "one-to-one"
  ) %>%
  filter(
    !is.na(V7),
    V7 != "BLANK"
  ) %>%
  mutate(
    SPECIES = gsub(" ", "_", toupper(V7)),
    DATE = as.Date(V11, format = "%Y-%m-%d"),
    COLONY_ID = paste(
      "NUNEZ",
      substr(SPECIES, 1, 3),
      gsub("_t[0-9]+$", "", SampleName),
      sep = "_"
    ),
    LATITUDE = 21.461,
    LONGITUDE = -157.823
  ) %>%
  filter(
    !is.na(DATE),
    !is.na(SampleName),
    !is.na(SPECIES)
  ) %>%
  select(
    ID, COLONY_ID, DATE,
    SPECIES, LATITUDE, LONGITUDE
  )

# Check metadata consistency within colonies
otu.meta %>%
  group_by(COLONY_ID) %>%
  summarise(
    n_species = n_distinct(SPECIES),
    n_lat = n_distinct(LATITUDE),
    n_lon = n_distinct(LONGITUDE),
    .groups = "drop"
  ) %>%
  filter(n_species > 1 | n_lat > 1 | n_lon > 1)

# Check for multiple samples per colony-date
otu.meta %>%
  count(COLONY_ID, DATE) %>%
  filter(n > 1)


# Keep samples with valid metadata
otus <- otus[
  rownames(otus) %in% otu.meta$ID,
  ,
  drop = FALSE
]

# Exclude samples with no remaining Symbiodiniaceae reads
otus <- otus[
  rowSums(otus) > 0,
  ,
  drop = FALSE
]

# Normalize to relative abundance
otus.rel <- sweep(otus, 1, rowSums(otus), "/")


# Each biological sample contributes equally,
# irrespective of sequencing depth.

otu.averaged <- as.data.frame(
  otus.rel,
  check.names = FALSE
) %>%
  rownames_to_column("ID") %>%
  left_join(
    otu.meta %>%
      select(ID, COLONY_ID, DATE),
    by = "ID",
    relationship = "one-to-one"
  ) %>%
  group_by(COLONY_ID, DATE) %>%
  summarise(
    across(all_of(colnames(otus.rel)), mean),
    .groups = "drop"
  ) %>%
  mutate(
    sample_name = paste(
      COLONY_ID,
      format(DATE, "%Y%m%d"),
      sep = "__"
    )
  )

# Rebuild sample metadata after averaging technical replicates.
relevant.meta <- otu.averaged %>%
  select(sample_name, COLONY_ID, DATE)

# Rebuild ASV matrix
otus.rel <- otu.averaged %>%
  select(sample_name, all_of(colnames(otus.rel))) %>%
  as.data.frame(check.names = FALSE) %>%
  column_to_rownames("sample_name") %>%
  as.matrix()

# Verify alignment and uniqueness
stopifnot(
  !anyDuplicated(relevant.meta$sample_name),
  !anyDuplicated(relevant.meta[c("COLONY_ID", "DATE")]),
  all(rownames(otus.rel) %in% relevant.meta$sample_name)
)


# Genus is still classified for Nunez Pons, but NO UniFrac is calculated.
# ASV identifiers are sequences rather than SymPortal DIV tip IDs.
nunez.taxonomy <- tax.assign %>%
  mutate(genus_letter = case_when(
    grepl("Symbiodinium", Genus, ignore.case = TRUE) ~ "A",
    grepl("Breviolum", Genus, ignore.case = TRUE) ~ "B",
    grepl("Cladocopium", Genus, ignore.case = TRUE) ~ "C",
    grepl("Durusdinium", Genus, ignore.case = TRUE) ~ "D",
    grepl("Effrenium", Genus, ignore.case = TRUE) ~ "E",
    grepl("Fugacium", Genus, ignore.case = TRUE) ~ "F",
    grepl("Gerakladium", Genus, ignore.case = TRUE) ~ "G",
    TRUE ~ NA_character_
  )) %>%
  select(ASV, genus_letter) %>%
  tibble::deframe()

nunez.genus <- get_last_genus(otus.rel, taxonomy = nunez.taxonomy)

# Calculate distances within each colony rather
# than creating a full study-wide distance matrix.

colony.groups <- relevant.meta %>%
  arrange(COLONY_ID, DATE) %>%
  group_split(COLONY_ID)

nunez.dists <- lapply(colony.groups, function(m) {
  
  # Only longitudinally sampled colonies
  if (nrow(m) < 2) return(NULL)
  
  x <- otus.rel[
    m$sample_name,
    ,
    drop = FALSE
  ]
  
  d <- as.matrix(
    vegan::vegdist(x, method = "bray")
  )

  g <- nunez.genus$last_genus[
    match(m$sample_name, nunez.genus$sample_name)
  ]

  tibble(
    COLONY_ID = m$COLONY_ID[-1],
    FROM_DATE = m$DATE[-nrow(m)],
    TO_DATE = m$DATE[-1],
    DIST = d[cbind(
      seq_len(nrow(m) - 1),
      seq_len(nrow(m) - 1) + 1
    )],
    DIST_INITIAL = d[cbind(
      rep(1, nrow(m) - 1),
      seq.int(2, nrow(m))
    )],
    UNIFRAC = NA_real_,
    UNIFRAC_INITIAL = NA_real_,
    last_genus = g[-length(g)],
    initial_genus = rep(g[1], length(g) - 1)
  )
  
}) %>%
  bind_rows()


colony.meta <- otu.meta %>%
  select(
    COLONY_ID, LATITUDE,
    LONGITUDE, SPECIES
  ) %>%
  distinct()

# Ensure one metadata record per colony
stopifnot(!anyDuplicated(colony.meta$COLONY_ID))

nunezcolonydists <- nunez.dists %>%
  left_join(
    colony.meta,
    by = "COLONY_ID",
    relationship = "many-to-one"
  ) %>%
  arrange(COLONY_ID, FROM_DATE) %>%
  select(
    COLONY_ID, FROM_DATE, TO_DATE,
    DIST, DIST_INITIAL,
    UNIFRAC, UNIFRAC_INITIAL,
    last_genus, initial_genus,
    LATITUDE, LONGITUDE, SPECIES
  )


nunezcolonydists %>%
  group_by(SPECIES) %>%
  summarise(
    n_colonies = n_distinct(COLONY_ID),
    n_comparisons = n(),
    mean_dist = mean(DIST),
    .groups = "drop"
  )


#########################################################
# COMBINE STUDIES, VALIDATE DISTANCES, AND SAVE OUTPUTS #
#########################################################
# Bind the study-level results, calibrate the hybrid metric across studies,
# run final coverage checks, and save the reusable interval-level objects.

# 1. Estimate global alpha from observed sequence-level UniFrac
alpha <- calibrate_unifrac_alpha()

print(alpha)


# 2. Combine study-specific dataframes
all.colonydists <- list(
  LEE = lee.colonydists,
  TERRANEO = terrcolonydists,
  CLAAR = claarcolonydists,
  BUZZONI = buzzcolonydists,
  VAN_NYNATTEN = vncolonydists,
  LEINBACH = lbachcolonydists,
  HOWELLS = howcolonydists,
  JUNG = jungcolonydists,
  ALESSI = alessicolonydists,
  GLASS = glasscolonydists,
  POWELL = powellcolonydists,
  HUFFMYER = huffcolonydists,
  CUNNING = cuncolonydists,
  DE_SOUZA = desouzacolonydists,
  STARKO = starkcolonydists,
  NUNEZ_PONS = nunezcolonydists
) %>%
  purrr::imap_dfr(
    ~ recalibrate_colony_distances(.x, alpha) %>%
      rename_with(
        ~ rep("COLONY_ID", length(.x)),
        any_of("COLONYID")
      ) %>%
      mutate(STUDY = .y)
  )


# 3. Verify exact UniFrac distances were preserved

# Sequential distances
valid.seq <- !is.na(all.colonydists$UF_EXACT)

stopifnot(
  all(
    abs(
      all.colonydists$UNIFRAC[valid.seq] -
        all.colonydists$UF_EXACT[valid.seq]
    ) < 1e-10
  )
)

# Initial-community distances
valid.initial <- !is.na(all.colonydists$UF_EXACT_INITIAL)

stopifnot(
  all(
    abs(
      all.colonydists$UNIFRAC_INITIAL[valid.initial] -
        all.colonydists$UF_EXACT_INITIAL[valid.initial]
    ) < 1e-10
  )
)


# 4. Retain final variables
all.colonydists <- all.colonydists %>%
  select(
    STUDY, COLONY_ID, SPECIES,
    FROM_DATE, TO_DATE,
    DIST, DIST_INITIAL,
    UNIFRAC, UNIFRAC_INITIAL,
    UNIFRAC_METHOD, UNIFRAC_INITIAL_METHOD,
    last_genus, initial_genus,
    LATITUDE, LONGITUDE
  ) %>%
  arrange(STUDY, COLONY_ID, FROM_DATE)


# 5. Inspect distance estimation methods
method.summary <- all.colonydists %>%
  count(STUDY, UNIFRAC_METHOD) %>%
  group_by(STUDY) %>%
  mutate(
    proportion = n / sum(n)
  ) %>%
  ungroup()

print(method.summary, n = Inf)


# 6. Save results
save(
  all.colonydists,
  alpha,
  method.summary,
  file = "data/ColonyDistances.Rdata"
)

ggplot(all.colonydists, aes(x=DIST, y=UNIFRAC, color= UNIFRAC_METHOD)) + geom_point()
