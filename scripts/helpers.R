###########################
# SHARED WORKFLOW HELPERS #
###########################
# Functions in this file are called directly by the active numbered analysis
# scripts. Historical and exploratory helpers are retained separately with the
# original full helper file.


#############################
# GEOGRAPHIC CLASSIFICATION #
#############################

#-------- Assign sampling coordinates to ocean basins ------------
# Normalize longitudes to -180 to 180 degrees, then classify the study region
# using the geographic boundaries applied throughout the analysis.
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


######################################
# POSTERIOR RANDOM-EFFECT EXTRACTION #
######################################

#-------- Extract one species-level random-effect block ------------
# Build the Stan column names for one random-effect block and coefficient,
# confirm that all expected columns are present, and return them as a matrix.
# This helper uses n_spec and draws_df created by the calling Bayesian script.
extract_re <- function(block, coef_slot) {
  cols <- paste0("r_", block, "[", seq_len(n_spec), ",", coef_slot, "]")
  missing <- setdiff(cols, names(draws_df))
  if (length(missing) > 0) {
    stop("Missing columns: ", paste(head(missing, 3), collapse = ", "), "...")
  }
  as.matrix(draws_df[, cols, drop = FALSE])
}
