  #########################################################################
# ENVIRONMENTAL DATA EXTRACTION & PCA FOR CORAL SYMBIONT SWITCHING DATASET
#########################################################################
# Revised extraction:
# - Query a consistent, wider neighborhood around each sampling location.
# - Choose the nearest valid ocean grid cells, rather than averaging every
#   returned pixel or changing the spatial window if the nearest is missing.
# - Calculate spatial means at each satellite time step, then summarize time.
# - Retain extraction metadata and avoid global-mean substitution before PCA.
# - Use NEW filenames so older saved extractions aren't silently reused.
#########################################################################

library(tidyverse)
library(lubridate)
library(future)
library(furrr)
library(rerddap)
library(geosphere)
source("scripts/helpers.R")

# ==============================================================================
# 1. SET UP SAMPLING DATES AND LOCATIONS
# ==============================================================================

all_datapoints = readxl::read_excel("data/Long_Symbiont_Manual.xlsx") %>%
  mutate(DATE = as.Date(DATE))

load("data/ColonyDistances.Rdata")

bray_datapoints = all.colonydists %>%
  select(FROM_DATE, TO_DATE, LONGITUDE, LATITUDE) %>%
  pivot_longer(FROM_DATE:TO_DATE, names_to = "whichdate", values_to = "DATE") %>%
  select(-whichdate) %>%
  mutate(DATE = as.Date(DATE)) %>%
  distinct()

get_lee_env = all.colonydists %>%
  filter(STUDY == "LEE") %>%
  select(FROM_DATE, TO_DATE, LONGITUDE, LATITUDE) %>%
  pivot_longer(FROM_DATE:TO_DATE, names_to = "whichdate", values_to = "DATE") %>%
  select(-whichdate) %>%
  mutate(
    DATE = as.Date(DATE),
    LONGITUDE = 103.06064931859011,
    LATITUDE = 5.610409515086883
  ) %>%
  distinct()

# A date recorded as the 1st indicates month-level sampling information.
# Use the actual calendar month; otherwise use the preceding 30 days.
# If the 1st was an actual sampling date for any study, distinguish that
# using an explicit date-precision field rather than this convention.
to_grab = bind_rows(
  all_datapoints %>% select(LONGITUDE, LATITUDE, DATE),
  bray_datapoints,
  get_lee_env
) %>%
  distinct() %>%
  filter(!is.na(LONGITUDE), !is.na(LATITUDE), !is.na(DATE)) %>%
  mutate(
    from_date = if_else(day(DATE) == 1, DATE, DATE - days(30)),
    to_date = if_else(day(DATE) == 1,
                      as.Date(ceiling_date(DATE, "month") - days(1)), DATE)
  )



# ---- Standardize Lee coordinates ----

lee.lat <- 5.610409515086883
lee.lon <- 103.06064931859011

lee.dates <- unique(as.Date(get_lee_env$DATE))

to_grab <- to_grab %>%
  mutate(
    LATITUDE = if_else(
      DATE %in% lee.dates &
        abs(LATITUDE - lee.lat) < 0.001 &
        abs(LONGITUDE - lee.lon) < 0.001,
      lee.lat, LATITUDE
    ),
    LONGITUDE = if_else(
      DATE %in% lee.dates &
        abs(LATITUDE - lee.lat) < 0.001 &
        abs(LONGITUDE - lee.lon) < 0.001,
      lee.lon, LONGITUDE
    )
  ) %>%
  distinct(LATITUDE, LONGITUDE, DATE, .keep_all = TRUE)

# Check that Lee has exactly one row per date
to_grab %>%
  filter(
    DATE %in% lee.dates,
    LATITUDE == lee.lat,
    LONGITUDE == lee.lon
  ) %>%
  count(DATE)

# ==============================================================================
# 2. FUNCTIONS FOR CONSISTENT SPATIAL AND TEMPORAL SUMMARIES
# ==============================================================================

# Half-width of the initial search area (degrees). This replaces the
# 0.01 -> 0.05 fallback, rather than eliminating nearby-pixel retrieval.
search_radius = 1 #let's make a large radius
n_nearest_cells = 3L

dir.create("data/envdata", showWarnings = FALSE, recursive = TRUE)

safe_mean = function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
safe_min  = function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_max  = function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
safe_sd   = function(x) if (sum(!is.na(x)) < 2) NA_real_ else sd(x, na.rm = TRUE)
safe_sum  = function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)

# Select the nearest k cells with at least one finite observation in the
# reference field across the requested interval. Keep the selected cell
# identities fixed for all dates within that extraction. This makes nearby
# sites likely to share pixels, while respecting coastal missingness.
# IMPORTANT: this is a spatial mean of the nearest valid cells, NOT a local
# interpolation of land pixels or a global-mean imputation.
get_nearby_series = function(tmp, lat, lon, fields,
                             reference_field = fields[1],
                             n_cells = n_nearest_cells) {

  required = c("latitude", "longitude", "time", fields)
  missing = setdiff(required, names(tmp))
  if (length(missing)) stop("Missing ERDDAP columns: ", paste(missing, collapse = ", "))

  tmp = tmp %>%
    mutate(across(all_of(fields), ~ {
      x = as.numeric(.)
      x[!is.finite(x)] = NA_real_
      x
    }))

  valid_cells = tmp %>%
    group_by(latitude, longitude) %>%
    summarise(valid = any(!is.na(.data[[reference_field]])), .groups = "drop") %>%
    filter(valid) %>%
    mutate(
      dist_m = geosphere::distHaversine(
        cbind(longitude, latitude), c(lon, lat)
      )
    ) %>%
    arrange(dist_m) %>%
    slice_head(n = n_cells)

  if (nrow(valid_cells) == 0) {
    empty = tibble(time = unique(tmp$time))
    for (v in fields) empty[[v]] = NA_real_
    return(list(
      series = empty,
      n_cells = 0L,
      nearest_dist_km = NA_real_,
      farthest_dist_km = NA_real_,
      n_times = n_distinct(tmp$time)
    ))
  }

  series = tmp %>%
    semi_join(valid_cells %>% select(latitude, longitude),
              by = c("latitude", "longitude")) %>%
    group_by(time) %>%
    summarise(across(all_of(fields), safe_mean), .groups = "drop")

  list(
    series = series,
    n_cells = nrow(valid_cells),
    nearest_dist_km = min(valid_cells$dist_m) / 1000,
    farthest_dist_km = max(valid_cells$dist_m) / 1000,
    n_times = nrow(series)
  )
}

make_id = function(i) {
  paste(to_grab$LATITUDE[i], to_grab$LONGITUDE[i], to_grab$DATE[i], sep = "_")
}

# Keep results from this revised method separate from original txt files.
# readr CSV output includes headers; don't combine with old headerless txt.

# ==============================================================================
# 3. NOAA SST / DEGREE HEATING WEEKS / SST ANOMALY
# ==============================================================================
#NOTE TO USER: erdapp often fails to recover data that exists if too many server requests are made in a short period
#if your final data has many NAs, it's worth reruning this script 2-3 times to see if they fill
#the script includes a check to make sure you only rerun it for lon/lat/date combos that do not yet exist in the data
#(this will make subsequent runs much faster)
future::plan(future::multisession, workers = 11)

furrr::future_map(seq_len(nrow(to_grab)), function(row) {
  site_id = make_id(row)
  outfile = paste0("data/envdata/sst_v2_", site_id, ".csv")
  if (file.exists(outfile) && file.size(outfile) > 0) {
    old = tryCatch(readr::read_csv(outfile, show_col_types = FALSE),
                   error = function(e) NULL)
    
    if (!is.null(old) &&
        "mean_SST" %in% names(old) &&
        any(is.finite(old$mean_SST))) {
      return(NULL)
    }
  }
  tryCatch({
    tmp = rerddap::griddap(
      "NOAA_DHW",
      time = c(as.character(to_grab$from_date[row]),
               as.character(to_grab$to_date[row])),
      latitude = to_grab$LATITUDE[row] + c(-search_radius, search_radius),
      longitude = to_grab$LONGITUDE[row] + c(-search_radius, search_radius),
      fields = c("CRW_SST", "CRW_DHW", "CRW_SSTANOMALY")
    )$data

    nearby = get_nearby_series(
      tmp, to_grab$LATITUDE[row], to_grab$LONGITUDE[row],
      fields = c("CRW_SST", "CRW_DHW", "CRW_SSTANOMALY"),
      reference_field = "CRW_SST"
    )
    x = nearby$series

    env_data = tibble(
      lat = to_grab$LATITUDE[row], lon = to_grab$LONGITUDE[row],
      date = to_grab$DATE[row],
      min_SST = safe_min(x$CRW_SST),
      max_SST = safe_max(x$CRW_SST),
      sd_SST = safe_sd(x$CRW_SST),
      mean_SST = safe_mean(x$CRW_SST),
      max_DHW = safe_max(x$CRW_DHW),
      min_DHW = safe_min(x$CRW_DHW),
      mean_DHW = safe_mean(x$CRW_DHW),
      mean_ANOM = safe_mean(x$CRW_SSTANOMALY),
      sd_ANOM = safe_sd(x$CRW_SSTANOMALY),
      min_ANOM = safe_min(x$CRW_SSTANOMALY),
      max_ANOM = safe_max(x$CRW_SSTANOMALY),
      # One spatially averaged anomaly per satellite time step.
      # Interpret as the sum of observed time steps; only a day-sum
      # if the product has one time step per day.
      cum_ANOM = safe_sum(x$CRW_SSTANOMALY),
      sst_n_cells = nearby$n_cells,
      sst_nearest_km = nearby$nearest_dist_km,
      sst_farthest_km = nearby$farthest_dist_km,
      sst_n_times = nearby$n_times
    )

    readr::write_csv(env_data, outfile, na = "NA")
  }, error = function(e) {
    message("SST row ", row, " failed: ", conditionMessage(e))
    NULL
  })
}, .options = furrr::furrr_options(seed = TRUE))

future::plan(future::sequential)

# ==============================================================================
# 4. OCEAN COLOR (CHLOROPHYLL AND LIGHT ATTENUATION)
# ==============================================================================

future::plan(future::multisession, workers = 11)

furrr::future_map(seq_len(nrow(to_grab)), function(row) {
  site_id = make_id(row)
  outfile = paste0("data/envdata/color_v2_", site_id, ".csv")
  if (file.exists(outfile) && file.size(outfile) > 0) {
    old = tryCatch(readr::read_csv(outfile, show_col_types = FALSE),
                   error = function(e) NULL)
    
    if (!is.null(old) &&
        "mean_chla" %in% names(old) &&
        any(is.finite(old$mean_chla))) {
      return(NULL)
    }
  }

  tryCatch({
    tmp = rerddap::griddap(
      "occci_V6_8day_4km",
      url = "https://comet.nefsc.noaa.gov/erddap/",
      time = c(as.character(to_grab$from_date[row]),
               as.character(to_grab$to_date[row])),
      latitude = to_grab$LATITUDE[row] + c(-search_radius, search_radius),
      longitude = to_grab$LONGITUDE[row] + c(-search_radius, search_radius),
      fields = c("chlor_a", "kd_490")
    )$data

    nearby = get_nearby_series(
      tmp, to_grab$LATITUDE[row], to_grab$LONGITUDE[row],
      fields = c("chlor_a", "kd_490"),
      reference_field = "chlor_a"
    )
    x = nearby$series

    env_data = tibble(
      lat = to_grab$LATITUDE[row], lon = to_grab$LONGITUDE[row],
      date = to_grab$DATE[row],
      min_chla = safe_min(x$chlor_a),
      max_chla = safe_max(x$chlor_a),
      sd_chla = safe_sd(x$chlor_a),
      mean_chla = safe_mean(x$chlor_a),
      min_turb = safe_min(x$kd_490),
      max_turb = safe_max(x$kd_490),
      sd_turb = safe_sd(x$kd_490),
      mean_turb = safe_mean(x$kd_490),
      color_n_cells = nearby$n_cells,
      color_nearest_km = nearby$nearest_dist_km,
      color_farthest_km = nearby$farthest_dist_km,
      color_n_times = nearby$n_times
    )

    readr::write_csv(env_data, outfile, na = "NA")
  }, error = function(e) {
    message("COLOR row ", row, " failed: ", conditionMessage(e))
    NULL
  })
}, .options = furrr::furrr_options(seed = TRUE))

future::plan(future::sequential)

# ==============================================================================
# 5. MONTHLY PRECIPITATION
# ==============================================================================

future::plan(future::multisession, workers = 11)

furrr::future_map(seq_len(nrow(to_grab)), function(row) {
  site_id = make_id(row)
  outfile = paste0("data/envdata/precip_v2_", site_id, ".csv")
  if (file.exists(outfile) && file.size(outfile) > 0) {
    old = tryCatch(readr::read_csv(outfile, show_col_types = FALSE),
                   error = function(e) NULL)
    
    if (!is.null(old) &&
        "mean_precip" %in% names(old) &&
        any(is.finite(old$mean_precip))) {
      return(NULL)
    }
  }

  tryCatch({
    tmp = rerddap::griddap(
      "IMERG_monthly_global_precip",
      url = "https://erddap.marine.usf.edu/erddap/",
      time = c(as.character(to_grab$from_date[row]),
               as.character(to_grab$to_date[row])),
      latitude = to_grab$LATITUDE[row] + c(-search_radius, search_radius),
      longitude = to_grab$LONGITUDE[row] + c(-search_radius, search_radius),
      fields = c("precipitation")
    )$data

    nearby = get_nearby_series(
      tmp, to_grab$LATITUDE[row], to_grab$LONGITUDE[row],
      fields = "precipitation"
    )
    x = nearby$series

    env_data = tibble(
      lat = to_grab$LATITUDE[row], lon = to_grab$LONGITUDE[row],
      date = to_grab$DATE[row],
      # Monthly data may yield only a single time step in a 30-day window.
      # Keep only mean precipitation as a PCA predictor below.
      min_precip = safe_min(x$precipitation),
      max_precip = safe_max(x$precipitation),
      sd_precip = safe_sd(x$precipitation),
      mean_precip = safe_mean(x$precipitation),
      precip_n_cells = nearby$n_cells,
      precip_nearest_km = nearby$nearest_dist_km,
      precip_farthest_km = nearby$farthest_dist_km,
      precip_n_times = nearby$n_times
    )

    readr::write_csv(env_data, outfile, na = "NA")
  }, error = function(e) {
    message("PRECIP row ", row, " failed: ", conditionMessage(e))
    NULL
  })
}, .options = furrr::furrr_options(seed = TRUE))

future::plan(future::sequential)

# ==============================================================================
# 6. READ REVISED ENVIRONMENTAL FILES AND GENERATE ENV PCs
# ==============================================================================

read_product = function(pattern) {
  files = list.files("data/envdata", full.names = TRUE, pattern = pattern)
  if (!length(files)) stop("No files found for: ", pattern)
  purrr::map_dfr(files, readr::read_csv, show_col_types = FALSE) %>%
    mutate(date = as.Date(date)) %>%
    distinct(lat, lon, date, .keep_all = TRUE)
}

sst.data = read_product("^sst_v2_.*\\.csv$")
precip.data = read_product("^precip_v2_.*\\.csv$")
color.data = read_product("^color_v2_.*\\.csv$")

# Start from all requested rows, so failures in SST are not silently dropped.
join.data = to_grab %>%
  transmute(lat = LATITUDE, lon = LONGITUDE, date = DATE) %>%
  left_join(sst.data, by = c("lat", "lon", "date")) %>%
  left_join(precip.data, by = c("lat", "lon", "date")) %>%
  left_join(color.data, by = c("lat", "lon", "date"))

# Check extraction coverage and whether faraway pixels were used.
coverage = join.data %>%
  summarise(
    n = n(),
    sst_missing = mean(is.na(mean_SST)),
    color_missing = mean(is.na(mean_chla)),
    precip_missing = mean(is.na(mean_precip)),
    max_sst_km = safe_max(sst_farthest_km),
    max_color_km = safe_max(color_farthest_km),
    max_precip_km = safe_max(precip_farthest_km)
  )
print(coverage)

##----- impute NAs based on lat/lon/month and year. There are not many NA values left
library(missForest)
library(lubridate)

env.vars = c(
  "mean_ANOM", "min_ANOM", "max_ANOM",
  "cum_ANOM", "mean_DHW", "sd_ANOM",
  "mean_SST", "min_SST", "max_SST", "sd_SST",
  "mean_precip", "min_precip", "max_precip",
  "mean_chla", "min_chla", "max_chla", "sd_chla",
  "mean_turb", "min_turb", "max_turb", "sd_turb"
)

# Prepare environmental and spatial predictors
impute.dat = join.data %>%
  mutate(
    date = as.Date(date),
    month = as.numeric(month(date)),
    year = year(date)
  ) %>%
  select(lat, lon, month, year, all_of(env.vars)) %>%
  mutate(
    across(where(is.numeric),
           ~replace(., !is.finite(.), NA_real_))
  )

# Check missingness before imputation
colSums(is.na(impute.dat))

# Run missForest
set.seed(42)
impute.dat = as.data.frame(impute.dat)
imp = missForest(
  impute.dat,
  maxiter = 5,
  ntree = 300,
  verbose = TRUE,
  variablewise = TRUE
)

# Extract imputed environmental data
means.only = imp$ximp %>%
  select(all_of(env.vars))

# Check missingness after imputation
colSums(is.na(means.only))

# Check imputation error
imp$OOBerror

# ---- PCA ----

env.pc = prcomp(
  means.only,
  center = TRUE,
  scale. = TRUE
)

screeplot(env.pc, bstick = T)
env.score = as.data.frame(env.pc$x)

env.data.final = bind_cols(
  join.data,
  env.score
) %>%
  mutate(DATE = as.Date(date)) %>%
  select(!date) %>%
  distinct()

save(env.data.final, file = "data/EnvData_PCs.Rdata")

# Optional loadings heatmap
mat = env.pc$rotation[, 1:3]
lim = max(abs(mat), na.rm = TRUE)
svg("results/pc_loadings_v2.svg", width = 4, height = 4)
pheatmap::pheatmap(
  mat,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-lim, lim, length.out = 101)
)
dev.off()

plot(env.data.final$sd_SST, env.data.final$PC3)

###Or, maybe more intuitive, correlations
# Correlations of environmental variables with first 3 PCs
pc_cor = cor(
  means.only,
  env.pc$x[, 1:3],
  use = "complete.obs"
)

lim = max(abs(pc_cor), na.rm = TRUE)

svg("results/pc_correlations_v2.svg", width = 4, height = 6)
pheatmap::pheatmap(
  pc_cor,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-lim, lim, length.out = 101),
  main = "Correlation of environmental variables with PCs"
)
dev.off()
#eh, it's fundamentally the same plot (good), but wouldn't help in intepretation



cbind(join.data[,c("lat", "lon", "date")], data.frame(env.score)) %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d")) %>%
  mutate(mon = month(date)) %>%
  ggplot(aes(x=PC1, y=PC2, col = mon)) + geom_point() +
  scale_color_viridis_c()



#### ---------- NICE BIPLOT OF ENV PCS ------------- 
#This code also annotates the Kiribati and FL Keys mass bleaching event with reference point to recovery conditions
score_scale <- 0.4  # adjust if points feel cramped or spread
vec_scale   <- 2    # adjust arrow length to taste

scores_df <- data.frame(
  PC1     = env.pc$x[, 1],
  PC2     = env.pc$x[, 2],
  lat = join.data$lat,
  lon = join.data$lon,
  date = join.data$date,
  region = assign_ocean_basin(join.data$lon, join.data$lat)
)
label_df <- scores_df %>%
  filter(!is.na(region)) %>%
  group_by(region) %>%
  slice(1) 

vecs_df <- data.frame(
  varnames = rownames(env.pc$rotation),
  PC1      =  env.pc$rotation[, 1] * vec_scale,
  PC2      = env.pc$rotation[, 2] * vec_scale
) %>%
  mutate(label = rownames(env.pc$rotation)) %>%
  filter(varnames %in% c("mean_SST", "mean_precip",
                         "sd_SST", "mean_turb"))


label_df <- scores_df %>%
  filter(PC1 > quantile(PC1, 0.95) | PC1 < quantile(PC1, 0.05) |
           PC2 > quantile(PC2, 0.95) | PC2 < quantile(PC2, 0.05)) %>%
  #filter(lat < -10) %>% 
  filter(date%in%c("2005-09-01", "2015-07-01", "2016-11-01", "2005-12-01")) %>% 
  group_by(region, date) %>%
  slice(1) 

label_df$plotlab = c("FL Keys:Sept 2005", 
                     "FL Keys:Dec 2005",
                     "Kiribati:July 2015",
                     "Kiribati:Nov 2016"
)

# Region colors
region_colors <- c(
  "Indo-Pacific" = "#2a7c7c",
  "Caribbean"    = "#c1692b",
  "ETP"          = "#6b4c9a",
  "Red Sea"      = "#b5282a"
)

# --- Plot ---
library(ggforce)
library(ggrepel)
ggplot() +
  xlim(c(-12.5, 7.5)) + 
  ylim(c(-7.5, 10))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray90", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray90", linewidth = 0.4) +
  geom_point(
    data = scores_df,
    aes(x = PC1, y = PC2),
    col = "gray75",
    size = 0.8, alpha = 0.8
  ) +
  geom_segment(
    data = vecs_df,
    aes(x = 0, y = 0, xend = PC1*10, yend = PC2*10),
    arrow     = arrow(length = unit(0.2, "cm"), type = "closed"),
    color     = "#c0392b", linewidth = 0.5,
    alpha = 1
  ) +
  geom_text(
    data = vecs_df,
    aes(x = ifelse(PC1 > 0, PC1 * 10 + 1, PC1*10 - 1) , y = ifelse(PC2 > 0, (PC2 * 10)+.5, (PC2*10)-.5), label = label),
    color    = "#c0392b",
    alpha = 1,
    size     = 2.5,
    fontface = "bold"
  ) +
  # Labeled sites: slightly larger points
  geom_point(
    data = label_df,
    aes(x = PC1 , y = PC2 ),
    color = "#2e5f74",
    size = 2.2, shape = 22
  ) +
  # Site labels
  geom_label(
    data = label_df,
    aes(x = PC1-.75, y = ifelse(PC2 > 0, PC2+1, PC2 -1),
        label = plotlab),
    color = "#2e5f74",
    fontface = "italic",
    size          = 1.8,
    alpha = .6
  ) +
  scale_color_viridis_d(option = "viridis") +
  # scale_color_manual(values = c("#5b8fa8","#c17d4a","#7a9e5f","#b05c7a"), name = "Region") +
  #scale_fill_manual(values  = c("#5b8fa8","#c17d4a","#7a9e5f","#b05c7a"), name = "Region") +
 # labs(
 #   x = paste0("PC1 (", round(varexp[1], 2)*100, "%)"),
 #   y = paste0("PC2 (", round(varexp[2], 2)*100, "%)")
 # ) +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(size = 8),
    axis.line   = element_line(color = "gray40"),
    axis.ticks  = element_line(color = "gray40"),
    plot.margin = margin(10, 15, 10, 10),
    legend.position = "right"
  )

ggsave("results/env.pcs.svg", width = 3.25, height = 3.25)#, units = "in", dpi = 1000)


# ==============================================================================
# VISUALIZE GEOGRAPHIC AND TEMPORAL EXTENT
# ==============================================================================
library(maps)

points_df <- distinct(to_grab[,c("LONGITUDE", "LATITUDE")])
world <- map_data("world")
points_density <- all_datapoints %>%
  select(LONGITUDE, LATITUDE) %>%
  mutate(
    lon_bin = floor(LONGITUDE),
    lat_bin = floor(LATITUDE)
  ) %>%
  group_by(lon_bin, lat_bin) %>%
  summarise(
    nobs    = n(),
    mean_lon = mean(LONGITUDE),
    mean_lat = mean(LATITUDE),
    .groups = "drop"
  ) 
plot_cdist = all.colonydists %>% select(LATITUDE, LONGITUDE) %>%
  mutate(
    lon_bin = floor(LONGITUDE),
    lat_bin = floor(LATITUDE)
  ) %>%
  group_by(lon_bin, lat_bin) %>%
  summarise(
    nobs = n(), 
    mean_lon = mean(LONGITUDE),
    mean_lat = mean(LATITUDE),
    .groups = "drop"
  ) 
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "gray80", color = "gray80", linewidth = 0.2) +
  geom_point(data = points_density, aes(x = mean_lon, y = mean_lat, size = nobs),
             color = "#c0392b", alpha = 0.5) +
  scale_size_continuous(
    #trans  = "sqrt",
    range  = c(1, 5),
    breaks = c(50, 100, 500, 1000),
    labels = scales::comma,
    name   = "N Obs:"
  ) +
  geom_point(data = plot_cdist, aes(y=mean_lat, x=mean_lon, size = nobs), pch = 1) +
  
  coord_fixed(1.3) + 
  theme_void() + 
  theme(legend.position = "bottom",
        text = element_text(size = 8))#, #xlim = c(-151, -149), ylim = c(-18, -17)) #+

ggsave("results/sample_map.svg", width = 4, height = 3)
#have a bad point entry


#### ---------- Viz Time -------------
library(lubridate)
by_year <- all_datapoints %>%
  mutate(year = year(DATE)) %>%
  group_by(year) %>%
  summarise(
    nobs      = n(),
    n_studies = n_distinct(STUDY),
    .groups   = "drop"
  )

by_year %>%
  ggplot(aes(x = year, y = nobs, fill = n_studies)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c(
    name   = "N Studies",
    breaks = scales::pretty_breaks(n = 5)
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", y = "N Observations") +
  theme_classic() + 
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.height   = unit(0.125, "in"),   # shrinks the color bar
        legend.key.width    = unit(0.125, "in"))
ggsave(file = "results/studiesperyear.svg", width = 2.5, height = 1.5)


### ----------- Viz Species ----------------

obs_per_spec = all_datapoints %>% group_by(SPECIES) %>%
  summarise(nobs      = n(),
            # n_studies = n_distinct(STUDY),
            .groups   = "drop"
  ) %>% arrange(desc(nobs)) 

meanobs = median(obs_per_spec$nobs)

p <- obs_per_spec %>% ggplot(aes(x=nobs)) +
  geom_histogram(bins = 20, fill = "gray80", col = "black") +
  theme_classic() + 
  xlab("N Observations") +
  ylab("N Species") + 
  geom_vline(xintercept = meanobs, col = "#c0392b") + 
  annotate("text", x=meanobs+200, y= 17, label = paste("Median =", meanobs), col = "#c0392b",
           size = 2) +
  theme(text = element_text(size = 8),
        plot.margin = margin(5, 15, 5, 5, "pt"))
p
ggsave(p, file = "results/Obs_per_species.svg", width = 2.5, height = 1.5)



### ----------- Viz Sites ----------------
all_datapoints %>% mutate(site = paste(LONGITUDE, LATITUDE)) %>%
  group_by(STUDY) %>% summarise(sites = n_distinct(site))

### ----------- Assign Reef Region to Unique Sites ----------------

# Build a table of unique lon/lat combinations from all_datapoints
site_coords <- all_datapoints %>%
  distinct(LONGITUDE, LATITUDE) %>%
  filter(!is.na(LONGITUDE), !is.na(LATITUDE))

aims_region_path <- "data/World_AIMS_Reef-regions-main/finaldata/World_AIMS_Reef-regions_V1.shp"
if (file.exists(aims_region_path)) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for the region join. Install it with install.packages('sf').")
  }
  aims_regions <- sf::st_read(aims_region_path, quiet = TRUE)
  pts <- sf::st_as_sf(site_coords, coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE)
  
  joined <- suppressWarnings(sf::st_join(pts, aims_regions[, "Region"]))
  missing <- is.na(joined$Region)
  if (any(missing)) {
    nearest_idx <- sf::st_nearest_feature(pts[missing, ], aims_regions)
    joined$Region[missing] <- aims_regions$Region[nearest_idx]
    fallback_dist_km <- as.numeric(sf::st_distance(pts[missing, ], aims_regions[nearest_idx, ], by_element = TRUE)) / 1000
    message(sum(missing), " site(s) fell outside all AIMS region polygons -- ",
            "assigned to nearest region instead (max fallback distance: ",
            round(max(fallback_dist_km), 1), " km). Worth spot-checking these.")
  }
  
  site_coords$site_region <- joined$Region
  message("Joined ", nrow(site_coords), " unique site(s) to AIMS reef regions.")
} else {
  message(aims_region_path, " not found -- site_region will be NA for all sites. ",
          "Clone https://github.com/eatlas/World_AIMS_Reef-regions and copy finaldata/* to ",
          "data/World_AIMS_Reef-regions/ to enable this.")
  site_coords$site_region <- NA_character_
}

# Merge the region assignment back onto every row of all_datapoints
all_datapoints <- all_datapoints %>%
  left_join(site_coords, by = c("LONGITUDE", "LATITUDE"))

region = all_datapoints %>% mutate(site = paste(LONGITUDE, LATITUDE), year = year(DATE)) %>%
  group_by(STUDY, site_region) %>% summarise(sites = n_distinct(site), year = median(year)) %>% arrange(year)

writexl::write_xlsx(region, path = "results/Sites_per_Study.xlsx")

