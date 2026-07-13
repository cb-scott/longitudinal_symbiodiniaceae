#########################################################################
# ENVIRONMENTAL DATA EXTRACTION & PCA FOR CORAL SYMBIONT SWITCHING DATASET
#########################################################################
# Pulls satellite env data (SST/DHW, chlorophyll/turbidity, precipitation)
# for each sampling point via rerddap, reduces variables via PCA, and
# generates summary figures (sample map, coverage over time, species obs).
# Intermediate env data files provided in data/envdata/ so the rerddap 
#queries (slow) can be skipped.
#########################################################################


library(tidyverse)
library(future)
source("scripts/helpers.R")
###NOTE: I highly recommend you use read_excel for this. For my own sanity,
#I collated data in excel, including the *dates* for each sampling point
#These dates are in excel format, which reads into R as a POSXIct
#(and gives time in seconds, which necessitates some of the funky arithmetic below)
#Converting dates to an as.Date object in R will fix this; I've done this in some later scripts

all_datapoints = readxl::read_excel("data/Long_Symbiont_Data_Manual.xlsx")

# ==============================================================================
# GET ENVIRONMENTAL DATA
# ==============================================================================
library(rerddap)

#get unique lon, lats dates
library(zoo)
to_grab = all_datapoints %>% select(LONGITUDE, LATITUDE, DATE) %>% distinct() %>% 
  mutate(to_date = 
           case_when(day(DATE) == 1 ~ DATE + (30*86400),
                     day(DATE) != 1 ~ DATE),
         from_date = case_when(day(DATE) != 1 ~ DATE - (30*86400),
                               day(DATE) == 1 ~ DATE) ) #the dates are in seconds, so convert
#For studies that only report the month, grab that month's mean. For studies that reported an exact sampling date, grab all data a month prior
#In Long_Symbiont_Data_Manual, if the exact date wasn't given or was reported as "early MONTH", I recorded it as the first of the month
#If an uncertain timerframe was given (e.g., mid-month or end-month) I recorded the date as the 15th or the 30th, respectively

#####---- Extract Data from gridded climate products------ 
###NOTE: This code takes a little while to run (~10-20 minutes per variable). 
#I've provided the intermediate products in the folder 'data/envdata' which can be read in below

library(furrr)
library(purrr)


##--------- NOAA SST-RELATED PRODUCTS --------------
plan(multisession, workers = 11)  # adjust worker count for your computer; I think that sometimes parallel code is funky on Windows too
#I ran on Mac Tahoe 26.4.1

future_map(1:nrow(to_grab), function(row) {
  outfile <- paste0("data/envdata/sst_", row, "env.txt")
  if (file.exists(outfile) && file.size(outfile) > 0) return(NULL)
  
  tryCatch({
    tmp <- rerddap::griddap(
      "NOAA_DHW",
      time      = c(as.character(to_grab$from_date[row]), as.character(to_grab$to_date[row])),
      latitude  = c(to_grab$LATITUDE[row]  - 0.01, to_grab$LATITUDE[row]  + 0.01),
      longitude = c(to_grab$LONGITUDE[row] - 0.01, to_grab$LONGITUDE[row] + 0.01),
      fields    = c("CRW_SST", "CRW_DHW", "CRW_SSTANOMALY")
    )
    tmp <- tmp$data
    
    if(length(is.na(tmp$CRW_SST)[is.na(tmp$CRW_SST) == T]) == length(tmp$CRW_SST)){
      tmp <- rerddap::griddap(
        "NOAA_DHW",
        time      = c(as.character(to_grab$from_date[row]), as.character(to_grab$to_date[row])),
        latitude  = c(to_grab$LATITUDE[row]  - 0.05, to_grab$LATITUDE[row]  + 0.05),
        longitude = c(to_grab$LONGITUDE[row] - 0.05, to_grab$LONGITUDE[row] + 0.05),
        fields    = c("CRW_SST", "CRW_DHW", "CRW_SSTANOMALY")
      )
      tmp <- tmp$data
    }
    
    
    env_data <- data.frame(
      lat       = to_grab$LATITUDE[row],
      lon       = to_grab$LONGITUDE[row],
      date      = to_grab$DATE[row],
      min_SST   = min(tmp$CRW_SST,         na.rm = TRUE),
      max_SST   = max(tmp$CRW_SST,         na.rm = TRUE),
      sd_SST    = sd(tmp$CRW_SST,          na.rm = TRUE),
      mean_SST  = mean(tmp$CRW_SST,        na.rm = TRUE),
      max_DHW   = max(tmp$CRW_DHW,         na.rm = TRUE),
      min_DHW   = min(tmp$CRW_DHW,         na.rm = TRUE),
      mean_DHW  = mean(tmp$CRW_DHW,        na.rm = TRUE),
      mean_ANOM = mean(tmp$CRW_SSTANOMALY, na.rm = TRUE),
      sd_ANOM   = sd(tmp$CRW_SSTANOMALY,   na.rm = TRUE),
      min_ANOM  = min(tmp$CRW_SSTANOMALY,  na.rm = TRUE),
      max_ANOM  = max(tmp$CRW_SSTANOMALY,  na.rm = TRUE),
      cum_ANOM = sum(tmp$CRW_SSTANOM, na.rm = T)
    )
    
    Sys.sleep(1)
    
    write.table(env_data,
                file      = paste0("data/envdata/sst_", row, "env.txt"),
                col.names = FALSE,
                row.names = FALSE)
    #env_data
    
  }, error = function(e) {
    message(sprintf("%s — row %d failed: %s", Sys.time(), row, conditionMessage(e)))
    NULL
  })
}, .options = furrr_options(seed = TRUE))

plan(sequential)

# ------------- OCEAN COLOR FOR TURBIDITY AND CHLA ------------

plan(multisession, workers = 11)  # adjust worker count
#ACTUALLY ONLY NEED OCEAN COLOR FOR BIOGEOCHEMICAL
future_map(1:nrow(to_grab), function(row) {
  outfile <- paste0("data/envdata/color_", row, "env.txt")
  if (file.exists(outfile) && file.size(outfile) > 0) return(NULL)
  
  tryCatch({
    tmp <- rerddap::griddap(
      "occci_V6_8day_4km",
      url = "https://comet.nefsc.noaa.gov/erddap/",
      time      = c(as.character(to_grab$from_date[row]), as.character(to_grab$to_date[row])),
      latitude  = c(to_grab$LATITUDE[row]  - 0.01, to_grab$LATITUDE[row]  + 0.01),
      longitude = c(to_grab$LONGITUDE[row] - 0.01, to_grab$LONGITUDE[row] + 0.01),
      fields    = c("chlor_a", "kd_490")
    )
    tmp <- tmp$data
    
    if(length(is.na(tmp$chlor_a)[is.na(tmp$chlor_a) == T]) == length(tmp$chlor_a)){
      tmp <- rerddap::griddap(
        "occci_V6_8day_4km",
        url = "https://comet.nefsc.noaa.gov/erddap/",
        time      = c(as.character(to_grab$from_date[row]), as.character(to_grab$to_date[row])),
        latitude  = c(to_grab$LATITUDE[row]  - 0.051, to_grab$LATITUDE[row]  + 0.051),
        longitude = c(to_grab$LONGITUDE[row] - 0.051, to_grab$LONGITUDE[row] + 0.051),
        fields    = c("chlor_a", "kd_490")
      )
      tmp <- tmp$data
    }
    
    env_data <- data.frame(
      lat       = to_grab$LATITUDE[row],
      lon       = to_grab$LONGITUDE[row],
      date      = to_grab$DATE[row],
      min_ChlA   = min(tmp$chlor_a,         na.rm = TRUE),
      max_ChlA   = max(tmp$chlor_a,         na.rm = TRUE),
      sd_ChlA    = sd(tmp$chlor_a,          na.rm = TRUE),
      mean_ChlA = mean(tmp$chlor_a, na.rm = T), 
      min_turb   = min(tmp$kd_490,         na.rm = TRUE),
      max_turb   = max(tmp$kd_490,         na.rm = TRUE),
      sd_turb    = sd(tmp$kd_490,          na.rm = TRUE),
      mean_turb = mean(tmp$kd_490, na.rm = T)
    )
    
    Sys.sleep(2)
    
    write.table(env_data,
                file      = paste0("data/envdata/color_", row, "env.txt"),
                col.names = FALSE,
                row.names = FALSE)
    #env_data
    
  }, error = function(e) {
    message(sprintf("%s — row %d failed: %s", Sys.time(), row, conditionMessage(e)))
    NULL
  })
}, .options = furrr_options(seed = TRUE))

plan(sequential)

#----------------- PRECIPITATION ---------------

plan(multisession, workers = 11)  # adjust worker count

future_map(1:nrow(to_grab), function(row) {
  outfile <- paste0("data/envdata/precip_", row, "env.txt")
  if (file.exists(outfile) && file.size(outfile) > 0) return(NULL)
  
  tryCatch({
    tmp <- rerddap::griddap(
      "IMERG_monthly_global_precip",
      url = "https://erddap.marine.usf.edu/erddap/",
      time      = c(as.character(to_grab$from_date[row]), as.character(to_grab$to_date[row])),
      latitude  = c(to_grab$LATITUDE[row]  - 0.01, to_grab$LATITUDE[row]  + 0.01),
      longitude = c(to_grab$LONGITUDE[row] - 0.01, to_grab$LONGITUDE[row] + 0.01),
      fields    = c("precipitation")
    )
    tmp <- tmp$data
    
    env_data <- data.frame(
      lat       = to_grab$LATITUDE[row],
      lon       = to_grab$LONGITUDE[row],
      date      = to_grab$DATE[row],
      min_precip   = min(tmp$precipitation,         na.rm = TRUE),
      max_precip   = max(tmp$precipitation,         na.rm = TRUE),
      sd_precip    = sd(tmp$precipitation,          na.rm = TRUE),
      mean_precip = mean(tmp$precipitation, na.rm = T)
    )
    
    Sys.sleep(2)
    
    write.table(env_data,
                file      = paste0("data/envdata/precip_", row, "env.txt"),
                col.names = FALSE,
                row.names = FALSE)
    #env_data
    
  }, error = function(e) {
    message(sprintf("%s — row %d failed: %s", Sys.time(), row, conditionMessage(e)))
    NULL
  })
}, .options = furrr_options(seed = TRUE))

plan(sequential)


# ==============================================================================
# GENERATE AND SAVE ENV PCs
# ==============================================================================

#read files in 
#NOTE YOU CAN START HERE WITH PROVIDED INTERMEDIATE DATA
library(vegan)

sst.files <- list.files("data/envdata", full.names = TRUE, pattern = "sst")

sst.data <- do.call(rbind, lapply(sst.files, read.table, header = FALSE, nrows = 1))
colnames(sst.data) = c("lat", "lon", "date", 
                       "min_SST", "max_SST", "sd_SST", "mean_SST",
                       "max_DHW", "min_DHW", "mean_DHW", 
                       "mean_ANOM", "sd_ANOM", "min_ANOM", "max_ANOM", "cum_ANOM")

precip.files <- list.files("data/envdata/", full.names = T, pattern = "precip")
precip.data <- do.call(rbind, lapply(precip.files, read.table, header = FALSE, nrows = 1))
colnames(precip.data) = c("lat", "lon", "date", 
                          "min_precip", "max_precip", "sd_precip", "mean_precip")

color.files <- list.files("data/envdata/", full.names = T, pattern = "color")
color.data <- do.call(rbind, lapply(color.files, read.table, header = FALSE, nrows = 1))

colnames(color.data) = c("lat", "lon", "date", 
                         "min_chla", "max_chla", "sd_chla", "mean_chla",
                         "min_turb", "max_turb", "sd_turb", "mean_turb")

join.data = left_join(left_join(sst.data, precip.data), color.data)

means.only = join.data %>% select( mean_ANOM, min_ANOM, max_ANOM,
                                   cum_ANOM, mean_DHW, sd_ANOM,
                                   mean_SST, min_SST, max_SST, sd_SST,
                                   mean_precip, min_precip, max_precip, sd_precip,
                                   mean_chla, min_chla, max_chla, sd_chla,
                                   mean_turb, min_turb, max_turb, sd_turb) %>%
  mutate(across(everything(), ~ if_else(is.infinite(.), NA, .))) %>%
  mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))

mean.stand = decostand(means.only, method = "standardize")
env.pc = prcomp(mean.stand, center = TRUE, scale = TRUE)
env.score = scores(env.pc)

dev.off()
screeplot(env.pc, bstick = T)
varexp=eigenvals(env.pc)/sum(eigenvals(env.pc))
varexp
plot(env.score[,1:2])
biplot(env.pc)
#make a nice biplot

cbind(join.data[,c("lat", "lon", "date")], data.frame(env.score)) %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d")) %>%
  mutate(mon = month(date)) %>%
  ggplot(aes(x=PC1, y=PC2, col = mon)) + geom_point() +
  scale_color_viridis_c()

env.data.final = cbind(join.data, data.frame(env.score)) %>%
  mutate(DATE = as.Date(date, format="%Y-%m-%d")) %>% select(!date)

plot(env.data.final$PC1, env.data.final$mean_ANOM)
#there's a problem with the data in moorea from 2011-2012. Changing the lon/lats slightly did not fix.
tiff("results/pc_loadings.tiff", width = 3.25, height = 3.25, units = "in", res = 600)
pheatmap::pheatmap(env.pc$rotation[,1:3], cluster_cols = F)
dev.off()
save(env.data.final, file = "data/EnvData_PCs.Rdata")


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
  labs(
    x = paste0("PC1 (", round(varexp[1], 2)*100, "%)"),
    y = paste0("PC2 (", round(varexp[2], 2)*100, "%)")
  ) +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(size = 8),
    axis.line   = element_line(color = "gray40"),
    axis.ticks  = element_line(color = "gray40"),
    plot.margin = margin(10, 15, 10, 10),
    legend.position = "right"
  )
ggsave("results/env.pcs.tiff", width = 3.25, height = 3.25, units = "in", dpi = 1000)


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
  coord_fixed(1.3) + 
  theme_void() + 
  theme(legend.position = "bottom",
        text = element_text(size = 8))#, #xlim = c(-151, -149), ylim = c(-18, -17)) #+
ggsave("results/sample_map.png", width = 4, height = 3, units = "in", dpi = 1000)

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
ggsave(file = "results/studiesperyear.png", width = 2.5, height = 1.5, units = "in", dpi = 1000)

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
ggsave(p, file = "results/Obs_per_species.tiff", width = 2.5, height = 1.5, units = "in", dpi = 1000)
