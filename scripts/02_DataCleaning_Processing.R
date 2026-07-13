#########################################################################
# DATASET SUMMARY STATS, SPECIES COMPLEX ASSIGNMENT, & SYMBIONT SWITCH CALLING
#########################################################################
# Loads env PCs + manual symbiont data, generates basic dataset summary
# stats and a by-study summary table (results/StudySummary.xlsx), collapses
# select species into recognized species complexes, joins env PC data by
# lat/lon/date, and calls genus-level symbiont switching events per colony
# (NONE/MINOR/MAJOR/SUBTYPE_SWITCH) based on sequential sampling within
# each colony.
#########################################################################


library(tidyverse)
library(phyloseq)

### ------------- Load in data -------------
load("data/EnvData_PCs.Rdata")
all_datapoints = readxl::read_excel("data/Long_Symbiont_Data_Manual.xlsx")
#basic stats
#number colonies represented
paste("Number Unique Colonies =", length(unique(all_datapoints$COLONY_ID)))
paste("Number Unique Sites =", length(unique(paste(all_datapoints$LONGITUDE, all_datapoints$LATITUDE))))
paste("Number Studies =", length(unique(all_datapoints$STUDY)))
paste("Number Species Represented =", length(unique(all_datapoints$SPECIES)))

# ------------ Get a by study summary- -----------------
study_sum = all_datapoints %>% group_by(STUDY) %>%
  summarise(lat_min = round(min(LATITUDE), 2),
            lat_max = round(max(LATITUDE), 2),
            lon_min = round(min(LONGITUDE), 2),
            lon_max = round(max(LONGITUDE), 2),
            n_specs = n_distinct(SPECIES),
            n_colonies = n_distinct(COLONY_ID),
            time_min = year(as.Date(min(DATE))), 
            time_max = year(as.Date(max(DATE))),
            n_obs = n()) %>%
  mutate(lat_range = ifelse(lat_min > 0, paste(lat_min, "N", "-", lat_max, "N"),
                                          paste(-1*lat_min, "S", "-", -1*lat_min, "S")), 
         lon_range = ifelse(lon_min > 0, paste(lon_min, "E -", lon_max, "E", sep = " "),
                                        paste(-1*lon_min, "W -", -1*lon_max, "W", sep = " ")),
         date_range= ifelse(time_min != time_max, paste0(time_min, "-", time_max), time_min))
writexl::write_xlsx(study_sum, path = "results/StudySummary.xlsx")

#### ----- ASSIGN SPECIES COMPLEXES ------------
#revise for species complexes
#only consider really well recognized species complexes
#Ultimately, I do not use most of these species complexes downstream, but rather let the small phylogenetic distances account for this

all_datapoints = all_datapoints %>% 
  mutate(HOST_SPEC_REVISED = case_when(
    SPECIES %in% c("PORITES_LUTEA", "PORITES_LOBATA", "PORITES_SP.") ~ "MASSIVE PORITES SPECIES COMPLEX",
    SPECIES %in% c("POCILLOPORA_ACUTA", "POCILLOPORA_DAMICORNIS", "POCILLOPORA_VERRUCOSA") ~ "POCILLOPORA SPECIES COMPLEX",
    .default = SPECIES
  )) 


all_datapoints = left_join(all_datapoints, env.data.final, by = c("DATE"= "DATE", "LONGITUDE"="lon", "LATITUDE" = "lat"))

#### ---- DETERMINE SWITCHING EVENTS ----------
genlevel <- all_datapoints %>%
  separate(GENUS, into = c("Gen1", "Gen2", "Gen3"), sep = "_", remove = FALSE) %>%
  filter(!is.na(GENUS)) %>%
  group_by(COLONY_ID) %>%
  arrange(DATE, .by_group = TRUE) %>%
  mutate(
    lag_GENUS   = lag(GENUS, 1),
    lag_SUBTYPE = lag(SUBTYPE, 1),
    
    # Check if any Gen component of current row appears in lagged GENUS
    is_minor = {
      curr_parts <- pmap(list(Gen1, Gen2, Gen3), function(g1, g2, g3) {
        na.omit(c(g1, g2, g3))
      })
      map2_lgl(curr_parts, lag_GENUS, function(parts, lag_g) {
        if (is.na(lag_g) || length(parts) == 0) return(FALSE)  # <-- key guard
        any(parts != "" & str_detect(lag_g, paste(parts, collapse = "|")))
      })
    },
    
    GenShift = case_when(
      is.na(lag_GENUS)   ~ NA_character_,   # first observation
      GENUS == lag_GENUS ~ "NONE",
      is_minor           ~ "MINOR",
      TRUE               ~ "MAJOR"
    ),
    
    GenShift = case_when(
      GenShift == "NONE" &
        !is.na(SUBTYPE) &
        !is.na(lag_SUBTYPE) &
        SUBTYPE != lag_SUBTYPE ~ "SUBTYPE_SWITCH",
      TRUE ~ GenShift
    )
  ) %>%
  select(-lag_GENUS, -lag_SUBTYPE, -is_minor)

save(genlevel, all_datapoints, file = "data/LongZooxInputData_wEnv.Rdata")



