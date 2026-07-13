# Data from: Effects of an extreme marine heatwave on *Siderastrea siderea* in Panama

Dataset includes files for data associated with manuscript by Glass et al.

## Data and file structure

### File 1: Health_status_data.csv

This file is for data pertaining to the health of S. siderea colonies over time.

This file has headers, which are:
* Year = year during which data were collected
* Site = geographical location (i.e., study site) where colony was located
* Coral_ID = unique identifying number for each coral colony
* Bleaching_status = visually assigned coral bleaching status (pigmented, partially bleached, or fully bleached)
* Percent_mortality_since_2022 = colony-level tissue loss expressed as percent since 2022, scored visually in 10% intervals
* Signs_of_disease = whether colony displayed dark purple tissue discoloration indicative of dark spot syndrome

### Files 2 and 3: PostMED_abs_abundance.csv and PostMED_relative_abundance.csv

This file is for data pertaining to the absolute (abs) and relative abundances of algal symbiont sequences in corals quantified through ITS2 amplicon sequencing.

This file has headers, which are:
* Sample_ID = unique sample ID
* Year = year during which data were collected
* Coral_ID = unique identifying number for each coral colony
* Site = geographical location (i.e., study site) where colony was located
* A3–D2u = names of each symbiont sequence

Other columns were not used in analyses.

### File 4: Site_temperature_data.csv

This file contains temperature data from the study sites recorded using in situ loggers.

This file has headers, which are:
* Site = the name of the study site from which data were collected
* Date_time = the date and time of the temperature measurement
* Temperature_C = seawater temperature in °C

Other columns were not used in analyses.

### File 5: 16S_metadata.csv

This file contains sample metadata associated with 16S rRNA gene sequencing data.

This file has headers, which are:
* Sample_ID = unique sample ID
* Year = year during which data were collected
* Coral_ID = unique identifying number for each coral colony
* Site = geographical location (i.e., study site) where colony was located

Other columns were not used in analyses.

### File 6: Bocas_weather_station_temperature_data.csv

This file contains temperature data from a weather station present in Bocas del Toro near the Smithsonian Tropical Research Institute.

This file has headers, which are: 
* Date_time = the date and time of the temperature measurement
* Temperature_C = seawater temperature in °C

## Sharing/Access information

Links to other publicly accessible locations of the data:

* NA

Data was derived from the following sources:

* Files 2 and 3 were derived from analysis of raw ITS2 amplicon sequencing data using SymPortal (Hume et al., 2019 Molecular Ecology)
* File 4 was derived from data compiled by Gantt et al. (2025) Coral Reefs
* File 6 was derived from data collected by the Smithsonian Tropical Research Institute and downloaded from the Institute website ([https://stri.si.edu/](https://stri.si.edu/)).

## Code/Software

### File 1: Panama_heatwave_symbioses.Rmd

This file contains code in R Markdown format for producing analyses and visualizations of the data.
