#########################################################################
# CLEAN SUPPLEMENTAL DATA FROM PUBLISHED STUDIES
#########################################################################
# NOTE: YOU SHOULD NOT HAVE TO RUN THIS TO REPLICATE RESULTS
# CLEANED DATA IS GIVEN IN Long_Symbiont_Data_Manual.xlsx

# Each section below pulls in a different published study's raw
# supplemental data (SymPortal profiles/counts, ITS2 OTU tables, or
# per-study metadata sheets) and reshapes it into the common schema used
# in Long_Symbiont_Data_Manual.xlsx:
#   SPECIES, COLONY_ID, STUDY, DATE, GENUS, SUBTYPE, BLEACHING_OBS, DEAD,
#   DEPTH_M, LATITUDE, LONGITUDE
# Each block ends by writing a cleaned .xlsx, which then gets manually
# pasted into Long_Symbiont_Data_Manual.xlsx.
# NOTE: only colonies with >1 timepoint (nobs > 1) are kept, since we
# need repeat sampling to detect a switching event.

#Not every study is contained here, as many reported data in figure form
#in the text, which was manually extracted.
#########################################################################

library(tidyverse)
library(writexl)
#Van nyatten 2025
vn.meta <- read.csv("data/van_nynatten_colony_metadata.csv")
vn.seqid <- read.csv("data/van_nynatten_sample_metadata.csv")
vn.sites <- read.csv("data/van_nynatten_survey_metadata.csv")
vn.exps <- read.csv("data/van_nynatten_expeditions.csv") %>% select(Expedition.code, Mean.date) %>% data.frame()
vn.exps <- vn.exps %>%
  mutate(Expedition.code = ifelse(Expedition.code == "2023", "2023b", Expedition.code))
vn.coms <- read.table("data/van_nynatten20250527T170756.profiles.absolute.abund_and_meta.txt", sep = '\t', header =F, skip = 7)
vn.head <- read.table("data/van_nynatten20250527T170756.profiles.absolute.abund_and_meta.txt", sep = '\t', header =F, skip = 2) %>% head(1)
colnames(vn.coms)=vn.head
colnames(vn.coms)[1:2] = c("sym_id", "sampleID")
vn.coms = janitor::clean_names(vn.coms)
vn.coms = vn.coms %>% dplyr::select(!sym_id)
vn.all = vn.coms %>% pivot_longer(!sample_id, names_to="zoox_type", values_to="count") %>%
  left_join(vn.seqid, by = c("sample_id"="library_ID")) %>% 
  left_join(vn.exps, by = c("expedition"="Expedition.code"))

dominance = vn.all %>% 
  filter(sample_id != "") %>% 
  group_by(sample_id) %>% 
  mutate(count = as.numeric(count)) %>%
  mutate(relabun = count/sum(count)) %>%
  filter(relabun > .1) %>% separate(zoox_type, into = c("maj.zoox"), sep = "_") %>%
  mutate(maj.zoox = toupper(gsub("[a-z]+\\b", "", maj.zoox))) %>%
  mutate(GENUS = substr(maj.zoox, 1, 1)) %>%
  mutate(COLONY_ID=paste0("VN_", substr(toupper(host_genus), 1, 1),
                          substr(toupper(host_species), 1, 3), "_",
                         coral_tag)) %>%
  mutate(SPECIES = paste(toupper(host_genus), toupper(host_species), sep = "_")) %>%
  left_join(vn.meta, by=c("coral_tag" = "coral_tag", "expedition" = "expedition")) %>%
  mutate(DEAD = ifelse(grepl("dead", health_status) == T, "Y", "N")) %>%
  mutate(BLEACHING_OBS = NA) %>%
  left_join(vn.sites, by = c("site_name" = "site_name")) %>%
  group_by(sample_id) %>%
  mutate(SUBTYPE = paste(unique(maj.zoox), collapse = "_")) %>%
  mutate(GENUS = paste(unique(GENUS), collapse = "_")) %>%
  ungroup() %>%
  select(-maj.zoox) %>%
  distinct() %>%
  mutate(DATE = Mean.date) %>%
  mutate(DEPTH_M = NA) %>%
  mutate(LATITUDE = lat, 
           LONGITUDE = lon, 
         STUDY = "Van Nynatten et al., 2025") %>%
  group_by(COLONY_ID) %>% mutate(counts = n()) %>%
  filter(counts > 1) %>% 
  select(SPECIES, COLONY_ID, STUDY, DATE, GENUS, SUBTYPE, BLEACHING_OBS, DEAD, DEPTH_M, LATITUDE, LONGITUDE)

write_xlsx(dominance, "data/cleaned_vn_2025.xlsx")


#Get Manzello Data -- just need teh bleaching ids
x=load("data/clones_meta_me6.RData")
write_xlsx(me6 %>% select(Tag.Number, bleach) %>% mutate(BLEACH_OBS=ifelse(bleach %in% c("1","2"), "Y", "N")),path = "data/Manzello_BleachScores.xlsx")


##clean up Claar 2020 data
meta_in = read.csv("data/Claar_Metadata.csv") %>% 
  mutate(new_date = as.Date(paste0("01 ", collection_date), format = "%d %Y %B"))
meta_in_dates1 = meta_in %>% filter(is.na(new_date) == F)
meta_in_dates2 = meta_in %>% filter(is.na(new_date) == T) %>% 
  mutate(new_date = as.Date(as.numeric(collection_date), origin = "1899-12-30"))
meta_in_all = rbind(meta_in_dates1, meta_in_dates2) %>%
  select(sample_name, host_genus, host_species, collection_latitude, collection_longitude, new_date, collection_depth, coral_tag)
zoox_in = read.csv("data/Claar_RawRead_matrix.csv") %>% dplyr::select(!X) %>% 
  pivot_longer(!sample_name, names_to = "Clade", values_to="count") %>%
  left_join(meta_in_all)

filtered = zoox_in %>% group_by(sample_name) %>%
  filter(Clade != "Sequence_count") %>%
  mutate(relabun = count/sum(count)) %>%
  filter(relabun > .1) %>%
  mutate(SUBTYPE=sub("^([^a-z.]*).*", "\\1", Clade)) %>% # a little dirty, maybe I want to keep the full symportal types? won't use them for most things anyways
  mutate(TYPE = substr(SUBTYPE, 1, 1)) %>%
  group_by(sample_name) %>%
  mutate(SUBTYPE = paste(unique(SUBTYPE), collapse = "_")) %>%
  mutate(GENUS = paste(unique(TYPE), collapse = "_")) %>%
  select(!c(TYPE, Clade, count, relabun)) %>%
  ungroup() %>%
  mutate(DEPTH = 10) %>%
  mutate(LONGITUDE = -1*as.numeric(gsub(" W", "", collection_longitude))) %>%
  mutate(LATITUDE = -1*as.numeric(gsub(" N", "", collection_latitude)),
         SPECIES = paste(toupper(host_genus), toupper(host_species), sep= "_"),
         COLONY_ID = paste("CLAAR", coral_tag, sep = "_")) %>%
 select(SPECIES, COLONY_ID, new_date, GENUS, SUBTYPE, DEPTH, LATITUDE, LONGITUDE ) %>%
  distinct() %>% group_by(COLONY_ID) %>%
  mutate(nobs = n()) %>% filter(nobs > 1) 

write_xlsx(filtered, "data/Claar_Cleaned.xlsx")


### starko
load("data/Starko_Porites_DIV.RData")
#full metadata is included, thank you!
meta = sample_data(Porites_DIV)
#and it has a proper OTU table at the symbiont level
p.otus = data.frame(otu_table(Porites_DIV)) #this is relative abundance
p.otus = p.otus %>% mutate(ID = rownames(p.otus)) %>% pivot_longer(!ID, names_to = "zoox", values_to = "relabun") %>%
  mutate(short.sub = sub("^([^a-z.]*).*", "\\1", zoox)) %>%
  filter(grepl("^X", zoox) == F) %>%
  group_by(ID, short.sub) %>% summarise(relabun = sum(relabun)) %>%
  filter(relabun > .1) %>% 
  mutate(GENUS = substr(short.sub, 1,1))
p.otus

#get metadata
meta.df = data.frame(meta) %>% mutate(ID = rownames(meta)) %>% select(ID, host_genus, host_species,
                                                                      collection_latitude, collection_longitude,
                                                                      collection_date, collection_depth, coral_tag) %>%
  mutate(DEPTH = 10) %>%
  mutate(DATE = as.Date(collection_date, format = "%Y%m%d"))

p.meta = left_join(p.otus, meta.df) %>%
  group_by(coral_tag) %>%
  mutate(SUBTYPE = paste(sort(unique(short.sub)), collapse = "_")) %>%
  mutate(GENUS = paste(sort(unique(GENUS)), collapse = "_")) %>%
  select(!c(short.sub, relabun)) %>%
  ungroup() %>% distinct()
p.meta$year = year(p.meta$DATE)

mortality_data2015 = read.csv("data/2015_Starko_Mortality_and_bleaching.csv") %>% select(Bleaching_McClanahan5, coral_tag) %>%
  mutate(year = 2015)
mortalitiy_data2016 = read.csv("data/2016_Starko_Mortality_and_bleaching.csv")%>% select(Bleaching_McClanahan5, coral_tag) %>%
  mutate(year = 2016)

mort.all=rbind(mortality_data2015, mortalitiy_data2016)
p.wmort = left_join(p.meta, mort.all, by = c("coral_tag" = "coral_tag", "year"="year")) %>%
  mutate(BLEACHING_OBS = ifelse(Bleaching_McClanahan5 >=3, "Y", "N"))

stark.out = p.wmort %>% select(host_genus, host_species, coral_tag, DATE, GENUS, SUBTYPE, BLEACHING_OBS,
                   DEPTH, collection_latitude, collection_longitude) %>%
  mutate(SPECIES = paste(toupper(host_genus), toupper(host_species), sep = "_")) %>%
  mutate(COLONY_ID = paste("STARK_PLOB", coral_tag, sep = "_")) %>%
  select(!c(host_genus, host_species, coral_tag)) %>% distinct() %>%
  group_by(COLONY_ID) %>% mutate(nobs = n()) %>% filter(nobs > 1)

write_xlsx(stark.out, path = "data/CleanedStarko2023.xlsx")


#### Grab Cunning et al., 2024
cun.sym = read.table("data/Cunning_Symportal_111_20200624_2020-06-24_03-40-50.238182.seqs.absolute.abund_and_meta.txt", sep = "\t", header = T) %>%
  select(sample_name,
         collection_latitude, collection_longitude,
         collection_depth, collection_date,
         host_genus, host_species,
         A4:G3b)

#peep his metadata
cun.meta = read.csv("data/Cunning_Meta_STJ_Metadata_2017_2019.csv") #redundant, all of it is contained in teh symportal metadata

cun.all = cun.sym %>% pivot_longer(A4:G3b, names_to = "zoox", values_to = "count") %>%
  group_by(sample_name) %>%
  mutate(short.sub = sub("^([^a-z.]*).*", "\\1", zoox)) %>%
  filter(grepl("^X", zoox) == F) %>%
  group_by(sample_name) %>%
  mutate(relabun = count/sum(count)) %>%
  filter(relabun > .1) %>% 
  mutate(GENUS = substr(short.sub, 1,1)) %>% 
  separate(sample_name, into = c("SITE", "SPEC", "ID", "YEAR"), sep = "_", remove = F) %>%
  mutate(COLONY_ID = paste("CUN", SITE, SPEC, ID, sep = "_")) %>%
  select(!c( zoox, count, relabun)) %>%
  group_by(sample_name) %>%
  mutate(SUBTYPE = paste(sort(unique(short.sub)), collapse = "_")) %>%
  mutate(GENUS = paste(sort(unique(GENUS)), collapse = "_")) %>%
  select(!c(short.sub)) %>%
  ungroup() %>% distinct() %>%
  mutate(DATE = as.Date(collection_date, format = "%Y.%m.%d"),
         LONGITUDE = -1*collection_longitude,
         SPECIES = paste(toupper(host_genus), toupper(host_species), sep="_")) %>%
  group_by(COLONY_ID) %>% mutate(nobs=n()) %>% filter(nobs > 1)

write_xlsx(cun.all, path = "data/CleanCunning2024.xlsx")


####### TERREANO 2024
meta = readxl::read_excel("data/Tarraneo_Meta.xlsx")
meta = meta[-1,]
meta = janitor::clean_names(meta)

ter.sym = readxl::read_excel("data/Tarreano_Symportal.xlsx") %>% pivot_longer(!ID, names_to = "zoox", values_to = "relabun") %>%
  mutate(short.sub = sub("^([^a-z.]*).*", "\\1", zoox)) %>% group_by(short.sub, ID) %>%
  mutate(relabun = as.numeric(relabun)) %>%
  summarise(relabun = sum(relabun)) %>% drop_na() %>%
  filter(relabun > .1)

ter.out = ter.sym %>% left_join(meta, by = c("ID" = "sample_number")) %>%
  mutate(DATE = as.Date(paste(sampling_time_in_2016, "01", "2016"), format = "%B %d %Y")) %>% 
  mutate(BLEACHING_OBS = ifelse(bleaching_prevalence == "NB", "N",
                                ifelse(bleaching_prevalence %in% c("PB", "B"), "Y", bleaching_prevalence))) %>%
  mutate(DEAD = ifelse(bleaching_prevalence == "DEAD", "Y", "N")) %>%
  mutate(GENUS = substr(short.sub, 1,1)) %>%
  group_by(ID) %>%
  mutate(SUBTYPE = paste(sort(unique(short.sub)), collapse = "_")) %>%
  mutate(GENUS = paste(sort(unique(GENUS)), collapse = "_")) %>%
  select(!c(short.sub, relabun)) %>%
  ungroup() %>% distinct() %>%
  mutate(COLONY_ID = paste("TERR", toupper(substr(GENUS, 1, 1)),
                           toupper(substr(x4, 1, 3)), tagged_colony, sep = "_")) %>%
  mutate(SPECIES = paste(toupper(genus), toupper(x4), sep = "_")) %>%
  group_by(COLONY_ID) %>% mutate(nobs=n()) %>% filter(nobs>1) %>% arrange(sampling_locality)

write_xlsx(ter.out, "Terraneo_Clean.xlsx")  

##### HUFFMYER ET AL 2025

meta = read.csv("data/Huffmyer_Meta_20230425_Cunning_Symportal_datasheet.csv") 
huff.sym = read.table("data/Huffmyer-Symportal.txt", sep = "\t", header = T)

huff.proc = huff.sym %>% pivot_longer(A1ee:G3ab, names_to = "zoox", values_to="count") %>%
  mutate(short.sub = sub("^([^a-z.]*).*", "\\1", zoox)) %>% 
  group_by(sample_name) %>%
  mutate(relabun = count/sum(count)) %>%
  select(!zoox) %>% distinct() %>%
  group_by(sample_name, short.sub) %>%
  mutate(relabun = sum(relabun)) %>%
  filter(relabun > .1) %>% 
  mutate(GENUS = substr(short.sub, 1,1)) %>%
  group_by(sample_name) %>%
  mutate(SUBTYPE = paste(sort(unique(short.sub)), collapse = "_")) %>%
  mutate(GENUS = paste(sort(unique(GENUS)), collapse = "_")) %>%
  select(!c(short.sub, count, relabun)) %>%
  ungroup() %>% distinct()%>%
  select(sample_name, 
         collection_longitude, collection_latitude,
         collection_date, 
         collection_depth,
         host_genus, host_species, GENUS,
         SUBTYPE) %>%
  filter(collection_longitude != 999) %>%
  mutate(DATE = as.Date(collection_date, format = "%Y%m%d")) %>%
  separate(sample_name, into = c("spec", "ID", "tp")) %>%
  mutate(COLONY_ID = paste("HUFF", spec, ID, sep = "_")) %>%
  mutate(SPECIES = paste(toupper(host_genus), toupper(host_species), sep = "_")) %>%
  filter(GENUS != "C_X") %>% ungroup() %>%
  group_by(COLONY_ID) %>% mutate(nobs = n()) %>% filter(nobs > 1)

write_xlsx(huff.proc, "HuffClean.xlsx")


###### LEINBACH et al., 2023
lbach = read.csv("data/Leinbach_profs.csv") %>% pivot_longer(A1:D6, names_to = "zoox", values_to = "count")
#May = May 2019
#August = August 2019
#October = October 2019
lbach.clean = lbach %>% group_by(Sample, month) %>%
  mutate(short.sub = sub("^([^a-z.]*).*", "\\1", zoox)) %>%
  group_by(Sample, depth_m, reef_zone, month, response, bleaching_status, Individual) %>%
  mutate(relabun = count/sum(count)) %>%
  group_by(Sample, depth_m, reef_zone, month, response, bleaching_status, Individual, short.sub) %>%
  summarise(relabun = sum(relabun)) %>% filter(relabun > .1) %>%
  mutate(DATE = as.Date(paste(month, "01", "2019"), format = "%B %d %Y")) %>%
  mutate(GENUS = substr(short.sub, 1,1)) %>%
  group_by(Sample) %>%
  mutate(SUBTYPE = paste(sort(unique(short.sub)), collapse = "_")) %>%
  mutate(GENUS = paste(sort(unique(GENUS)), collapse = "_")) %>%
  select(!c(short.sub, relabun)) %>%
  ungroup() %>% distinct() %>%
  mutate(COLONY_ID = paste("LBACH_AHYC", Individual, sep = "_")) %>%
  mutate(DEAD = case_when(response == "died"~"Y",
                       response == "recovered"~"N",
                       response == "resistant"~"N",
                       response == "unknown"~NA)) %>%
  mutate(BLEACHING_OBS = ifelse(bleaching_status == "bleached" | bleaching_status == "recovering", "Y", "N")) %>%
  group_by(Individual) %>% mutate(nobs = n()) %>% filter(nobs > 1) %>%
  mutate(LATITUDE = case_when(
    reef_zone == "crest" ~ -17.4741,
    reef_zone == "forereef" ~ 17.4731
  ),
  LONGITUDE = case_when(
    reef_zone == "crest" ~ -149.8173,
    reef_zone == "forereef" ~ -149.8176
  )) 

write_xlsx(lbach.clean, path = "data/Lbach_clean.xlsx")

###Remotely ran Nunez et al. 2023 #######
### bring inputs here

load("data/get.taxa.nunez.Rdata")
load("data/Nunez_OTUtable.ITS2.Rdata")

otus = seqtab.nochim
tax.assign = get.taxa
#keep only samples with 5000 reads
otus = otus[rowSums(otus) >=5000,]
#keep only OTUs with at least 10 reads across samples
otus = otus[,colSums(otus)>=10]
dim(otus)

tax.assign = data.frame(tax.assign) %>%
  mutate(ASV = rownames(tax.assign))

#get the metadata
srr_2_samn = read.csv("data/PRJNA794042_meta.csv") %>% select(Run, BioSample, SampleName)
nunez.meta = read.table("data/NunezMeta.tsv", sep = "\t") # this is what we need

prep.otus = data.frame(otus) %>% mutate(ID = rownames(otus)) %>%
  pivot_longer(!ID, names_to = "ASV", values_to = "count") %>%
  left_join(tax.assign) %>% group_by(ID, Genus, Species) %>%
  filter(is.na(Genus) == F) %>%
  summarise(count = sum(count)) %>% ungroup() %>%
  group_by(ID) %>%
  mutate(relabun = count/sum(count)) %>% filter(relabun >= .1)

filt.otus = prep.otus %>% mutate(Species=gsub( "s__", "", Species),
                     Genus = gsub("g__", "", Genus)) %>%
  mutate(Genus = substr(Genus, 1, 1)) %>%
  group_by(ID) %>%
  mutate(SUBTYPE = paste(sort(unique(Species)), collapse = "_")) %>%
  mutate(GENUS = paste(sort(unique(Genus)), collapse = "_")) %>%
  select(!c(Species, Genus,relabun)) %>%
  mutate(SUBTYPE = ifelse(SUBTYPE=="", NA, SUBTYPE))

otu.meta = filt.otus %>% left_join(srr_2_samn, by = c("ID" = "Run")) %>% select(!count) %>%
  distinct() %>% left_join(nunez.meta, by = c("BioSample" = "V1")) %>%
  filter(V7 != "BLANK") %>%
  ungroup() %>%
  mutate(SPECIES = gsub(" ", "_", toupper(V7)),
         DATE = as.Date(V11, format = "%Y-%m-%d"),
         DEPTH_M = as.numeric(gsub(" m", "", V21)),
         LATLON = V27,
         BLEACHING_OBS = V37,
         COLONY_ID = paste(substr(SPECIES, 1, 3), gsub("_t[0-9]", "", SampleName), sep = "_"),
         LATITUDE = 21.461,
         LONGITUDE = -157.823
  ) %>%
  select(SPECIES, COLONY_ID, DATE, GENUS, SUBTYPE, DEPTH_M, BLEACHING_OBS, LATITUDE, LONGITUDE) %>%
  group_by(COLONY_ID) %>%
  mutate(nobs = n()) %>%
  filter(nobs > 1)

writexl::write_xlsx(otu.meta, path = "data/NunezClean2023.xlsx")


##### Process Turnham et al., 2025
turnin = read.csv("data/Palau_symbiont_diversity_data.csv") %>%
  separate(Coordinates, into = c("lat", "lon"), sep = "N,|N ") %>%
  mutate(LATITUDE = parse_lat(lat),
         LONGITUDE = parse_lon(lon)) %>%
  mutate(COLONY_ID=paste("TURN", Host.species, Sample.number))


### -------------- de Souza ---------------
#Assume that "early 2018" means like Jan 1? A lot are labeled as 2017.
#2019 = october of 2019
in.deSouza = readxl::read_excel("data/deSouza_Symportal.xlsx") %>% select(sample_name, 
                                                             host_genus, 
                                                             host_species,
                                                             B1:F3.1
                                                          ) %>%
  separate(sample_name, into = c("Date", "ColonyID"), sep = "_")

deSouza_filt = in.deSouza %>% pivot_longer(B1:F3.1, names_to = "zoox", values_to = "count") %>%
  mutate(short.sub = sub("^([^a-z.]*).*", "\\1", zoox)) %>%
  group_by(Date, ColonyID, host_genus, host_species) %>%
  mutate(relabun = count/sum(count)) %>%
  group_by(Date, ColonyID, host_genus, host_species, short.sub) %>%
  summarise(relabun = sum(relabun)) %>% filter(relabun > .1) %>%
  mutate(DATE = if_else(Date == "Clo2017", as.Date("01-01-2018", format = "%m-%d-%Y"),
                       as.Date("10-01-2019", format = "%m-%d-%Y"))) %>%
   mutate(GENUS = substr(short.sub, 1,1)) %>%
   group_by(Date, ColonyID) %>%
   mutate(SUBTYPE = paste(sort(unique(short.sub)), collapse = "_")) %>%
   mutate(GENUS = paste(sort(unique(GENUS)), collapse = "_")) %>%
   select(!c(short.sub, relabun)) %>% distinct()

deSouza_filt = deSouza_filt %>%
  mutate(COLONY_ID = paste("deSouz_MCAP_", ColonyID, sep = "")) %>%
  mutate(SPECIES = paste(toupper(host_genus), toupper(host_species), sep = "_")) %>%
  mutate(LATITUDE = 21.47098,
         LONGITUDE = -157.82333,
         DEPTH = 3) %>% ungroup() %>%
  group_by(COLONY_ID) %>% filter(n() > 1)

writexl::write_xlsx(deSouza_filt, path = "data/deSouzaClean.xlsx")

### ------------ POWELL 2026 ------------------
powell26 = read.csv(file = "data/Powell_OT_ITS2.dominanttype.CW_2020.csv")
powell.filt = powell26 %>% mutate(DATE = as.Date(collection_date, format = c("%m/%d/%y"))) %>%
  mutate(DATE = if_else(is.na(DATE)==T, as.Date(paste(collection_month, "01", collection_year), format = "%B %d %Y"), DATE)) %>%
  mutate(COLONY_ID = paste("POWELL_",
                           toupper(substr(host_genus, 1, 1)),
                           toupper(substr(host_species, 1,3)),
                           "_",
                           number, sep = ""),
         DEPTH = collection_depth,
         SPECIES = ifelse(host_genus == "Porites",
                               "PORITES_FURCATA",
                               paste(toupper(host_genus), toupper(host_species), sep = "_")
                               ),
         LONGITUDE = -1*collection_longitude) %>%
  pivot_longer(A4_sum:D1_sum, names_to = "zoox",
               values_to = "prop") %>%
  select(COLONY_ID, DEPTH, SPECIES, LONGITUDE, collection_latitude, DATE, zoox, prop) %>%
  filter(prop >= .1) %>%
  mutate(GENUS = substr(zoox, 1,1),
         short.sub = gsub("_sum", "", zoox)) %>%
  group_by(DATE, COLONY_ID) %>%
  mutate(SUBTYPE = paste(sort(unique(short.sub)), collapse = "_")) %>%
  mutate(GENUS = paste(sort(unique(GENUS)), collapse = "_")) %>%
  select(!c(short.sub, prop, zoox)) %>% distinct() %>%
  group_by(COLONY_ID) %>% mutate(nobs = n()) %>% filter(nobs > 1)
  
writexl::write_xlsx(powell.filt, path = "data/PowellClean.xlsx")

#### ---------------  GLASS et al. 2026 ------------------------
#https://onlinelibrary-wiley-com.libproxy.lib.unc.edu/doi/10.1111/gcb.70853

glass26 = read.csv(file = "data/Glass_2026/PostMED_relative_abundance.csv") %>%
  mutate(DATE = case_when(Year == 2022 ~ as.Date("05/01/2022", format = "%m/%d/%Y"),
                          Year == 2023 ~ as.Date("08/15/2023", format = "%m/%d/%Y"),
                          Year == 2024 ~ as.Date("04/28/2024", format = "%m/%d/%Y"))) %>%
  mutate(longitude = case_when(Site == "Cristobal Island" ~ -82.242280,
                               Site == "Hospital Point" ~ -82.217650,
                               Site == "Punta Donato" ~ -82.36804,
                               Site == "STRI Point" ~ -82.2627),
         latitude = case_when(Site == "Cristobal Island" ~ 9.268270,
                              Site == "Hospital Point" ~ 9.333190,
                              Site == "Punta Donato" ~ 9.35845,
                              Site == "STRI Point" ~ 9.34868),
         DEPTH = case_when(Site == "Cristobal Island" ~ 6.5,
                           Site == "Hospital Point" ~ 6,
                           Site == "Punta Donato" ~ 3,
                           Site == "STRI Point" ~ 4)) %>%
  mutate(COLONY_ID = paste("GLASS_SSID",
                           Coral_ID, sep = "_"),
         SPECIES = toupper("Siderastrea_siderea")
        ) %>%
  select(COLONY_ID, DEPTH, SPECIES, longitude, latitude, DATE, A3:D2u) %>%
  pivot_longer(A3:D2u, names_to = "zoox",
               values_to = "prop") %>%
  filter(prop >= .1) %>%
  mutate(GENUS = substr(zoox, 1,1),
         short.sub = sub("^([^a-z.]*).*", "\\1", zoox)) %>%
  group_by(DATE, COLONY_ID) %>%
  mutate(SUBTYPE = paste(sort(unique(short.sub)), collapse = "_")) %>%
  mutate(GENUS = paste(sort(unique(GENUS)), collapse = "_")) %>%
  select(!c(short.sub, prop, zoox)) %>% distinct() %>%
  group_by(COLONY_ID) %>% mutate(nobs = n()) %>% filter(nobs > 1)

writexl::write_xlsx(glass26, path = "data/GlassClean.xlsx")



######## CLEAN JUNG ET AL 2021 #####
injung = readxl::read_excel("data/Jung_2021_Symportal_post_med_absolut abundances_meta.xlsx") %>%
  pivot_longer(`11843_A`:`11852_C`, names_to = "zoox", values_to = "abund") %>%
  group_by(`Sample ID`) %>% mutate(relabun = abund/sum(abund)) %>%
  select(`Sample ID`, sample_name, zoox, relabun) %>% separate(`Sample ID`, into = c("Mon", "Tidal", "Bleached", "sample_name")) %>%
  filter(relabun >= .1) %>%
  mutate(GENUS = substr(zoox, 1,1)) %>%
  #mutate(SUBTYPE = ifelse(GENUS %in% LETTERS == T, 0, NA)) %>%
  mutate(GENUS = ifelse(GENUS %in% LETTERS == F, substr(zoox, nchar(zoox), nchar(zoox)), GENUS)) %>% 
  mutate(DATE = case_when(Mon == "Apr" ~ as.Date("04/01/2016", format = "%m/%d/%Y"), 
                          Mon == "Nov" ~ as.Date("11/01/2016", format = "%m/%d/%Y"))) %>%
  mutate(COLONY_ID = sample_name) %>% 
  mutate(zoox = ifelse(grepl("_", zoox) == T, NA, zoox)) %>% 
  group_by(DATE, COLONY_ID) %>%
  mutate(SUBTYPE = paste(sort(na.omit(unique(zoox))), collapse = "_"))  %>% select(!c(zoox, Mon, relabun)) %>% distinct() %>%
  group_by(COLONY_ID) %>% mutate(nobs = n()) %>% filter(nobs > 1)

outjung = injung %>% mutate(lat = -16.479, lon = 123.0435, 
                  BLEACHED = ifelse(Bleached == "NB", "N", "Y"), 
                  Notes = ifelse(Tidal == "IT", "Intertidal", "Subtidal"),
                  COLONY_ID = paste("Jung_AASP_", COLONY_ID, sep = ""),
                  SPECIES = "ACROPORA_ASPERA"
)

writexl::write_xlsx(outjung, path = "data/JungClean.xlsx")

######## CLEAN ALESSI ET AL 2024 #####
inalessi <- read.csv("data/Alessi_ITS2_results.csv") %>% 
  mutate(BLEACHED = ifelse(grepl("Z", Treat) == T, "Y", "N"),
         COLONY_ID = paste("ALESSI", "AMUR", ID, sep = "_")) %>%
  pivot_longer(C1:C21bh, names_to = "zoox", values_to = "relabun") %>% 
 mutate(DATE = ifelse(Time == "T1", "01/18/2021",
              "05/05/2021")) %>%
  filter(relabun >= .1) %>% 
  mutate(GENUS = substr(zoox, 1,1)) %>%
  mutate(GENUS = ifelse(GENUS %in% LETTERS == F, substr(zoox, nchar(zoox), nchar(zoox)), GENUS)) %>% 
  mutate(COLONY_ID = paste("ALESSI_AMUR_", ID, sep = "")) %>% 
  mutate(zoox = ifelse(grepl("_", zoox) == T, NA, zoox)) %>% 
  group_by(DATE, COLONY_ID) %>%
  mutate(SUBTYPE = paste(sort(na.omit(unique(gsub("[a-z]", "", zoox)))), collapse = "_"))  %>%
  select(!c(zoox, ID, Time, Treat, others, relabun)) %>% distinct() %>%
  group_by(COLONY_ID) %>% mutate(nobs = n()) %>% filter(nobs > 1) %>%
  mutate(LAT = ifelse(Site == "Bouraké", 165.990, 165.945),
         LON = ifelse(Site == "Bouraké", -21.952, -21.965)) %>% arrange(COLONY_ID)

writexl::write_xlsx(outjung, path = "data/JungClean.xlsx")

