# R Thornley
# 17/03/2025
# For each focal site species 
# Look at the effect on cumulative cover of picking a certain number of species most abundant
# do this at the level of the treatment per year sampled

library(tidyverse)

################################################################################
# Load DRAGNet master data and tidy
################################################################################

# dragnet cover data load
drag <- read.csv("data/full-cover-2025-09-11.csv")
unique(drag$year_trt)

# before filtering the dragnet data set has this number of sites and species
drag %>% summarize(unique_values = n_distinct(site_name)) # 58 sites
drag %>% summarize(unique_values = n_distinct(Taxon)) # 1752 species
# BUT these numbers include subspecies / non plant inputs etc..

# list these non taxon entries to filter out in the following pipe
none_taxa <- c("Fungi", "Other_litter", "Other_standing_water", "Ground", "Bryophyte", 
               "Other_animal_diggings", "Other_woody_overstory", "Lichen",
               "Other_animal_digging", "Other_animal_droppings", "Other_rock")

# mutate the Taxon list in the Drag data so it is formatted well
# filter for the non taxonomic entries
# select the relevant cols
drag <- drag %>%
  mutate(New_taxon = str_to_sentence(Taxon)) %>%
  mutate(New_taxon = str_replace_all(New_taxon, " ", "_")) %>%
  filter(!str_detect(New_taxon, ".sp")) %>% # get rid of entries not to taxon level
  filter(!str_detect(New_taxon, "_x_")) %>% # get rid of any hybrids
  filter(!str_detect(New_taxon, "Unknown")) %>% # get rid of unknown species
  filter(!New_taxon %in% none_taxa) %>% # get rid of non taxa entries 
  mutate(New_taxon = case_when(New_taxon == "Helianthemum_nummularium_var._Grandiflorum" ~ "Helianthemum_nummularium",
                               New_taxon == "Mimosa_quadrivalvis_var._Platycarpa" ~ "Mimosa_quadrivalvis",
                               New_taxon == "Sebaea_sedoides_var._Schoenlandii" ~ "Sebaea_sedoides", 
                               TRUE ~ New_taxon)) %>%
  arrange(New_taxon) %>%
  select(site_name, year_trt, trt, block, plot, max_cover, New_taxon)

#################################################################################
# Get the cumulative cover we would obtain per site under different numbers of
# species sample to get an idea of required sampling effort
################################################################################

# For each year in each site
# get the percentage of the total plant cover of all species

# First get the total cover per treatment, per yearm per, site of each species present at the site
taxa_cover <- drag %>% group_by(site_name, year_trt, New_taxon, trt) %>% nest() %>% 
  mutate(total_cover_plot = map(data, ~ .x %>% summarise(total = sum(max_cover)))) %>%
  select(site_name:New_taxon, total_cover_plot) %>%
  unnest(cols = c("total_cover_plot"))

# now use this data to calculate if we were to sample the top 3 most abundant species 
# what percentage of the total plant cover we would achieve
all_sites_5 <-
  taxa_cover %>% group_by(site_name, year_trt, trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>% # divide total amount by the total veg over recorded
  select(site_name, year_trt, trt, New_taxon, prop) %>%
  group_by(site_name, year_trt, trt) %>% slice_max(order_by = prop, n = 5) %>% 
  summarise(top_05_sum = sum(prop)) 

all_sites_8 <-
  taxa_cover %>% group_by(site_name, year_trt, trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>% # divide total amount by the total veg over recorded
  select(site_name, year_trt, trt, New_taxon, prop) %>%
  group_by(site_name, year_trt, trt) %>% slice_max(order_by = prop, n = 8) %>% 
  summarise(top_08_sum = sum(prop)) 

all_sites_10 <-
  taxa_cover %>% group_by(site_name, year_trt, trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>% # divide total amount by the total veg over recorded
  select(site_name, year_trt, trt, New_taxon, prop) %>%
  group_by(site_name, year_trt, trt) %>% slice_max(order_by = prop, n = 10) %>% 
  summarise(top_10_sum = sum(prop)) 

all_sites_15 <-
  taxa_cover %>% group_by(site_name, year_trt, trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>%
  select(site_name, year_trt, New_taxon, trt, prop) %>%
  group_by(site_name, year_trt, trt) %>% slice_max(order_by = prop, n = 15) %>% 
  summarise(top_15_sum = sum(prop)) 

all_sites_20 <-
  taxa_cover %>% group_by(site_name, year_trt, trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>%
  select(site_name, year_trt, New_taxon, trt, prop) %>%
  group_by(site_name, year_trt, trt) %>% slice_max(order_by = prop, n = 20) %>% 
  summarise(top_20_sum = sum(prop)) 

all_sites_30 <-
  taxa_cover %>% group_by(site_name, year_trt, trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>%
  select(site_name, year_trt, New_taxon, trt, prop) %>%
  group_by(site_name, year_trt, trt) %>% slice_max(order_by = prop, n = 30) %>% 
  summarise(top_30_sum = sum(prop)) 

dfs <- list(all_sites_5, all_sites_8, all_sites_10, all_sites_15, all_sites_20, all_sites_30)
dat <- reduce(dfs, left_join)
names(dat)
dat <- dat %>% pivot_longer(cols = top_05_sum:top_30_sum)

write_csv(dat, "results/cumulative_cover_number_species_sites.csv")


