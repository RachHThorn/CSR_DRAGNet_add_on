# R Thornley
# 21/02/2025
# Find out how many species we need to sample per DRAGNet site for leaf trait collection

library(tidyverse)
library(vegan)

################################################################################
# Load DRAGNet master data and tidy
################################################################################

# dragnet cover data load
# may need to update the latest master data set
drag <- read.csv("data/full-cover-2025-09-11.csv")
unique(drag$year_trt)

# before filtering the dragnet data set has
drag %>% summarize(unique_values = n_distinct(site_name)) # 47 sites
drag %>% summarize(unique_values = n_distinct(Taxon)) # 1293 species
# BUT these numbers include subspecies etc..

# list some of these non taxon entries to filter out in the following pipe
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

################################################################################
# create richness and evenness measures for each of our sites
################################################################################

# create barchart of site level  richness
rich <- drag %>% group_by(site_name) %>% summarise(richness = n_distinct(New_taxon))
range(rich$richness)
ggplot(rich, aes(reorder(site_name, richness), richness)) + geom_col() + 
  theme(axis.text.x = element_text(angle = 90)) + xlab("site name")

even <- drag %>% 
  select(max_cover, New_taxon, site_name) %>%
  pivot_wider(names_from = New_taxon, values_from = max_cover, values_fn = ~ mean(.x, na.rm = TRUE)) %>%
  replace(is.na(.), 0) %>%
  arrange(site_name) %>%
  select(!site_name) %>%
  mutate(shannon = vegan::diversity(., "shannon", MARGIN = 1)) %>%
  mutate(simpson = vegan::diversity(., "simpson", MARGIN = 1)) %>%
  mutate(invsimp = vegan::diversity(., "invsimpson", MARGIN = 1)) %>%
  select(shannon, simpson, invsimp)

sites <- drag %>% arrange(site_name) %>% pull(site_name) %>% unique() 
even$site_name <- sites

div <- even %>% left_join(rich) 

# create Pielou evenness from shannon index and put all metrics into long format
div <- div %>% 
  mutate(pielou = shannon/log(richness)) %>%
  pivot_longer(cols = c(shannon, simpson, invsimp, pielou, richness))

# this data can now be filtered by the site to see values
div$site_name
div %>% filter(site_name == "Piedmont Prairie") # richness 24, pielou 0.730
div %>% filter(site_name == "Agroscope Changins") # richness 43, pielou 0.827
div %>% filter(site_name == "Wytham Woods") # richness 84, pielou 0.909

# save file for later plotting etc
write_csv(div, "results/site_level_diversity_metrics.csv")
