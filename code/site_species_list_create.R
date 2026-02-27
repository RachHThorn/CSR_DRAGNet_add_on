# R Thornley
# 27/02/2026
# Create site species lists for top 20 species

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

################################################################################
# get top ten 
################################################################################

# For each year in each site
# get the percentage of the total plant cover of each species

# look at the total cover per site of each species
taxa_cover <- drag %>% group_by(site_name, year_trt, New_taxon) %>% nest() %>% 
  mutate(total_cover_plot = map(data, ~ .x %>% summarise(total = sum(max_cover)))) %>%
  select(site_name, total_cover_plot) %>%
  unnest(cols = c("site_name", "total_cover_plot"))

# get the top 20 taxa for each site and time - this is the absolute max. we would ever ask people to collect
taxon_lists <-
  taxa_cover %>% 
  group_by(site_name, year_trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>% # divide total amount by the total veg over recorded
  select(site_name, year_trt, New_taxon, prop) %>%
  group_by(site_name, year_trt) %>% slice_max(order_by = prop, n = 20)

write_csv(taxon_lists, "results/taxon_lists_top_20_all_sites.csv")

# visualise for the different sites

taxon_lists %>% filter(site_name == "Hazelrigg") %>% 
  mutate(New_taxon = fct_reorder(New_taxon, prop, .desc = TRUE)) %>%
  ggplot(aes(year_trt, prop, fill = New_taxon)) + geom_col() + ylim(0, 100) + theme_bw() +ggtitle("Hazelrigg")
# ggsave("Hazelrigg prop by species; n= 10.jpeg", height = 5, width = 6)

# we can now filter this per site 
wytham_lists <- taxon_lists %>% filter(site_name == "Wytham Woods")

# we can now filter this per site 
wytham_lists <- 
  taxon_lists %>% 
  filter(site_name == "Wytham Woods") %>% 
  group_by(site_name, New_taxon) %>% 
  summarise(total_prop = sum(prop)) %>% 
  slice_max(order_by = total_prop, n = 20) %>% 
  mutate(sample_order = rank(-total_prop))





################################################################################
# we can now filter this per site 
wytham_lists <- 
  taxon_lists %>% 
  filter(site_name == "Wytham Woods") %>% 
  group_by(site_name, New_taxon) %>% 
  summarise(total_prop = sum(prop)) %>% 
  slice_max(order_by = total_prop, n = 20) %>% 
  mutate(sample_order = rank(-total_prop))


