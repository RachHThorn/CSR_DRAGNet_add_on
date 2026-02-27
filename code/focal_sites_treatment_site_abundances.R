# R Thornley
# 25/02/2025
# For each focal site species 
# 1) Look at the effect on cumulative cover of picking a certain number of species most abundant

library(tidyverse)

################################################################################

# dragnet cover data load
drag <- read.csv("data/full-cover-2024-07-27.csv")
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

# For each year in each site
# get the percentage of the total plant cover of one species

# look at the total cover per site of each species
taxa_cover <- drag %>% group_by(site_name, year_trt, New_taxon) %>% nest() %>% 
  mutate(total_cover_plot = map(data, ~ .x %>% summarise(total = sum(max_cover)))) %>%
  select(site_name, total_cover_plot) %>%
  unnest(cols = c("site_name", "total_cover_plot"))

all_sites_5 <-
  taxa_cover %>% group_by(site_name, year_trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>% # divide total amount by the total veg over recorded
  select(site_name, year_trt, New_taxon, prop) %>%
  group_by(site_name, year_trt) %>% slice_max(order_by = prop, n = 5) %>% 
  summarise(top_05_sum = sum(prop)) 

all_sites_8 <-
  taxa_cover %>% group_by(site_name, year_trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>% # divide total amount by the total veg over recorded
  select(site_name, year_trt, New_taxon, prop) %>%
  group_by(site_name, year_trt) %>% slice_max(order_by = prop, n = 8) %>% 
  summarise(top_08_sum = sum(prop)) 

all_sites_10 <-
  taxa_cover %>% group_by(site_name, year_trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>% # divide total amount by the total veg over recorded
  select(site_name, year_trt, New_taxon, prop) %>%
  group_by(site_name, year_trt) %>% slice_max(order_by = prop, n = 10) %>% 
  summarise(top_10_sum = sum(prop)) 

all_sites_15 <-
  taxa_cover %>% group_by(site_name, year_trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>%
  select(site_name, year_trt, New_taxon, prop) %>%
  group_by(site_name, year_trt) %>% slice_max(order_by = prop, n = 15) %>% 
  summarise(top_15_sum = sum(prop)) 

all_sites_20 <-
  taxa_cover %>% group_by(site_name, year_trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>%
  select(site_name, year_trt, New_taxon, prop) %>%
  group_by(site_name, year_trt) %>% slice_max(order_by = prop, n = 20) %>% 
  summarise(top_20_sum = sum(prop)) 

all_sites_30 <-
  taxa_cover %>% group_by(site_name, year_trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>%
  select(site_name, year_trt, New_taxon, prop) %>%
  group_by(site_name, year_trt) %>% slice_max(order_by = prop, n = 30) %>% 
  summarise(top_30_sum = sum(prop)) 

dfs <- list(all_sites_5, all_sites_8, all_sites_10, all_sites_15, all_sites_20, all_sites_30)
dat <- reduce(dfs, left_join)
names(dat)
dat <- dat %>% pivot_longer(cols = top_05_sum:top_30_sum)

unique(dat$site_name)
site_wanted <- c("Wytham Woods", "Ainsdale Dune Slacks", "Hazelrigg", "Algaida", "Caracoles")
dat %>% filter(site_name %in% site_wanted) %>%
  ggplot(aes(year_trt, value, group = name, colour = name)) + 
  geom_point()+geom_line() + theme_bw()+
  facet_wrap(~site_name)+
  xlim(-1, 2)+
  ylim(40, 105)+
  ylab("Percentage of site level cover capture")+
  ggtitle("Site level cover capture by year of experiment and by number of species sampled")
# ggsave("Traits_subnetwork_sites_site_level_cover_capture_by_nos_species_smapled.jpeg", height = 6, width = 8)

################################################################################
# get the top ten taxa for each site and time
taxon_lists <-
  taxa_cover %>% 
  filter(site_name %in% site_wanted) %>%
  group_by(site_name, year_trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>% # divide total amount by the total veg over recorded
  select(site_name, year_trt, New_taxon, prop) %>%
  group_by(site_name, year_trt) %>% slice_max(order_by = prop, n = 10)

# ggplot(taxon_lists, aes(year_trt, prop, fill = New_taxon)) + geom_col() +facet_wrap(~site_name)

taxon_lists %>% filter(site_name == "Hazelrigg") %>% 
  mutate(New_taxon = fct_reorder(New_taxon, prop, .desc = TRUE)) %>%
  ggplot(aes(year_trt, prop, fill = New_taxon)) + geom_col() + ylim(0, 100) + theme_bw() +ggtitle("Hazelrigg")
# ggsave("Hazelrigg prop by species; n= 10.jpeg", height = 5, width = 6)

taxon_lists %>% filter(site_name == "Ainsdale Dune Slacks") %>% 
  mutate(New_taxon = fct_reorder(New_taxon, prop, .desc = TRUE)) %>%
  ggplot(aes(year_trt, prop, fill = New_taxon)) + geom_col() + ylim(0, 100) + theme_bw() +ggtitle("Ainsdale")
ggsave("Ainsdale prop by species; n= 10.jpeg", height = 5, width = 6)

taxon_lists %>% filter(site_name == "Wytham Woods") %>% 
  mutate(New_taxon = fct_reorder(New_taxon, prop, .desc = TRUE)) %>%
  ggplot(aes(year_trt, prop, fill = New_taxon)) + geom_col() + ylim(0, 100) + theme_bw()+ggtitle("Wytham")
ggsave("Wytham prop by species n = 10.jpeg", height = 5, width = 6)

taxon_lists %>% filter(site_name == "Caracoles") %>% 
  mutate(New_taxon = fct_reorder(New_taxon, prop, .desc = TRUE)) %>%
  ggplot(aes(year_trt, prop, fill = New_taxon)) + geom_col() + ylim(0, 100) + theme_bw()+ggtitle("Caracoles")
ggsave("Caracoles prop by species n = 10.jpeg", height = 5, width = 6)

taxon_lists %>% filter(site_name == "Algaida") %>% 
  mutate(New_taxon = fct_reorder(New_taxon, prop, .desc = TRUE)) %>%
  ggplot(aes(year_trt, prop, fill = New_taxon)) + geom_col() + ylim(0, 100) + theme_bw()+ggtitle("Algaida")
ggsave("Algaida prop by species n = 10.jpeg", height = 6, width = 8)


################################################################################

taxon_lists <-
  taxa_cover %>% 
  filter(site_name %in% site_wanted) %>%
  group_by(site_name, year_trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>% # divide total amount by the total veg over recorded
  select(site_name, year_trt, New_taxon, prop) %>%
  group_by(site_name, year_trt) %>% slice_max(order_by = prop, n = 20)

taxon_lists %>% filter(site_name == "Wytham Woods") %>% 
  mutate(New_taxon = fct_reorder(New_taxon, prop, .desc = TRUE)) %>%
  ggplot(aes(year_trt, prop, fill = New_taxon)) + geom_col() + ylim(0, 100) + theme_bw()+ggtitle("Wytham")
ggsave("Wytham prop by species n = 20.jpeg", height = 5, width = 10)

################################################################################

taxa_cover_trt <- drag %>% group_by(site_name, year_trt, trt, New_taxon) %>% nest() %>% 
  mutate(total_cover_plot = map(data, ~ .x %>% summarise(total = sum(max_cover)))) %>%
  select(site_name, total_cover_plot) %>%
  unnest(cols = c("site_name", "total_cover_plot"))

site_wanted <- c("Wytham Woods", "Ainsdale Dune Slacks", "Hazelrigg", "Algaida", "Caracoles")
all_sites_5 <-
  taxa_cover_trt %>%
  filter(site_name %in% site_wanted) %>%
  group_by(site_name, year_trt, trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>% # divide total amount by the total veg over recorded
  select(site_name, year_trt, trt, New_taxon, prop) %>%
  group_by(site_name, year_trt, trt) %>% slice_max(order_by = prop, n = 5) %>% 
  summarise(top_05_sum = sum(prop)) 

all_sites_10 <-
  taxa_cover_trt %>%
  filter(site_name %in% site_wanted) %>%
  group_by(site_name, year_trt, trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>% # divide total amount by the total veg over recorded
  select(site_name, year_trt, trt, New_taxon, prop) %>%
  group_by(site_name, year_trt, trt) %>% slice_max(order_by = prop, n = 10) %>% 
  summarise(top_10_sum = sum(prop)) 

all_sites_20 <-
  taxa_cover_trt %>%
  filter(site_name %in% site_wanted) %>%
  group_by(site_name, year_trt, trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>% # divide total amount by the total veg over recorded
  select(site_name, year_trt, trt, New_taxon, prop) %>%
  group_by(site_name, year_trt, trt) %>% slice_max(order_by = prop, n = 20) %>% 
  summarise(top_20_sum = sum(prop)) 


dfs <- list(all_sites_5, all_sites_10, all_sites_20)
dat <- reduce(dfs, left_join)
names(dat)
dat <- dat %>% pivot_longer(cols = top_05_sum:top_20_sum)

unique(dat$site_name)
site_wanted <- c("Wytham Woods", "Ainsdale Dune Slacks", "Hazelrigg", "Algaida", "Caracoles")
dat %>% filter(site_name %in% site_wanted) %>%
  ggplot(aes(year_trt, value, group = name, colour = name)) + 
  geom_point()+geom_line() + theme_bw()+
  facet_grid(trt~site_name)+
  xlim(-1, 2)+
  ylim(30, 105)+
  ylab("Percentage of site level cover capture")+
  ggtitle("Site level cover capture by treatment, year of experiment and by number of species sampled")

taxon_lists <-
  taxa_cover_trt %>% 
  filter(site_name %in% site_wanted) %>%
  group_by(site_name, year_trt, trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>% # divide total amount by the total veg over recorded
  select(site_name, year_trt, trt, New_taxon, prop) %>%
  group_by(site_name, year_trt, trt) %>% slice_max(order_by = prop, n = 20)

names(taxon_lists)
taxon_lists %>% filter(site_name == "Wytham Woods") %>% 
  mutate(New_taxon = fct_reorder(New_taxon, prop, .desc = TRUE)) %>%
  ggplot(aes(year_trt, prop, fill = New_taxon)) + geom_col() + facet_wrap(~trt) + 
  ylim(0, 100) + theme_bw()+ggtitle("Wytham")
# ggsave("Wytham prop by species and trt n = 20.jpeg", height = 5, width = 10)

taxon_lists <-
  taxa_cover_trt %>% 
  filter(site_name %in% site_wanted) %>%
  group_by(site_name, year_trt, trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>% # divide total amount by the total veg over recorded
  select(site_name, year_trt, trt, New_taxon, prop) %>% 
  group_by(site_name, year_trt, trt) %>% 
  mutate(sp_rank = rank(-prop))
range(taxon_lists$sp_rank)

taxon_lists %>% 
  filter(site_name == "Wytham Woods") %>% 
  filter(sp_rank == "1") %>%
  ggplot(aes(year_trt, prop, group = New_taxon, colour = New_taxon)) +
  geom_point()+
  geom_line()+
  theme_bw()+
  facet_wrap(~trt)

ranks_wanted <- c("1", "2", "3", "4", "5")
taxon_lists %>% 
  filter(site_name == "Ainsdale Dune Slacks") %>% 
  filter(sp_rank %in% ranks_wanted) %>%
  ggplot(aes(year_trt, prop, group = New_taxon, colour = New_taxon)) +
  geom_point()+
  geom_line()+
  theme_bw()+
  facet_wrap(~trt)

################################################################################
# so now do this for the 

# For each year in each site
# get the percentage of the total plant cover of one species

# look at the total cover per site of each species
taxa_cover <- drag %>% group_by(site_name, year_trt, New_taxon) %>% nest() %>% 
  mutate(total_cover_plot = map(data, ~ .x %>% summarise(total = sum(max_cover)))) %>%
  select(site_name, total_cover_plot) %>%
  unnest(cols = c("site_name", "total_cover_plot"))

all_sites_5 <-
  taxa_cover %>% group_by(site_name, year_trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>% # divide total amount by the total veg over recorded
  select(site_name, year_trt, New_taxon, prop) %>%
  group_by(site_name, year_trt) %>% slice_max(order_by = prop, n = 5) %>% 
  summarise(top_05_sum = sum(prop)) 

all_sites_8 <-
  taxa_cover %>% group_by(site_name, year_trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>% # divide total amount by the total veg over recorded
  select(site_name, year_trt, New_taxon, prop) %>%
  group_by(site_name, year_trt) %>% slice_max(order_by = prop, n = 8) %>% 
  summarise(top_08_sum = sum(prop)) 

all_sites_10 <-
  taxa_cover %>% group_by(site_name, year_trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>% # divide total amount by the total veg over recorded
  select(site_name, year_trt, New_taxon, prop) %>%
  group_by(site_name, year_trt) %>% slice_max(order_by = prop, n = 10) %>% 
  summarise(top_10_sum = sum(prop)) 

all_sites_15 <-
  taxa_cover %>% group_by(site_name, year_trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>%
  select(site_name, year_trt, New_taxon, prop) %>%
  group_by(site_name, year_trt) %>% slice_max(order_by = prop, n = 15) %>% 
  summarise(top_15_sum = sum(prop)) 

all_sites_20 <-
  taxa_cover %>% group_by(site_name, year_trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>%
  select(site_name, year_trt, New_taxon, prop) %>%
  group_by(site_name, year_trt) %>% slice_max(order_by = prop, n = 20) %>% 
  summarise(top_20_sum = sum(prop)) 

all_sites_30 <-
  taxa_cover %>% group_by(site_name, year_trt) %>% mutate(sum = sum(total)) %>% 
  mutate(prop = total/sum *100) %>%
  select(site_name, year_trt, New_taxon, prop) %>%
  group_by(site_name, year_trt) %>% slice_max(order_by = prop, n = 30) %>% 
  summarise(top_30_sum = sum(prop)) 

dfs <- list(all_sites_5, all_sites_8, all_sites_10, all_sites_15, all_sites_20, all_sites_30)
dat <- reduce(dfs, left_join)
names(dat)
dat <- dat %>% pivot_longer(cols = top_05_sum:top_30_sum)

################################################################################
# Some code that visualises the culmulative cover under different number of species
################################################################################


unique(dat$site_name)
site_wanted <- c("Wytham Woods", "Ainsdale Dune Slacks", "Hazelrigg", "Algaida", "Caracoles")
dat %>% filter(site_name %in% site_wanted) %>%
  ggplot(aes(year_trt, value, group = name, colour = name)) + 
  geom_point()+geom_line() + theme_bw()+
  facet_wrap(~site_name)+
  xlim(-1, 2)+
  ylim(40, 105)+
  ylab("Percentage of site level cover capture")+
  ggtitle("Site level cover capture by year of experiment and by number of species sampled")
# ggsave("Traits_subnetwork_sites_site_level_cover_capture_by_nos_species_smapled.jpeg", height = 6, width = 8)
