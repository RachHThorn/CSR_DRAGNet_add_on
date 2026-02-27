# R Thornley
# 27/02/206
# useful plots for justifying the approach / methods

library(tidyverse)
library(viridis) # for the plotting colour palette

################################################################################
# Site level DIVERSITY visualise
###############################################################################

# read in diversity data
div <- read_csv("results/site_level_diversity_metrics.csv")

names(div)
wanted <- c("invsimp", "richness")
div %>% filter(name == "richness") %>% 
  ggplot(aes(reorder(site_name, value), value, fill = site_name)) + geom_col() +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.position = "none")+
  xlab("DRAGNet site")+
  ylab("Species richness")
ggsave("figures/all_site_richness_barchart.jpeg", height = 6, width = 6)

div %>% filter(name == "invsimp") %>% 
  ggplot(aes(reorder(site_name, value), value, fill = site_name)) + geom_col() +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.position = "none")+
  xlab("DRAGNet site")+
  ylab("Inverse Simpson's Index")
ggsave("figures/all_site_invsimp_barchart.jpeg", height = 6, width = 6)

div %>% filter(name == "pielou") %>% 
  drop_na() %>%
  ggplot(aes(reorder(site_name, value), value, fill = site_name)) + geom_col() +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.position = "none")+
  xlab("DRAGNet site")+
  ylab("Pielou Evenness Index")
ggsave("figures/all_site_pielou_barchart.jpeg", height = 6, width = 6)

div %>% 
  pivot_wider() %>%
  ggplot(aes(richness, pielou, colour = site_name)) + geom_point() + 
  theme_classic() + geom_text(aes(label = site_name)) +
  theme(legend.position = "none")+
  xlab("Species Richness")+
  ylab("Pielou Evenness Index")
ggsave("figures/all_site_richness_by_evenness.jpeg", height = 6, width = 7)

#################################################################################
# Visualise the cumulative cover we would obtain for a few representative sites
################################################################################

# read in cumulative cover data
cover <- read_csv("results/cumulative_cover_number_species_sites.csv")

# select a list of sites we want to see data for
site_wanted <- c("Wytham Woods", "Piedmont Prairie", "Agroscope Changins")
# and the levels of samling effort (i.e. number of species)
levels_wanted <- c("top_05_sum", "top_10_sum", "top_08_sum", "top_15_sum")

# filter data / tidy names and plot out as a grid
cover %>% 
  filter(site_name %in% site_wanted) %>%
  filter(name %in% levels_wanted) %>%
  mutate(plot_label = case_when(site_name == "Piedmont Prairie" ~ "Piedmont: rich 24, even 0.7",
                                site_name == "Agroscope Changins" ~ "Agroscope: rich 43, even 0.8",
                                site_name == "Wytham Woods" ~ "Wytham: rich 84, even 0.9")) %>%
  mutate(plot_label = factor(plot_label, 
                             levels = c("Piedmont: rich 24, even 0.7", 
                                        "Agroscope: rich 43, even 0.8",
                                        "Wytham: rich 84, even 0.9"))) %>%
  mutate(new_names = case_when(name == "top_05_sum" ~ "5 species sampled",
                               name == "top_08_sum" ~ "8 species sampled",
                               name == "top_10_sum" ~ "10 species sampled",
                               name == "top_15_sum" ~ "15 species sampled")) %>%
  mutate(new_names = factor(new_names, levels = c("5 species sampled", "8 species sampled",
                                                  "10 species sampled", "15 species sampled"))) %>%
  ggplot(aes(year_trt, value, group = new_names, colour = new_names)) + 
  geom_point() + geom_line() + theme_bw()+
  facet_grid(trt ~ plot_label)+
  xlim(0, 3)+
  ylim(30, 105)+
  ylab("Percentage of site level cover capture")+
  geom_hline(yintercept = 80, colour = "black", linetype = "dashed")+
  theme(legend.title = element_blank())+
  xlab("Year of experiment")+
  theme(strip.background = element_rect(fill = "white", color = "black"))+
  scale_y_continuous(breaks = seq(30, 100, by = 10))+
  scale_color_viridis_d(option = "D")+
  theme(strip.text = element_text(size = 10))+
  theme(legend.text = element_text(size = 12))

# save plot
ggsave("Figures/Number_of_species_needed_representative_3_sites.jpeg", height = 8, width = 9)

################################################################################
#
################################################################################

# visualise cumulative quadrat cover captured for one site
names(cover)
Jena <- cover %>% filter(site_name == "Jena") 
Jena %>% 
  ggplot(aes(year_trt, value, group = New_taxon, colour = New_taxon)) + 
  geom_point() + geom_line()
Jena %>% group_by(year_trt) %>%  
  ggplot(aes(year_trt, total, group = New_taxon, colour = New_taxon)) + 
  geom_line() + geom_text(data = Jena[Jena$year_trt == 2, ],  
                          aes(label = New_taxon),
                          vjust = -1, hjust = 1.1, color = "grey10", size = 4) +
  theme_bw()+
  theme(legend.position = "none")

# ggsave("Jena_all_species.jpeg", height = 10, width = 3)
names(Jena)
Jena <- Jena %>% group_by(year_trt) %>% mutate(sum = sum(total)) %>% mutate(prop = total/sum *100)
Jena %>% 
  ggplot(aes(year_trt, prop, group = New_taxon, colour = New_taxon)) + 
  geom_line() + geom_text(data = Jena[Jena$year_trt == 2, ],  
                          aes(label = New_taxon),
                          vjust = -1, hjust = 1.1, color = "grey10", size = 4) +
  theme_bw()+
  theme(legend.position = "none")+
  ylab("Proportion of all species cover at site")
ggsave("Jena_all_species_prop.jpeg", height = 8, width = 3)

Jena %>% group_by(year_trt) %>% slice_max(order_by = prop, n = 10) %>% 
  ggplot(aes(year_trt, prop, group = New_taxon, colour = New_taxon)) + 
  geom_line() + geom_text(data = Jena[Jena$year_trt == 2, ],  
                          aes(label = New_taxon),
                          vjust = -1, hjust = 1.1, color = "grey10", size = 4) +
  theme_bw()+
  theme(legend.position = "none")+
  ylab("Proportion of all species cover at site")

dat_1 <- 
  Jena %>% group_by(year_trt) %>% slice_max(order_by = prop, n = 10) %>% 
  distinct(New_taxon, .keep_all = "TRUE") %>%
  group_by(year_trt) %>% summarise(top_10_sum = sum(prop)) 

dat_2 <- 
  Jena %>% group_by(year_trt) %>% slice_max(order_by = prop, n = 20) %>% 
  distinct(New_taxon, .keep_all = "TRUE") %>%
  group_by(year_trt) %>% summarise(top_20_sum = sum(prop))

dat_3 <- 
  Jena %>% group_by(year_trt) %>% slice_max(order_by = prop, n = 30) %>% 
  distinct(New_taxon, .keep_all = "TRUE") %>%
  group_by(year_trt) %>% summarise(top_30_sum = sum(prop))

dat_4 <- 
  Jena %>% group_by(year_trt) %>% slice_max(order_by = prop, n = 40) %>% 
  distinct(New_taxon, .keep_all = "TRUE") %>%
  group_by(year_trt) %>% summarise(top_40_sum = sum(prop))

dat_5 <- 
  Jena %>% group_by(year_trt) %>% slice_max(order_by = prop, n = 50) %>% 
  distinct(New_taxon, .keep_all = "TRUE") %>%
  group_by(year_trt) %>% summarise(top_50_sum = sum(prop))

dat_6 <- 
  Jena %>% group_by(year_trt) %>% slice_max(order_by = prop, n = 60) %>% 
  distinct(New_taxon, .keep_all = "TRUE") %>%
  group_by(year_trt) %>% summarise(top_60_sum = sum(prop))

dat_7 <- 
  Jena %>% group_by(year_trt) %>% slice_max(order_by = prop, n = 70) %>% 
  distinct(New_taxon, .keep_all = "TRUE") %>%
  group_by(year_trt) %>% summarise(top_70_sum = sum(prop))


dfs <- list(dat_1, dat_2, dat_3, dat_4, dat_5, dat_6, dat_7)
dat <- reduce(dfs, left_join)
names(dat)
dat <- dat %>% pivot_longer(cols = top_10_sum:top_70_sum)

dat %>%
  ggplot(aes(year_trt, value, group = name, colour = name)) + 
  geom_point()+geom_line() + theme_bw()+
  xlim(-1, 2)+
  ylim(0, 100)
# ggsave("Jena_sp_cumulative_cover.jpeg", height = 3, width = 3)

###############################################################################
# calculate and visualise plot level turnover between disturbance and control in year 1
##############################################################################

# 

treatments <- c("Control", "Disturbance")
nested_dfs <- drag %>% 
  filter(year_trt == 1) %>%
  filter(trt %in% treatments) %>%
  group_by(site_name, year_trt, trt, New_taxon) %>% 
  nest() %>% 
  mutate(total_cover_plot = map(data, ~ .x %>% summarise(total_cover_plot = sum(max_cover)))) %>%
  select(site_name, total_cover_plot) %>%
  unnest(cols = c("site_name", "total_cover_plot")) %>% 
  group_by(site_name, year_trt, trt) %>% 
  mutate(sum = sum(total_cover_plot)) %>% 
  mutate(prop = total_cover_plot/sum *100) %>%
  select(prop, New_taxon, site_name, trt) %>%
  group_by(site_name) %>%
  nest()

data <- nested_dfs %>% pull(data) %>% pluck(1)

get_dist <- function(data) {
  data <- 
    data %>% 
    select(!year_trt) %>%
    pivot_wider(names_from = New_taxon, values_from = prop) %>%
    replace(is.na(.), 0) %>%
    column_to_rownames(var = "trt")
  data <- vegan::vegdist(data)
  return(data)
}

dissim <- nested_dfs %>%
  mutate(dissim = map(data, get_dist)) %>%
  select(site_name, dissim) %>%
  mutate(dissim = map(dissim, ~ paste(.x, collapse = ", "))) %>%
  unnest(cols = c(dissim)) %>%
  filter(!dissim == "") %>%
  mutate(dissim = as.numeric(dissim))

range(dissim$dissim)
# values range from 0 to 0.98
# 0 identical and close to one - most different

dissim %>% 
  ggplot(aes(reorder(site_name, dissim), dissim, fill = site_name)) + 
  geom_col() +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.position = "none")+
  xlab("DRAGNet site")+
  ylab("Sorensens Dissimilarity between DIST and CONTROL plots")
ggsave("sorensens_barchart.jpeg", height = 7, width = 6)

################################################################################

# 17/03/2025
# do the analysis so we compare also between the treatments for the cumulative cover analysis
# look at the total cover per site of each species
# use this to normalise the cover values at the site, year, trt levels
taxa_cover <- drag %>% group_by(site_name, year_trt, trt ,New_taxon) %>% nest() %>% 
  mutate(total_cover_plot = map(data, ~ .x %>% summarise(total = sum(max_cover)))) %>%
  select(site_name:New_taxon, total_cover_plot) %>%
  unnest(cols = c("total_cover_plot"))

