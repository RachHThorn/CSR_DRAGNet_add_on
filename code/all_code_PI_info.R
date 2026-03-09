# rough code that needs to be functionalised

# name the focal site
focal_site <- "Wytham Woods"

# read in diversity data
div <- read_csv("results/site_level_diversity_metrics.csv")
p1 <- div %>% 
  pivot_wider() %>%
  mutate(focal_site_flag = if_else(site_name == focal_site, "yes", "no")) %>% 
  ggplot(aes(richness, pielou, colour = focal_site_flag)) + geom_point() + 
  theme_classic() + geom_text(aes(label = site_name)) +
  theme(legend.position = "none")+
  xlab("Species Richness")+
  ylab("Pielou Evenness Index")+
  scale_color_manual(values = c(c("yes" = "red", 
                                  "no" = "black")))
p1
ggsave("figures/all_site_richness_by_evenness.jpeg", p1, height = 6, width = 7)

# read in cumulative cover data
cover <- read_csv("results/cumulative_cover_number_species_sites.csv")
cover$name
p2 <- cover %>% 
  filter(site_name == focal_site) %>%
  mutate(new_names = case_when(name == "top_05_sum" ~ "5 species sampled",
                               name == "top_08_sum" ~ "8 species sampled",
                               name == "top_10_sum" ~ "10 species sampled",
                               name == "top_15_sum" ~ "15 species sampled",
                               name == "top_20_sum" ~ "20 species sampled",
                               name == "top_30_sum" ~ "30 species sampled")) %>%
  mutate(new_names = factor(new_names, levels = c("5 species sampled", "8 species sampled",
                                                  "10 species sampled", "15 species sampled",
                                                  "20 species sampled", "30 species sampled"))) %>%
  ggplot(aes(year_trt, value, group = new_names, colour = new_names)) + 
  geom_point() + geom_line() + theme_bw()+
  facet_wrap(~trt)+
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
p2
ggsave("Leaf_traits_CSR_number_of_species_needed_3 sites.jpeg", p2, height = 8, width =9)

# read in taxon lists
taxa <- read_csv("results/taxon_lists_top_20_all_sites.csv")
# visualise for the different sites

p3 <- taxon_lists %>% filter(site_name == focal_site) %>% 
  mutate(New_taxon = fct_reorder(New_taxon, prop, .desc = TRUE)) %>%
  ggplot(aes(year_trt, prop, fill = New_taxon)) + geom_col() + ylim(0, 100) + theme_bw() + ggtitle(paste0(focal_site))
ggsave("MOst_important_species_top_20_treatment_year.jpeg", p3, height = 8, width =9)


# we can now filter this per site 
table <- 
  taxon_lists %>% 
  filter(site_name == focal_site) %>% 
  group_by(site_name, New_taxon) %>% 
  summarise(total_prop = sum(prop)) %>% 
  slice_max(order_by = total_prop, n = 20) %>% 
  mutate(sample_order = rank(-total_prop)) %>% 
  dplyr::select(site_name, New_taxon, sample_order) %>%
  rename(Site = site_name, Species = New_taxon, Species_sampling_priority = sample_order)
table


