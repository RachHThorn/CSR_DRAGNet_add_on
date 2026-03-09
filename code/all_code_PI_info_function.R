# R Thornley
# 09/03/2026
# Function that prod

library(tidyverse)
library(viridis)

see <- read_csv("results/site_level_diversity_metrics.csv")
unique(see$site_name)

make_site_report <- function(
    focal_site,
    div_file   = "results/site_level_diversity_metrics.csv",
    cover_file = "results/cumulative_cover_number_species_sites.csv",
    taxa_file  = "results/taxon_lists_top_20_all_sites.csv",
    output_dir = "results/site_sampling_reports"
) {
  
  # make safe folder name
  safe_site_name <- focal_site %>%
    str_replace_all("[^[:alnum:] ]", "") %>%
    str_squish() %>%
    str_replace_all(" ", "_")
  
  site_dir <- file.path(output_dir, safe_site_name)
  dir.create(site_dir, recursive = TRUE, showWarnings = FALSE)
  
  #-----------------------------#
  # 1. Diversity scatter plot
  #-----------------------------#
  div <- read_csv(div_file, show_col_types = FALSE)
  
  p1 <- div %>%
    pivot_wider() %>%
    mutate(focal_site_flag = if_else(site_name == focal_site, "yes", "no")) %>%
    ggplot(aes(richness, pielou, colour = focal_site_flag)) +
    geom_point() +
    geom_text(aes(label = site_name), hjust = 0, nudge_x = 0.05) +
    theme_classic() +
    theme(legend.position = "none") +
    xlab("Species Richness") +
    ylab("Pielou Evenness Index") +
    scale_color_manual(values = c("yes" = "red", "no" = "black"))
  
  ggsave(
    filename = file.path(site_dir, paste0(safe_site_name, "_richness_by_evenness.jpeg")),
    plot = p1,
    height = 6,
    width = 7
  )
  
  #-----------------------------#
  # 2. Cumulative cover plot
  #-----------------------------#
  cover <- read_csv(cover_file, show_col_types = FALSE)
  
  p2 <- cover %>%
    filter(site_name == focal_site) %>%
    mutate(
      new_names = case_when(
        name == "top_05_sum" ~ "5 species sampled",
        name == "top_08_sum" ~ "8 species sampled",
        name == "top_10_sum" ~ "10 species sampled",
        name == "top_15_sum" ~ "15 species sampled",
        name == "top_20_sum" ~ "20 species sampled",
        name == "top_30_sum" ~ "30 species sampled",
        TRUE ~ NA_character_
      ),
      new_names = factor(
        new_names,
        levels = c(
          "5 species sampled", "8 species sampled", "10 species sampled",
          "15 species sampled", "20 species sampled", "30 species sampled"
        )
      )
    ) %>%
    ggplot(aes(year_trt, value, group = new_names, colour = new_names)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    facet_wrap(~trt) +
    xlim(0, 3) +
    ylim(30, 105) +
    ylab("Percentage of site level cover capture") +
    geom_hline(yintercept = 80, colour = "black", linetype = "dashed") +
    theme(
      legend.title = element_blank(),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(size = 10),
      legend.text = element_text(size = 12)
    ) +
    xlab("Year of experiment") +
    scale_y_continuous(breaks = seq(30, 100, by = 10)) +
    scale_color_viridis_d(option = "D")
  
  ggsave(
    filename = file.path(site_dir, paste0(safe_site_name, "_cover_capture.jpeg")),
    plot = p2,
    height = 8,
    width = 9
  )
  
  #-----------------------------#
  # 3. Taxon composition plot
  #-----------------------------#
  taxa <- read_csv(taxa_file, show_col_types = FALSE)
  
  p3 <- taxa %>%
    filter(site_name == focal_site) %>%
    mutate(New_taxon = fct_reorder(New_taxon, prop, .desc = TRUE)) %>%
    ggplot(aes(year_trt, prop, fill = New_taxon)) +
    geom_col() +
    ylim(0, 100) +
    theme_bw() +
    ggtitle(focal_site) +
    xlab("Year / treatment") +
    ylab("Proportion")
  
  ggsave(
    filename = file.path(site_dir, paste0(safe_site_name, "_top20_taxa.jpeg")),
    plot = p3,
    height = 8,
    width = 9
  )
  
  #-----------------------------#
  # 4. Priority table
  #-----------------------------#
  priority_table <- taxa %>%
    filter(site_name == focal_site) %>%
    group_by(site_name, New_taxon) %>%
    summarise(total_prop = sum(prop, na.rm = TRUE), .groups = "drop") %>%
    slice_max(order_by = total_prop, n = 20) %>%
    mutate(sample_order = rank(-total_prop, ties.method = "first")) %>%
    arrange(sample_order) %>%
    dplyr::select(site_name, New_taxon, sample_order) %>%
    rename(
      Site = site_name,
      Species = New_taxon,
      Species_sampling_priority = sample_order
    )
  
  write_csv(
    priority_table,
    file.path(site_dir, paste0(safe_site_name, "_species_priority_table.csv"))
  )
  
  # return useful outputs
  return(list(
    site_folder = site_dir,
    richness_evenness_plot = p1,
    cover_plot = p2,
    taxa_plot = p3,
    priority_table = priority_table
  ))
}


make_site_report(focal_site = "Wytham Woods")

make_site_report(focal_site = "Millville")
make_site_report(focal_site = "Jena")
make_site_report(focal_site = "CEREEP - Ecotron IDF")
make_site_report(focal_site = "Algaida")
make_site_report(focal_site = "Agroscope Changins")
