# R Thornley
# 27/02/2026
# Create site species lists for top 20 species

library(tidyverse)

################################################################################
# FUNCTION
################################################################################

make_site_report <- function(
    focal_site,
    div_file   = "results/site_level_diversity_metrics.csv",
    cover_file = "results/cumulative_cover_number_species_sites.csv",
    taxa_file  = "results/taxon_lists_top_20_all_sites.csv",
    figures_dir = "figures"
) {
  # packages used
  require(dplyr)
  require(readr)
  require(ggplot2)
  require(forcats)
  require(viridis)
  
  # helper to make a safe folder name
  safe_name <- function(x) {
    x %>%
      gsub("[/\\\\:*?\"<>|]", "", .) %>%  # remove illegal path chars
      gsub("\\s+", "_", .) %>%            # spaces -> underscores
      trimws()
  }
  
  site_folder <- file.path(figures_dir, safe_name(focal_site))
  dir.create(site_folder, showWarnings = FALSE, recursive = TRUE)
  
  # -------------------------
  # PLOT 1: richness vs evenness
  # -------------------------
  div <- read_csv(div_file, show_col_types = FALSE)
  
  p1 <- div %>%
    pivot_wider() %>%   # <-- this was missing
    mutate(focal_site_flag = if_else(site_name == focal_site, "yes", "no")) %>%
    ggplot(aes(richness, pielou, colour = focal_site_flag)) +
    geom_point() +
    geom_text(aes(label = site_name), vjust = -0.4, size = 3) +
    theme_classic() +
    theme(legend.position = "none") +
    xlab("Species Richness") +
    ylab("Pielou Evenness Index") +
    scale_colour_manual(values = c("yes" = "red", "no" = "black"))
  
  ggsave(
    filename = file.path(site_folder, "all_site_richness_by_evenness.jpeg"),
    plot = p1, height = 6, width = 7, dpi = 300
  )
  
  # -------------------------
  # PLOT 2: cumulative cover
  # -------------------------
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
        TRUE ~ as.character(name)
      ),
      new_names = factor(
        new_names,
        levels = c("5 species sampled", "8 species sampled",
                   "10 species sampled", "15 species sampled",
                   "20 species sampled", "30 species sampled")
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
    scale_colour_viridis_d(option = "D")
  
  ggsave(
    filename = file.path(site_folder, "cumulative_cover_capture_by_species_count.jpeg"),
    plot = p2, height = 8, width = 9, dpi = 300
  )
  
  # -------------------------
  # PLOT 3: top-20 taxa stacked bars
  # -------------------------
  taxon_lists <- read_csv(taxa_file, show_col_types = FALSE)
  
  p3 <- taxon_lists %>%
    filter(site_name == focal_site) %>%
    mutate(New_taxon = fct_reorder(New_taxon, prop, .desc = TRUE)) %>%
    ggplot(aes(year_trt, prop, fill = New_taxon)) +
    geom_col() +
    ylim(0, 100) +
    theme_bw() +
    ggtitle(focal_site) +
    xlab("Year of experiment") +
    ylab("Proportion (%)")
  
  ggsave(
    filename = file.path(site_folder, "top_20_species_importance_by_treatment_year.jpeg"),
    plot = p3, height = 8, width = 9, dpi = 300
  )
  
  # -------------------------
  # TABLE: sampling priority top-20
  # -------------------------
  table_out <- taxon_lists %>%
    filter(site_name == focal_site) %>%
    group_by(site_name, New_taxon) %>%
    summarise(total_prop = sum(prop), .groups = "drop") %>%
    slice_max(order_by = total_prop, n = 20, with_ties = FALSE) %>%
    mutate(sample_order = rank(-total_prop)) %>%
    select(site_name, New_taxon, sample_order) %>%
    rename(
      Site = site_name,
      Species = New_taxon,
      Species_sampling_priority = sample_order
    ) %>%
    arrange(Species_sampling_priority)
  
  write_csv(table_out, file.path(site_folder, "top_20_species_sampling_priority.csv"))
  
  invisible(list(
    folder = site_folder,
    p1 = p1, p2 = p2, p3 = p3,
    table = table_out
  ))
}

################################################################################
# RUN per site
################################################################################

# example site Wytahm Woods
make_site_report(focal_site = "Wytham Woods", 
                 div_file   = "results/site_level_diversity_metrics.csv",
                 cover_file = "results/cumulative_cover_number_species_sites.csv",
                 taxa_file  = "results/taxon_lists_top_20_all_sites.csv",
                 figures_dir = "figures")

