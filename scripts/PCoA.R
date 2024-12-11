library(pacman)
p_load(tidyverse, magrittr, purrr, patchwork)
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/community_functions.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/tax_glom2.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/rarefy_even_depth2.R")

urbanbio.path <- '~/Desktop/ip34/urbanBio'

ps.ls <- read_rds(file.path(urbanbio.path, 'data/ps.list.rds'))
cities <- ps.ls$BACT@sam_data$city %>% unique
seasons <- c('Spring' = 'springgreen3', 'Summer' = 'skyblue3', 'Fall' = 'violetred3')
barcodes <- c('BACT', 'FUNG', 'PLAN')

ps_rare.ls <- lapply(ps.ls, function(ps) {
  set.seed(1234)
  prune_samples(sample_sums(ps) >= 2000, ps) %>% 
    rarefy_even_depth2(ncores = 7) 
})

pcoa_bray.ls <- lapply(ps.ls, function(ps) {
  compute_pcoa(ps, dist = 'bray')
})

# Split dataset by city
ps_byCity.ls <- list()
ps_byCity.ls <- lapply(cities, function(city_name) {
  lapply(ps.ls, function(ps) {
    # Subset samples by city
    selected_samples <- rownames(sample_data(ps))[sample_data(ps)$city == city_name]
    
    # Subset the phyloseq object and remove zero-count taxa
    prune_samples(selected_samples, ps) %>%
      prune_taxa(taxa_sums(.) > 0, .)
  })
}); names(ps_byCity.ls) <- cities

# Compute pcoa, save in same structure list
pcoa_bray_byCity.ls <- imap(
  ps_byCity.ls, function(ps.ls, city) {
    imap(
      ps.ls, function(ps, barcode) {
        compute_pcoa(ps, 'bray')
  })
})

# FUNCTION PLOT PCOA
plot_pcoa <- function(pcoa.ls, ellipse) {
  # extract pcoa eigenvalues
  eig <- (100*pcoa.ls$eig[1:2]/sum(pcoa.ls$eig))  %>% round(1)
  
  #Plot 
  pcoa.ls$metadata %>% 
    mutate(time = factor(time, levels = names(seasons))) %>% 
    ggplot(aes(x = PCo1, y = PCo2, colour = !!sym(ellipse))) +
    geom_point(size = 2) +
    stat_ellipse(level = 0.95, geom = 'polygon', 
                 alpha = 0.2, aes(fill = !!sym(ellipse))) +
    theme_minimal() +
    scale_colour_manual(values = seasons) +
    scale_fill_manual(values = seasons) +
    labs(x = paste0('PCo1 (',eig[1],'% )'),
         y = paste0('PCo2 (',eig[2],'% )'))
}

# Example usage:
plot_pcoa(pcoa_bray.ls$BACT, "city")

# Plot pcoas, save plots in list
pcoa_bray_byCity.plot <- imap(
  pcoa_bray_byCity.ls, function(city.ls, city) {
    imap(
      city.ls, function(pcoa.ls, barcode) {
        plot_pcoa(pcoa.ls, 'time')
      } 
    )
  }
)


column_headers <- map(cities, ~ ggplot() + 
                        ggtitle(.x) + 
                        theme_void() + 
                        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))

# Add row headers (lvl2 names)
row_labels <- map(barcodes, ~ ggplot() + 
                    ggtitle(.x) + 
                    theme_void() + 
                    theme(plot.title = element_text(angle = 90, hjust = 0.5, size = 14, face = "bold")))


# Create plots for the main grid
# plot_rows <- map(pcoa_bray_byCity.plot, ~wrap_plots(.x)) # Combine plots for each lvl1
# 
# # Combine row labels, main grid, and column headers
# main_grid <- wrap_plots(
#   c(
#     list(wrap_plots(column_headers, ncol = length(column_headers))),  # Top row with column headers
#     map2(row_labels, plot_rows, ~wrap_plots(list(.x, .y), ncol = 1))     # Add row labels to each row
#   ),
#   ncol = 1
# ) +
#   plot_layout(heights = c(0.1, rep(1, length(row_labels)))) +  # Adjust height ratios
#   plot_annotation(title = "Grid of Plots")
# 
# # Display the grid
# main_grid
# 
# 
# 
# # Combine plots within each city
# patchwork_grid <- map(pcoa_bray_byCity.plot, ~wrap_plots(.x)) %>%  # Create grids for each lvl1
#   wrap_plots() +                                          # Combine all grids
#   plot_annotation(title = "Grid of Plots", tag_levels = "A")
# 



# Example for a single city
ps.ls$FUNG %>% 
  subset_samples(city == "Québec") %>% 
  prune_taxa(taxa_sums(.) >0, .) %>% 
  compute_pcoa(dist = 'bray') %$% metadata %>% 
  ggplot(aes(x = PCo1, y = PCo2, colour = time)) +
  geom_point(size = 2) +
  stat_ellipse(level = 0.95, geom = 'polygon', alpha = 0.2, aes(fill = time)) +
  theme_minimal()

# tax_glom to genus, or class level?

# Plot all samples coloured by city
pcoa_bray.ls$FUNG$metadata %>%
  ggplot(aes(x = PCo1, y = PCo2, colour = city)) +
  geom_point(size = 2) +
  stat_ellipse(level = 0.95, geom = 'polygon', alpha = 0.2, aes(fill = city)) +
  theme_minimal()
  # add eig% to axes
