library(pacman)
p_load(tidyverse, magrittr, purrr, patchwork, grid, MetBrewer, parallel)
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/community_functions.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/tax_glom2.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/rarefy_even_depth2.R")
source("https://raw.githubusercontent.com/knights-lab/LMdist/refs/heads/main/lib/lmdist_source.r")
urbanbio.path <- '~/Desktop/ip34/urbanBio'

ps.ls <- read_rds(file.path(urbanbio.path, 'data/ps.ls.rds'))
ps_rare.ls <- read_rds(file.path(urbanbio.path, 'data/ps_rare.ls.rds'))

cities <- ps.ls$BACT@sam_data$city %>% unique
seasons <- c('Spring' = 'springgreen3', 'Summer' = 'skyblue3', 'Fall' = 'orange3')
barcodes <- c('BACT' = 'Bacteria', 'FUNG' = 'Fungi', 'PLAN' = 'Plants')
dist <- 'bray'

# Compute pcoas for each dataset separately
# pcoa_bray.ls <- mclapply(ps.ls, function(ps) {
#   compute_pcoa(ps, dist = dist)
# })

pcoa_bray_rare.ls <- mclapply(ps_rare.ls, function(ps) {
  compute_pcoa(ps, dist = dist)
})

# Split every dataset by city
ps_byCity.ls <- lapply(cities, function(city_name) { # 1st level: city
  lapply(ps_rare.ls, function(ps) { # 2nd level: barcode
    
    # Subset samples by city
    selected_samples <- rownames(sample_data(ps))[sample_data(ps)$city == city_name]
    
    # Subset the phyloseq object and remove zero-count taxa
    prune_samples(selected_samples, ps) %>%
      prune_taxa(taxa_sums(.) > 0, .)
  })
}); names(ps_byCity.ls) <- cities

# Compute pcoa for each city dataset
pcoa_bray_byCity.ls <- imap(
  ps_byCity.ls, function(ps.ls, city) {
    imap(
      ps.ls, function(ps, barcode) {
        compute_pcoa(ps, dist)
  })
})

# Mega dataframe with all results
# Imap iterates over lists and creates a variable with the list name
pcoa.df <- imap(pcoa_bray_byCity.ls, function(pcoa.ls, city) {
  imap(pcoa.ls, function(pcoa, barcode) {
    pcoa$metadata %>% 
      rownames_to_column('Sample') %>% 
      select(Sample, city, time, PCo1, PCo2) %>% 
      mutate(barcode = barcode) # add barcode name for each iteration
  }) %>% list_rbind
}) %>% list_rbind %>% 
  mutate(time = factor(time, levels = names(time)), 
         barcode = recode(barcode, !!!barcodes)) # for plots

# Eigenvalues dataframe to annotate the plots
eig.df <- imap(pcoa_bray_byCity.ls, function(pcoa.ls, city) {
  imap(pcoa.ls, function(pcoa, barcode) {
    pcoa$eig %>% 
      data.frame(Eig = .) %>% 
      rownames_to_column('MDS') %>% 
      tibble %>% 
      mutate(barcode = barcode,
             city = city) %>% 
      group_by(barcode, city) %>% 
      mutate(Eig = 100*Eig/sum(Eig)) %>% # Compute %eig
      filter(MDS %in% c('MDS1', 'MDS2')) # keep the 1st two
    }) %>% list_rbind
}) %>% list_rbind %>%
  # Data formatting for the plot
  group_by(barcode, city, MDS) %>%
  mutate(barcode = recode(barcode, !!!barcodes),
         MDS = case_when(MDS == 'MDS1' ~ 'PCo1',
                         MDS == 'MDS2' ~ 'PCo2')) %>% 
  summarize(Eig = paste0(MDS, ": ", round(Eig, 1), "%"), .groups = "drop") %>%
  pivot_wider(names_from = MDS, values_from = Eig)

# Plot !
pcoa.df %>% 
  ggplot(aes(x = PCo1, y = PCo2, colour = time)) +
  geom_point(size = 1) +
  stat_ellipse(level = 0.95, geom = 'polygon', 
               alpha = 0.2, aes(fill = time)) +
  theme_light() +
  facet_grid(barcode~city) +
  scale_colour_manual(values = time) +
  scale_fill_manual(values = time) +
  geom_text(data = eig.df, 
            aes(x = -1.8, y = -2.3, label = paste(PCo1, PCo2)),
            inherit.aes = FALSE, hjust = 0, vjust = 0) +
  labs(colour = 'Season', fill = 'Season')

ggsave('~/Desktop/ip34/urbanBio/out/pcoa_season_rare.pdf',
       bg = 'white', width = 2200, height = 2000, 
       units = 'px', dpi = 220)



# Testing models
test_perm <- adonis2(
  pcoa_bray.ls$BACT$dist.mx ~ median_income*vegetation_index_NDVI_landsat + city + mean_temperature + precip + mean_relative_humidity,
  data = pcoa_bray.ls$BACT$metadata,
  by = 'terms'); test_perm

pcoa_bray_rare.ls$BACT$metadata %>% 
  ggplot(aes(x = PCo1, y = PCo2, colour = temp_moy)) + 
  geom_point() 
# 
# test.pcoa <- pcoa_bray_byCity.ls$Montréal$FUNG
# 
# res <- adonis2(formula = test.pcoa$dist.mx ~ time + mean_relative_humidity +mean_temperature + median_income_bracket. + veg_index_bracket. + mean_wind_speed + concDNA, 
#                permutations = 1000,
#                data = test.pcoa$metadata,
#                by = 'margin',
#                na.action = na.exclude,
#                parallel = 8)
# res

###################################################################
############## SANDBOX **##########################################

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
    labs(x = paste0('PCo1 (',eig[1],'%)'),
         y = paste0('PCo2 (',eig[2],'%)'))
}

# Example usage:
plot_pcoa(pcoa_bray.ls$BACT, "time") +
  scale_colour_manual(values = seasons) +
  scale_fill_manual(values = seasons)

pcoa_bray.ls$PLAN$metadata %>% 
  #filter(seqDepth<20000) %>% 
  ggplot(aes(x = PCo1, y = PCo2, colour = seqDepth)) +
  geom_point(size = 2) +
  theme_minimal()

adonis2(pcoa_bray_rare.ls$PLAN$dist.mx~seqDepth,
        permutations = 10000,
        data = pcoa_bray_rare.ls$PLAN$metadata,
        by = 'terms',
        parallel = 8)


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

# Example for a single city
ps.ls$FUNG %>% 
  subset_samples(city == "Québec") %>% 
  prune_taxa(taxa_sums(.) >0, .) %>% 
  compute_pcoa(dist = dist) %$% metadata %>% 
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
