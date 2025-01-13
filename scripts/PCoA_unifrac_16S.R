library(pacman)
p_load(tidyverse, magrittr, purrr, patchwork, grid, unifrac)
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/community_functions.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/tax_glom2.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/rarefy_even_depth2.R")

urbanbio.path <- '~/Desktop/ip34/urbanBio'

ps <- read_rds(file.path(urbanbio.path, 'data/ps_mASV.rds')) # Masked ASV with tree
cities <- ps@sam_data$city %>% unique
seasons <- c('Spring' = 'springgreen3', 'Summer' = 'skyblue3', 'Fall' = 'orange3')
dist <- 'unifrac.w'

ps_rare <- prune_samples(sample_sums(ps) >= 2000, ps) %>% 
  rarefy_even_depth2(ncores = 7)

pcoa_uf.ls <- compute_pcoa(ps_rare, dist)

# Split every dataset by city
ps_byCity.ls <- lapply(cities, function(city_name) { # 1st level: city
  selected_samples <- rownames(sample_data(ps))[sample_data(ps)$city == city_name]
  prune_samples(selected_samples, ps) %>%
    prune_taxa(taxa_sums(.) > 0, .)
}); names(ps_byCity.ls) <- cities

# Compute pcoa for each city dataset
pcoa_byCity.ls <- imap(
  ps_byCity.ls, function(ps.ls, city) {
    compute_pcoa(ps, dist)
  })

# Mega dataframe with all results
# Imap iterates over lists and creates a variable with the list name
pcoa.df <- imap(pcoa_byCity.ls, function(pcoa, city) {
    pcoa$metadata %>% 
      rownames_to_column('Sample') %>% 
      select(Sample, city, time, PCo1, PCo2) 
  }) %>% list_rbind %>% 
  mutate(time = factor(time, levels = names(seasons))) # for plots


# Eigenvalues dataframe to annotate the plots
eig.df <- imap(pcoa_byCity.ls, function(pcoa, city) {
    pcoa$eig %>% 
      data.frame(Eig = .) %>% 
      rownames_to_column('MDS') %>% 
      tibble %>% 
      mutate(city = city) %>% 
      group_by(city) %>% 
      mutate(Eig = 100*Eig/sum(Eig)) %>% # Compute %eig
      filter(MDS %in% c('MDS1', 'MDS2')) # keep the 1st two
}) %>% list_rbind %>%
  # Data formatting for the plot
  group_by(city, MDS) %>%
  mutate(MDS = case_when(MDS == 'MDS1' ~ 'PCo1',
                         MDS == 'MDS2' ~ 'PCo2')) %>% 
  summarize(Eig = paste0(MDS, ": ", round(Eig, 1), "%"), .groups = "drop") %>%
  pivot_wider(names_from = MDS, values_from = Eig)

# PLOT 
pcoa.df %>% 
  ggplot(aes(x = PCo1, y = PCo2, colour = time)) +
  geom_point(size = 1) +
  stat_ellipse(level = 0.95, geom = 'polygon', 
               alpha = 0.2, aes(fill = time)) +
  theme_light() +
  facet_grid(~city, scales = 'free', space = 'free') +
  scale_colour_manual(values = seasons) +
  scale_fill_manual(values = seasons) +
  geom_text(data = eig.df, 
            aes(x = -1, y = -1, label = paste(PCo1, PCo2)),
            inherit.aes = FALSE, hjust = 0, vjust = 0) +
  labs(colour = 'Season', fill = 'Season')

ggsave('~/Desktop/ip34/urbanBio/out/pcoa_unifrac_rare.pdf',
       bg = 'white', width = 2200, height = 2000, 
       units = 'px', dpi = 220)
