# Diversity metrics for visualisation and statistics

#############
# X. SETUP ###
###############
library(pacman)
p_load(tidyverse, phyloseq)

# Diversity functions:
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/community_functions.R")
source('scripts/0_config.R') # Variable naming and such

ps_rare.ls <- read_rds('data/ps_rare.ls.rds')

#######################
#=== ALPHA DIVERSITY ###
#########################

# A bunch of indices: 
diversity.df <- imap(ps_rare.ls, function(ps, barcode) {
  map(c('Richness', 'Shannon', 'Simpson', 'Tail'), function(index){
    estimate_diversity(ps, index) %>% 
      as.data.frame() %>% 
      setNames(index)
    }) %>% list_cbind %>% 
    rownames_to_column('Sample') %>% 
    mutate(Hill_1 = exp(Shannon),
           Hill_2 = 1/Simpson) %>% 
    left_join(samdat_as_tibble(ps), # add samdat
              by = 'Sample') %>% 
    mutate(Barcode = barcode,
           Barcode = recode(Barcode, !!!kingdoms),
           time = factor(time, periods)) %>% 
    tibble() %>% 
    dplyr::rename(City = city) #plyr doesnt work
}) %>% list_rbind 

write_rds(diversity.df, 'data/diversity/alpha_diversity.rds')

######################
#=== BETA DIVERSITY ###
########################

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

# Iterate permanova for each city, barcode, distance
pcoa.ls <- imap(ps_byCity.ls, function(ps.ls, city) {
  imap(ps.ls, function(ps, barcode){
    out <- list()
    out[['robust.aitchison']] <- compute_pcoa(ps, dist = 'robust.aitchison')
    out[['bray']] <- compute_pcoa(ps, dist = 'bray', vst = TRUE)
    out
  })
})

# Extract first two eigenvalues and create strings for plots
eig.df <- imap(pcoa.ls, function(pcoa_barcode.ls, city) {
  imap(pcoa_barcode.ls, function(pcoa_dist.ls, barcode) {
    imap(pcoa_dist.ls, function(pcoa, dist) {
      pcoa$eig %>% 
        data.frame(Eig = .) %>% 
        rownames_to_column('MDS') %>% 
        tibble %>% 
        mutate(Barcode = barcode,
               City = city,
               Dist = dist) %>% 
        group_by(Barcode, City, Dist) %>% 
        mutate(Eig = 100*Eig/sum(Eig)) %>% # Compute %eig
        filter(MDS %in% c('MDS1', 'MDS2')) # keep the 1st two
    }) %>% list_rbind
  }) %>% list_rbind
}) %>% list_rbind %>% 
  group_by(Barcode, City, Dist, MDS) %>%
  mutate(Barcode = recode(Barcode, !!!kingdoms),
         MDS = case_when(MDS == 'MDS1' ~ 'PCo1',
                         MDS == 'MDS2' ~ 'PCo2')) %>% 
  summarize(Eig = paste0(MDS, ": ", round(Eig, 1), "%"), .groups = "drop") %>%
  pivot_wider(names_from = MDS, values_from = Eig) 

# Compile pcoa data points for plots
plot.df <- imap(pcoa.ls, function(pcoa_barcode.ls, city) {
  imap(pcoa_barcode.ls, function(pcoa_dist.ls, barcode) {
    imap(pcoa_dist.ls, function(pcoa, dist) {
      pcoa$metadata %>% 
        rownames_to_column('Sample') %>% 
        select(Sample, time, PCo1, PCo2) %>% 
        mutate(Barcode = barcode,
               Dist = dist,
               City = city) # add barcode name for each iteration
    }) %>% list_rbind
  }) %>% list_rbind
}) %>% list_rbind %>% tibble %>% 
  mutate(Barcode = recode(Barcode, !!!kingdoms),
         time = factor(time, levels = periods),
         # time = recode_factor(time, !!!c(
         #   'Spring' = '1',
         #   'Summer' = '2',
         #   'Fall' = '3')),
         # time = factor(time, levels = c('1', '2', '3'))
         )

pcoa.ls$plot.df <- plot.df
pcoa.ls$eig.df <- eig.df

write_rds(pcoa.ls,  file.path(urbanbio.path,'data/diversity/beta_diversity.ls.rds'))

# Sanity check plot : 
pcoa.ls$plot.df %>% 
  filter(Dist == 'bray') %>% 
  ggplot(aes(x = PCo1, y = PCo2, colour = time)) +
  geom_point() +
  stat_ellipse(level = 0.95, geom = 'polygon', 
               alpha = 0.2, aes(fill = time)) +
  facet_grid(City ~ Barcode)








