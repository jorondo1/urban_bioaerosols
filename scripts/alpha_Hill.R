library(pacman)
p_load(tidyverse, magrittr, purrr, patchwork, grid)
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/tax_glom2.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/rarefy_even_depth2.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/community_functions.R")

urbanbio.path <- '~/Desktop/ip34/urbanBio'
source(file.path(urbanbio.path, 'scripts/myFunctions.R'))


ps.ls <- read_rds(file.path(urbanbio.path, 'data/ps.ls.rds'))

# Rarefy phyloseq tables
ps_rare.ls <- lapply(ps.ls, function(ps) {
  set.seed(1234)
  prune_samples(sample_sums(ps) >= 2000, ps) %>% 
    rarefy_even_depth2(ncores = 7) 
})

###################
# Hill Numbers #####
###################
div_Hill <- imap(ps_rare.ls, function(ps, barcode) {
  div.fun(ps,  c(0,1,2)) %>% 
    data.frame %>% 
    rownames_to_column('Sample') %>% 
    left_join(samdat_as_tibble(ps), 
              by = 'Sample') %>% 
    mutate(barcode = barcode,
           time = factor(time, levels = c('Spring', 'Summer', 'Fall')))
}) %>% list_rbind %>% 
  pivot_longer(cols = c('H_0', 'H_1', 'H_2', 'Tail'), 
               names_to = 'Index',
               values_to = 'Effective_taxa') 

for (bc in barcode_mapping) {
  p <- div_Hill %>% 
    filter(barcode == bc) %>% 
    ggplot(aes(x = time, y = Effective_taxa, colour = time)) +
    geom_boxplot() +
    facet_grid(city~Index, scales = 'free_y') +
    theme_light() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = 'bottom') +
    scale_colour_brewer(palette = 'Dark2') +
    labs(y = 'Effective number of taxa')
  
  ggsave(paste0('~/Desktop/ip34/urbanBio/out/alpha_', bc,'_Hill.pdf'),
         plot = p, bg = 'white', width = 2000, height = 2000, 
         units = 'px', dpi = 220)
}
########################
# Classic div indices ###
########################

div <- imap(ps_rare.ls, function(ps, barcode) {
  div_estimate <- list()
  for (div in c('Richness', 'Shannon', 'Simpson')) {
    div_estimate[[div]] <- estimate_diversity(ps, index = div)  
  }
  div_estimate %>% 
    data.frame %>% 
    rownames_to_column('Sample') %>% 
    left_join(samdat_as_tibble(ps), 
              by = 'Sample') %>% 
    mutate(barcode = barcode,
           time = factor(time, levels = c('Spring', 'Summer', 'Fall')))
}) %>% list_rbind %>% 
  pivot_longer(cols = c('Richness', 'Simpson', 'Shannon'), 
               names_to = 'Index',
               values_to = 'Effective_taxa') 

div %>% 
  filter(barcode == "BACT") %>% 
  ggplot(aes(x = city, y = Effective_taxa, colour = city)) +
  geom_boxplot() +
  facet_grid(Index~time, scales = 'free_y') +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'bottom') +
  scale_colour_brewer(palette = 'Dark2') +
  labs(y = 'Effective number of taxa')





