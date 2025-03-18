library(pacman)
p_load(tidyverse, magrittr, purrr, patchwork, grid, phyloseq,
       ggdist, gghalves, tidyquant, ggthemes)
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/tax_glom2.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/rarefy_even_depth2.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/community_functions.R")

urbanbio.path <- '~/Desktop/ip34/urbanBio'
source(file.path(urbanbio.path, 'scripts/myFunctions.R'))

ps.ls <- read_rds(file.path(urbanbio.path, 'data/ps.ls.rds'))
ps_rare.ls <- read_rds(file.path(urbanbio.path, 'data/ps_rare.ls.rds'))


ps.ls$BACT@sam_data %>% 
  as_tibble %>% 
  group_by(city) %>% 
  ggplot(aes(x = date, groups = city, fill = city)) +
  geom_dotplot(binwidth = 6, method = 'histodot') 

# 
# for (bc in barcode_mapping) {
#   p <- div_Hill %>% 
#     filter(barcode == bc) %>% 
#     group_by(date, city, barcode, Index) %>% 
#     summarise(meanEffTax = mean(Effective_taxa),
#               sdEffTax = sd(Effective_taxa)) %>%
#     ggplot(aes(x = date, y = meanEffTax, colour = city)) +
#     geom_point()+ geom_line()+
#     facet_grid(Index ~ ., scales = 'free')
#   
#   ggsave(paste0('~/Desktop/ip34/urbanBio/out/alpha_', bc,'_Hill_lines.pdf'),
#          plot = p, bg = 'white', width = 2000, height = 2000, 
#          units = 'px', dpi = 220)
# }

# 
# div_Hill %>% 
#   filter(barcode == "BACT" & Index != 'Tail') %>% 
#   ggplot(aes(x = time, y = Effective_taxa, fill = time)) +
#   stat_halfeye(
#     adjust = 0.6,
#    # justification = -0.2,
#     normalize = 'panels',
#     .width = 0,
#     point_colour = NA
#   ) +
#   stat_dots(
#     # ploting on left side
#     side = "left",
#     # adjusting position
#    # justification = 1.2
#   ) +
#   geom_boxplot(
#     width = 0.2,
#     outlier.colour = NA,
#     alpha = 0
#   ) +
#   facet_grid(Index ~ city, scales = 'free') +  # Faceting
#   theme_light() +  # Optional: Use a minimal theme
#   labs(x = "City", y = "Effective Taxa", colour = "Season", fill = "Season")+  # Customize labels 
#   scale_fill_brewer(palette = 'Set2') 
#   

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





