# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 
#   /$$$$$$$ /$$       /$$$$$$ /$$$$$$$$/$$$$$$         /$$$$$$ 
#  | $$__  $| $$      /$$__  $|__  $$__/$$__  $$       /$$__  $$ 
#  | $$  \ $| $$     | $$  \ $$  | $$ | $$  \__/      |__/  \ $$ 
#  | $$$$$$$| $$     | $$  | $$  | $$ |  $$$$$$         /$$$$$$/
#  | $$____/| $$     | $$  | $$  | $$  \____  $$       /$$____|_
#  | $$     | $$     | $$  | $$  | $$  /$$  \ $$      | $$       
#  | $$     | $$$$$$$|  $$$$$$/  | $$ |  $$$$$$/      | $$$$$$$$ 
#  |__/     |________/\______/   |__/  \______/       |________/ 
#                                                                              
#
#            Beta-diversity at all cities level
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#############
# X. SETUP ###
###############
library(pacman)
p_load(tidyverse, RColorBrewer, phyloseq, patchwork, magrittr)

source('scripts/myFunctions.R')
source('scripts/0_config.R') # Variable naming and such

theme_set(theme_light())

#############################
# 1. PCoA ###
###########################

# From scripts/3_metrics.R
betadiv.ls <- read_rds('data/diversity/beta_diversity_full.ls.rds')


bdiv.plot.ls <- map(kingdoms, function(barcode) {
  
  PCo <- betadiv.ls$eig.df %>% 
    filter(Dist == 'bray' & Barcode == barcode) %>% 
    select(PCo1, PCo2) %>% as.vector
  
  betadiv.ls$plot.df %>% 
    filter(Dist == 'bray' & Barcode == barcode) %>%
    ggplot(aes(x = PCo1, y = PCo2, colour = city)) +
    geom_point(size = 2, colour = 'black', shape = 21, aes(fill = city)) +
    stat_ellipse(level = 0.95, geom = 'polygon', 
                 alpha = 0.2, aes(fill = city)) +
  #  scale_fill_manual(values = period_colours) +
    labs(colour = 'Sampling\nperiod', fill = 'Sampling\nperiod', 
         x = PCo$PCo1, y = PCo$PCo2) +
    theme(axis.title.y = element_text(margin = margin(r = -2, l=0)),
          axis.text = element_blank()) %>% 
    return() 
})


bdiv.plot.ls$BACT + bdiv.plot.ls$FUNG + bdiv.plot.ls$PLAN +
  plot_layout(guides = 'collect')

# EXPORT
ggsave(paste0('out/_MAIN/community_allCities.pdf'),
       bg = 'white', width = 2800, height = 1000,
       units = 'px', dpi = 200)


