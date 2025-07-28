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

betadiv.ls$plot.df$time<-c(betadiv.ls$BACT$bray$metadata$time,betadiv.ls$BACT$bray$metadata$time,
                           betadiv.ls$FUNG$bray$metadata$time,betadiv.ls$FUNG$bray$metadata$time,
                           betadiv.ls$PLAN$bray$metadata$time,betadiv.ls$PLAN$bray$metadata$time)
betadiv.ls$plot.df$time<-factor(betadiv.ls$plot.df$time,levels=levels(betadiv.ls$plot.df$time)[c(2,3,1)])
bdiv.plot.ls <- map(kingdoms, function(barcode) {
  
  PCo <- betadiv.ls$eig.df %>% 
    filter(Dist == 'bray' & Barcode == barcode) %>% 
    select(PCo1, PCo2) %>% as.vector
  
  betadiv.ls$plot.df %>% 
    filter(Dist == 'bray' & Barcode == barcode) %>%
    ggplot(aes(x = PCo1, y = PCo2, colour = city)) +
    geom_point(size = 0.5, colour = 'black', aes(fill = city)) +
    stat_ellipse(level = 0.95, geom = 'polygon', 
                 alpha = 0.2, aes(fill = city,colour = city)) +
    scale_color_manual(values = c("black", "black", "black")) +
    scale_fill_manual(values = c("#C60C30", "#1D4E89", "#2E8B57")) +
    labs(colour = 'City', fill = 'City', shape ='Sampling\nperiod',
         x = PCo$PCo1, y = PCo$PCo2) +
    theme(axis.title.y = element_text(margin = margin(r = -2, l=0)),
          axis.text = element_blank())+
    geom_point(size = 3, colour = 'black', aes(fill = city,shape = time))+
    scale_shape_manual(values=c(21,22,23))+
    guides(fill="none",shape="none",colour="none")%>% 
    return() 
})


bdiv.plot.ls$BACT + bdiv.plot.ls$FUNG + bdiv.plot.ls$PLAN +
  plot_annotation(tag_levels = "A")

# EXPORT
ggsave(paste0('out/_MAIN/community_allCities.pdf'),
       bg = 'white', width = 3000, height = 1000,
       units = 'px', dpi = 200)


