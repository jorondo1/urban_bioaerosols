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
p_load(mgx.tools, tidyverse, RColorBrewer, phyloseq, patchwork, magrittr)

source('scripts/0_config.R') # Variable naming and such

theme_set(theme_light())

#############################
# 1. PCoA ###
###########################

# From scripts/3_metrics.R
betadiv.ls <- read_rds('data/diversity/beta_diversity_full.ls.rds')

betadiv.ls$plot.df$time<-factor(betadiv.ls$plot.df$time,levels=levels(betadiv.ls$plot.df$time)[c(2,3,1)])
bdiv.plot.ls <- map(kingdoms, function(barcode) {
  
  PCo <- betadiv.ls$eig.df %>% 
    filter(Dist == 'bray' & Barcode == barcode) %>% 
    select(PCo1, PCo2) %>% as.vector
  
  betadiv.ls$plot.df %>% 
    filter(Dist == 'bray' & Barcode == barcode) %>%
    ggplot(aes(x = PCo1, y = PCo2, colour = city)) +
    stat_ellipse(level = 0.95, geom = 'polygon', linewidth = 0,
                 alpha = 0.2, aes(fill = city,colour = city)) +
    geom_point(size = 3, colour = 'black', aes(fill = city, shape = time))+
    scale_color_manual(values = c("black", "black", "black")) +
    scale_fill_manual(values = c("#C60C30", "#1D4E89", "#2E8B57")) +
    scale_shape_manual(values=c(21,22,23))+
    
    labs(
      colour = 'City', 
      fill = 'City', 
      shape ='Sampling period',
      x = PCo$PCo1, y = PCo$PCo2) +
    #guides(fill="none",shape="none",colour="none")+
    theme(
      axis.title.y = element_text(margin = margin(r = -2, l=0)),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = 'none') +
    guides(colour = "none")
  
})

# extract legend
# Extract legend with overridden keys to remove black outline
p_with_legend <- bdiv.plot.ls$BACT +
  guides(
    fill = guide_legend(override.aes = list(colour = NA, shape = 22)),
    shape = guide_legend(override.aes = list(colour = "black"))
  ) +
  theme(
    legend.position = "right",
    legend.box = "horizontal",
    legend.text = element_text(size = 12),      # Increase legend text size
    legend.title = element_text(size = 15)      # Increase legend title size
    
    ) 

blank_with_legend <- ggplot() + 
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 17, "pt"), 
    axis.title.y = element_text(margin = margin(r = -2, l=0)) ,
  ) +
  annotation_custom(
    cowplot::get_legend(p_with_legend),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )

(bdiv.plot.ls$BACT + ggtitle("A. Bacteria") | 
    bdiv.plot.ls$FUNG + ggtitle("B. Fungi")) / 
  (bdiv.plot.ls$PLAN + ggtitle("C. Plant particles") | 
     blank_with_legend) +
  plot_layout(widths = c(1, 1), heights = c(1, 1))


# EXPORT
ggsave(paste0('out/_MAIN/community_allCities.pdf'),
       bg = 'white', width = 2000, height = 2000,
       units = 'px', dpi = 280)


