
library(pacman)
p_load(mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
       tidyverse, phyloseq, patchwork, ggh4x, ggpubr, #facet_nested
       magrittr, readxl, RColorBrewer)

source('scripts/0_config.R') # Variable naming and such

tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
ps.ls <- readRDS('data/ps.ls.rds')

# From scripts/3_metrics.R
alphadiv.df <- read_rds('data/diversity/alpha_diversity.rds')
theme_set(theme_light())

######################################################################
# S1 #################################################################
# Site values for NDVI and Income gradients for each City ############
######################################################################

plot.gradients <- 
  alphadiv.df %>%
  filter(Barcode == "Bacteria",time=="Spring") %>% 
  ggplot(aes(x = as.factor(veg_index_bracket.), y = as.factor(median_income_bracket.), color=City, fill=City)) + 
  geom_jitter(size=3,width = 0.1,height = 0.1,shape=21,alpha=0.9)+
  xlab("Vegetation density")+
  ylab("Median Household Income")+
  scale_color_manual(values = c("black", "black", "black")) +
  scale_fill_manual(values = city_colours) +
  facet_grid(.~City)+
  theme(strip.text = element_text(color = "black",size = 14,face = "bold"),
        legend.position = c(0.95, 0.16),  # Adjust position inside plot area
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(face = "bold")
  );plot.gradients

ggsave(paste0('out/_SUPP/S1_gradients.pdf'),
       plot = plot.gradients, bg = 'white', width = 2000, height = 800,
       units = 'px', dpi = 160)

# Not used, smoothing graphs rather than boxplots
# map(cities, function(ci) {
#   adiv.plot.smooth <- 
#     alphadiv.df %>%
#     filter(City == ci) %>% 
#     ggplot(aes(x = veg_index_bracket., y = Shannon, fill=as.factor(median_income_bracket.))) +
#     geom_point()+
#     geom_smooth()+
#     facet_grid(.~Barcode)+
#     labs(fill = "Income")+
#     theme(strip.text = element_text(color = "black",size = 14,face = "bold"),
#           legend.position = c(0.06, 0.16),  # Adjust position inside plot area
#           legend.background = element_rect(fill = "white", color = "black"),
#           legend.title = element_text(face = "bold"))
#   adiv.plot.smooth2 <- 
#     alphadiv.df %>%
#     filter(City == ci) %>% 
#     ggplot(aes(x = median_income_bracket., y = Shannon, fill=as.factor(veg_index_bracket.))) +
#     geom_point()+
#     geom_smooth()+
#     facet_grid(.~Barcode)+
#     labs(fill = "Vegetation")+
#     theme(strip.text = element_text(color = "black",size = 14,face = "bold"),
#           legend.position = c(0.06, 0.16),  # Adjust position inside plot area
#           legend.background = element_rect(fill = "white", color = "black"),
#           legend.title = element_text(face = "bold"))
# })


###################################
# Community composition overview ###
#####################################

which_taxrank <- 'Family'

# Created melted table by barcode
melted_glom.ls <- imap(ps.ls, function(ps, barcode){
  psflashmelt(ps) %>%
    # Sum abundance by taxrank for each sample
    filter(Abundance != 0) %>% 
    group_by(Sample, date, city, time, site_id, !!sym(which_taxrank)) %>% 
    summarise(Abundance = sum(Abundance), .groups = 'drop')  %>% 
    group_by(Sample, date, city, time, site_id) %>% 
    mutate(relAb = Abundance/sum(Abundance),
           site_date = fct_inorder(paste0(date,'_', site_id)),
           time = factor(time, levels=c('Spring', 'Summer', 'Fall'))) %>% 
    ungroup() %>% 
    arrange(date) %>% 
    mutate(site_date = fct_inorder(paste0(date,'_', site_id))) 
})

nTaxa <- 16
expanded_palette <- colorRampPalette(brewer.pal(12, 'Paired'))(nTaxa+2) 

# Build relab plots
comm_plot_data.ls <- imap(melted_glom.ls, function(melted.df, barcode){
  
  # Compute top taxa and create "Others" category
  (top_taxa <- topTaxa(melted.df, which_taxrank, nTaxa))
  
  top_taxa_lvls <- top_taxa %>% 
    group_by(aggTaxo) %>% 
    aggregate(relAb ~ aggTaxo, data = ., FUN = sum) %>% 
    arrange(relAb) %$% aggTaxo %>% 
    as.character() %>% # Others first:
    setdiff(., c('Others', 'Unclassified')) %>% c('Others', 'Unclassified', .)
  
  melted.df %>% 
    left_join(top_taxa %>% select(-relAb), by = which_taxrank) %>%
    filter(!is.na(time)) %>% 
    mutate(aggTaxo = factor(aggTaxo, levels = top_taxa_lvls)) %>%
    group_by(Sample, aggTaxo, city, time, site_date) %>% 
    summarise(relAb = sum(relAb), .groups = 'drop')
})

# Plot !
imap(comm_plot_data.ls, function(comm_plot_data.df, barcode){
  
  # PLOT + 
  comm_plot_data.df %>% 
    ggplot(aes(x = site_date, y = relAb, fill = aggTaxo)) +
    geom_col() +
    #theme_light() +
    facet_nested(cols=vars(city,time), scales = 'free', space = 'free',
                 strip = strip_nested(
                   clip = "off",  # Prevents text clipping
                   size = "variable",  # Adjusts text size
                   bleed = TRUE  # Merges strip spaces
                 )) +
    scale_fill_manual(values = expanded_palette) +
    labs(fill = which_taxrank, 
         x = 'Samples ordered by sampling date',
         y = paste('Relative abundance of', kingdoms[[barcode]], 'ASVs')) +
    theme(strip.text = element_text(color = "black",size = 13,face = "bold"),
          panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 14)
          ) +
    guides(fill = guide_legend(ncol = 1)) +
    # Entirely remove padding around cols : 
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0))
  
  # Dynamic figure naming!
  figNum <- 4 + which(names(barcodes) == barcode)
  
  ggsave(paste0('out/_SUPP/S',figNum,'_composition_',which_taxrank,'_',barcode,'.pdf'),
         bg = 'white', width = 3600, height = 2000, 
         units = 'px', dpi = 260)
})


###########################################################
# Figures S8-S12 codes by Isabelle
###########################################################

# Plot for Figure S8
# Bioaerosol alpha diversity across season
# All data from all cities together

#S8A All cities a-div ~ Sampling period
adiv.plot.all<-alphadiv.df %>%
  ggplot(aes(x = time, y = Shannon, fill=time)) + 
  scale_fill_manual(values = period_colours) +
  geom_violin(width = 1, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA)+  # Removes points from boxplot
  geom_pwc(method = "wilcox_test", label = "p.adj.signif")+
  facet_grid(.~Barcode) +
  labs(fill = 'Sampling\nperiod',tag="A")+
  xlab("Sampling period");adiv.plot.all

#S8B All cities a-div ~ NDVI
adiv.plot.ndvi.all <- alphadiv.df %>%
  ggplot(aes(x = as.factor(veg_index_bracket.), y = Shannon, fill=as.factor(veg_index_bracket.))) + 
  scale_fill_manual(values = nvdi_colours) +
  geom_violin(width = 1, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA)+  
  #geom_pwc(method = "wilcox_test", label = "p.adj.signif")+ #NS
  facet_grid(.~Barcode) +
  labs(fill = 'Vegetation\nindex',tag="B")+
  xlab("Vegetation index");adiv.plot.ndvi.all

#S8C All cities a-div ~ Median Income
adiv.plot.income.all <- alphadiv.df %>%
  ggplot(aes(x = as.factor(median_income_bracket.), y = Shannon, fill=as.factor(median_income_bracket.))) + 
  scale_fill_manual(values = c("#6A51A3", "#807DBA", "#9ECAE1", "#FFFFB2")) +
  geom_violin(width = 1, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA)+  
  #geom_pwc(method = "wilcox_test", label = "p.adj.signif")+
  facet_grid(.~Barcode) +
  labs(fill = 'Median\nHousehold\nIncome',tag="C")+
  xlab("Median Household Income"); adiv.plot.income.all

plot_adiv.all <- adiv.plot.all/ adiv.plot.ndvi.all / adiv.plot.income.all &
  theme(strip.text = element_text(color = "black",size = 14,face = "bold"),
        legend.position = 'none') &
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) &
  ylab("Shannon Index") ; plot_adiv.all

ggsave(paste0('out/_SUPP/S8_alphadiv_allcities.pdf'),plot = plot_adiv.all,
       bg = 'white', width = 1500, height = 1500,
       units = 'px', dpi = 160)

# Plot for Figure S9
# Bioaerosol alpha diversity across season
# All data from all cities together

#S9A All cities in Spring a-div ~ City
alphadiv.dfSpring<-alphadiv.df %>%
  filter(time=="Spring") %>%
  ggplot(aes(x = City, y = Shannon, fill=City)) + 
  geom_violin(width = 1, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA)+
  facet_wrap(.~Barcode) +
  theme(axis.title.x = element_blank())+
  labs(tag="A")

#S9B All cities in Summer a-div ~ City
alphadiv.dfSummer<-alphadiv.df %>%
  filter(time=="Summer" & City != "Quebec") %>%
  ggplot(aes(x = City, y = Shannon, fill=City)) + 
  geom_violin(width = 1, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA)+
  facet_wrap(.~Barcode) +
  theme(axis.title.x = element_blank())+
  labs(tag="B")

#S9C All cities in Fall a-div ~ City
alphadiv.dfFall<-alphadiv.df %>%
  filter(time=="Fall") %>%
  ggplot(aes(x = City, y = Shannon, fill=City)) + 
  geom_violin(width = 1, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA)+
  facet_wrap(.~Barcode) +
  labs(tag="C")

# Combine S9ABC
plot_adiv <- 
  alphadiv.dfSpring/ alphadiv.dfSummer / alphadiv.dfFall &
  scale_fill_manual(values = city_colours) &
  geom_pwc(method = "wilcox_test", label = "p.adj.signif") &
  theme(legend.position = "none",
        strip.text = element_text(color = "black",size = 14,face = "bold")
  ) &
  ylab("Shannon Index") & ylim(0,9); plot_adiv

# Export
ggsave(paste0('out/_SUPP/S9_adiversityseasonallcities.pdf'),
       plot = plot_adiv, bg = 'white', width = 1500, height = 1500,
       units = 'px', dpi = 160)

# Plots for Figures S10-11-12
# Bioaerosol alpha diversity across NDVI and INCOME
# For each city separated

map(cities, function(ci) {
  
  #NDVI
  adiv.plot.ndvi <- 
    alphadiv.df %>%
    filter(City == ci) %>% 
    ggplot(aes(x = as.factor(veg_index_bracket.), y = Shannon, fill=as.factor(veg_index_bracket.))) + 
    scale_fill_manual(values = nvdi_colours) +
    geom_violin(width = 1, alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA)+
    facet_grid(.~Barcode) +
    labs(fill = 'Vegetation\nindex',tag="A",
         x = "Vegetation Index") 
  
  # Income 
  adiv.plot.income <- 
    alphadiv.df %>%
    filter(City == ci) %>% 
    ggplot(aes(x = as.factor(median_income_bracket.), y = Shannon, fill=as.factor(median_income_bracket.))) + 
    scale_fill_manual(values = c("#6A51A3", "#807DBA", "#9ECAE1", "#FFFFB2")) +
    geom_violin(width = 1, alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA)+
    facet_grid(.~Barcode) +
    labs(fill = 'Median\nHousehold\nIncome',tag="B", x= "Median Household Income")
  
  # Income & NDVI
  adiv.plot.income.by.veg <- 
    alphadiv.df %>%
    filter(City == ci) %>% 
    ggplot(aes(x = as.factor(veg_index_bracket.), y = Shannon, fill=as.factor(median_income_bracket.))) + 
    scale_fill_manual(values = c("#6A51A3", "#807DBA", "#9ECAE1", "#FFFFB2")) +
    geom_boxplot(width = 1, outlier.shape = NA)+
    facet_grid(.~Barcode) +
    labs(fill = 'Median\nHousehold\nIncome',tag="C", x = "Vegetation Index")
  
  # Pop density (not used)
  adiv.plot.popdensity <- 
    alphadiv.df %>%
    filter(City == ci) %>% 
    ggplot(aes(x = as.factor(population_density_index), y = Shannon, fill=as.factor(population_density_index))) + 
    geom_boxplot()+  # Removes points from boxplot
    geom_pwc(method = "wilcox_test", label = "p.adj.signif")+
    facet_grid(.~Barcode) 
    labs(fill = 'Population\ndensity',tag="D", x = "Population density")
    
  ### BUILD MEGAPLOT (not using pop density)
  plot <- adiv.plot.ndvi / adiv.plot.income / adiv.plot.income.by.veg &
    ylab("Shannon Index") &
#    geom_pwc(method = "wilcox_test", label = "p.adj.signif") &
    theme(strip.text = element_text(color = "black",size = 14,face = "bold"),
          legend.position = 'none') &
    ylim(0,NA)
  
  # Dynamic figure number naming!
  figNum <- 9 + which(cities == ci)
  
  # EXPORT
  ggsave(paste0('out/_SUPP/S',figNum,'_adiversity', ci,'.pdf'),
         plot = plot, bg = 'white', width = 1500, height = 1950,
         units = 'px', dpi = 160)
})

