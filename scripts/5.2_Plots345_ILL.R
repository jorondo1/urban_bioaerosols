# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 
#   /$$$$$$$ /$$       /$$$$$$ /$$$$$$$$/$$$$$$         /$$$$$$       /$$$$$$$ 
#  | $$__  $| $$      /$$__  $|__  $$__/$$__  $$       /$$__  $$     | $$____/ 
#  | $$  \ $| $$     | $$  \ $$  | $$ | $$  \__/      |__/  \ $$     | $$      
#  | $$$$$$$| $$     | $$  | $$  | $$ |  $$$$$$         /$$$$$$/$$$$$| $$$$$$$ 
#  | $$____/| $$     | $$  | $$  | $$  \____  $$       /$$____|______|_____  $$
#  | $$     | $$     | $$  | $$  | $$  /$$  \ $$      | $$            /$$  \ $$
#  | $$     | $$$$$$$|  $$$$$$/  | $$ |  $$$$$$/      | $$$$$$$$     |  $$$$$$/
#  |__/     |________/\______/   |__/  \______/       |________/      \______/ 
#                                                                              
#
#            Community composition barcharts
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#############
# X. SETUP ###
###############
library(pacman)
p_load(tidyverse, RColorBrewer, phyloseq, patchwork, magrittr,ggrain)

source('scripts/myFunctions.R')
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/psflashmelt.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/community_functions.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/myFunctions.R")
source('scripts/0_config.R') # Variable naming and such

ps_rare.ls <- read_rds('data/ps_rare.ls.rds')
theme_set(theme_light())
# Dataframe with taxrank mean relative abundance by (taxrank, city, date)
which_taxrank <- 'Family'

# Number of taxa to display by barcode
nTax_by_barcode <- c(
  'BACT' = 14,
  'FUNG' = 14,
  'PLAN' = 8
)

##############################
# 1. Prepare plot dataframe ###
################################

# relAb by taxrank by sample with relevant metadata
melted <- imap(ps_rare.ls, function(ps, barcode) {
  
  psflashmelt(ps) %>%
    # Sum abundance by taxrank for each sample
    group_by(Sample, date, city, time, !!sym(which_taxrank)) %>% 
    summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
    
    # Then for all samples within a date/city/time combination
    group_by(date, city, time, !!sym(which_taxrank)) %>% 
    summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
    
    # Compute relab
    group_by(date, time, city) %>% 
    mutate(relAb = Abundance/sum(Abundance)) %>% 
    ungroup() %>% 
    
    # Select and format variables
    select(!!sym(which_taxrank), relAb, city, time, date) %>% 
    mutate(Barcode = barcode, 
           # To prevent gap in date x axis, we recode date. 
           date_label = paste(month(date, label = TRUE, abbr=TRUE), 
                              day(date)),
           # Recode date as character factor 
           date = as.character(date))
}) %>% list_rbind

# Count number of sample by unique (taxrank, city, date) for each barcode
sample_counts <- imap(ps_rare.ls, function(ps, barcode){
  samdat_as_tibble(ps) %>% 
    group_by(date, city, time) %>% 
    summarise(n_samples = n(), .groups = 'drop') %>% 
    mutate(Barcode = barcode,
           date_label = paste(month(date, label = TRUE, abbr=TRUE), 
                              day(date)),
           date = as.character(date))
}) %>% list_rbind

###############################
# 2. Community plots by city ###
#################################

# Compile community plots for each city 
community_plots.ls <- map(cities, function(ci) {
  
  # Each timeslot needs to be created separately, to better manage the x (time) axis
  time_plots.ls <- map(names(period_colours), function(time_period) {
    
    # Create custom (text) date labels for x axis with no date spacing
    date_labels <- melted %>% 
      filter(city == ci & time == time_period) %>% 
      select(date, date_label) %>% distinct %>% 
      deframe()
    # This is required for vertical x axis consistency !
    
    if (ci == 'Quebec' & time_period == 'Summer') {
      return(NULL)  # Skip this iteration
    }
    
    # Each city plot contains a subplot for each barcode
    imap(kingdoms, function(barcode_name, barcode) {
      
      # First for all periods in the barcode, define top taxa
      # (redundant for each time iteration!)
      melted_filt <- melted %>% 
        filter(city == ci & Barcode == barcode) %>% 
        mutate(Barcode = recode_factor(Barcode, !!!kingdoms),
               date = factor(date, levels = names(date_labels)))
      
      # How many taxa to display:
      nTaxa <- nTax_by_barcode[barcode] 
      top_taxa <- topTaxa(melted_filt, taxLvl = which_taxrank, topN = nTaxa)
      
      # Define top taxa levels
      top_taxa_lvls <- top_taxa %>%
        group_by(aggTaxo) %>% 
        aggregate(relAb ~ aggTaxo, data = ., FUN = sum) %>% 
        arrange(relAb) %$% aggTaxo %>% 
        as.character() %>% # Set Others and Unclassified first:
        setdiff(., c('Others', 'Unclassified')) %>% c('Others', 'Unclassified', .)
      
      # Plot data
      plot.df <- melted_filt %>% 
        left_join(top_taxa %>% select(-relAb), by = which_taxrank) %>%
        group_by(date, time, Barcode, aggTaxo) %>% 
        summarise(relAb = sum(relAb), .groups = 'drop') %>% 
        mutate(aggTaxo = factor(aggTaxo, levels = top_taxa_lvls)
        ) %>% # Filter at the end, to keep all taxa levels in the df
        filter(time == time_period)
      
      # Count number of samples by sample groups
      sample_counts_filt <- sample_counts %>% 
        filter(city == ci & Barcode == barcode & time == time_period) %>% 
        select(date, n_samples)
      
      # MAIN PLOT
      plot.df %>%
        ggplot(aes(x = date, y = relAb, fill = aggTaxo)) +
        geom_col() +
        # Sample count labels:
        geom_text(data = sample_counts_filt,
                  aes(x = date,
                      y = 1,  # Place at top of column
                      label = paste0("n = ", n_samples)), #small n
                  vjust = -1,  size = 2,
                  inherit.aes = FALSE  # Ignore fill aesthetic
        ) +
        coord_cartesian(clip = "off") +  # Disable clipping
        facet_grid(Barcode~.) +
        scale_fill_manual(values = palettes[[barcode]]) +
        scale_x_discrete(labels = date_labels, drop = FALSE) + # 
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # expand space between top of cols and panel border
        labs(fill = which_taxrank, y = '') +
        theme(panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid = element_blank(),
              strip.text = element_text(size = 14))
    })
  }); names(time_plots.ls) <- names(period_colours)
  
  # Remove empty lists (e.g. Québec summer)
  time_plots.ls %<>% compact()
  
  # === COMPOSITION PLOTS
  # BACTERIA (Top)
  p_bact <- list(
    time_plots.ls[['Spring']][['BACT']] + 
      theme(strip.text.y = element_blank()) +
      labs(tag = "A"), # Tags
    
    # Conditionally include Summer
    if (ci != 'Quebec') {
      time_plots.ls[['Summer']][['BACT']] +
        theme(strip.text.y = element_blank(),
              axis.text.y = element_blank())
    } else { NULL },
    
    time_plots.ls[['Fall']][['BACT']] +
      theme(axis.text.y = element_blank())
  ) %>% compact()
  
  p_bact <- purrr::reduce(p_bact, `|`) + 
    plot_layout(guides = "collect") & 
    theme(axis.text.x = element_blank(),
          axis.title = element_blank())
  
  # FUNGI (Middle)
  p_fung <- list(
    time_plots.ls[['Spring']][['FUNG']] + 
      theme(strip.text.y = element_blank()) +
      labs(y = 'Mean relative abundance of amplicons', tag = "B"), # Tags
    
    if (ci != 'Quebec') {
      time_plots.ls[['Summer']][['FUNG']] +
        theme(strip.text.y = element_blank(),
              axis.text.y = element_blank(),
              axis.title = element_blank()) 
    } else { NULL },
    
    time_plots.ls[['Fall']][['FUNG']] +
      theme(axis.text.y = element_blank(),
            axis.title = element_blank()) 
  ) %>% compact()
  
  p_fung <- purrr::reduce(p_fung, `|`) + 
    plot_layout(guides = "collect") & 
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank())
  
  # POLLEN (Bottom)
  p_plan <- list(
    time_plots.ls[['Spring']][['PLAN']] + 
      theme(strip.text.y = element_blank())+
      labs(tag = "C"), # Tags
    
    if (ci != 'Quebec') {
      time_plots.ls[['Summer']][['PLAN']] +
        theme(strip.text.y = element_blank(),
              axis.text.y = element_blank()
        ) 
    } else { NULL },
    
    time_plots.ls[['Fall']][['PLAN']] +
      theme(axis.text.y = element_blank())
  ) %>% compact
  
  p_plan <- purrr::reduce(p_plan, `|`) + 
    plot_layout(guides = "collect") & 
    theme(axis.text.x = element_text(angle = 90),
          axis.title = element_blank())
  
  # ==== SUPER PATCHWORK 
  p_bact / p_fung / p_plan &
    theme(
      legend.text = element_text(margin = margin(l = 5),
                                 size = 10),
      legend.key.size = unit(12,"pt")) 
});community_plots.ls

#############################
# 3. Add Diversity plots ###
###########################

# From scripts/3_metrics.R
alphadiv.df <- read_rds('data/diversity/alpha_diversity.rds')
betadiv.ls <- read_rds('data/diversity/beta_diversity_byCity.ls.rds')

map(cities, function(ci) {
  
  adiv.plot <- 
    alphadiv.df %>%
    filter(City == ci) %>% 
    ggplot(aes(x = time, y = Shannon, fill=time)) + 
    scale_fill_manual(values = period_colours) +
    geom_violin(width = 1, alpha = 0.5) + #Add violin & boxplot next line
    geom_boxplot(width = 0.1, outlier.shape = NA)+  # Removes points from boxplot
    geom_pwc(method = "dunn_test", label = "p.adj.signif")+ #Add stats
    facet_grid(.~Barcode) +
    labs(fill = 'Sampling\nperiod', tag = "D")+ # Tags
    ylim(0,8.5)
  
  bdiv.plot.ls <- map(kingdoms, function(barcode) {
    
    PCo <- betadiv.ls$eig.df %>% 
      filter(City == ci & Dist == 'bray' & Barcode == barcode) %>% 
      select(PCo1, PCo2) %>% as.vector
    
    betadiv.ls$plot.df %>% 
      filter(City == ci & Dist == 'bray' & Barcode == barcode) %>%
      ggplot(aes(x = PCo1, y = PCo2, colour = time)) +
      geom_point(size = 1) +
      stat_ellipse(level = 0.95, geom = 'polygon', 
                   alpha = 0.2, aes(fill = time)) +
      scale_colour_manual(values = period_colours) +
      scale_fill_manual(values = period_colours) +
      labs(colour = 'Sampling\nperiod', fill = 'Sampling\nperiod', 
           x = PCo$PCo1, y = PCo$PCo2) +
      theme(axis.title.y = element_text(margin = margin(r = -2, l=0)),
            axis.text = element_blank()) %>% 
      return() 
  })
  
  bdiv.plot <- 
    (bdiv.plot.ls$BACT+labs(tag="E")) + bdiv.plot.ls$FUNG + bdiv.plot.ls$PLAN + # Tags
    plot_layout(guides = 'collect')
  
  ### BUILD MEGAPLOT
  plot <- 
    community_plots.ls[[ci]] / adiv.plot / bdiv.plot &
    theme(legend.position = "right",
          legend.justification = "left",
          legend.box.spacing = unit(0, "cm"),
    )
  
  # EXPORT
  ggsave(paste0('out/_MAIN/community', ci,'.pdf'),
         plot = plot, bg = 'white', width = 2000, height = 2600,
         units = 'px', dpi = 160)
  
})

# Plot for Figure S8
# Bioaerosol alpha diversity across season
# All data from all cities together

#S8A All cities a-div ~ Sampling period
adiv.plot.all<-alphadiv.df %>%
  ggplot(aes(x = time, y = Shannon, fill=time)) + 
  scale_fill_manual(values = period_colours) +
  geom_violin(width = 1, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA)+  # Removes points from boxplot
  geom_pwc(method = "dunn_test", label = "p.adj.signif")+
  facet_grid(.~Barcode) +
  labs(fill = 'Sampling\nperiod',tag="A")+
  ylim(0,NA)+xlab("Sampling period")+ylab("Shannon Index")+
  theme(strip.text = element_text(color = "black",size = 14,face = "bold"),
        legend.position = c(0.04, 0.18),  # Adjust position inside plot area
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(face = "bold"));adiv.plot.all

#S8B All cities a-div ~ NDVI
adiv.plot.ndvi.all <- alphadiv.df %>%
  ggplot(aes(x = as.factor(veg_index_bracket.), y = Shannon, fill=as.factor(veg_index_bracket.))) + 
  scale_fill_manual(values = nvdi_colours) +
  geom_violin(width = 1, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA)+  # Removes points from boxplot
  #geom_pwc(method = "dunn_test", label = "p.adj.signif")+
  facet_grid(.~Barcode) +
  labs(fill = 'Vegetation\nindex',tag="B")+
  ylim(0,NA)+xlab("Vegetation index")+ylab("Shannon Index")+
  theme(strip.text = element_text(color = "black",size = 14,face = "bold"),
        legend.position = c(0.04, 0.18),  # Adjust position inside plot area
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(face = "bold"));adiv.plot.ndvi.all

#S8C All cities a-div ~ Median Income
adiv.plot.income.all <- alphadiv.df %>%
  ggplot(aes(x = as.factor(median_income_bracket.), y = Shannon, fill=as.factor(median_income_bracket.))) + 
  scale_fill_manual(values = c("#6A51A3", "#807DBA", "#9ECAE1", "#FFFFB2")) +
  geom_violin(width = 1, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA)+  # Removes points from boxplot
  #geom_pwc(method = "dunn_test", label = "p.adj.signif")+
  facet_grid(.~Barcode) +
  labs(fill = 'Median\nHousehold\nIncome',tag="C")+
  ylim(0,NA)+xlab("Median Household Income")+ylab("Shannon Index")+
  theme(strip.text = element_text(color = "black",size = 14,face = "bold"),
        legend.position = c(0.04, 0.23),  # Adjust position inside plot area
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(face = "bold")); adiv.plot.income.all
plot_adiv.all <- adiv.plot.all/ adiv.plot.ndvi.all / adiv.plot.income.all
ggsave(paste0('out/_MAIN/alphadiv_allcities.pdf'),plot = plot_adiv.all,
       bg = 'white', width = 3000, height = 3000,
       units = 'px', dpi = 200)

# Plot for Figure S9
# Bioaerosol alpha diversity across season
# All data from all cities together

#S9A All cities in Spring a-div ~ City
alphadiv.dfSpring<-alphadiv.df %>%
  filter(time=="Spring") %>%
  ggplot(aes(x = City, y = Shannon, fill=City)) + 
  scale_fill_manual(values = c("#C60C30", "#1D4E89", "#2E8B57")) +
  geom_violin(width = 1, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA)+
  geom_pwc(method = "dunn_test", label = "p.adj.signif")+
  facet_wrap(.~Barcode) +
  theme(strip.text = element_text(color = "black",size = 14,face = "bold"),
        axis.title.x = element_blank())+
  guides(fill="none")+ylab("Shannon Index")+
  labs(tag="A")+
  #theme(strip.text.x = element_blank())+
  ylim(0,9)

#S9B All cities in Summer a-div ~ City
alphadiv.dfSummer<-alphadiv.df %>%
  filter(time=="Summer" & City != "Quebec") %>%
  ggplot(aes(x = City, y = Shannon, fill=City)) + 
  scale_fill_manual(values = c("#C60C30", "#1D4E89", "#2E8B57")) +
  geom_violin(width = 1, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA)+
  geom_pwc(method = "dunn_test", label = "p.adj.signif")+
  facet_wrap(.~Barcode) +
  theme(strip.text = element_text(color = "black",size = 14,face = "bold"),
        axis.title.x = element_blank())+
  guides(fill="none")+ylab("Shannon Index")+
  labs(tag="B")+
  #theme(strip.text.x = element_blank())+
  ylim(0,9)

#S9C All cities in Fall a-div ~ City
alphadiv.dfFall<-alphadiv.df %>%
  filter(time=="Fall") %>%
  ggplot(aes(x = City, y = Shannon, fill=City)) + 
  scale_fill_manual(values = c("#C60C30", "#1D4E89", "#2E8B57")) +
  geom_violin(width = 1, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA)+
  geom_pwc(method = "dunn_test", label = "p.adj.signif")+
  facet_wrap(.~Barcode) +
  guides(fill="none")+ylab("Shannon Index")+
  labs(tag="C")+
  theme(strip.text = element_text(color = "black",size = 14,face = "bold"))+
  ylim(0,9)

# Combine S9ABC
plot_adiv <- 
  alphadiv.dfSpring/ alphadiv.dfSummer / alphadiv.dfFall &
  theme(legend.position = "right",
        legend.justification = "left",
        legend.box.spacing = unit(0, "cm"),
  )

# Export
ggsave(paste0('out/_MAIN/adiversityseasonallcities.pdf'),
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
    labs(fill = 'Vegetation\nindex',tag="A")+
    xlab("Vegetation Index")+ylab("Shannon Index")+
    theme(strip.text = element_text(color = "black",size = 14,face = "bold"),
          legend.position = c(0.045, 0.16),  # Adjust position inside plot area
          legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_text(face = "bold"))+
    ylim(0,NA)

  # Income 
  adiv.plot.income <- 
    alphadiv.df %>%
    filter(City == ci) %>% 
    ggplot(aes(x = as.factor(median_income_bracket.), y = Shannon, fill=as.factor(median_income_bracket.))) + 
    scale_fill_manual(values = c("#6A51A3", "#807DBA", "#9ECAE1", "#FFFFB2")) +
    geom_violin(width = 1, alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA)+
    facet_grid(.~Barcode) +
    xlab("Median Household Income")+ylab("Shannon Index")+
    theme(strip.text = element_text(color = "black",size = 14,face = "bold"),
          legend.position = c(0.045, 0.20),  # Adjust position inside plot area
          legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_text(face = "bold"))+
    labs(fill = 'Median\nHousehold\nIncome',tag="B")+
    ylim(0,NA)
  
  # Income & NDVI
  adiv.plot.income.by.veg <- 
    alphadiv.df %>%
    filter(City == ci) %>% 
    ggplot(aes(x = as.factor(veg_index_bracket.), y = Shannon, fill=as.factor(median_income_bracket.))) + 
    scale_fill_manual(values = c("#6A51A3", "#807DBA", "#9ECAE1", "#FFFFB2")) +
    geom_boxplot(width = 1, outlier.shape = NA)+
    facet_grid(.~Barcode) +
    xlab("Vegetation Index")+ylab("Shannon Index")+
    theme(strip.text = element_text(color = "black",size = 14,face = "bold"),
          legend.position = c(0.046, 0.20),  # Adjust position inside plot area
          legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_text(face = "bold"))+
    labs(fill = 'Median\nHousehold\nIncome',tag="C")+
    ylim(0,NA)
  
  # Pop density (not used)
  adiv.plot.popdensity <- 
    alphadiv.df %>%
    filter(City == ci) %>% 
    ggplot(aes(x = as.factor(population_density_index), y = Shannon, fill=as.factor(population_density_index))) + 
    geom_boxplot()+  # Removes points from boxplot
    geom_pwc(method = "dunn_test", label = "p.adj.signif")+
    facet_grid(.~Barcode) +
    theme(strip.text = element_text(color = "black",size = 14,face = "bold"),
          legend.position = c(0.04, 0.23),  # Adjust position inside plot area
          legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_text(face = "bold"))+
    labs(fill = 'Population\ndensity',tag="D")+
    xlab("Population density")+ylab("Shannon Index")+
    ylim(0,NA)

### BUILD MEGAPLOT (not using pop density)
plot <- adiv.plot.ndvi / adiv.plot.income / adiv.plot.income.by.veg

# EXPORT
ggsave(paste0('out/_MAIN/adiversity', ci,'.pdf'),
       plot = plot, bg = 'white', width = 2000, height = 2600,
       units = 'px', dpi = 160)
})

# Plots for Figures S1
# Site values for NDVI and Income gradients for each City
  plot.gradients <- 
    alphadiv.df %>%
    filter(Barcode == "Bacteria",time=="Spring") %>% 
    ggplot(aes(x = as.factor(veg_index_bracket.), y = as.factor(median_income_bracket.), color=City, fill=City)) + 
    geom_jitter(size=3,width = 0.1,height = 0.1,shape=21,alpha=0.9)+
    xlab("Vegedation density")+
    ylab("Median Household Income")+
    scale_color_manual(values = c("black", "black", "black")) +
    scale_fill_manual(values = c("#C60C30", "#1D4E89", "#2E8B57")) +
    facet_grid(.~City)+
    theme(strip.text = element_text(color = "black",size = 14,face = "bold"),
          legend.position = c(0.95, 0.16),  # Adjust position inside plot area
          legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_text(face = "bold")
    );plot.gradients
ggsave(paste0('out/_MAIN/gradients.pdf'),
         plot = plot.gradients, bg = 'white', width = 2000, height = 800,
         units = 'px', dpi = 160)

# Not used, smoothing graphs rather than boxplots
map(cities, function(ci) {
  adiv.plot.smooth <- 
    alphadiv.df %>%
    filter(City == ci) %>% 
    ggplot(aes(x = veg_index_bracket., y = Shannon, fill=as.factor(median_income_bracket.))) +
    geom_point()+
    geom_smooth()+
    facet_grid(.~Barcode)+
    labs(fill = "Income")+
    theme(strip.text = element_text(color = "black",size = 14,face = "bold"),
          legend.position = c(0.06, 0.16),  # Adjust position inside plot area
          legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_text(face = "bold"))
  adiv.plot.smooth2 <- 
    alphadiv.df %>%
    filter(City == ci) %>% 
    ggplot(aes(x = median_income_bracket., y = Shannon, fill=as.factor(veg_index_bracket.))) +
    geom_point()+
    geom_smooth()+
    facet_grid(.~Barcode)+
    labs(fill = "Vegetation")+
    theme(strip.text = element_text(color = "black",size = 14,face = "bold"),
          legend.position = c(0.06, 0.16),  # Adjust position inside plot area
          legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_text(face = "bold"))
})
