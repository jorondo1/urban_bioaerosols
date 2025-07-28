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
p_load(tidyverse, RColorBrewer, phyloseq, patchwork, magrittr)

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
                      label = paste0("N = ", n_samples)),
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
      labs(title = ci),
    
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
      labs(y = 'Mean relative abundance of amplicons'),
    
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
      theme(strip.text.y = element_blank()),
    
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
})

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
    ggplot(aes(x = time, y = Hill_1, fill = factor(veg_index_bracket.))) + 
    scale_fill_manual(values = nvdi_colours) +
    geom_boxplot(linewidth = 0.2) +
    facet_grid(.~Barcode) +
    labs(fill = 'Vegetation index')+
    ylim(0,NA)
  
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
    bdiv.plot.ls$BACT + bdiv.plot.ls$FUNG + bdiv.plot.ls$PLAN +
    plot_layout(guides = 'collect')
  
  ### BUILD MEGAPLOT
  plot <- 
    community_plots.ls[[ci]] / adiv.plot / bdiv.plot &
    theme(legend.position = "right",
          legend.justification = "left",
          legend.box.spacing = unit(0, "cm"),
    )
  plot
  # EXPORT
  # ggsave(paste0('out/_MAIN/community', ci,'.pdf'),
  #        plot = plot, bg = 'white', width = 2000, height = 2600,
  #        units = 'px', dpi = 160)
  # 
})

