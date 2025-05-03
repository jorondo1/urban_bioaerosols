# PLOT 2-5 : Community composition + alpha + beta div  (one per city)

#############
# X. SETUP ###
###############
library(pacman)
p_load(tidyverse, RColorBrewer, phyloseq, patchwork)
urbanbio.path <- '~/Desktop/ip34/urbanBio'

source(file.path(urbanbio.path, 'scripts/myFunctions.R'))
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/psflashmelt.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/myFunctions.R")

ps_rare.ls <- read_rds(file.path(urbanbio.path, 'data/ps_rare.ls.rds'))

barcodes <- c('BACT' = "Bacteria",
              'FUNG' = "Fungi",
              'PLAN' = "Pollen")

Hill_indices <- c('H_0' = 'Richness',
                  'H_1' = 'exp^(Shannon)',
                  'H_2' = 'Inverse Simpson')

cities <- c('Montreal' = 'Montreal', 
            'Quebec' = 'Quebec', 
            'Sherbrooke' = 'Sherbrooke')

# ==== COLOUR PALETTES
nTax_by_barcode <- c(
  'BACT' = 14,
  'FUNG' = 10,
  'PLAN' = 8
)

period_colours <- c('Spring' = 'springgreen3',
                    'Summer' = 'skyblue3', 
                    'Fall' = 'orange3')
nvdi_colours <- c("#E1BE6A", "#40B0A6", "#E66100")
# Define the bacteria palette with named colors
palettes <- list() 
palettes$BACT <- c(
  Others = "#C6C2C2",
  Unclassified = "#E6E5C1",
  Beijerinckiaceae = "#FFFFFF",
  `67-14` = "hotpink1",
  Streptomycetaceae = "darkmagenta",
  Acetobacteraceae = "brown2",
  Oxalobacteraceae = "sienna2",
  Sphingomonadaceae = "#FFA500",
  Paracoccaceae = "tomato",
  Pseudonocardiaceae = "#FFD700",
  Kineosporiaceae = "lightslateblue",
  Intrasporangiaceae = "mediumvioletred",
  Hymenobacteraceae = "#D700D7",
  Deinococcaceae = "#FF6B6B",
  Microbacteriaceae = "#B2006B", 
  Geodermatophilaceae = "#8A008A",
  Micrococcaceae = "violet",
  Nocardioidaceae = "#6A008A")

palettes$FUNG <- c(
  Others="#C6C2C2",
  Phanerochaetaceae="#FFFFFF",
  Mycosphaerellaceae="paleturquoise1",
  Rickenellaceae="royalblue4",
  Omphalotaceae="lightsteelblue",
  Cortinariaceae="seagreen2",
  Discinellaceae="darkgreen",
  Gloeophyllaceae="palegreen4",
  Atheliaceae="seagreen",
  Strophariaceae="darkslategrey",
  Cerrenaceae="cyan1",
  Incrustoporiaceae="#00CED1",
  Ganodermataceae="#40E0D0",
  Peniophoraceae="#20B2AA",
  Irpicaceae="#008B8B",
  Pleosporaceae="#5F9EA0",
  Fomitopsidaceae="#4682B4",
  Polyporaceae="#1E90FF",
  Cladosporiaceae="deepskyblue4")

palettes$PLAN <- c(
  Others="#C6C2C2",
  Unclassified = "#E6E5C1",
  Poaceae="#32CD32",
  Cupressaceae="#00FF00",
  Betulaceae="#ADFF2F",
  Equisetaceae="#74C476",
  Oleaceae="#238B45",
  Sapindaceae="#41AB5D",
  Asteraceae="#006D2C",
  Pinaceae="#00441B")

#period_colour = c('Sampling time_period 1' = 'springgreen4', 
#                  'Sampling time_period 2' = 'skyblue3', 
#                  'Sampling time_period 3' = 'orange3')


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

# Dataframe with taxrank mean relative abundance by (taxrank, city, date)
which_taxrank <- 'Family'

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
    imap(barcodes, function(barcode_name, barcode) {
      
      # First for all periods in the barcode, define top taxa
      # (redundant for each time iteration!)
      melted_filt <- melted %>% 
        filter(city == ci & Barcode == barcode) %>% 
        mutate(Barcode = recode_factor(Barcode, !!!barcodes),
               date = factor(date, levels = names(date_labels)))
        
      # How many taxa to display:
      nTaxa <- nTax_by_barcode[barcode] 
      top_taxa <- topTaxa(melted_filt, taxLvl = which_taxrank, topN = nTaxa)
      
      # Define top taxa levels
      top_taxa_lvls <- top_taxa %>%
        group_by(aggTaxo) %>% 
        aggregate(relAb ~ aggTaxo, data = ., FUN = sum) %>% 
        arrange(relAb) %$% aggTaxo %>% 
        as.character %>% # Set Others and Unclassified first:
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
              strip.text = element_text(size = 13))
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

#=== ADD DIVERSITIES TO COMMUNITY PLOT
alphadiv.df <- read_rds(file.path(urbanbio.path, 'data/diversity/alpha_diversity.rds'))
betadiv.ls <- read_rds(file.path(urbanbio.path, 'data/diversity/beta_diversity.ls.rds'))

map(cities, function(ci) {
  
  adiv.plot <- alphadiv.df %>%
    filter(city == ci) %>% 
    ggplot(aes(x = time, y = Shannon, fill = factor(veg_index_bracket.))) + 
    scale_fill_manual(values = nvdi_colours) +
    geom_boxplot(linewidth = 0.2) +
    facet_grid(.~Barcode) +
    labs(fill = 'Vegetation index') 
  
  bdiv.plot.ls <- map(barcodes, function(barcode) {
    
    PCo <- betadiv.ls$eig.df %>% 
      filter(City == ci & Dist == 'bray' & Barcode == barcode) %>% 
      select(PCo1, PCo2) %>% as.vector
    
    betadiv.ls$plot.df %>% 
      filter(City == ci & Dist == 'bray' & Barcode == barcode) %>%
      ggplot(aes(x = PCo1, y = PCo2, colour = time)) +
      geom_point(size = 1) +
      stat_ellipse(level = 0.95, geom = 'polygon', 
                   alpha = 0.2, aes(fill = time)) +
      theme_light() +
      scale_colour_manual(values = period_colours) +
      scale_fill_manual(values = period_colours) +
      labs(colour = 'Sampling\nperiod', fill = 'Sampling\nperiod', 
           x = PCo$PCo1, y = PCo$PCo2) +
      theme(axis.title.y = element_text(margin = margin(r = -2, l=0)),
            axis.text = element_blank())
  })
  
  bdiv.plot <- bdiv.plot.ls$BACT + bdiv.plot.ls$FUNG + bdiv.plot.ls$PLAN +
    plot_layout(guides = 'collect')
  
  plot <- 
    community_plots.ls[[ci]] / adiv.plot / bdiv.plot &
    theme(legend.position = "right",
          legend.justification = "left",
          legend.box.spacing = unit(0, "cm"),
    )
  
   ggsave(paste0('~/Desktop/ip34/urbanBio/out/community', ci,'.pdf'),
          plot = plot, bg = 'white', width = 2000, height = 2600,
          units = 'px', dpi = 160)
  
})

# X4 Community composition all samples 
# Same as 2., but order by date within facets, time time_period nested in barcode


# X2 Median income vs. vegetation index design (dotplot)
