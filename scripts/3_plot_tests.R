library(pacman)
p_load(tidyverse, RColorBrewer, phyloseq, patchwork)
urbanbio.path <- '~/Desktop/ip34/urbanBio'

source(file.path(urbanbio.path, 'scripts/myFunctions.R'))
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/community_functions.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/psflashmelt.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/myFunctions.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/tax_glom2.R")

#############
# X. SETUP ###
###############

# Agglomerate
ps_rare.ls <- read_rds(file.path(urbanbio.path, 'data/ps_rare.ls.rds'))

barcodes <- c('BACT' = "Bacteria",
              'FUNG' = "Fungi",
              'PLAN' = "Pollen")

Hill_indices <- c('H_0' = 'Richness',
                  'H_1' = 'exp^(Shannon)',
                  'H_2' = 'Inverse Simpson')

periods <- c('Spring', 'Summer', 'Fall')

#period_colour = c('Sampling time_period 1' = 'springgreen4', 
#                  'Sampling time_period 2' = 'skyblue3', 
#                  'Sampling time_period 3' = 'orange3')


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
#  /$$$$$$$  /$$        /$$$$$$  /$$$$$$$$        /$$$$$$ 
# | $$__  $$| $$       /$$__  $$|__  $$__/       /$$__  $$
# | $$  \ $$| $$      | $$  \ $$   | $$         |__/  \ $$
# | $$$$$$$/| $$      | $$  | $$   | $$           /$$$$$$/
# | $$____/ | $$      | $$  | $$   | $$          /$$____/ 
# | $$      | $$      | $$  | $$   | $$         | $$      
# | $$      | $$$$$$$$|  $$$$$$/   | $$         | $$$$$$$$
# |__/      |________/ \______/    |__/         |________/
# 
#            Community composition barcharts
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Dataframe with taxrank mean relative abundance by (taxrank, city, date)
which_taxrank <- 'Family'
melted <- imap(ps_rare.ls, function(ps, barcode) {
  
  psflashmelt(ps) %>%
    # filter(Abundance != 0) %>% # KEEP them !
    # Sum abundance by taxrank for each sample
    group_by(Sample, date, city, time, !!sym(which_taxrank)) %>% 
    summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
    # Then for all samples within a date/city/time combination
    group_by(date, city, time, !!sym(which_taxrank)) %>% 
    summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
    # Compute relab
    group_by(date, time, city) %>% 
    mutate(relAb = Abundance/sum(Abundance)) %>% 
    select(!!sym(which_taxrank), relAb, city, time, date) %>% 
    ungroup() %>% 
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

nTax_by_barcode <- c(
  'BACT' = 14,
  'FUNG' = 10,
  'PLAN' = 8
)

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

cities <- c('Montreal' = 'Montreal', 
            'Quebec' = 'Quebec', 
          'Sherbrooke' = 'Sherbrooke')

# Create all plots 
all_plots.ls <- map(cities, function(ci) {
  
  # Each timeslot needs to be created separately, to better manage the x (time) axis
  time_plots.ls <- map(periods, function(time_period) {
    
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
      
      # Is Unclassified part of the topN? 
      if(('Unclassified' %in% top_taxa$aggTaxo)) {
        fixed_colours <- c("#C6C2C2", "#E6E5C1")
      } else {
        fixed_colours <- c("#C6C2C2")
      }
      
      # Define palette
      expanded_palette <- colorRampPalette(brewer.pal(12, 'Set3'))(nTaxa) %>% 
        setdiff(., fixed_colours) %>% # colours reserved
        c(fixed_colours, .)
      
      # 
      plot.df <- melted_filt %>% 
        left_join(top_taxa %>% select(-relAb), by = which_taxrank) %>%
        group_by(date, time, Barcode, aggTaxo) %>% 
        summarise(relAb = sum(relAb), .groups = 'drop') %>% 
        mutate(aggTaxo = factor(aggTaxo, levels = top_taxa_lvls)
        ) %>% # Filter at the end, to keep all taxa levels in the df
        filter(time == time_period)
      
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
        facet_grid(Barcode~.) +
        scale_fill_manual(values = palettes[[barcode]]) +
        scale_x_discrete(labels = date_labels, drop = FALSE) +
        labs(fill = which_taxrank, y = 'Mean relative amplicon sequence abundance') +
        theme(panel.grid.minor = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid = element_blank(),
              legend.justification = "left",
              legend.box.spacing = unit(0, "cm"),
              strip.text = element_text(size = 13))
    })
  }); names(time_plots.ls) <- periods
  
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
  
  p_bact <- reduce(p_bact, `|`) + 
    plot_layout(guides = "collect") & 
    theme(axis.text.x = element_blank(),
          axis.title = element_blank())
  
  # FUNGI (Middle)
  p_fung <- list(
    time_plots.ls[['Spring']][['FUNG']] + 
      theme(strip.text.y = element_blank()),
    
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
  
  p_fung <- reduce(p_fung, `|`) + 
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
  
  p_plan <- reduce(p_plan, `|`) + 
    plot_layout(guides = "collect") & 
    theme(axis.text.x = element_text(angle = 90),
          axis.title = element_blank())
  
  # ==== SUPER PATCHWORK
  p_bact / p_fung / p_plan &
    theme(
      legend.position = "right",
      legend.justification = "left",
      legend.box.spacing = unit(0, "cm"),
      legend.text = element_text(margin = margin(l = 5)))
  
})

map(cities, function(ci) {
  map(barcodes, function(barcode) {
    
  })
})

imap(all_plots.ls, function(a_plot, ci) {
  ggsave(paste0('~/Desktop/ip34/urbanBio/out/community', ci,'.pdf'),
         plot = a_plot, bg = 'white', width = 2200, height = 2500,
         units = 'px', dpi = 200)
})

 melted.ls$PLAN %>% 
  left_join(top_taxa %>% select(-relAb), by = 'OTU') %>%
  mutate(aggTaxo = factor(aggTaxo, levels = top_taxa_lvls)) %>% 
  ggplot(aes(x = date, y = relAb, fill = aggTaxo)) +
  geom_col() +
  facet_grid(city~time, scales = 'free') +
  scale_fill_manual(values = expanded_palette) +
  labs(fill = 'OTU') +
   theme(panel.grid.minor = element_blank(),
        #axis.text.x = element_blank(),
       axis.ticks.x = element_blank())

# BFP grid / wrap with nested time_period, MQS rows, mean per time time_period
# Time time_period represented by ~median date ? distinct by city


### From here on, individual plots will be patched together
### per city, to contain both diversity analyses as well as
### some differential testing and maybe nestedness? 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
#  /$$$$$$$  /$$        /$$$$$$  /$$$$$$$$        /$$$$$$ 
# | $$__  $$| $$       /$$__  $$|__  $$__/       /$$__  $$
# | $$  \ $$| $$      | $$  \ $$   | $$         |__/  \ $$
# | $$$$$$$/| $$      | $$  | $$   | $$            /$$$$$/
# | $$____/ | $$      | $$  | $$   | $$           |___  $$
# | $$      | $$      | $$  | $$   | $$          /$$  \ $$
# | $$      | $$$$$$$$|  $$$$$$/   | $$         |  $$$$$$/
# |__/      |________/ \______/    |__/          \______/ 
#
#                  Alpha diversity 
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
 
# ps_rare_Genus.ls <- lapply(ps_rare.ls, function(ps){
#   tax_glom2(ps, taxrank = "Genus")
# })

# Compute diversity 
div_Hill <- imap(ps_rare.ls, function(ps, barcode) {
  div.fun(ps, c(0,1,2)) %>% 
    data.frame %>% 
    rownames_to_column('Sample') %>% 
    left_join(samdat_as_tibble(ps), # add samdat
              by = 'Sample') %>% 
    mutate(barcode = barcode) # Add barcode
}) %>% list_rbind %>% # long format:
  pivot_longer(cols = c('H_0', 'H_1', 'H_2', 'Tail'), 
               names_to = 'Index',
               values_to = 'Effective_taxa')

# Recode factors
div_Hill %<>%
  mutate_time %>% 
  mutate(barcode = recode_factor(barcode, !!!barcodes),
  Index = recode_factor(Index, !!!Hill_indices))

# One plot per city
for (ci in cities) {
  p <- list()
  # one subplot per barcode, because different scales prevent using facet_grid
  for (bc in barcodes) {
    p[[bc]] <- div_Hill %>% 
      filter(city == ci & Index != 'Tail' & barcode == bc) %>% 
      ggplot(aes(x = time, y = Effective_taxa, fill = time)) +
      geom_boxplot(outliers = TRUE) +
      facet_grid(Index~barcode, scales = 'free') +
      theme_light() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            legend.position = 'bottom') +
     # scale_fill_manual(values = period_colour) +
      labs(y = 'Effective number of taxa') +
      ylim(0,NA) 
  }
  
  # Some patchwork gymnastics to stitch them back together
  # and make it look like a facet_grip with column-flexible y scale
  p[[barcodes[1]]] <- p[[barcodes[1]]] +
    theme(strip.text.y = element_blank())
  
  p[[barcodes[2]]] <- p[[barcodes[2]]] +
    theme(strip.text.y = element_blank(),
          axis.title.y = element_blank())
  
  p[[barcodes[3]]] <- p[[barcodes[3]]] +
    theme(axis.title.y = element_blank())
  
  p2 <- p[[barcodes[1]]] + p[[barcodes[2]]] + p[[barcodes[3]]] +
    plot_layout(guides = "collect") &
    theme(legend.position = 'bottom')
  
  ggsave(paste0('~/Desktop/ip34/urbanBio/out/alpha_', ci,'_Hill.pdf'),
         plot = p2, bg = 'white', width = 2000, height = 2000,
         units = 'px', dpi = 220)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
#  /$$$$$$$  /$$        /$$$$$$  /$$$$$$$$       /$$   /$$
# | $$__  $$| $$       /$$__  $$|__  $$__/      | $$  | $$
# | $$  \ $$| $$      | $$  \ $$   | $$         | $$  | $$
# | $$$$$$$/| $$      | $$  | $$   | $$         | $$$$$$$$
# | $$____/ | $$      | $$  | $$   | $$         |_____  $$
# | $$      | $$      | $$  | $$   | $$               | $$
# | $$      | $$$$$$$$|  $$$$$$/   | $$               | $$
# |__/      |________/ \______/    |__/               |__/
#
#                   Beta diversity
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Betadiv en BFP
# 2 rows : BC + Aitchison? Ou UniFrac pour les bact?
# 3rd row : Turnover/nestedness
# Include variance partitioning?

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
    out[['bray']] <- compute_pcoa(ps, dist = 'bray')
    out
  })
})

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
  mutate(Barcode = recode(Barcode, !!!barcodes),
         MDS = case_when(MDS == 'MDS1' ~ 'PCo1',
                         MDS == 'MDS2' ~ 'PCo2')) %>% 
  summarize(Eig = paste0(MDS, ": ", round(Eig, 1), "%"), .groups = "drop") %>%
  pivot_wider(names_from = MDS, values_from = Eig) %>% unnest(PCo1, PCo2)

# Compile pcoa data 
pcoa.df <- imap(pcoa.ls, function(pcoa_barcode.ls, city) {
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
}) %>% list_rbind %>% 
  mutate(Barcode = recode(Barcode, !!!barcodes),
         time = recode_factor(time, !!!c(
             'Spring' = '1',
             'Summer' = '2',
             'Fall' = '3')),
         time = factor(time, levels = c('1', '2', '3')))

period_colour_num <- c(
  '1' = unname(period_colour[1]),
  '2' = unname(period_colour[2]),
  '3' = unname(period_colour[3])
)

for (ci in cities) {
  
  p.ls <- map(barcodes, function(barcode) {
    
    PCo <- eig.df %>% 
      filter(City == ci & Dist == 'bray' & Barcode == barcode) %>% 
      select(PCo1, PCo2) %>% as.vector
    
    pcoa.df %>% 
      filter(City == ci & Dist == 'bray' & Barcode == barcode) %>%
      ggplot(aes(x = PCo1, y = PCo2, colour = time)) +
      geom_point(size = 1) +
      stat_ellipse(level = 0.95, geom = 'polygon', 
                   alpha = 0.2, aes(fill = time)) +
      theme_light() +
      facet_grid(.~Barcode) +
      scale_colour_manual(values = period_colour_num) +
      scale_fill_manual(values = period_colour_num) +
      labs(colour = 'Sampling\nperiod', fill = 'Sampling\nperiod', 
           x = PCo$PCo1, y = PCo$PCo2) +
      theme(axis.title.y = element_text(margin = margin(r = -5, l=5)))
  })
  
  p <- p.ls$BACT + p.ls$FUNG + p.ls$PLAN +
    plot_layout(guides = 'collect')
  
  ggsave(paste0('~/Desktop/ip34/urbanBio/out/pcoa_', ci,'.pdf'),
         plot = p, bg = 'white', width = 2400, height = 1200,
         units = 'px', dpi = 200)
}
 
# X4 Community composition all samples 
# Same as 2., but order by date within facets, time time_period nested in barcode


# X2 Median income vs. vegetation index design (dotplot)
