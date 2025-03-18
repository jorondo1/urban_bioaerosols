library(pacman)
p_load(tidyverse, RColorBrewer, phyloseq, patchwork)
urbanbio.path <- '~/Desktop/ip34/urbanBio'

source(file.path(urbanbio.path, 'scripts/myFunctions.R'))
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/community_functions.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/tax_glom2.R")

# Agglomerate
ps_rare_Genus.ls <- lapply(ps_rare.ls, function(ps){
  tax_glom2(ps, taxrank = "Genus")
})

#######################################
# X3. Taxonomic classification rates ###
#######################################
# Informs following plots on which level to focus

#########################################
# 2. Community composition barcharts #####
###########################################
# Which level ??
# BFP grid / wrap with nested period, MQS rows, mean per time period
# Time period represented by ~median date ? distinct by city


### From here on, individual plots will be patched together
### per city, to contain both diversity analyses as well as
### some differential testing and maybe nestedness? 

#########################
# 3. Alpha diversity #####
###########################

ps_rare.ls <- read_rds(file.path(urbanbio.path, 'data/ps_rare.ls.rds'))

cities <- c('Montréal' = 'Montreal', 'Québec' = 'Quebec', 'Sherbrooke'='Sherbrooke')

month_colours = c('May' = 'springgreen3', 
                  'June' = 'skyblue3', 
                  'September' = 'orange3', 
                  'October' = 'darkred')

barcodes <- c('BACT' = "Bacteria",
              'FUNG' = "Fungi",
              'PLAN' = "Plants")

Hill_indices <- c('H_0' = 'Richness',
                  'H_1' = 'exp^(Shannon)',
                  'H_2' = 'Inverse Simpson')

mutate_time <- function(data) {
  data %>% 
    mutate(time = factor(case_when(
    date <= "2022-05-19" ~ 'May',
    date > "2022-05-30" & date <= "2022-06-30" ~ 'June',
    date > "2022-08-30" & date <= "2022-09-10" ~ 'September',
    date > "2022-09-10" ~ 'October',
    TRUE ~ NA
    ), levels = names(month_colours)),
    city = recode_factor(city, !!!cities))
}

# Compute diversity 
div_Hill <- imap(ps_rare_Genus.ls, function(ps, barcode) {
  
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
      scale_fill_manual(values = month_colours) +
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

#########################
# 4. Beta diversity ######
###########################
# Betadiv en BFP
# 2 rows : BC + Aitchison? Ou UniFrac pour les bact?
# 3rd row : Turnover/nestedness
# Include variance partitioning?

# Split every dataset by city
ps_byCity.ls <- lapply(cities, function(city_name) { # 1st level: city
  lapply(ps_rare_Genus.ls, function(ps) { # 2nd level: barcode
    
    # Subset samples by city
    selected_samples <- rownames(sample_data(ps))[sample_data(ps)$city == city_name]
    
    # Subset the phyloseq object and remove zero-count taxa
    prune_samples(selected_samples, ps) %>%
      prune_taxa(taxa_sums(.) > 0, .)
  })
}); names(ps_byCity.ls) <- cities

# Iterate permanova for each city, barcode, distance
pcoa_genus.ls <- imap(ps_byCity.ls, function(ps.ls, city) {
  imap(ps.ls, function(ps, barcode){
    out <- list()
    out[['robust.aitchison']] <- compute_pcoa(ps, dist = 'robust.aitchison')
    out[['bray']] <- compute_pcoa(ps, dist = 'bray')
    out
  })
})

eig.df <- imap(pcoa_genus.ls, function(pcoa_barcode.ls, city) {
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

pcoa.df <- imap(pcoa_genus.ls, function(pcoa_barcode.ls, city) {
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
         time = factor(time, levels = names(month_colours)))
 
for (ci in cities) {

  p <- pcoa.df %>% 
    filter(City == ci & Dist == 'bray') %>%
    ggplot(aes(x = PCo1, y = PCo2, colour = time)) +
    geom_point(size = 1) +
    stat_ellipse(level = 0.95, geom = 'polygon', 
                 alpha = 0.2, aes(fill = time)) +
    theme_light() +
    facet_grid(.~Barcode) +
    scale_colour_manual(values = month_colours) +
    scale_fill_manual(values = month_colours) +
    geom_text(data = eig.df %>% 
                filter(Dist == 'bray' & City == ci),
              aes(x = -1.2, y = -2, label = paste(PCo1, PCo2)),
              inherit.aes = FALSE, 
              hjust = 0, vjust = 0) +
    labs(colour = 'Season', fill = 'Season')
  
  ggsave(paste0('~/Desktop/ip34/urbanBio/out/pcoa_', ci,'.pdf'),
         plot = p, bg = 'white', width = 2200, height = 1200,
         units = 'px', dpi = 220)
}

# X4 Community composition all samples 
# Same as 2., but order by date within facets, time period nested in barcode


# X2 Median income vs. vegetation index design (dotplot)