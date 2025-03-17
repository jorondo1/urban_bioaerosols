library(pacman)
p_load(tidyverse, RColorBrewer, phyloseq, patchwork)
urbanbio.path <- '~/Desktop/ip34/urbanBio'

source(file.path(urbanbio.path, 'scripts/myFunctions.R'))
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/community_functions.R")

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
# One plot grid per city
# BFP columns, Idx rows

ps_rare.ls <- read_rds(file.path(urbanbio.path, 'data/ps_rare.ls.rds'))

###################
# Hill Numbers #####
###################
cities <- c('Montréal' =, 'Québec', 'Sherbrooke')
month_colours = c('May' = 'springgreen3', 
                  'June' = 'skyblue3', 
                  'September' = 'orange3', 
                  'October' = 'mediumpurple3')

barcodes <- c('BACT' = "Bacteria",
              'FUNG' = "Fungi",
              'PLAN' = "Plants")

Hill_indices <- c('H_0' = 'Richness',
                  'H_1' = 'exp^(Shannon)',
                  'H_2' = 'Inverse Simpson')

div_Hill <- imap(ps_rare.ls, function(ps, barcode) {
  div.fun(ps,  c(0,1,2)) %>% 
    data.frame %>% 
    rownames_to_column('Sample') %>% 
    left_join(samdat_as_tibble(ps), 
              by = 'Sample') %>% 
    mutate(barcode = barcode)
}) %>% list_rbind %>% 
  pivot_longer(cols = c('H_0', 'H_1', 'H_2', 'Tail'), 
               names_to = 'Index',
               values_to = 'Effective_taxa')

div_Hill %<>%
  mutate(time = factor(case_when(
    date <= "2022-05-19" ~ 'May',
    date > "2022-05-30" & date <= "2022-06-30" ~ 'June',
    date > "2022-08-30" & date <= "2022-09-10" ~ 'September',
    date > "2022-09-10" ~ 'October',
    TRUE ~ NA
  ), levels = names(month_colours)),
  barcode = recode_factor(barcode, !!!barcodes),
  Index = recode_factor(Index, !!!Hill_indices))


for (ci in cities) {
  p <- list()
  
  for (bc in barcodes) {
    p[[bc]] <- div_Hill %>% 
      filter(city == ci & Index != 'Tail' & barcode == bc) %>% 
      ggplot(aes(x = time, y = Effective_taxa, fill = time)) +
      geom_boxplot(outliers = FALSE) +
      facet_grid(Index~barcode, scales = 'free') +
      theme_light() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            legend.position = 'bottom') +
      scale_fill_manual(values = month_colours) +
      labs(y = 'Effective number of taxa') +
      ylim(0,NA) 
  }
  
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


# X4 Community composition all samples 
# Same as 2., but order by date within facets, time period nested in barcode


# X2 Median income vs. vegetation index design (dotplot)