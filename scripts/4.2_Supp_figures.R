
library(pacman)
p_load(tidyverse, phyloseq, patchwork, ggh4x, #facet_nested
       magrittr, readxl, RColorBrewer)

source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/tax_glom2.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/rarefy_even_depth2.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/psflashmelt.R')
source('scripts/myFunctions.R')

tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
ps.ls <- readRDS('data/ps.ls.rds')

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
    theme_light() +
    facet_nested(cols=vars(city,time), scales = 'free', space = 'free') +
    scale_fill_manual(values = expanded_palette) +
    labs(fill = which_taxrank, 
         x = 'Samples ordered by sampling date',
         y = paste('Relative abundance of', barcodes[[barcode]], 'sequences')) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_blank()) +
    guides(fill = guide_legend(ncol = 1))
  
  ggsave(paste0('out/composition_',which_taxrank,'_',barcode,'.pdf'),
         bg = 'white', width = 3600, height = 2000, 
         units = 'px', dpi = 220)
})

##########
# qPCR quantification
##########################


# Add qpcr data
parse_bacterial_load <- function(input_tibble) {
  require(magrittr, dplyr)
  input_tibble %>% 
    filter(is.na(note)) %>% 
    mutate(copy_number = as.numeric(copy_number))
}

bact_load <- read_xlsx('data/metadata/load_bacteria_16S.xlsx') %>% 
  set_names(c('Sample', 'replicate', 'copy_number', 'note')) %>% 
  parse_bacterial_load() %>% 
  group_by(Sample) %>% 
  summarise(mean_copy_number = mean(copy_number),
            sd_copy_number = sd(copy_number),
            n = n())

bact_load %>%  filter(n == 1) # at least 2 each, ok
bact_load %<>% # remove insanely high sds
  filter(! sd_copy_number > 0.25*mean_copy_number )

fung_load <- read_xlsx('data/metadata/load_fungi_ITS.xlsx') %>% 
  set_names(c('Sample', 'replicate', 'copy_number_16S','copy_number', 'note')) %>%
  parse_bacterial_load()

fung_load%>%  filter(n == 1) # many, not ideal

# Let's start with 16S... 

bact_qpcr.df <- comm_plot_data.ls$BACT %>% 
  inner_join(bact_load, by = 'Sample') %>% 
  mutate(qpcr_ab = relAb * mean_copy_number) 

bact_qpcr.df %>% 
  ggplot(aes(x = site_date, y = qpcr_ab, fill = aggTaxo)) +
  geom_col() +
  theme_light() +
  facet_nested(cols=vars(city,time), scales = 'free', space = 'free') +
  scale_fill_manual(values = expanded_palette) +
  labs(fill = which_taxrank, 
       x = 'Samples ordered by sampling date',
       y = paste('Relative abundance of', barcodes[['BACT']], 'sequences')) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 1))

ggsave(paste0('out/composition_QPCR_',which_taxrank,'_BACT.pdf'),
       bg = 'white', width = 3600, height = 2000, 
       units = 'px', dpi = 220)
