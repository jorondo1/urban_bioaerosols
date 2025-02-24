library(pacman)
p_load(dada2, tidyverse, magrittr, RColorBrewer, ggdist, tidyquant,
       phyloseq, ggh4x, Biostrings, ggridges, optparse)

source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/tax_glom2.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/rarefy_even_depth2.R')
source('~/Desktop/ip34/urbanBio/scripts/myFunctions.R')

tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
ps.ls <- readRDS('~/Desktop/ip34/urbanBio/data/ps.ls.rds')
ps_ctrl.ls <- readRDS('~/Desktop/ip34/urbanBio/data/ps_ctrl.ls.rds')

###################################
# Community composition overview ###
#####################################

which_taxrank <- 'Family'
melted <- ps.ls$PLAN %>% 
  tax_glom(taxrank = which_taxrank) %>% 
  psmelt %>%
  filter(Abundance != 0) %>% 
  group_by(Sample) %>% 
  mutate(relAb = Abundance/sum(Abundance),
         time = factor(time, levels=c('Spring', 'Summer', 'Fall'))) %>% 
  ungroup

# Compute top taxa and create "Others" category
nTaxa <- 10
(top_taxa <- topTaxa(melted, which_taxrank, nTaxa))
top_taxa_lvls <- top_taxa %>% 
  group_by(aggTaxo) %>% 
  aggregate(relAb ~ aggTaxo, data = ., FUN = sum) %>% 
  arrange(relAb) %$% aggTaxo %>% 
  as.character %>% # Others first:
  setdiff(., c('Others', 'Unclassified')) %>% c('Others', 'Unclassified', .)

# Plot !
expanded_palette <- colorRampPalette(brewer.pal(12, 'Set3'))(nTaxa+2) 

melted %>% 
  left_join(top_taxa %>% select(-relAb), by = which_taxrank) %>%
  filter(!is.na(time)) %>% 
  mutate(aggTaxo = factor(aggTaxo, levels = top_taxa_lvls)) %>% 
  ggplot(aes(x = Sample, y = relAb, fill = aggTaxo)) +
  geom_col() +
  theme_light() +
  facet_nested(cols=vars(city,time), scales = 'free', space = 'free') +
  scale_fill_manual(values = expanded_palette) +
  labs(fill = which_taxrank) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank())

ggsave(paste0('~/Desktop/ip34/urbanBio/out/composition_',which_taxrank,'.pdf'),
       bg = 'white', width = 3600, height = 2000, 
       units = 'px', dpi = 220)

###################################
# Controls ###
#####################################

which_taxrank <- 'Genus'

control.melt <- imap(ps_ctrl.ls, function(ps, barcode) {
  ps %>% 
    tax_glom(taxrank = which_taxrank) %>% 
    psmelt %>% 
    select(Sample, Abundance, site_id, !!sym(which_taxrank)) %>% 
    filter(Abundance != 0) %>% 
    group_by(Sample) %>% 
    mutate(relAb = Abundance, # NOT RELATIVE ABUNDANCES HERE, keeping the var cause it's hardcoded in topTaxa
           barcode = barcode)
}) %>% list_rbind %>% 
  filter(!Sample %in% c('D8-3046', 'D8-3044', 'Couloir-D8-3225'))

# Compute top taxa and create "Others" category
nTaxa <- 16
(top_taxa <- topTaxa(control.melt, which_taxrank, nTaxa))
top_taxa_lvls <- top_taxa %>% 
  group_by(aggTaxo) %>% 
  aggregate(relAb ~ aggTaxo, data = ., FUN = sum) %>% 
  arrange(relAb) %$% aggTaxo %>% 
  as.character %>% # Others first:
  setdiff(., c('Others', 'Unclassified')) %>% c('Others', 'Unclassified', .)

# Plot !
expanded_palette <- colorRampPalette(brewer.pal(12, 'Paired'))(nTaxa+2) 

control.melt %>% 
  left_join(top_taxa %>% select(-relAb), by = which_taxrank) %>%
  mutate(aggTaxo = factor(aggTaxo, levels = top_taxa_lvls),
         barcode = recode(barcode, !!!barcodes)) %>% 
  ggplot(aes(y = Sample, x = Abundance, fill = aggTaxo)) +
  geom_col() +
  theme_light() +
  facet_grid(~barcode, scales = 'free') +
  scale_fill_manual(values = expanded_palette) +
  labs(fill = which_taxrank) +
  theme(panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank()) +
  labs(x = 'Total number of sequences')

ggsave(paste0('~/Desktop/ip34/urbanBio/out/composition_',which_taxrank,'_controls.pdf'),
       bg = 'white', width = 3600, height = 2000, 
       units = 'px', dpi = 220)
