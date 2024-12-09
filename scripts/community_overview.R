library(pacman)
p_load(dada2, tidyverse, magrittr, RColorBrewer, phyloseq, ggh4x, Biostrings, ggridges, optparse)

source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/tax_glom2.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/rarefy_even_depth2.R')
source('scripts/myFunctions.R')

which_taxrank <- 'Family'
melted <- read_rds('data/trnL/5_out/ps_genus.RDS') %>% 
  tax_glom2(taxrank = which_taxrank)
  psmelt %>%
  filter(Abundance != 0) %>% 
  group_by(Sample) %>% 
  mutate(relAb = Abundance/sum(Abundance),
         time = factor(time, levels=c('Spring', 'Summer', 'Fall'))) %>% 
  ungroup

nTaxa <- 16
top_taxa <- topTaxa(melted, which_taxrank, nTaxa)
top_taxa_lvls <- top_taxa %>% 
  group_by(aggTaxo) %>% 
  aggregate(relAb ~ aggTaxo, data = ., FUN = sum) %>% 
  arrange(relAb) %$% aggTaxo %>% 
  as.character %>% 
  setdiff(., c('Others', 'Unclassified')) %>% c('Others', 'Unclassified', .)

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
