library(pacman)
p_load(tidyverse, RColorBrewer, phyloseq, patchwork, magrittr)
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/psflashmelt.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/rarefy_even_depth2.R")

ps_pas.ls <- read_rds('data/ps_passive.rds')

melted <- ps_pas.ls %>% 
  tax_glom(taxrank = 'Family') %>% 
  rarefy_even_depth2(rngseed = 123) %>% 
  psflashmelt() %>% 
  filter(Abundance > 0) %>% 
  group_by(Sample) %>% 
  mutate(relAb = Abundance / sum(Abundance))

which_taxrank <- 'Order'

nTaxa <- 16
expanded_palette <- colorRampPalette(brewer.pal(12, 'Paired'))(nTaxa+2) 

top_taxa <- topTaxa(melted, taxLvl = which_taxrank, topN = nTaxa)

# Define top taxa levels
top_taxa_lvls <- top_taxa %>%
  group_by(aggTaxo) %>% 
  aggregate(relAb ~ aggTaxo, data = ., FUN = sum) %>% 
  arrange(relAb) %$% aggTaxo %>% 
  as.character() %>% # Set Others and Unclassified first:
  setdiff(., c('Others', 'Unclassified')) %>% c('Others', 'Unclassified', .)

plot.df <- melted %>% 
  left_join(top_taxa %>% select(-relAb), by = which_taxrank) %>%
  group_by(Sample, Trap, `T`, City, aggTaxo) %>% 
  summarise(relAb = sum(relAb), .groups = 'drop') %>% 
  mutate(aggTaxo = factor(aggTaxo, levels = top_taxa_lvls))

plot.df %>% 
  filter(City == 'MTL') %>% 
  ggplot(aes(x = Sample, y = relAb, fill = aggTaxo)) +
  geom_col() +
  labs(x = '', y= 'Relative abundance') +
  facet_grid(.~Trap, scales = 'free', space = 'free')+
  scale_fill_manual(values = expanded_palette) +
  theme_minimal()+
  theme(axis.text.x = element_blank()) +
  labs(fill = which_taxrank) 

ggsave(filename = 'out/passives/composition_Order.pdf',
       bg = 'white', width = 2600, height = 1400, 
       units = 'px', dpi = 220)
