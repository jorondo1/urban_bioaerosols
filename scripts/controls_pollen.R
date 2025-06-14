library(pacman)
p_load(tidyverse, RColorBrewer, phyloseq, patchwork, magrittr, forcats)
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/psflashmelt.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/rarefy_even_depth2.R")
source("scripts/0_config.R")

ps_pas <- read_rds('data/ps_passive.rds')
ps_act <- read_rds('data/ps.ls.rds')[['PLAN']]

which_taxrank <- 'Order'
keep_taxranks <- all_of(taxRanks[1:which(taxRanks==which_taxrank)])

# Passive traps
melted_pas <- ps_pas %>% 
  tax_glom(taxrank = which_taxrank) %>% 
  rarefy_even_depth2(rngseed = 123) %>% 
  psflashmelt() %>% 
  filter(City == 'MTL') %>% 
  # define time periods
  mutate(time = case_when(
    TX %in% c('T2', 'T3') ~ 'Spring',
    TX == 'T5' ~ 'Summer',
    TX %in% c('T11', 'T12') ~ 'Fall')) %>% 
  select(Sample, Abundance, Start_date, site_id, time, all_of(keep_taxranks)) %>% 
  mutate(Type = 'P',
         date = Start_date,
         keep = 'unused')

# Active traps 
melted_act <- ps_act %>% 
  tax_glom(taxrank = which_taxrank) %>% 
  rarefy_even_depth2(rngseed = 123) %>% 
  psflashmelt() %>% 
  filter(city == 'Montreal') %>% 
  select(Sample, Abundance, site_id, date, time, all_of(keep_taxranks)) %>% 
  mutate(Type = 'A')

melted_all <- rbind(melted_pas, melted_act) %>% 
  filter(Abundance > 0) %>% 
  group_by(Sample) %>% 
  mutate(relAb = Abundance / sum(Abundance)) %>% 
  # Ensure samples are in chronological order
  mutate(Sample = fct_reorder(Sample, date, .fun = min))

# Check if site_id x time x Type is unique
plot.df %>% 
  select(Sample, site_id, time, Type) %>% 
  distinct() %>% 
  group_by(site_id, time, Type) %>% 
  summarise(n = n()) 

# Number of top taxa
nTaxa <- 16
expanded_palette <- colorRampPalette(brewer.pal(12, 'Paired'))(nTaxa+2) 

top_taxa <- topTaxa(melted_all, taxLvl = which_taxrank, topN = nTaxa)

# Define top taxa levels
top_taxa_lvls <- top_taxa %>%
  group_by(aggTaxo) %>% 
  aggregate(relAb ~ aggTaxo, data = ., FUN = sum) %>% 
  arrange(relAb) %$% aggTaxo %>% 
  as.character() %>% # Set Others and Unclassified first:
  setdiff(., c('Others', 'Unclassified')) %>% c('Others', 'Unclassified', .)

# Add top tax levels
plot.df <- melted_all %>% 
  left_join(top_taxa %>% select(-relAb), by = which_taxrank) %>%
  group_by(site_id, Type, time, date, aggTaxo) %>% 
  summarise(relAb = sum(relAb), .groups = 'drop') %>% 
  # Average rare cases of > samples by time x Type x site_id
  group_by(site_id, Type, time) %>% 
  mutate(relAb = relAb/sum(relAb),
         # set factors
         aggTaxo = factor(aggTaxo, levels = top_taxa_lvls),
         time = factor(time, levels = periods)) 

# Plot !
plot.df %>% 
  ggplot(aes(x = Type, y = relAb, fill = aggTaxo)) +
  geom_col() +
  labs(x = '', y= 'Relative abundance') +
  facet_grid(time~site_id,  space = 'free')+
  scale_fill_manual(values = expanded_palette) +
  theme_light() +
 # theme(axis.text.x = element_blank()) +
  labs(fill = which_taxrank, x = 'Trap type') 

ggsave(filename = 'out/passives/composition_Order.pdf',
       bg = 'white', width = 2600, height = 1400, 
       units = 'px', dpi = 220)
