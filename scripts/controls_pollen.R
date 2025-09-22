library(pacman)
p_load(tidyverse, RColorBrewer, phyloseq, patchwork, magrittr, forcats)
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/psflashmelt.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/rarefy_even_depth2.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/tax_glom2.R")
source("scripts/0_config.R")


# PASSIVE vs. ACTIVE TRAPS ------------------------------------------------

ps_pas <- read_rds('data/ps_passive.rds')
ps_act <- read_rds('data/ps.ls.rds')[['PLAN']]

which_taxrank <- 'Family'
keep_taxranks <- taxRanks[1:which(taxRanks==which_taxrank)]


# Passive traps
melted_pas <- ps_pas %>% 
  tax_glom2(taxrank = which_taxrank) %>% 
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
         .keep = 'unused')

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
  labs(fill = which_taxrank, x = 'Trap type') +
  guides(fill = guide_legend(ncol = 1))

ggsave(filename = 'out/passives/composition_Family.pdf',
       bg = 'white', width = 2600, height = 1400, 
       units = 'px', dpi = 220)


# POSITIVE CONTROLS -------------------------------------------------------

trnL_seqtab <- read_rds('~/Desktop/ip34/urbanBio/data/trnL/4_taxonomy_E22_100_trunc275/seqtab.RDS')
trnL_tax <-  read_rds('~/Desktop/ip34/urbanBio/data/trnL/4_taxonomy_E22_100_trunc275/taxonomy.RDS')

pos_ctrl <- trnL_seqtab['ctrl-PCR-pos',]
pos_ctrl <- c(pos_ctrl[which(pos_ctrl>0)])
tax_ctrl <- data.frame(trnL_tax[names(pos_ctrl),])
tax_ctrl$seq <- pos_ctrl[row.names(tax_ctrl)]

tax_ctrl %>% 
  tibble() %>% 
  group_by(Family, Genus) %>% 
  summarise(Abundance = sum(seq), .groups = 'drop') %>%  
  arrange(Family) %>% 
  mutate(relAb = Abundance/sum(Abundance))

tax_ctrl_relab <- tax_ctrl %>% 
  tibble() %>% 
  group_by(Family) %>% 
  summarise(Abundance = sum(seq)) %>%  
  arrange(Family) %>% 
  mutate(`Observed abundance` = Abundance/sum(Abundance))

syncom <- readxl::read_xlsx('~/Desktop/ip34/urbanBio/data/syncom.xlsx')

# Use the Angiosperm Phylogeny Group, as in NCBI taxonomy
syncom_corr_G <- syncom %>% 
  mutate(Family = case_when(
    Family == 'Aceraceae' ~ 'Sapindaceae',
    TRUE ~ Family
  )) %>% 
  group_by(Family) %>% 
  summarise(
    dna_conc = sum(dna_conc) ,
    .groups = 'drop'
  ) %>% 
  mutate(`Expected abundance` = dna_conc/sum(dna_conc))


relAb_comparison <- full_join(
  tax_ctrl_relab,
  syncom_corr_G,
  by = 'Family'
) %>% mutate(across(where(is.numeric), ~coalesce(., 0)))

relAb.plot <- relAb_comparison %>% 
  select(-Abundance, -dna_conc) %>% 
  pivot_longer(cols = -Family) %>% 
  ggplot(aes(x = name, y = value, fill = Family)) +
  geom_col() +
  theme_light() +
  facet_grid(~name, scales = 'free') +
  scale_fill_manual(values = c(
    "#78c52a",
    "#0052a0",
    "#f7ba11",
    "#3ca5ff",
    "#c6ce3c",
    "#ff8ad7",
    "#01c55d",
    "#a30a4a",
    "#01c996",
    "#ff7f4d",
    "#02adae",
    "#903518",
    "#78d8b0",
    "#803a63",
    "#8e6800",
    "#ecbe8d",
    "#5d5309"
  )) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  labs(y = 'Relative abundance'); relAb.plot 

ggsave('~/Repos/urbanBio/out/syncom/relAb.pdf', 
       bg = 'white', width = 1200, height = 900, 
       units = 'px', dpi = 180)

relAb.plot +
  labs(y = 'Abondance relative',
       fill = 'Famille') +
  facet_grid(
    ~name, scales = 'free',
    labeller = labeller(
      name = c("Expected abundance" = "Abondances attendues",
               "Observed abundance" = "Abondances observées")
  ))

ggsave('~/Library/CloudStorage/OneDrive-FreigegebeneBibliotheken–USherbrooke/Isabelle Laforest-Lapointe - RONDEAU_LECLAIRE_Jonathan/Memoire/Main/figures/relAb.pdf', 
       bg = 'white', width = 1200, height = 900, 
       units = 'px', dpi = 180)

ratios <- relAb_comparison %>% 
  filter(dna_conc>0 & Abundance > 0) %>% 
  mutate(Sequences = log10(Abundance) - mean(log10(Abundance)),
         DNA = log10(dna_conc) - mean(log10(dna_conc))) %>% 
  select(Family, Sequences, DNA) %>% 
  #mutate(across(where(is.numeric), ~case_when(.x==NaN~0, TRUE~.x))) %>% 
  pivot_longer(cols = -Family) 

ratios %>% 
  ggplot(aes(x = Family, y = value, fill = name)) +
  geom_col(position = 'dodge') +
  theme_light() +
  facet_grid(~Family, scales = 'free') +
  labs(y = 'Log10-fold increase relative to mean') +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.7,0.3),
        legend.background = element_rect(
          fill = "white", color = "black",
          linewidth = 0.2)
  )
ggsave('~/Repos/urbanBio/out/syncom/clr.pdf', 
       bg = 'white', width = 1700, height = 900, 
       units = 'px', dpi = 180)






