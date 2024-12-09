library(pacman)
p_load(dada2, tidyverse, magrittr, RColorBrewer, ggdist, tidyquant,
       phyloseq, ggh4x, Biostrings, ggridges, optparse)

source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/tax_glom2.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/rarefy_even_depth2.R')
source('scripts/myFunctions.R')

tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
# 
# ### PHYLOSEQ objects list 
# #ps.list <- list()
# ps.list[[paste0('trnL_Jo_',suffix)]] <- phyloseq(tax_table(taxa.subset),
#          otu_table(seqtab.subset, taxa_are_rows = FALSE),
#          sample_data(empty_df)) 
# write_rds(ps.list, 'merge_parameters_tests.ps.rds')
# 
# ps.list[['ITS']] <- read_rds('data/sarah/phyloITS.rds') %>% 
#   rarefy_even_depth2()
# 
# ps.list[['16S']] <- read_rds('data/sarah/phylo16S.rds') %>% 
#   rarefy_even_depth2()
# 
# ps.list <- read_rds('data/merge_parameters_tests.ps.rds')
# ps_sarah <- read_rds('data/metaITS.rds')


# Community composition overview
which_taxrank <- 'Family'
melted <- ps_trnL %>% 
  tax_glom(taxrank = which_taxrank) %>% 
  psmelt %>%
  filter(Abundance != 0) %>% 
  group_by(Sample) %>% 
  mutate(relAb = Abundance/sum(Abundance),
         time = factor(time, levels=c('Spring', 'Summer', 'Fall'))) %>% 
  ungroup
  
nTaxa <- 16
(top_taxa <- topTaxa(melted, which_taxrank, nTaxa))
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
  
ggsave('~/Desktop/ip34/urbanBio/out/composition_Family.pdf', bg = 'white', width = 3600, height = 2000, 
       units = 'px', dpi = 220)

ps.list <- list()
ps.list[['test']] <- ps_trnL

### Relationship between ASV length and classifiability?
tax_tibble <- ps_trnL@tax_table %>%
  data.frame() %>% 
  rownames_to_column('OTU') %>% tibble

length_classification <- tax_tibble %>% rowwise %>% 
  mutate(classification_resolution = {
    # Get the vector of column names for non-"unclassified" values
    valid_cols <- tax_ranks[unlist(across(all_of(tax_ranks), ~ .x != "Unclassified"))]
    # Pick the last column name or NA if none exists
    if (length(valid_cols) > 0) tail(valid_cols, 1) else 'No_classification'
  }) %>%
  ungroup() %>% 
  mutate(asv_len = nchar(OTU),
         classification_resolution = factor(classification_resolution,
                                  levels = c('No_classification',tax_ranks))) %>% 
  select(classification_resolution, asv_len)

# classification resolution of our ASVs along their length.
# This shows the length distribution at the highest classification 
# (so each colour represents a distinct set of ASVs)

# First, count ASVs per highest classification
count_labels <- length_classification %>% 
  dplyr::count(classification_resolution) %>% 
  mutate(label = paste0(classification_resolution, " (n = ", n, ")")) %>%
  select(-n) %>% 
  deframe()

length_classification %>% 
  ggplot(aes(x = asv_len, y = classification_resolution)) +
  geom_boxplot(position = position_nudge(y = 0.2), width = 0.4) +
  geom_density_ridges2(
    aes(fill = classification_resolution ),
    alpha = 0.6, show.legend = FALSE,
    scale = 3
    ) +
  scale_y_discrete(expand = c(0.01, 0),
                   labels = count_labels) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(x = 'ASV length', y = 'Highest classification resolution')+
  theme_classic() 
  
### What proportion of sequences are unclassified at each taxlev ?
### using all 3 datasets


# Full melt + dataset name as "marker" variable
# melted_all <- imap(ps.list, function(ps, dataset) {
#   psmelt(ps) %>% 
#     filter(Abundance != 0) %>% 
#     select(all_of(tax_ranks), Abundance, Sample, time, city) %>% 
#     mutate(marker = dataset,
#            across(all_of(tax_ranks),~ replace_na(., "Unclassified")))
# }) %>% list_rbind

melted_all <- psmelt(ps_trnL) %>% 
  select(all_of(tax_ranks), Abundance, Sample, time, city) %>% 
      mutate(across(all_of(tax_ranks),~ replace_na(., "Unclassified")))
  
# # Extract proportion of unclassified ASV for each rank per sample/approach
test <- map_dfr(tax_ranks, function(taxRank) {
  melted_all %>%
    group_by(Sample, !!sym(taxRank), time) %>%
    summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
    group_by(Sample) %>%
    mutate(relAb = Abundance/sum(Abundance)) %>%
    filter(!!sym(taxRank) == "Unclassified") %>%
    mutate(Rank = factor(taxRank, levels = tax_ranks))
})


# Boxplot
test %>% 
  filter(#marker %in% tests_list &
           time != 'NA') %>% 
  ggplot(aes(x = time, y = relAb, fill = time)) +
  geom_boxplot(linewidth = 0.2) +
  labs(y = 'Proportion of unclassified ASVs') +
  facet_grid(.~Rank) +
  theme_light() +
  labs(x = '') + 
  theme(
    axis.text.x = element_blank()
  ) +
  scale_fill_brewer(palette = "Set2", na.value = "beige")

ggsave('out/unclassified_Jo.pdf', bg = 'white', width = 1600, height = 2400, 
       units = 'px', dpi = 180)

# or ridge plot 
test %>%
  filter(Rank != 'Kingdom' &
           marker %in% tests_list &
           time != 'NA') %>% 
  ggplot(aes(x = relAb, y = Rank, fill = time)) +
  geom_density_ridges(
    linewidth = 0.2,
    scale = 0.9, 
    alpha = 0.7
  ) +
  facet_grid(marker~ ., scales = 'free_x') +
 theme_light() +
  labs(
    x = 'Proportion of unclassified ASVs',
    y = 'Taxonomic Rank'
  ) +
  theme(
    axis.title.y = element_text(vjust = 1.5),
    axis.title.x = element_text(vjust = -1)
    ) +
  scale_fill_brewer(palette = "Set2", na.value = "beige")

ggsave('out/unclassified_ridge.pdf', bg = 'white', width = 1600, height = 1200, 
       units = 'px', dpi = 180)

# Or even 'survival' lines?
test %>% 
  filter(marker %in% tests_list &
           time != 'NA') %>% 
  group_by(Rank, marker) %>% 
  summarise(
    mean_unclassified = mean(relAb),
    sd_unclassified = sd(relAb), 
    .groups = 'drop'
  ) %>% 
  ggplot(aes(x = Rank, y = mean_unclassified, group = marker, colour = marker)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean_unclassified - sd_unclassified,
                  ymax = mean_unclassified + sd_unclassified,
                  fill = marker), linetype = 'dotted',
              alpha = 0.1) +
  theme_minimal() +
  labs(
    y = 'Proportion of unclassified ASVs',
    x = 'Taxonomic Rank'
  )
  
# Look at sequence length in unclassified data :
melted_length <- ps.list$trnL_Jo_rescued50_pooled_EE42 %>% 
  psmelt %>% 
  filter(Abundance != 0) %>% 
  mutate(seqlen = str_length(OTU))

classified_bool <- map_dfr(tax_ranks, function(taxRank) {
  melted_length %>%
    mutate(classified = case_when(!!sym(taxRank) == 'Unclassified' ~ FALSE,
                                  TRUE~ TRUE),
           Rank = factor(taxRank, levels = tax_ranks))
  })
  
  
classified_bool %>% 
  ggplot(aes(x = seqlen, y = Rank, fill = classified)) +
    geom_density_ridges(
      aes(height = sqrt(after_stat(count))),
      linewidth = 0.1,
      stat = 'binline',
      binwidth = 2,
      scale=0.9, alpha = 0.7
      ) + theme_minimal() +
  labs(y = 'square root of read count')

ggsave('out/classification_seqlen.pdf', 
         bg = 'white', width = 1200, height = 1200, 
         units = 'px', dpi = 180)
  






















