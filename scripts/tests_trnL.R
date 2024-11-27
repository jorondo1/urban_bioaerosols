library(pacman)
p_load(dada2, tidyverse, magrittr, RColorBrewer, phyloseq, ggh4x, Biostrings, ggridges, optparse)

source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/tax_glom2.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/rarefy_even_depth2.R')
source('scripts/myFunctions.R')

##############################
# FILTER REFERENCE DATABASE ###
################################
trnl.ref <- 'data/trnL_hits.lineage.filtered.Genus.fa'

# check sequence length distribution
trnl.seq <- readDNAStringSet(trnl.ref, format = 'fasta')
seqlen <- width(trnl.seq)
ggplot(data.frame(length = seqlen), aes(x = length)) +
  geom_histogram(binwidth = 10, fill = "blue", color = "black") +
  labs(title = "Distribution of Sequence Lengths", x = "Sequence Length", y = "Frequency") +
  theme_minimal()

# filter within expected limits for trnL (???)
filtered_sequences <- trnl.seq[seqlen >= 300 & 
                                 width(trnl.seq) <= 800]

filtered_sequences %>% width %>% sort %>% hist
writeXStringSet(filtered_sequences, "/filtered_trnl_ref.fa")

######################
### ASSIGN TAXONOMY ###
########################
dirpath <- 'data'
suffix <- 'rescued50_pooled_EE42'
raw_sequences <- read_rds(paste0(dirpath,'/seqtab_',suffix,'.rds'))

# Filter out sequences smaller than 200bp
keep <- nchar(colnames(raw_sequences)) >= 200
raw_sequences200 <- raw_sequences[, keep, drop = FALSE]

# Assign taxonomy
taxa <- assignTaxonomy(raw_sequences200, 'data/filtered_trnl_ref.fa', multithread = TRUE, tryRC = TRUE)
taxa[taxa == ""] <- NA # some are NA, others ""
taxa[is.na(taxa)] <- 'Unclassified' # otherwise Tax_glom flushes out the NAs .
write_rds(taxa, paste0(dirpath,'/taxonomy_',suffix,'.rds'))
taxa <- read_rds(paste0(dirpath,'/taxonomy_',suffix,'.rds'))

#meta <- read_rds('data/metaTrnL.rds') #not usable yet

### TEMPORARY METADATA
  # Some samples are nearly empty, choose a rowsums cutoff:
  big_enough_samples <- rowSums(raw_sequences200)>1000
  seqtab.subset <- raw_sequences200[big_enough_samples,]
  
  detectable_taxa <- colSums(seqtab.subset) >100 # at least 100 sequences overall
  taxa.subset <- taxa[detectable_taxa,]
  
    # dummy samdata
  empty_df <- data.frame(Sample = rownames(seqtab.subset)) %>% 
    mutate(city = case_when(str_detect(Sample, 'MTL') ~ 'Montreal',
                            str_detect(Sample, 'SHER')~ 'Sherbrooke',
                            str_detect(Sample, 'QC') ~ 'Québec'),
           time = case_when(str_detect(Sample, '-A-') ~ "Fall",
                              str_detect(Sample, '-S-') ~ "Spring",
                              str_detect(Sample, '-E-') ~ "Summer")) %>% 
    column_to_rownames('Sample')
### /TEMPORARY

### PHYLOSEQ objects list 
#ps.list <- list()
ps.list[[paste0('trnL_Jo_',suffix)]] <- phyloseq(tax_table(taxa.subset),
         otu_table(seqtab.subset, taxa_are_rows = FALSE),
         sample_data(empty_df)) 
write_rds(ps.list, 'merge_parameters_tests.ps.rds')

ps.list[['ITS']] <- read_rds('data/sarah/phyloITS.rds') %>% 
  rarefy_even_depth2()

ps.list[['16S']] <- read_rds('data/sarah/phylo16S.rds') %>% 
  rarefy_even_depth2()



# Community composition overview
which_taxrank <- 'Family'
melted <- ps.list$trnL_Jo_rescued50_pooled_EE42 %>% 
  tax_glom(taxrank = which_taxrank) %>% 
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
  
ggsave('out/composition_Class_200.pdf', bg = 'white', width = 3600, height = 2000, 
       units = 'px', dpi = 220)



### What proportion of sequences are unclassified at each taxlev ?
### using all 3 datasets

tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# Full melt + dataset name as "marker" variable
melted_all <- imap(ps.list, function(ps, dataset) {
  psmelt(ps) %>% 
    filter(Abundance != 0) %>% 
    select(all_of(tax_ranks), Abundance, Sample, time, city) %>% 
    mutate(marker = dataset,
           across(all_of(tax_ranks),~ replace_na(., "Unclassified")))
}) %>% list_rbind

tests_list <- c('trnL_Sarah', 'trnL_Jo_concat',  'trnL_Jo_fwd','trnL_Jo_rescued', 'trnL_Jo_rescued50_pooled_EE42')

# Extract proportion of unclassified ASV for each rank per sample/approach
test <- map_dfr(tax_ranks, function(taxRank) {
  melted_all %>%
    filter(marker %in% tests_list) %>% 
    group_by(Sample, !!sym(taxRank), time, marker) %>% 
    summarise(Abundance = sum(Abundance), .groups = 'drop') %>% 
    group_by(Sample, marker) %>% 
    mutate(relAb = Abundance/sum(Abundance)) %>% 
    filter(!!sym(taxRank) == "Unclassified") %>%
  #  summarise(unclassified = mean(relAb, na.rm = TRUE)) %>% 
    mutate(Rank = factor(taxRank, levels = tax_ranks))
}) %>% 
  mutate(marker = factor(marker, levels = tests_list))

# Boxplot
test %>% 
  filter(marker %in% tests_list &
           time != 'NA') %>% 
  ggplot(aes(x = time, y = relAb, fill = time)) +
  geom_boxplot(linewidth = 0.2) +
  labs(y = 'Proportion of unclassified ASVs') +
  facet_grid(marker~Rank) +
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
  






















