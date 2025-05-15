# Migrate here from ps.ls creation
# Before rarefaction, impossible trnLs 
# Produce tables 
library(pacman)
p_load(phyloseq, tidyverse, kableExtra, ggridges)

# Import

source('scripts/myFunctions.R')
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/rarefy_even_depth2.R")
theme_set(theme_light())
ps.ls <- read_rds('data/ps.ls.rds')

# MANUAL CHECKS
# see current "test.trnL.R" for plant
# import decontam codes?

#############################
# ASV stats table #######
###############################
p_load(knitr, kableExtra, webshot2)

# Prep data
ps.stats <- imap(ps.ls, function(ps, barcode) {
  asv <- ps %>% otu_table 
  seq_per_sam <- rowSums(asv)
  asv_per_sam <- rowSums(asv>0)
  asv_prevalence <- colSums(asv>0)
  num_sam <- nrow(asv)
  
  # Table data:
  tibble(
    Dataset = barcode,
    Seq = sum(asv),
    ASVs = ncol(asv),
    N = num_sam,
    Mean_seq = mean(seq_per_sam),
    SD_seq = sd(seq_per_sam),
    Min_seq = min(seq_per_sam),
    Max_seq = max(seq_per_sam),
    Mean_asv = mean(asv_per_sam),
    SD_asv = sd(asv_per_sam),
    Min_asv = min(asv_per_sam),
    Max_asv = max(asv_per_sam),
    Mean_prev = mean(asv_prevalence),
    SD_prev = sd(asv_prevalence),
    Min_prev = min(asv_prevalence),
    Max_prev = max(asv_prevalence)
  )
}) %>% list_rbind

ps.stats %<>%
  mutate(across(where(is.numeric), ~ format(round(., 0),big.mark=',')))

ps.stats.k <- kable(ps.stats, "html", align = "l") %>%
  kable_styling(full_width = FALSE) %>%
  add_header_above(c(
    "Dataset" = 1, # specifies how many columns are covered
    "Sequences" = 1,
    "ASVs" = 1,
    "Samples" = 1,
    "Mean ± SD" = 2, 
    "[Min, Max]" = 2, 
    "Mean ± SD" = 2, 
    "[Min, Max]" = 2, 
    "Mean ± SD" = 2, 
    "[Min, Max]" = 2
  )) %>%  
  add_header_above(c(
    " " = 4, # no header for the first 4 columns
    "Sequences per sample" = 4, 
    "ASVs per sample" = 4, 
    "ASV prevalence" = 4
  )) %>%
  row_spec(0, extra_css = "display: none;") ; ps.stats.k # Hide the original column names

html_file <- "out/asv_processing/asv_summary.html"
save_kable(ps.stats.k, file = html_file)

###################################
# Taxonomic classification rates ###
#####################################

seqtab.ls <- read_rds('data/seqtab.ls.rds')
taxtab.ls <- read_rds('data/taxtab.ls.rds')

# Long df with sequence counts per ASV per sample, with taxonomy
merge_seq_tax <- function(seqtab, taxtab) {
  
  # Long seq table
  seqtab.long <- seqtab %>% 
    as.data.frame %>% 
    rownames_to_column('Sample') %>% 
    pivot_longer(cols = -Sample,
                 names_to = "ASV",
                 values_to = "Abundance") %>% 
    filter(Abundance>0) %>% 
    group_by(Sample) %>% 
    mutate(relAb = Abundance / sum(Abundance))
  
  # Tax table
  taxtab.df <- taxtab %>% 
    as.data.frame %>% 
    rownames_to_column('ASV')
  
  # Merge taxonomy
  left_join(
    seqtab.long,
    taxtab.df, 
    by = 'ASV'
  ) %>% return # a grouped tibble
}

# LOOP over datasets
classification.df <- map(c('BACT', 'FUNG', 'PLAN'), function(barcode){
  
  # Apply merging function:
  dataset <- merge_seq_tax(seqtab.ls[[barcode]], 
                       taxtab.ls[[barcode]])
  
  ranks <- c('Class', 'Order', 'Family', 'Genus', 'Species')
  # LOOP over taxranks
  map(ranks, function(rank) {
    if (!rank %in% colnames(dataset)) {
      return(tibble())  # Return an empty tibble if the variable is absent
    }
    
    dataset %>%
      select(Sample, !!sym(rank), relAb) %>% 
      mutate(classified = case_when(!!sym(rank)=='Unclassified' ~ 0, TRUE ~ 1)) %>% 
      summarise(  
        asv_prop = sum(classified)/n(), # proportion of classified asvs
        relAb_prop = sum(classified*relAb) # abundance-weighted prop of classified asvs
      ) %>%  
      pivot_longer(cols = c('relAb_prop','asv_prop'), 
                   names_to = 'proportion_type') %>% 
      mutate(taxRank = factor(rank, levels = ranks)) # add taxrank variable
  }) %>% list_rbind %>% 
    mutate(barcode = barcode) # add barcode variable
}) %>% list_rbind 

classification.df %>% 
  mutate(proportion_type = case_when(
    proportion_type == 'asv_prop' ~ 'Proportion of assigned ASVs',
    proportion_type == 'relAb_prop' ~ 'Proportion of assigned ASV reads'
    )) %>% 
  ggplot(aes(y = value, x = taxRank, colour = taxRank)) +
  geom_boxplot() +
  ylim(0,NA)+
  facet_grid(barcode~proportion_type) +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) 

ggsave('out/asv_processing/classification_rates.pdf', 
       bg = 'white', width = 1400, height = 2000, 
       units = 'px', dpi = 220)

############################
# Rarefy phyloseq objects ###
##############################

ps_rare.ls <- lapply(ps.ls, function(ps) {
  set.seed(1234)
  prune_samples(sample_sums(ps) >= 2000, ps) %>% 
    rarefy_even_depth2(ncores = 7) 
})

write_rds(ps_rare.ls, 'data/ps_rare.ls.rds',
          compress = 'gz')

###################################
# Most abundant taxa by barcode ###
#####################################


