# Defining the core microbiome

#############
# X. SETUP ###
###############
library(pacman)
p_load(tidyverse, phyloseq)

source('scripts/0_config.R') # Variable naming and such
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/community_functions.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/psflashmelt.R')

ps.ls <- read_rds('data/ps.ls.rds')

# Finding core ASVs and returning a table with taxonomy and abundance/prevalence stats
list_core_ASV <- function(ps, minPrev, minAbund) {  
  
  if(taxa_are_rows(ps)) {Margin <- c(1,2)} else {Margin <- c(2,1) }
  
  # Compute prevalence 
  prev <- apply(otu_table(ps), Margin[1], function(x) sum(x > 0)) / nsamples(ps)
  prevalent_ASVs <- names(prev[prev>minPrev])
  
  # Convert to relative abundance
  rel_abund <- apply(otu_table(ps), Margin[2], function(x) x / sum(x))
  
  # Subset to minPrev ASVs 
  rel_abund_prev <- rel_abund[prevalent_ASVs, , drop = FALSE]
  
  # Mean of non-zero relAb
  ASV_relab_means <- apply(rel_abund_prev, 1, function(row) {
    non_zero <- row[row != 0]
    if (length(non_zero) == 0) NA else mean(non_zero)
  })
  
  core_ASVs <- names(which(ASV_relab_means > minAbund))
  
  core_prev <- prev[core_ASVs]
  core_prev <- as.data.frame(core_prev) %>% 
    rownames_to_column('OTU')

  # Melt ps object  
  melt <- psflashmelt(ps)
  
  # If species is available, we want it
  tax_ranks_to_include <- intersect(colnames(melt), c(taxRanks[1:6], 'Species'))
  
  # Extract relative abundance from complete dataset
  melt %>% 
  dplyr::select(Sample, Abundance, all_of(tax_ranks_to_include), OTU) %>% # exclude Species
    group_by(Sample) %>% 
    mutate(relAb = Abundance/sum(Abundance), .keep = 'unused') %>% 
    ungroup() %>% 
    filter(OTU %in% core_ASVs) %>% 
    group_by(OTU, !!!syms(tax_ranks_to_include)) %>% 
    summarise(mean_relAb = mean(relAb),
              sd_relAb = sd(relAb), .groups = 'drop') %>% 
    left_join(core_prev, by = 'OTU') %>% 
    select(-OTU, -Kingdom)
}

# Execute on all and save using kable
core_ASVs.ls <- imap(ps.ls, function(ps, barcode) {
  core_table <- list_core_ASV(ps, minPrev = 0.95, minAbund = 0.01)
  core_table %>% 
    mutate(across(where(is.numeric), ~ round(., 3))) %>% 
    kable("html",
          caption = '') %>%
    kable_styling(full_width = FALSE) %>% 
    save_kable(file = paste0('out/core_microbiome/core_',barcode,'.html'))
})
