# Richness by taxonomic level

#############
# X. SETUP ###
###############
library(pacman)
p_load(tidyverse, phyloseq)
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/psflashmelt.R")
source('scripts/0_config.R') # Variable naming and such

ps_rare.ls <- read_rds('data/ps_rare.ls.rds')

# Melt all three datasets and subset for montreal
psmelt.ls <- map(ps_rare.ls, function(ps){
  ps %>% 
    psflashmelt() %>% 
    filter(city == 'Montreal')
})

# Richness by taxrank
for (taxRank in c('ASV', 'Genus', 'Family')) {
  glom_rank <- if(taxRank == 'ASV') {'OTU'} else {taxRank}
  
  richness.table <- imap(psmelt.ls, function(ps,barcode ){
    
    raw_table <- ps %>% 
      filter(Abundance != 0) %>% 
      select(all_of(glom_rank), Sample, Abundance, site_id, time)
    
    raw_table %>% 
      group_by(site_id, !!sym(glom_rank)) %>% 
      summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
      group_by(site_id) %>% 
      summarise(Richness = n()) %>% 
      mutate(Barcode = barcode)
    
  }) %>% list_rbind() %>% 
    mutate(Barcode = recode(Barcode, !!!kingdoms)) %>% 
    pivot_wider(names_from = Barcode, values_from = Richness) 
  
  write_tsv(richness.table, paste0('out/stats/',taxRank,'_richness_by_site.tsv'))
  
  # Taxonomic assignment rates 
  imap(psmelt.ls, function(ps,barcode ){
    
    classification_plot.df <- ps %>% 
      select(OTU, Sample, site_id, Abundance, !!sym(glom_rank)) %>% 
      filter(Abundance > 0) %>% 
      mutate(classified = case_when(!!sym(glom_rank) != 'Unclassified' ~ 'Classified',
                                    TRUE ~ 'Unclassified')) %>% 
      group_by(site_id, classified) %>% 
      summarise(Abundance = sum(Abundance), .groups = 'drop') %>% 
      group_by(site_id) %>% 
      mutate(relAb = Abundance/sum(Abundance)) 
    
    unclassified <- classification_plot.df %>% 
      group_by(classified) %>% 
      summarise(mean_relab = mean(relAb)) 
    
    message(paste0('Proportion unclassified for ', barcode, ' at ', glom_rank, 'level :', unclassified[2,2]))
    
    class_plot <- classification_plot.df %>% 
      ggplot(aes(x = site_id, y = relAb, fill = classified)) +
      geom_col() +
      labs(y = "Relative abundance of amplicons",
           x = 'Site ID',
           fill = paste('Classification of ',barcode, 'at the', glom_rank, 'level')) +
      theme(legend.position = 'bottom') 
    
    ggsave(plot = class_plot, paste0('Out/alainPaquette/', barcode, '_classification_',glom_rank,'.pdf'), 
           bg = 'white', width = 2000, height = 1500, 
           units = 'px', dpi = 200)
    
    
  })
}
