# Compute taxon prevalence, top 5 genera
library(pacman)
p_load(tidyverse, magrittr)
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/psflashmelt.R")

urbanbio.path <- '~/Desktop/ip34/urbanBio'
ps.ls <- readRDS(file.path(urbanbio.path, 'data/ps.ls.rds'))

melted  <- lapply(ps.ls, psflashmelt) %>% 
  list_rbind

# Genera of most prevalent ASVs

kingdom_names <- c('Bacteria' = 'Bacteria',
                   'Fungi' = 'Fungi',
                   'Pollen' = 'Viridiplantae')

imap(kingdom_names, function(k, barcode){
  
  # Per kingdom
  melted_k <- melted %>%
    filter(Kingdom == k)
  
  # Prevalence table
  prevalence_table <- melted_k %>% 
    # Compute prevalence
    mutate(Barcode = barcode) %>% 
    group_by(OTU, Barcode, Family, Genus, city) %>% 
    summarise(Prev = n(), .groups = 'drop') %>%
    # Sort by prevalence
    arrange(desc(Prev)) %>% 
    select(Barcode, Family, Genus, Prev, OTU, city) %>% 
    group_by(Barcode, city) %>% 
    # Create a unique ASV label
    mutate(ASV_rank = row_number()) %>% 
    rename(ASV = OTU) %>% 
    ungroup %>% 
    select(-Barcode)
  
  prev_SHB <- prevalence_table %>% filter(city == 'Sherbrooke') %>% select(-city)
  prev_MTL <- prevalence_table %>% filter(city == 'Montreal')%>% select(-city)
  prev_QBC <- prevalence_table %>% filter(city == 'Quebec')%>% select(-city)
  
  # We need to do a triple full_join! 
  # List the common keys:
  common_keys <- c('ASV', 'Genus', 'Family')
  
  # Create a list of tables to join and rename all columns with a custom suffix
  # except for the common keys
  tables <- list(
    prevalence_table %>% 
      filter(city == 'Montreal') %>% select(-city) %>% 
      rename_with(~ paste0(.x, '.MTL'), -all_of(common_keys)),
    
    prevalence_table %>% 
      filter(city == 'Quebec') %>% select(-city) %>% 
      rename_with(~ paste0(.x, '.QBC'), -all_of(common_keys)),
    
    prevalence_table %>% 
      filter(city == 'Sherbrooke') %>% select(-city) %>% 
      rename_with(~ paste0(.x, '.SHB'), -all_of(common_keys))
  )
  
  # Reduce the list using the full_join function
  joined_table <- reduce(tables, ~ full_join(.x, .y, by = common_keys)) %>% 
    mutate(ASV_label = paste0("ASV_", barcode, row_number()))
  

  ## sample count:
  nsam_SHB <- melted_k %>% filter(city == 'Sherbrooke') %>% 
    pull(Sample) %>% unique %>% length
  
  nsam_MTL <- melted_k %>% filter(city == 'Montreal') %>% 
    pull(Sample) %>% unique %>% length
  
  nsam_QBC <- melted_k %>% filter(city == 'Quebec') %>% 
    pull(Sample) %>% unique %>% length
  
  # Format table
  top_headers_names <- setNames(
    c(2, 2, 2, 2, 1),
    c(" ",
      paste0("Montreal (n=", nsam_MTL, ")"),
      paste0("Quebec (n=", nsam_QBC, ")"),
      paste0("Sherbrooke (n=", nsam_SHB, ")" ),
      " ")
  )

  # Table :
  joined_table %>%
    mutate(compound_score = ASV_rank.MTL + ASV_rank.SHB + ASV_rank.QBC,
           Prev.MTL = paste0(Prev.MTL, " (",round(100*Prev.MTL/nsam_MTL,0), " %)"),
           Prev.QBC = paste0(Prev.QBC, " (",round(100*Prev.QBC/nsam_QBC,0), " %)"),
           Prev.SHB = paste0(Prev.SHB, " (",round(100*Prev.SHB/nsam_SHB,0), " %)")) %>%
    arrange(compound_score) %>%
    select(-compound_score, -ASV) %>% 
     knitr::kable(col.names = c(
       'Family',
       'Genus',
       'Prev',
       'Rank',
       'Prev',
       'Rank',
       'Prev',
       'Rank',
       'ASV label'
       ), align = 'lccccr') %>%
  kableExtra::add_header_above(top_headers_names) %>% # Top header
  kableExtra::kable_styling(bootstrap_options = c('striped', 'hover'),
                            full_width = FALSE) %>%
  kableExtra::save_kable(paste0(urbanbio.path, "/out/ASV_prev_", barcode,".html"),
             self_contained = TRUE)
})




# Taxonomy table
joined_table %>% 
  select(ASV_label, ASV) %>% 
  left_join(melted %>% 
              select(OTU, Genus, Family, Order, Class, Phylum) %>% 
              unique,
            join_by(ASV == OTU)) %>% 
  knitr::kable() %>% 
  kableExtra::kable_styling('striped') %>% 
  save_kable("2023/out/ASV_labels.html",self_contained = TRUE)


