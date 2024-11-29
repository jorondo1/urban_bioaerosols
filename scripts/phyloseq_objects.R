library(pacman)
p_load(tidyverse, phyloseq, magrittr)
meta <- read_delim("data/metadata_2022_samples_final.csv")
meta %<>% mutate(sample_id = case_when(sample_id == "Contrôle-blanc-12h00-PM" ~ "Controle-blanc-12h00-PM",
                                      TRUE ~ sample_id))


# 16S #########################
# Somehow assignSpecies drops all ranks excetp genus... fix 
taxa <- read_rds("data/16S/4_taxonomy/taxa.RDS")
taxa_species <- read_rds("data/16S/4_taxonomy/taxa_species.RDS")
Species <- taxa_species[,2]
names(Species) <- rownames(taxa_species)
taxa_complete <- cbind(taxa, Species)

# Subset only samples for which we have metadata
seqtab.nochim <- read_rds("data/16S/4_taxonomy/seqtab.nochim.RDS")
seqtab.samples <- seqtab.nochim[meta$sample_id,]

phyloseq(tax_table(taxa_complete),
           otu_table(seqtab.samples, taxa_are_rows = FALSE),
           sample_data(meta %>% column_to_rownames('sample_id')))
