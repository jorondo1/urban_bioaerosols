library(pacman)
p_load(tidyverse, phyloseq, magrittr)

source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/tax_glom2.R")

meta <- read_delim("data/metadata_2022_samples_final.csv")
#meta %<>% mutate(sample_id = paste0("X", sample_id))
meta_controls <- meta %>% 
  filter(time =="None") %>% 
  mutate(sample_id = case_when(sample_id == "Contrôle-blanc-12h00-PM" ~ "Controle-blanc-12h00-PM",
                                      TRUE ~ sample_id))

meta_samples <- meta %>% 
  filter(time != "None")
  

#########
# 16S ####
###########
# Somehow assignSpecies drops all ranks excetp genus... fix 
taxa_16S <- read_rds("data/16S/4_taxonomy/taxa.RDS")
taxa_16S_species <- read_rds("data/16S/4_taxonomy/taxa_species.RDS")
Species_16S <- taxa_16S_species[,2]
names(Species_16S) <- rownames(taxa_16S_species)
taxa_complete_16S <- cbind(taxa, Species_16S)

# Subset only samples for which we have metadata
seqtab_16S <- read_rds("data/16S/4_taxonomy/seqtab.nochim.RDS")
seqtab_16S_samples <- seqtab_16S[meta_samples$sample_id,] # subset to samples with metadata

# Phyloseq object
ps_16S <- phyloseq(tax_table(taxa_complete_16S),
           otu_table(seqtab_16S_samples, taxa_are_rows = FALSE),
           sample_data(meta %>% column_to_rownames('sample_id')))
ps_16S_genus <- tax_glom2(ps_16S, 'Genus')

write_rds(ps_16S, "data/16S/5_out/ps.RDS")

##########
# trnL ####
############
taxa_trnL <- read_rds("data/trnL/4_taxonomy_E22_100/taxonomy.RDS")
seqtab_trnL <- read_rds("data/trnL/4_taxonomy_E22_100/seqtab.RDS")
 #rownames(seqtab_trnL) <- sub("-[^-]*$", "", rownames(seqtab_trnL))

# Keep samples with metadata info
trnL_samples <- intersect(rownames(seqtab_trnL),meta_samples$sample_id)
seqtab_trnL_samples <- seqtab_trnL[trnL_samples,]
dim(seqtab_trnL_samples)

# Remove ASVs with no hits
seqtab_trnL_samples <- seqtab_trnL_samples[, colSums(seqtab_trnL_samples) > 0]
dim(seqtab_trnL_samples)

# And ASVs with no taxonomy at all 
unclassified_asvs <- taxa_trnL[taxa_trnL[, "Kingdom"] != "Unclassified", ]
asv_subset <- intersect(rownames(unclassified_asvs), colnames(seqtab_trnL_samples))
taxa_trnL_present <- taxa_trnL[asv_subset,]
seqtab_trnL_present <- seqtab_trnL_samples[,asv_subset]

# Remove samples with fewer than 5 sequences
seqtab_trnL_final <- seqtab_trnL_present[rowSums(seqtab_trnL_present) > 5,]

dim(taxa_trnL_present); dim(seqtab_trnL_final)

# Phyloseq object
ps_trnL <- phyloseq(tax_table(taxa_trnL_present),
         otu_table(seqtab_trnL_present, taxa_are_rows = FALSE),
         sample_data(meta %>% column_to_rownames('sample_id')))
ps_trnL_genus <- tax_glom(ps_trnL, taxrank = 'Genus')

write_rds(ps_trnL, "data/trnL/5_out/ps.RDS")
write_rds(ps_trnL_genus, "data/trnL/5_out/ps_genus.RDS")
