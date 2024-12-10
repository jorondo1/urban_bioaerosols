library(pacman)
p_load(tidyverse, phyloseq, magrittr)

source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/tax_glom2.R")

urbanbio.path <- '~/Desktop/ip34/urbanBio'
meta <- read_delim(file.path(urbanbio.path,"data/metadata_2022_samples_final.csv"))
#meta %<>% mutate(sample_id = paste0("X", sample_id))
meta_controls <- meta %>% 
  filter(time =="None") %>% 
  mutate(sample_id = case_when(sample_id == "Contrôle-blanc-12h00-PM" ~ "Controle-blanc-12h00-PM",
                                      TRUE ~ sample_id))
meta_samples <- meta %>% 
  filter(time != "None")
  
sample.names <- meta_samples$sample_id


###############
# Functions ####
#################

# Keep samples with metadata
subset_samples <- function(seqtab, samples) {
  seqtab[rownames(seqtab) %in% samples, ] %>% # subset
    .[, colSums(.) > 0] # Remove ASVs with no hits
}

# ASVs classified at the kingdom level and present in seqtab
subset_asvs <- function(taxonomy, seqtab) {
  if (!is.data.frame(taxonomy)) {
    taxonomy <- as.data.frame(taxonomy)
  }
  asvs <- subset(taxonomy, Kingdom != "Unclassified") %>% # subset needs the input to be a df
    rownames %>% 
    intersect(colnames(seqtab)) # only keep asvs still present in seqtab
  taxonomy[asvs, ] %>% as.matrix()
}

# Remove samples with fewer than n sequences once taxa removed
remove_ultra_rare <- function(seqtab, taxonomy, n) {
  seqtab[,rownames(taxonomy)] %>% .[rowSums(.) > 10, ]
}

#########
# 16S ####
###########
path_16S <- file.path(urbanbio.path,'data/16S')
taxa_16S_genus <- read_rds(file.path(path_16S, '4_taxonomy/taxonomy.RDS'))
taxa_16S_species <- read_rds(file.path(path_16S, '4_taxonomy/taxonomy_species.RDS'))
seqtab_16S <- read_rds(file.path(path_16S, '4_taxonomy/seqtab.RDS'))

# Somehow assignSpecies drops all ranks except genus... fix 
Species_16S <- taxa_16S_species[,2]
names(Species_16S) <- rownames(taxa_16S_species)
taxa_16S <- cbind(taxa_16S_genus, Species_16S) %>% data.frame

# Keep samples with metadata info
seqtab_16S_samples <- subset_samples(seqtab_16S, sample.names)

# Subset ASVs
taxa_16S_present <- subset_asvs(taxa_16S, seqtab_16S_samples) 

# Remove near-empty samples
seqtab_16S_present <- remove_ultra_rare(seqtab_16S_samples, taxa_16S_present, 10)
dim(seqtab_16S_samples); dim(seqtab_16S_present); dim(taxa_16S_present)

# Phyloseq object
ps_16S <- phyloseq(
  tax_table(taxa_16S_present),
  otu_table(seqtab_16S_present, taxa_are_rows = FALSE),
  sample_data(meta %>% column_to_rownames('sample_id'))
)

#########
# ITS ####
###########
path_ITS <- file.path(urbanbio.path,'data/ITS')
taxa_ITS <- read_rds(file.path(path_ITS, '4_taxonomy/taxonomy.RDS'))
seqtab_ITS <- read_rds(file.path(path_ITS, '4_taxonomy/seqtab.RDS'))

# Keep samples with metadata info
seqtab_ITS_samples <- subset_samples(seqtab_ITS, sample.names)

# Subset ASVs
taxa_ITS_present <- subset_asvs(taxa_ITS, seqtab_ITS_samples)

# Remove near-empty samples
seqtab_ITS_present <- remove_ultra_rare(seqtab_ITS_samples, taxa_ITS_present, 10)
dim(seqtab_ITS_samples); dim(seqtab_ITS_present); dim(taxa_ITS_present)

# Phyloseq object
ps_ITS <- phyloseq(
  tax_table(taxa_ITS_present),
  otu_table(seqtab_ITS_present, taxa_are_rows = FALSE),
  sample_data(meta %>% column_to_rownames('sample_id'))
  )

##########
# trnL ####
############
path_trnL <- file.path(urbanbio.path,'data/trnL')
taxa_trnL <- read_rds(file.path(path_trnL, '4_taxonomy/taxonomy.RDS'))
seqtab_trnL <- read_rds(file.path(path_trnL, '4_taxonomy/seqtab.RDS'))

# Keep samples with metadata info
seqtab_trnL_samples <- subset_samples(seqtab_trnL, sample.names)

# Subset ASVs
taxa_trnL_present <- subset_asvs(taxa_trnL, seqtab_trnL_samples)

# Remove near-empty samples
seqtab_trnL_present <- remove_ultra_rare(seqtab_trnL_samples, taxa_trnL_present, 10)
dim(seqtab_trnL_samples); dim(seqtab_trnL_present); dim(taxa_trnL_present)

# Phyloseq object
ps_trnL <- phyloseq(
  tax_table(taxa_trnL_present),
  otu_table(seqtab_trnL_present, taxa_are_rows = FALSE),
  sample_data(meta %>% column_to_rownames('sample_id'))
)
ps.list <- list()
ps.list[["16S"]] <- ps_16S
ps.list[["ITS"]] <- ps_ITS
ps.list[["trnL"]] <- ps_trnL

saveRDS(ps.list, file.path(urbanbio.path,'data/ps.list1.rds'))
