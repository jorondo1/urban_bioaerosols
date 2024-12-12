library(pacman)
p_load(tidyverse, phyloseq, magrittr, decontam)
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/tax_glom2.R")

###############
# Functions ####
#################

# Keep samples with metadata
subset_samples <- function(seqtab, samples) {
  seqtab[rownames(seqtab) %in% samples, ] %>% # subset
    .[, colSums(.) > 0] # Remove ASVs with no hits
}

# ASVs classified at the kingdom level and present in seqtab
subset_asvs <- function(taxonomy, seqtab, min_seq) {
  if (!is.data.frame(taxonomy)) {
    taxonomy <- as.data.frame(taxonomy)
  }
  
  asvs <- subset(taxonomy, Kingdom != "Unclassified") %>% # subset needs the input to be a df
    rownames %>% 
    intersect(
      colnames(seqtab)[colSums(seqtab) >= min_seq]
    ) # only keep asvs still present in seqtab
  taxonomy[asvs, ] %>% as.matrix()
}

# Remove samples with fewer than n sequences once taxa removed
remove_ultra_rare <- function(seqtab, taxonomy, n) {
  seqtab[,rownames(taxonomy)] %>% .[rowSums(.) > 10, ]
}

# Parse DNA concentration xlsx files
parse_CERMO_xlsx <- function(files) {
  require(readxl, dplyr)
  files %>% 
    map(read_xlsx, range = "C13:D120", # CERMO-files specific
        col_names = c('sample_id', 'concDNA')) %>% 
    list_rbind %>% 
    filter(!is.na(sample_id) 
           & !concDNA %in% c('n/a', 'Too Low')) %>% 
    mutate(sample_id = sub("_[^_]*$", "", sample_id))
}

# Add sequencing depth and dna concentration to sample metadata
add_seq_depth <- function(seqtab, meta, dnaConc) {
  require(tibble)
  meta_subset <- meta %>% 
    left_join(dnaConc, by = 'sample_id') %>% # adds seq depth column
    column_to_rownames('sample_id') %>% 
    .[rownames(seqtab),]
  
  seqtab %>% rowSums %>% 
    data.frame(seqDepth = .) %>% 
    cbind(meta_subset, .) 
}

##########
# Setup ###
############

urbanbio.path <- '~/Desktop/ip34/urbanBio'
dna.path <- file.path(urbanbio.path,"data/metadata")

# Metadata
meta <- read_delim(file.path(urbanbio.path,"data/metadata/metadata_2022_samples_final.csv"))

meta_ctrl <- meta %>% 
  filter(time =="None" | control=='TRUE') %>% 
  mutate(sample_id = case_when(
    sample_id == "Contrôle-blanc-12h00-PM" ~ "Controle-blanc-12h00-PM",
    TRUE ~ sample_id))

meta_samples <- meta %>% 
  filter(time != "None" & control == 'FALSE')
  
sample.names <- meta_samples$sample_id
ctrl.names <- meta_ctrl$sample_id

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
seqtab_16S_sam <- subset_samples(seqtab_16S, sample.names)
seqtab_16S_ctrl <- subset_samples(seqtab_16S, ctrl.names)

# Subset ASVs
taxa_16S_sam <- subset_asvs(taxa_16S, seqtab_16S_sam, 100) 
taxa_16S_ctrl <- subset_asvs(taxa_16S, seqtab_16S_ctrl, 100) 

# Remove near-empty samples
seqtab_16S_sam_filt <- remove_ultra_rare(seqtab_16S_sam, taxa_16S_sam, 10)
seqtab_16S_ctrl_filt <- remove_ultra_rare(seqtab_16S_ctrl, taxa_16S_ctrl, 10)
dim(seqtab_16S_sam); dim(seqtab_16S_sam_filt); dim(taxa_16S_sam)
dim(seqtab_16S_ctrl); dim(seqtab_16S_ctrl_filt); dim(taxa_16S_ctrl)

# Add sequencing effort and dna concentration to metadata
dna_16S <- Sys.glob(file.path(dna.path,'CERMO_*16s*.xlsx')) %>% parse_CERMO_xlsx()
meta_samples_16S <- add_seq_depth(seqtab_16S_sam_filt, meta_samples, dna_16S)
meta_ctrl_16S <- add_seq_depth(seqtab_16S_ctrl_filt, meta_controls, dna_16S)

# Phyloseq object
ps_16S <- phyloseq(
  tax_table(taxa_16S_sam),
  otu_table(seqtab_16S_sam_filt, taxa_are_rows = FALSE),
  sample_data(meta_samples_16S)
)

ps_16S_ctrl <- phyloseq(
  tax_table(taxa_16S_ctrl),
  otu_table(seqtab_16S_ctrl_filt, taxa_are_rows = FALSE),
  sample_data(meta_ctrl_16S)
)

#########
# ITS ####
###########

path_ITS <- file.path(urbanbio.path,'data/ITS')
taxa_ITS <- read_rds(file.path(path_ITS, '4_taxonomy/taxonomy.RDS'))
seqtab_ITS <- read_rds(file.path(path_ITS, '4_taxonomy/seqtab.RDS'))

# Keep samples with metadata info
seqtab_ITS_sam <- subset_samples(seqtab_ITS, sample.names)
seqtab_ITS_ctrl <- subset_samples(seqtab_ITS, ctrl.names)

# Subset ASVs
taxa_ITS_sam <- subset_asvs(taxa_ITS, seqtab_ITS_sam, 100)
taxa_ITS_ctrl <- subset_asvs(taxa_ITS, seqtab_ITS_ctrl, 100)

# Remove near-empty samples
seqtab_ITS_sam_filt <- remove_ultra_rare(seqtab_ITS_sam, taxa_ITS_sam, 10)
seqtab_ITS_ctrl_filt <- remove_ultra_rare(seqtab_ITS_ctrl, taxa_ITS_ctrl, 10)

dim(seqtab_ITS_sam); dim(seqtab_ITS_sam_filt); dim(taxa_ITS_sam)
dim(seqtab_ITS_ctrl); dim(seqtab_ITS_ctrl_filt); dim(taxa_ITS_ctrl)

# Add sequencing effort and dna concentration to metadata
dna_ITS <- Sys.glob(file.path(dna.path,'CERMO_*ITS*.xlsx')) %>% parse_CERMO_xlsx
meta_samples_ITS <- add_seq_depth(seqtab_ITS_sam_filt, meta_samples, dna_ITS)
meta_ctrl_ITS <- add_seq_depth(seqtab_ITS_ctrl_filt, meta_controls, dna_ITS)

# Phyloseq objects
ps_ITS <- phyloseq(
  tax_table(taxa_ITS_sam),
  otu_table(seqtab_ITS_sam_filt, taxa_are_rows = FALSE),
  sample_data(meta_samples_ITS)
  )

ps_ITS_ctrl <- phyloseq(
  tax_table(taxa_ITS_ctrl),
  otu_table(seqtab_ITS_ctrl_filt, taxa_are_rows = FALSE),
  sample_data(meta_ctrl_ITS)
)

##########
# trnL ####
############

path_trnL <- file.path(urbanbio.path,'data/trnL')
taxa_trnL <- read_rds(file.path(path_trnL, '4_taxonomy/taxonomy.RDS'))
seqtab_trnL <- read_rds(file.path(path_trnL, '4_taxonomy/seqtab.RDS'))

# Keep samples with metadata info
seqtab_trnL_sam <- subset_samples(seqtab_trnL, sample.names)
seqtab_trnL_ctrl <- subset_samples(seqtab_trnL, ctrl.names)

# Subset ASVs
taxa_trnL_sam <- subset_asvs(taxa_trnL, seqtab_trnL_sam, 100)
taxa_trnL_ctrl <- subset_asvs(taxa_trnL, seqtab_trnL_ctrl, 100)

# Remove near-empty samples
seqtab_trnL_sam_filt <- remove_ultra_rare(seqtab_trnL_sam, taxa_trnL_sam, 10)
seqtab_trnL_ctrl_filt <- remove_ultra_rare(seqtab_trnL_ctrl, taxa_trnL_ctrl, 10)

dim(seqtab_trnL_sam); dim(seqtab_trnL_sam_filt); dim(taxa_trnL_sam)
dim(seqtab_trnL_ctrl); dim(seqtab_trnL_ctrl_filt); dim(taxa_trnL_ctrl)

# Add sequencing effort and dna concentration to metadata
dna_trnL <- Sys.glob(file.path(dna.path,'CERMO_*trnL*.xlsx')) %>% parse_CERMO_xlsx
meta_samples_trnL <- add_seq_depth(seqtab_trnL_sam_filt, meta_samples, dna_trnL)
meta_ctrl_trnL <- add_seq_depth(seqtab_trnL_ctrl_filt, meta_controls, dna_trnL)

# Phyloseq objects
ps_trnL <- phyloseq(
  tax_table(taxa_trnL_sam),
  otu_table(seqtab_trnL_sam_filt, taxa_are_rows = FALSE),
  sample_data(meta_samples_trnL)
)

ps_trnL_ctrl <- phyloseq(
  tax_table(taxa_trnL_ctrl),
  otu_table(seqtab_trnL_ctrl_filt, taxa_are_rows = FALSE),
  sample_data(meta_ctrl_trnL)
)

ps.ls <- list()
ps.ls[["BACT"]] <- ps_16S
ps.ls[["FUNG"]] <- ps_ITS
ps.ls[["PLAN"]] <- ps_trnL
saveRDS(ps.ls, file.path(urbanbio.path,'data/ps.ls.rds'))

ps_ctrl.ls <- list()
ps_ctrl.ls[["BACT"]] <- ps_16S_ctrl
ps_ctrl.ls[["FUNG"]] <- ps_ITS_ctrl
ps_ctrl.ls[["PLAN"]] <- ps_trnL_ctrl
saveRDS(ps_ctrl.ls, file.path(urbanbio.path,'data/ps_ctrl.ls.rds'))

# //DEV
# Decontamination DECONTAM
# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
ps_16S_load <- ps_16S %>%
  prune_samples(sample_data(.)$bacterial_load>0, .)
contam_freq<- ps_16S_load %>% 
  isContaminant(method = 'frequency', conc = 'bacterial_load')

table(contam_freq$contaminant)
head(which(contam_freq$contaminant))

ps_16S_load %>% 
  plot_frequency(., taxa_names(.)[c(516,1128)], conc="bacterial_load") + 
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

