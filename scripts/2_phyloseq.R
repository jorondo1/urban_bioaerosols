library(pacman)
p_load(tidyverse, phyloseq, magrittr, decontam, Biostrings,
       readxl, decontam)

###############
# Functions ####
#################

# Tax glom modified
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/tax_glom2.R")

# Phyloseq prep functions
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/phyloseq_functions.R")

# Parse DNA concentration xlsx files
parse_CERMO_xlsx <- function(files) {
  require(readxl, dplyr)
  files %>% 
    map(read_xlsx, range = "C13:D120", # CERMO-files specific
        col_names = c('sample_id', 'concDNA')) %>% 
    list_rbind %>% 
    filter(!is.na(sample_id) 
           & !concDNA %in% c('n/a', 'Too Low')) %>% 
    mutate(sample_id = sub("_[^_]*$", "", sample_id),
           concDNA = as.numeric(concDNA))
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

### Decontamination DECONTAM
# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

# 1. Find contaminants using the frequency method (based on DNA concentration)
decontaminate <- function(seqtab, samdata, var, p = 0.1) {
  require(decontam)
  
  idx <- which(samdata[,var] >0)
  conc <- samdata[idx, var]
  asv <- seqtab[idx,]
  
  out <- list()
  out$decontam <- isContaminant(seqtab = asv, method = 'frequency', conc = conc, threshold = p)
  how_many <- table(out$decontam$contaminant) # How many ASVs are contaminants ?
  asv_rank <- which(out$decontam$contaminant) # What is the abundance ranks of contaminants
  
  # Short report on decontam's findings; produces a plot
  out$p <- plot_frequency(asv, colnames(asv)[head(which(out$decontam$contaminant), n = 20)], conc=conc) + 
    xlab("DNA Concentration")
  
  message(paste(how_many[2], 'ASVs identified as contaminants, out of', sum(how_many)))
  message(paste('Top abundance ranks of contaminants: ', paste(head(asv_rank),collapse = ', '))) 
        # bigger numbers = lowest abundance ranks
  
  return(out)
}

# 2. Remove contaminants from ps object
prune_contam <- function(ps, decontam_table) {
  require(dplyr, phyloseq)
  
  # Save non-contaminant ASVs
  not_cont <- which(!decontam_table$contaminant) %>% 
    decontam_table[.,] %>% rownames
  
  #Remove them
  ps_clean <- prune_taxa(not_cont, ps)
  
  # Proportion of counts lost
  message(paste(round(100*(1-sum(ps_clean@otu_table)/sum(ps@otu_table)),2), 
                '% of reads were lost to decontamination.'))
  
  return(ps_clean)
}

##########
# Setup ###
############

urbanbio.path <- '~/Desktop/ip34/urbanBio'
dna.path <- file.path(urbanbio.path,"data/metadata")

# Metadata
meta <- read_delim(file.path(urbanbio.path,"data/metadata/metadata_2022_samples_final.csv")) 
meteo <- Sys.glob(file.path(urbanbio.path,"data/metadata/meteo_*.csv")) %>% 
  map_dfr(read_delim, locale = locale(decimal_mark = ","),show_col_types = FALSE) %>% 
  mutate(city = case_when(
    `Nom de la Station` == 'BEAUPORT' ~ 'Québec',
    `Nom de la Station` == 'SHERBROOKE' ~ 'Sherbrooke',
    `Nom de la Station` == 'MCTAVISH' ~ 'Montréal'
  )) %>% 
  transmute(city = city,
            temp_moy = `Temp moy.(°C)`,
            precip = `Précip. tot. (mm)`,
            date = `Date/Heure`)
meta %<>%
  mutate(across(where(is.character), as.factor)) %>% 
  select(-bacterial_load, -log_bacterial_load, fungal_load, log_fungal_load) %>% 
  left_join(meteo, by = c('date', 'city')) # Add precipitations & temperature data

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
taxa_16S_sam <- subset_asvs(taxa_16S, seqtab_16S_sam, 10) 
taxa_16S_ctrl <- subset_asvs(taxa_16S, seqtab_16S_ctrl, 10) 

seqtab_16S_sam %>% rowSums %>% sort %>% plot

## Check distribution of sample depth
hist(rowSums(seqtab_16S_sam), breaks = 100, xlab = "sample size", xaxt = "n")
axis(1, at = pretty(rowSums(seqtab_16S_sam), n = 40))  # adding ticks 

# Remove near-empty samples
seqtab_16S_sam_filt <- remove_ultra_rare(seqtab_16S_sam, taxa_16S_sam, 10000)
seqtab_16S_ctrl_filt <- remove_ultra_rare(seqtab_16S_ctrl, taxa_16S_ctrl, 10000)
dim(seqtab_16S_sam); dim(seqtab_16S_sam_filt); dim(taxa_16S_sam)
dim(seqtab_16S_ctrl); dim(seqtab_16S_ctrl_filt); dim(taxa_16S_ctrl)

# Add sequencing effort and dna concentration to metadata
dna_16S <- Sys.glob(file.path(dna.path,'CERMO_*16s*.xlsx')) %>% parse_CERMO_xlsx()
meta_samples_16S <- add_seq_depth(seqtab_16S_sam_filt, meta_samples, dna_16S)
meta_ctrl_16S <- add_seq_depth(seqtab_16S_ctrl_filt, meta_ctrl, dna_16S)

# Find contaminants
contam_freq_16S <- decontaminate(seqtab_16S_sam_filt, meta_samples_16S, 'concDNA')
contam_freq_16S$p # look at contaminants correlation with concDNA

# Phyloseq object
ps_16S <- phyloseq(
  tax_table(taxa_16S_sam),
  otu_table(seqtab_16S_sam_filt, taxa_are_rows = FALSE),
  sample_data(meta_samples_16S)
) %>% prune_contam(contam_freq_16S$decontam) # decontam

ps_16S_ctrl <- phyloseq(
  tax_table(taxa_16S_ctrl),
  otu_table(seqtab_16S_ctrl_filt, taxa_are_rows = FALSE),
  sample_data(meta_ctrl_16S)
)

# Export asvs as fasta
asv_to_fasta(seqtab_16S_sam_filt, file.path(path_16S, '4_taxonomy/asv.fa'))

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
taxa_ITS_sam <- subset_asvs(taxa_ITS, seqtab_ITS_sam, 10)
taxa_ITS_ctrl <- subset_asvs(taxa_ITS, seqtab_ITS_ctrl, 10)

seqtab_ITS_sam %>% rowSums %>% sort %>% data.frame(Y=.) %>% 
  ggplot(aes(x = Y)) + geom_histogram(bins=100)

## Check distribution of sample depth
hist(rowSums(seqtab_ITS_sam), breaks = 100, xlab = "sample size", xaxt = "n")
axis(1, at = pretty(rowSums(seqtab_ITS_sam), n = 40))  # adding ticks 

# Remove near-empty samples
seqtab_ITS_sam_filt <- remove_ultra_rare(seqtab_ITS_sam, taxa_ITS_sam, 1200)
seqtab_ITS_ctrl_filt <- remove_ultra_rare(seqtab_ITS_ctrl, taxa_ITS_ctrl, 1200)

dim(seqtab_ITS_sam); dim(seqtab_ITS_sam_filt); dim(taxa_ITS_sam)
dim(seqtab_ITS_ctrl); dim(seqtab_ITS_ctrl_filt); dim(taxa_ITS_ctrl)

# Add sequencing effort and dna concentration to metadata
dna_ITS <- Sys.glob(file.path(dna.path,'CERMO_*ITS*.xlsx')) %>% parse_CERMO_xlsx()
meta_samples_ITS <- add_seq_depth(seqtab_ITS_sam_filt, meta_samples, dna_ITS)
meta_ctrl_ITS <- add_seq_depth(seqtab_ITS_ctrl_filt, meta_ctrl, dna_ITS)

# Find contaminants
contam_freq_ITS <- decontaminate(seqtab_ITS_sam_filt, meta_samples_ITS, 'concDNA')
contam_freq_ITS$p # look at contaminants correlation with concDNA

# Phyloseq objects
ps_ITS <- phyloseq(
  tax_table(taxa_ITS_sam),
  otu_table(seqtab_ITS_sam_filt, taxa_are_rows = FALSE),
  sample_data(meta_samples_ITS)
  ) # %>% prune_contam(contam_freq_ITS$decontam) # a highly prevalent plant pathogen ...

ps_ITS_ctrl <- phyloseq(
  tax_table(taxa_ITS_ctrl),
  otu_table(seqtab_ITS_ctrl_filt, taxa_are_rows = FALSE),
  sample_data(meta_ctrl_ITS)
)

# Export asvs as fasta
asv_to_fasta(seqtab_ITS_sam_filt, file.path(path_ITS, '4_taxonomy/asv.fa'))

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
taxa_trnL_sam <- subset_asvs(taxa_trnL, seqtab_trnL_sam, 10)
taxa_trnL_ctrl <- subset_asvs(taxa_trnL, seqtab_trnL_ctrl, 10)

seqtab_trnL_sam %>% rowSums %>% sort %>% data.frame(Y=.) %>% 
  ggplot(aes(x = Y)) + geom_histogram(bins=100)

# Remove near-empty samples
seqtab_trnL_sam_filt <- remove_ultra_rare(seqtab_trnL_sam, taxa_trnL_sam, 2000)
seqtab_trnL_ctrl_filt <- remove_ultra_rare(seqtab_trnL_ctrl, taxa_trnL_ctrl, 2000)

dim(seqtab_trnL_sam); dim(seqtab_trnL_sam_filt); dim(taxa_trnL_sam)
dim(seqtab_trnL_ctrl); dim(seqtab_trnL_ctrl_filt); dim(taxa_trnL_ctrl)

# Add sequencing effort and dna concentration to metadata
dna_trnL <- Sys.glob(file.path(dna.path,'CERMO_*trnL*.xlsx')) %>% parse_CERMO_xlsx()
meta_samples_trnL <- add_seq_depth(seqtab_trnL_sam_filt, meta_samples, dna_trnL)
meta_ctrl_trnL <- add_seq_depth(seqtab_trnL_ctrl_filt, meta_ctrl, dna_trnL)

# Find contaminants
contam_freq_trnL <- decontaminate(seqtab_trnL_sam_filt, meta_samples_trnL, 'concDNA')
contam_freq_trnL$p 

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
) # %>% prune_contam(contam_freq_trnL$decontam) # birch?!

# Export asvs as fasta
asv_to_fasta(seqtab_trnL_sam_filt, file.path(path_trnL, '4_taxonomy/asv.fa'))

ps.ls <- list()
ps.ls[["BACT"]] <- ps_16S
ps.ls[["FUNG"]] <- ps_ITS
ps.ls[["PLAN"]] <- ps_trnL
saveRDS(ps.ls, file.path(urbanbio.path,'data/ps.ls.rds'))

# Rarefy phyloseq tables
ps_rare.ls <- lapply(ps.ls, function(ps) {
  set.seed(1234)
  prune_samples(sample_sums(ps) >= 2000, ps) %>% 
    rarefy_even_depth2(ncores = 7) 
})
write_rds(ps_rare.ls, file.path(urbanbio.path, 'data/ps_rare.ls.rds'),
          compress = 'gz')

ps_ctrl.ls <- list()
ps_ctrl.ls[["BACT"]] <- ps_16S_ctrl
ps_ctrl.ls[["FUNG"]] <- ps_ITS_ctrl
ps_ctrl.ls[["PLAN"]] <- ps_trnL_ctrl
saveRDS(ps_ctrl.ls, file.path(urbanbio.path,'data/ps_ctrl.ls.rds'))

# Stats table

imap(ps.ls, function(ps, barcode) {
  asv <- ps %>% otu_table
  seq_per_sam <- rowSums(asv)
  asv_per_sam <- colSums(asv>0)
  tibble(
    Dataset = barcode,
    Seq = sum(asv),
    ASVs = ncol(asv),
    N = nrow(asv),
    Mean_seq = mean(seq_per_sam),
    Min_seq = min(seq_per_sam),
    Max_seq = max(seq_per_sam),
    SD_seq = sd(seq_per_sam),
    Mean_asv = mean(asv_per_sam),
    Min_asv = min(asv_per_sam),
    Max_asv = max(asv_per_sam),
    SD_asv = sd(asv_per_sam)
  )
}) %>% list_rbind
## Not sure why min is 0 ??
# //DEV





# 
# # 
# # Only check samples with positive DNA concententration values
# # 
# ctrl_samples <- c('22-A-QC-Q15-B', '22-E-MTL-23A-B', '22-S-SHER-SH07-B')
# neg_ctrl_sam <- ps_16S_ctrl@otu_table %>% .[ctrl_samples,]
# 
# contam_freq %>%
#   rownames_to_column('ASV') %>%
#   left_join(
#     t(neg_ctrl_sam) %>% as.data.frame %>% rownames_to_column('ASV'),
#     by = 'ASV'
#   ) %>% 
#   filter()
# 
# 
# # Subset OTU table to samples with concDNA info
# asv_table <- ps_16S@otu_table %>%
#   data.frame %>%
#   .[names(concDNA),]
# 
