# Author : Jonathan Rondeau-Leclaire 2025

##############
# Contents ####
################

# 0.1. Load functions
# 0.2. Setup / metadata parsing
# 1.1. Process 16S samples
# 1.2. Process ITS samples
# 1.3. Process trnL samples
# 2.1. Export ps objects 
# 3.1. ASV stats table
# 3.2. Taxonomic classification rate

###############################
# 0.1. Functions & packages ####
#################################

library(pacman)
p_load(tidyverse, phyloseq, magrittr, decontam, Biostrings,
       readxl, decontam)

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

######################
# 0.2. Parse metadata ###
########################

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

####################
# 1.1. 16S samples ####
######################

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

# Seq count distribution
viz_seqdepth(seqtab_16S_sam)

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

####################
# 1.2. ITS samples ####
######################

path_ITS <- file.path(urbanbio.path,'data/ITS')
taxa_ITS <- read_rds(file.path(path_ITS, '4_taxonomy/taxonomy.RDS'))
seqtab_ITS <- read_rds(file.path(path_ITS, '4_taxonomy/seqtab.RDS'))

# Keep samples with metadata info
seqtab_ITS_sam <- subset_samples(seqtab_ITS, sample.names)
seqtab_ITS_ctrl <- subset_samples(seqtab_ITS, ctrl.names)

# Subset ASVs
taxa_ITS_sam <- subset_asvs(taxa_ITS, seqtab_ITS_sam, 10)
taxa_ITS_ctrl <- subset_asvs(taxa_ITS, seqtab_ITS_ctrl, 10)

# Seq count distribution
viz_seqdepth(seqtab_ITS_sam)

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
# contam_freq_ITS <- decontaminate(seqtab_ITS_sam_filt, meta_samples_ITS, 'concDNA')
# contam_freq_ITS$p # look at contaminants correlation with concDNA

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

#######################
# 1.3. trnL samples ####
#########################

path_trnL <- file.path(urbanbio.path,'data/trnL')
taxa_trnL <- read_rds(file.path(path_trnL, '4_taxonomy/taxonomy.RDS'))
seqtab_trnL <- read_rds(file.path(path_trnL, '4_taxonomy/seqtab.RDS'))

# Keep samples with metadata info
seqtab_trnL_sam <- subset_samples(seqtab_trnL, sample.names)
seqtab_trnL_ctrl <- subset_samples(seqtab_trnL, ctrl.names)

# Subset ASVs
taxa_trnL_sam <- subset_asvs(taxa_trnL, seqtab_trnL_sam, 10)
taxa_trnL_ctrl <- subset_asvs(taxa_trnL, seqtab_trnL_ctrl, 10)

# Seq count distribution
viz_seqdepth(seqtab_16S_sam)

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
# contam_freq_trnL <- decontaminate(seqtab_trnL_sam_filt, meta_samples_trnL, 'concDNA')
# contam_freq_trnL$p 

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

###########################################
# 2.1. Export ps object lists and rarefy ###
#############################################

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

#############################
# 3.1. ASV stats table #######
###############################
p_load(knitr, kableExtra, webshot2)

# Prep data
ps.stats <- imap(ps.ls, function(ps, barcode) {
  asv <- ps %>% otu_table 
  seq_per_sam <- rowSums(asv)
  asv_per_sam <- rowSums(asv>0)
  asv_prevalence <- colSums(asv>0)
  num_sam <- nrow(asv)
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

html_file <- file.path(urbanbio.path, "out/output_table.html")
save_kable(ps.stats.k, file = html_file)

# Use webshot to convert the HTML file to PDF
webshot(html_file, 
        file.path(urbanbio.path, "out/output_table.pdf"),
        cliprect = "viewport")

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



# Set browser
# Sys.setenv(CHROMOTE_CHROME = "/Applications/Brave Browser.app/Contents/MacOS/Brave Browser")
