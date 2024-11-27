### 16S STANDARD DADA2 PIPELINE

library(pacman)
p_load(dada2, tidyverse, Biostrings, ShortRead, parallel)
source('scripts/myFunctions.R')

# List raw files 
barcode <- '16S'
path_dbio <- paste0('/Volumes/DBio_Rech_Data/LaforestLapointeI/PROJECTS/SARAH_POIRIER/DATA/', barcode)
path_dbio_jo <- paste0(path_dbio,'/',barcode,'_JRL')
path_raw <- paste0(path_dbio,'/sass_samples_2022_',barcode,'/raw')

fnFs <- sort(list.files(path_raw, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path_raw, pattern="_R2_001.fastq", full.names = TRUE))

# List and write out sample names
(sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1))
write_delim(data.frame(sample.names), 'data/sample_names_',barcode,'.tsv')

########################
# 1. N-FILTERING ########
##########################

# Inspect qc
plotQualityProfile(fnFs[5:10])
plotQualityProfile(fnRs[5:10])

### PRE-FILTER (no qc, just remove Ns)
fnFs.filtN <- file.path(path_dbio_jo, "1_filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path_dbio_jo, "1_filtN", basename(fnRs))
out.N <- filterAndTrim(fnFs, fnFs.filtN,
                         fnRs, fnRs.filtN, 
                         rm.lowcomplex = TRUE, # added because of https://github.com/benjjneb/dada2/issues/2045#issuecomment-2452299127
                         maxN = 0, 
                         multithread = TRUE)

###########################
# 2. PRIMER REMOVAL ########
#############################
# Identify primers
FWD <- "AACMGGATTAGATACCCKG"
REV <- "AGGGTTGCGCTCGTTG"

# Analyse primer occurence
primer_occurence(fnFs.filtN, fnRs.filtN, FWD, REV)

### CUTADAPT
cutadapt <- "/Users/jorondo/miniconda3/envs/cutadapt/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path_dbio_jo, "2_cutadapt")

if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# Run Cutadapt (custom function defined in myFunctions.R)
mclapply(seq_along(fnFs), run_cutadapt, mc.cores = 8)

# Check if it worked?
primer_occurence(fnFs.cut, fnRs.cut, FWD, REV)

##############################
# 3. QUALITY FILTERING ########
################################

### Filter and trim for real
# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Check quality
plotQualityProfile(cutFs[10:21])
plotQualityProfile(cutRs[10:21])

# Filter samples; define out files
filtFs <- file.path(path_dbio_jo, "3_filtered_EE42", basename(cutFs))
filtRs <- file.path(path_dbio_jo, "3_filtered_EE42", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, 
                     maxN = 0, maxEE = c(4, 4), truncQ = 2,
                     minLen = 100, #RERUN AT 100 
                     rm.phix = TRUE, 
                     compress = TRUE, multithread = TRUE) 
plotQualityProfile(filtFs[10:21])
plotQualityProfile(filtRs[10:21])

####################################
# 4. ERROR MODEL AND CLEANUP ########
######################################

# Filtering with minlen may yield empty samples (e.g. neg. controls);
# list files that did survive filtering:
filtFs_survived <- filtFs[file.exists(filtFs)]
filtRs_survived <- filtRs[file.exists(filtRs)]

# Learn errors from the data
errF <- learnErrors(filtFs_survived, multithread = TRUE)
errR <- learnErrors(filtRs_survived, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# Infer sample composition
dadaFs <- dada(filtFs_survived, err = errF, pool = TRUE, multithread = TRUE)
dadaRs <- dada(filtRs_survived, err = errR, pool = TRUE, multithread = TRUE)

# Modified version of mergePairs that rescues non-merged reads by concatenation
# source('scripts/mergePairsRescue.R')
mergers_pooled <- mergePairsRescue(dadaFs, filtFs_survived, dadaRs, filtRs_survived,
                 returnRejects = TRUE,
                 minOverlap = 12,
                 maxMismatch = 0,
                 rescueUnmerged = TRUE
)

# Intersect the merge and concat; allows merge to fail when overlap is mismatched,
# but recovers non-overlapping pairs by concatenating them. 
# Motivated by https://github.com/benjjneb/dada2/issues/537#issuecomment-412530338

seqtab <- makeSequenceTable(mergers_pooled)
dim(seqtab) # 220 25606
table(nchar(getSequences(seqtab))) %>% sort %>% plot # distrib of seq len

# Clean chimeras
seqtab.nochim <- removeBimeraDenovo(
  seqtab, method="consensus", multithread=TRUE, verbose=TRUE
  )
table(nchar(getSequences(seqtab.nochim))) %>% sort %>% plot # distrib of seq len

# Writeout
write_rds(seqtab.nochim, paste0('data/seqtab_', barcode,'.rds'))

### TRACK PIPELINE READS
track_change <- track_dada(mergers = mergers_pooled)

# samples lost at merge stage
track_change %>%  filter(lost_merged>0.1) %>% dim 

# Check chimera proportion
track_change %>% 
  pivot_longer(names_to = 'variable', values_to = 'values') %>% 
  ggplot(aes(x = variable, y = values)) +
  geom_boxplot() + theme_minimal() # overall very low chimeric rate

### ASSIGN TAXONOMY
taxaITS <- assignTaxonomy(seqtab.nochim, 
                          paste0(path_dbio, '/sass_samples_2022_ITS/sh_general_release_dynamic_04.04.2024.fasta'), 
                          multithread=TRUE, tryRC = TRUE, verbose = TRUE)

taxaITS %>% as.data.frame() %>%
  mutate(across(Kingdom:Species, ~ gsub("^[a-z]__", "", gsub("_", " ", .))))




