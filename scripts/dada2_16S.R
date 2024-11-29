### 16S STANDARD DADA2 PIPELINE

library(pacman)
p_load(dada2, tidyverse, Biostrings, ShortRead, parallel, phyloseq)
source('scripts/myFunctions.R')
 
# CONFIG 
barcode <- '16S'
FWD <- "AACMGGATTAGATACCCKG"
REV <- "AGGGTTGCGCTCGTTG"

ncores <- 72
path_data <- paste0('data/',barcode)
path_raw <- paste0(path_data, '/0_raw')

fnFs <- sort(list.files(path_raw, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path_raw, pattern="_R2_001.fastq", full.names = TRUE))

(sample.names <- sapply(fnFs, get.sample.name, USE.NAMES = FALSE))
write_delim(data.frame(sample.names), paste0('data/sample_names_',barcode,'.tsv'))

########################
# 1. N-FILTERING ########
##########################

# Inspect qc
plotQualityProfile(fnFs[10:21])
plotQualityProfile(fnRs[10:21])

### PRE-FILTER (no qc, just remove Ns)
fnFs.filtN <- file.path(path_data, "1_filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path_data, "1_filtN", basename(fnRs))
out.N <- filterAndTrim(fnFs, fnFs.filtN,
                         fnRs, fnRs.filtN, 
                         rm.lowcomplex = TRUE, # added because of https://github.com/benjjneb/dada2/issues/2045#issuecomment-2452299127
                         maxN = 0, 
                         multithread = ncores)

###########################
# 2. PRIMER REMOVAL ########
#############################

# Identify primers
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Analyse primer occurence
primer_occurence(fnFs.filtN, fnRs.filtN, FWD, REV)

### CUTADAPT
cutadapt <- "/Users/jorondo/miniconda3/envs/cutadapt/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path_data, "2_cutadapt")

if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

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

# Check quality
plotQualityProfile(cutFs[10:21])
plotQualityProfile(cutRs[10:21])

# Filter samples; define out files
filtFs <- file.path(path_data, "3_filtered", basename(cutFs))
filtRs <- file.path(path_data, "3_filtered", basename(cutRs))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, 
                     maxN = 0, maxEE = c(1, 1), truncQ = 2,
                     truncLen = c(260, 150),
                     rm.phix = TRUE, 
                     compress = TRUE, multithread = ncores) 

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
errF <- learnErrors(filtFs_survived, multithread = ncores)
errR <- learnErrors(filtRs_survived, multithread = ncores)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# Infer sample composition
dadaFs <- dada(filtFs_survived, err = errF, pool = 'pseudo', multithread = 36, verbose = TRUE)
dadaRs <- dada(filtRs_survived, err = errR, pool = 'pseudo', multithread = 36, verbose = TRUE)

# Merge pairs
mergers <- mergePairs(dadaFs, filtFs_survived, dadaRs, filtRs_survived, verbose=TRUE)
write_rds(mergers, paste0('data/mergers_',barcode,'.RDS'))

# Sequence table
seqtab <- makeSequenceTable(mergers)
rownames(seqtab) <- sample.names
dim(seqtab) # 174 27683
table(nchar(getSequences(seqtab))) %>% sort %>% plot # distrib of seq len

# Clean chimeras
seqtab.nochim <- removeBimeraDenovo(
  seqtab, method="consensus", multithread = ncores, verbose=TRUE
  )
table(nchar(getSequences(seqtab.nochim))) %>% sort %>% plot # distrib of seq len

# Writeout
write_rds(seqtab.nochim, paste0('data/',barcode,'/seqtab_', barcode,'.RDS'))

### TRACK PIPELINE READS
track_change <- track_dada(out.N = out.N, out = out,
                           dadaFs = dadaFs, dadaRs = dadaRs,
                           mergers = mergers_pooled,
                           seqtab.nochim = seqtab.nochim)
# samples lost at merge stage
track_change %>%  filter(lost_merged>0.1) %>% dim 

# Check chimera proportion
track_change %>% 
  pivot_longer(names_to = 'variable', values_to = 'values') %>% 
  ggplot(aes(x = variable, y = values)) +
  geom_boxplot() + theme_minimal() # overall very low chimeric rate

### ASSIGN TAXONOMY
taxa <- assignTaxonomy(seqtab.nochim, 
                          'data/16S/4_taxonomy/silva_nr99_v138.2_toGenus_trainset.fa.gz', 
                          multithread = ncores, tryRC = TRUE, verbose = TRUE)

taxa_species <- assignSpecies(taxaITS, "data/16S/4_taxonomy/silva_v138.2_assignSpecies.fa.gz", verbose = TRUE, tryRC = TRUE, n = 5000)




