# Reproducing Sarah's script

library(pacman)
p_load(dada2, tidyverse, Biostrings, ShortRead, parallel)
source('scripts/myFunctions.R')
source('scripts/mergePairsRescue.R')

# List raw files 
barcode <- 'ITS'
path_dbio <- paste0('/Volumes/DBio_Rech_Data/LaforestLapointeI/PROJECTS/SARAH_POIRIER/DATA/',barcode)
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
head(out.N)

###########################
# 2. PRIMER REMOVAL ########
#############################
# Identify primers
FWD <- "CTTGGTCATTTAGAGGAAGTAA"
REV <- "GCTGCGTTCTTCATCGATGC"
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# Analyse primer occurence
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]])) 

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

# Run Cutadapt multicore (each sample is run single-core)

mclapply(seq_along(fnFs), run_cutadapt, mc.cores = 8)

# Check if it worked?
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

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

# WRITE OUT
write_rds(seqtab.nochim, paste0('data/seqtab_', barcode,'.rds'))

### TRACK PIPELINE READS
getN <- function(x) sum(getUniques(x))
track <- cbind(out.N, out[,2], 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers_pooled, getN),
               rowSums(seqtab.nochim))

colnames(track) <- c("input", "removeNs", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

track_change <- track %>% data.frame %>% 
  rownames_to_column('Sample') %>% 
  tibble %>% 
  filter(filtered>10) %>% 
  mutate(lost_Ns = (input-removeNs)/input,
         lost_filt = (removeNs-filtered)/removeNs,
         lost_noise = (filtered-denoisedR)/filtered,
         lost_merged = (denoisedR-merged)/denoisedR, # Proportion of reads lost to merging
         prop_chimera = (merged-nonchim)/merged) # proportion of chimeras
track_change
# Number of samples lost at merge stage
track_change %>%  filter(lost_merged>0.1) %>% dim 

# Check chimera proportion
track_change %>% 
  #filter(prop_chimera>=0) %>% 
  ggplot(aes(x = NA, y = prop_chimera))+
  geom_boxplot() + theme_minimal() # overall very low chimeric rate


