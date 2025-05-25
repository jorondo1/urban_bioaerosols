### 16S STANDARD DADA2 PIPELINE

library(pacman)
p_load(dada2, tidyverse, Biostrings, ShortRead, parallel, phyloseq)
source('scripts/myFunctions.R')
 
# CONFIG 
barcode <- '16S'
FWD <- "AACMGGATTAGATACCCKG"
REV <- "AGGGTTGCGCTCGTTG"

ncores <- 80
path_data <- paste0('data/',barcode)
path_raw <- paste0(path_data, '/0_raw')
if(!dir.exists(path_raw)) dir.create(path_raw)

fnFs <- sort(list.files(path_raw, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path_raw, pattern="_R2_001.fastq", full.names = TRUE))

(sample.names <- sapply(fnFs, get.sample.name, USE.NAMES = FALSE))
write_delim(data.frame(sample.names), file.path(path_data, paste0('sample_names_',barcode,'.tsv')))

########################
# 1. N-FILTERING ########
##########################

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

# Analyse primer occurence
primer_occurence(fnFs.filtN, fnRs.filtN, FWD, REV)

### CUTADAPT
# cutadapt <- "/Users/jorondo/miniconda3/envs/cutadapt/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
cutadapt <- '/cvmfs/soft.mugqic/CentOS6/software/cutadapt/cutadapt-2.10/bin/cutadapt'
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path_data, "2_cutadapt")

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
mclapply(seq_along(fnFs), run_cutadapt, mc.cores = ncores)

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
plotQualityProfile(cutFs[10:15])
plotQualityProfile(cutRs[10:15])

# Filter samples; define out files
filtFs <- file.path(path_data, "3_filtered", basename(cutFs))
filtRs <- file.path(path_data, "3_filtered", basename(cutRs))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, 
                     maxN = 0, 
                     maxEE = c(2, 2), 
                     truncQ = 2,
                     truncLen = c(230, 120),
                     rm.phix = TRUE,
                     compress = TRUE, multithread = ncores) 

plotQualityProfile(filtFs[10:15])
plotQualityProfile(filtRs[10:15])

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
dadaFs <- dada(filtFs_survived, errF, pool = 'pseudo', 
               multithread = min(ncores,36), verbose = TRUE)
dadaRs <- dada(filtRs_survived, errR, pool = 'pseudo', 
               multithread = min(ncores,36), verbose = TRUE)

# Merge pairs
mergers <- mergePairs(
  dadaFs, filtFs_survived, 
  dadaRs, filtRs_survived,
  verbose=TRUE)

# Sequence table
path.tax <- file.path(path_data, "4_taxonomy")
if(!dir.exists(path.tax)) dir.create(path.tax)

seqtab <- makeSequenceTable(mergers)

# Clean chimeras
seqtab.nochim <- removeBimeraDenovo(
  seqtab, method="consensus", multithread = ncores, verbose=TRUE
  )

dim(seqtab); dim(seqtab.nochim) # 174 27683
#table(nchar(getSequences(seqtab))) %>% sort %>% plot # distrib of seq len
#table(nchar(getSequences(seqtab.nochim))) %>% sort %>% plot # distrib of seq len

# Writeout
write_rds(seqtab.nochim, file.path(path.tax,'seqtab.RDS'))

### TRACK PIPELINE READS
track_change <- track_dada(out.N = out.N, out = out,
                           dadaFs = dadaFs, dadaRs = dadaRs,
                           mergers = mergers,
                           seqtab.nochim = seqtab.nochim)
track_change %>% 
  filter(values>=0) %>% 
  plot_track_change() %>% 
  ggsave(paste0('out/change_',barcode,'.pdf'), plot = ., 
         bg = 'white', width = 1600, height = 1200, 
         units = 'px', dpi = 180)

### ASSIGN TAXONOMY
taxa <- assignTaxonomy(
  seqtab.nochim,
  file.path(path_data, 'silva_nr99_v138.2_toGenus_trainset.fa.gz'), 
                          multithread = ncores, tryRC = TRUE, verbose = TRUE)

taxa_species <- assignSpecies(
  taxa, 
  file.path(path_data, "silva_v138.2_assignSpecies.fa.gz"), 
  verbose = TRUE, tryRC = TRUE, n = 20000)

write_rds(taxa, file.path(path.tax, 'taxonomy.RDS'))
write_rds(taxa_species, file.path(path.tax, 'taxonomy_species.RDS'))


