# trnL
library(pacman)
p_load(dada2, tidyverse, Biostrings, ShortRead, parallel)
source('scripts/myFunctions.R')
source('scripts/mergePairsRescue.R')

# CONFIG
barcode <- 'trnL'
FWD <- "CGAAATCGGTAGACGCTACG"
REV <- "GGGGATAGAGGGACTTGAAC"

ncores <- 72
path_data <- paste0('data/',barcode)
path_raw <- paste0(path_data, '/0_raw')

fnFs <- sort(list.files(path_raw, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path_raw, pattern="_R2_001.fastq", full.names = TRUE))

<<<<<<< HEAD
=======
(sample.names <- sapply(fnFs, get.sample.name, USE.NAMES = FALSE))
write_delim(data.frame(sample.names), paste0('data/sample_names_',barcode,'.tsv'))

>>>>>>> dc237bc3c5115199b89feae0c1c1bde6c34b6b1e
########################
# 1. N-FILTERING ########
##########################
# Inspect qc
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

### PRE-FILTER (no qc, just remove Ns)
fnFs.filtN <- file.path(path_data, "1_filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path_data, "1_filtN", basename(fnRs))
out.N <- filterAndTrim(fnFs, fnFs.filtN,
                         fnRs, fnRs.filtN, 
                         rm.lowcomplex = TRUE, # added because of https://github.com/benjjneb/dada2/issues/2045#issuecomment-2452299127
                         maxN = 0, 
                         multithread = ncores)
head(out.N)

###########################
# 2. PRIMER REMOVAL ########
#############################

# Analyse primer occurence
primer_occurence(fnFs.filtN, fnRs.filtN, FWD, REV)

### CUTADAPT
cutadapt_path <- "/Users/jorondo/miniconda3/envs/cutadapt/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt_path, args = "--version") # Run shell commands from R

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

# Run Cutadapt multicore (each sample is run single-core)
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

# Check quality
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[10:20])

# Filter samples; define out files
filtFs <- file.path(path_data, "3_filtered", basename(cutFs))
filtRs <- file.path(path_data, "3_filtered", basename(cutRs))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, 
                     maxN = 0, maxEE = c(2, 4), truncQ = 2,
                     minLen = 100, 
                     rm.phix = TRUE, 
                     compress = TRUE, multithread = ncores)  # on windows, set multithread = FALSE

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
<<<<<<< HEAD
dadaFs <- dada(filtFs_survived, err = errF, pool = TRUE, multithread = TRUE, verbose = TRUE)
dadaRs <- dada(filtRs_survived, err = errR, pool = TRUE, multithread = TRUE, verbose = TRUE)
=======
dadaFs <- dada(filtFs_survived, err = errF, pool = TRUE, multithread = ncores)
dadaRs <- dada(filtRs_survived, err = errR, pool = TRUE, multithread = ncores)
>>>>>>> dc237bc3c5115199b89feae0c1c1bde6c34b6b1e

#####################
### SEQUENCE TABLE ###
#######################
# Modified version of mergePairs that rescues non-merged reads by concatenation
# source('scripts/mergePairsRescue.R')
mergers_pooled <- mergePairsRescue(
  dadaFs, filtFs_survived, 
  dadaRs, filtRs_survived,
  returnRejects = TRUE,
  minOverlap = 12,
  maxMismatch = 0,
  rescueUnmerged = TRUE
)

# Intersect the merge and concat; allows merge to fail when overlap is mismatched,
# but recovers non-overlapping pairs by concatenating them. 
# Motivated by https://github.com/benjjneb/dada2/issues/537#issuecomment-412530338

seqtab <- makeSequenceTable(mergers_pooled)
dim(seqtab) # 218 18579
table(nchar(getSequences(seqtab))) %>% sort %>% plot # distrib of seq len

# Clean chimeras
seqtab.nochim <- removeBimeraDenovo(
  seqtab, method="consensus", multithread=TRUE, verbose=TRUE
  )

# WRITE OUT
write_rds(seqtab.nochim, paste0('data/seqtab_', barcode,'.rds'))

### TRACK PIPELINE READS
track_change <- track_dada(out.N = out.N, out = out,
                           dadaFs = dadaFs, dadaRs = dadaRs,
                           mergers = mergers_pooled,
                           seqtab.nochim = seqtab.nochim)

track_change %>% 
  filter(values>=0) %>% 
  plot_track_change()

###################################
# TRNL: CUSTOM TAXONOMY DATABASE ###
#####################################
trnl.ref <- 'data/trnL_hits.lineage.filtered.Genus.fa'

# check sequence length distribution
trnl.seq <- readDNAStringSet(trnl.ref, format = 'fasta')
seqlen <- width(trnl.seq)
ggplot(data.frame(length = seqlen), aes(x = length)) +
  geom_histogram(binwidth = 10, fill = "blue", color = "black") +
  labs(title = "Distribution of Sequence Lengths", x = "Sequence Length", y = "Frequency") +
  theme_minimal()

# filter within expected limits for trnL (???)
filtered_sequences <- trnl.seq[seqlen >= 300 & 
                                 width(trnl.seq) <= 800]

filtered_sequences %>% width %>% sort %>% hist
writeXStringSet(filtered_sequences, "data/filtered_trnl_ref.fa")

######################
### ASSIGN TAXONOMY ###
########################
seqtab_trnL <- read_rds(paste0('data/seqtab_',suffix,'.rds'))

# Filter out sequences smaller than 200bp
keep <- nchar(colnames(seqtab_trnL)) >= 200
seqtab200 <- seqtab_trnL[, keep, drop = FALSE]

# Assign taxonomy
taxa <- assignTaxonomy(seqtab200, 'data/filtered_trnl_ref.fa', multithread = ncores, tryRC = TRUE)
taxa[taxa == ""] <- NA # some are NA, others ""
taxa[is.na(taxa)] <- 'Unclassified' # otherwise Tax_glom flushes out the NAs .
write_rds(taxa, paste0(dirpath,'/taxonomy_',suffix,'.rds'))
taxa <- read_rds(paste0(dirpath,'/taxonomy_',suffix,'.rds'))
taxa %>% View

# Phyloseq Object
meta_trnL <- read_rds('data/meta_trnL.rds')
phyloseq(tax_table(taxa),
         otu_table(seqtab200, taxa_are_rows = FALSE),
         sample_data(meta_trnL))

