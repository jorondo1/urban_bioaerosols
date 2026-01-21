# This pipeline is optimised for using only the forward read, as trnL genes
# have extremely variable length and often cannot be merged using 300x300 reads.
# The most important step is deciding at what length to trim, where there is a 
# tradeoff between choosing a lower truncLen to keep more reads but ending up
# with shorter reads and thus lower resolution for taxonomic assignemnt.
# Visualisation is recommended (included in step 3).

# Several steps were retained both Fwd and Rev (#commented) in this pipeline, as 
# the decision not to use merged read was based on comparing different approaches.

#  ml StdEnv/2023 r/4.4.0 mugqic/cutadapt/2.10
### 16S STANDARD DADA2 PIPELINE

library(pacman)
p_load(mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
       dada2, tidyverse, Biostrings, ShortRead, parallel, phyloseq)

# CONFIG
barcode <- 'trnL'
FWD <- "CGAAATCGGTAGACGCTACG"
REV <- "GGGGATAGAGGGACTTGAAC"

ncores <- 60
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
                       # fnRs, fnRs.filtN, 
                       #trimLeft = c(nchar(FWD),nchar(REV)),
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

# Run Cutadapt multicore (each sample is run single-core)
mclapply(seq_along(fnFs), run_cutadapt, mc.cores = ncores)

# Check if it worked?
primer_occurence(fnFs.cut, fnRs.cut, FWD, REV)

##################################
# 3.1. QUALITY FILTERING ########
################################

### Filter and trim for real
# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
#cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Check quality
quality_plot_before <- plotQualityProfile(cutFs[10:20])
ggsave(plot = quality_plot_before,
       filename = file.path(path_data, 'qualPlot_raw.pdf'),
       bg = 'white', width = 2400, height = 2400, 
       units = 'px', dpi = 220)

# Filter samples; define out files
filtFs <- file.path(path_data, 
                    "3_filtered_E22_100", 
                    str_remove(basename(cutFs), '.gz'))
#filtRs <- file.path(path_data, "3_filtered_E22_100", basename(cutRs))

names(filtFs) <- sample.names
#names(filtRs) <- sample.names

out <- filterAndTrim(cutFs, filtFs, #cutRs, filtRs, 
                     maxEE=c(2#,2
                             ),
                     truncQ = 2,
                     minLen = 100, 
                     rm.phix = TRUE, 
                     compress = FALSE, 
                     multithread = ncores)  

quality_plot_after <- plotQualityProfile(filtFs[10:20])
ggsave(plot = quality_plot_after,
       filename = file.path(path_data, 'qualPlot_filtered.pdf'),
       bg = 'white', width = 2400, height = 2400, 
       units = 'px', dpi = 220)

##################################
# 3.2. LENGTH TRUNCATION ########
################################

# Number of sequences per ASV 
all_lengths <- integer(0)
for (file in filtFs) {
  # Read the FASTA file
  sequences <- readDNAStringSet(file, format = 'fastq')  # Use AAStringSet for protein sequences
  
  # Get sequence lengths
  seq_lengths <- width(sequences)
  
  # Append to the main vector
  all_lengths <- c(all_lengths, seq_lengths)
}

# Number of sequences per ASV length (summarise)
asv_len_count <- tibble(all_lengths) %>% 
  group_by(all_lengths) %>% 
  summarise(n = n()) 

write_rds(asv_len_count, 'data/trnL/asv_length_summary.rds')

# Visualize the cumulative number of sequences at various lengths.
# Remember that truncLen will not only truncate the reads at the
# specified length; it will also discards any shorter reads.
asv_len_count %>% 
  arrange(all_lengths, n) %>% 
  mutate(cum_count = cumsum(n)) %>% 
  ggplot(aes(x = all_lengths, y = cum_count)) +
  geom_line() +
  labs(x = 'ASV length', y = 'Cumulative sum of sequences at that length')

asv_len_count %>% 
  mutate(group = case_when(all_lengths >= 254 ~ 'Keep',
                           TRUE ~ 'Drop')) %>% 
  group_by(group) %>% 
  summarise(seq_sum = sum(n)) %>% 
  mutate(prop = seq_sum/sum(seq_sum))

# This is specific to variable-length amplicons where gene length is expected
# to exceed the sequencing range, e.g. 300x300 can theoretically merge amplicons
# of length (300+300) - 2*primer_len - minOverlap = 548 bp, but the real number
# is lower because reads will have been quality trimmed. trnL can go up to 800
# so it is advised to use only the forward reads to prevent introducing a 
# detection bias against species with long trnL.

# Filter samples; define out files
filtFs_survived_trunc <- file.path(
  path_data,
  "3.2_filtered_truncated_275",
  str_replace(basename(filtFs[file.exists(filtFs)]), 'fastq', 'fastq.gz')
) %>% 
  setNames(names(filtFs)[file.exists(filtFs)]) 

out_trunc <- filterAndTrim(filtFs, filtFs_survived_trunc, #cutRs, filtRs, 
                     truncLen = 275,
                     compress = TRUE, 
                     multithread = ncores)  

####################################
# 4. ERROR MODEL AND CLEANUP ########
######################################

# Filtering with minlen may yield empty samples (e.g. neg. controls);
# list files that did survive filtering:
filtFs_survived <- filtFs_survived_trunc[file.exists(filtFs_survived_trunc)]
# filtRs_survived <- filtRs[file.exists(filtRs)]

# Learn errors from the data
errF <- learnErrors(filtFs_survived, multithread = ncores)
# errR <- learnErrors(filtRs_survived, multithread = ncores)

plotErrors(errF, nominalQ = TRUE)
# plotErrors(errR, nominalQ = TRUE)

# Infer sample composition
dadaFs <- dada(filtFs_survived, err = errF, pool = 'pseudo', 
               multithread = min(ncores, 36))

# dadaRs <- dada(filtRs_survived, err = errR, pool = 'pseudo', 
#                multithread = min(ncores, 36))

#####################
### SEQUENCE TABLE ###
#######################
# Modified version of mergePairs that rescues non-merged reads by concatenation
# source('scripts/mergePairsRescue.R'). 
# Intersect the merge and concat; allows merge to fail when overlap is mismatched,
# but recovers non-overlapping pairs by concatenating them. 
# Motivated by https://github.com/benjjneb/dada2/issues/537#issuecomment-412530338
# mergers_pooled <- mergePairsRescue(
#   dadaFs, filtFs_survived, 
#   dadaRs, filtRs_survived,
#   returnRejects = TRUE,
#   minOverlap = 12,
#   maxMismatch = 0,
#   rescueUnmerged = TRUE
# )
### THIS EXPERIMENTAL METHOD seems less effective in increasing the 
### taxonomic assignment rate than just using the forwards. 
### So we stick to the forwards, as recommended here https://github.com/benjjneb/dada2/issues/2091

path.tax <- file.path(path_data, "4_taxonomy_E22_100_trunc275")
if(!dir.exists(path.tax)) dir.create(path.tax)

seqtab <- # makeSequenceTable(mergers_pooled) 
  makeSequenceTable(dadaFs) ## to use FWD READS ONLY
dim(seqtab)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(
  seqtab, method="consensus", multithread = ncores, verbose = TRUE
)

dim(seqtab); dim(seqtab.nochim) # from 17410 to 9024
# distrib of seq len:
# table(nchar(getSequences(seqtab.nochim))) %>% sort %>% plot 
# table(nchar(getSequences(seqtab))) %>% sort %>% plot 

# WRITE OUT
write_rds(seqtab.nochim, paste0(path.tax,'/seqtab.RDS'))

### TRACK PIPELINE READS
track_change <- track_dada(out.N = out.N, out = out,
                           dadaFs = dadaFs, dadaRs = dadaRs,
                           mergers = mergers_pooled,
                           seqtab.nochim = seqtab.nochim)
track_change %>% 
  filter(values>=0) %>% 
  plot_track_change() %>% 
  ggsave(paste0('out/change_',barcode,'.pdf'), plot = ., 
         bg = 'white', width = 1600, height = 1200, 
         units = 'px', dpi = 180)


###################################
# TRNL: CUSTOM TAXONOMY DATABASE ###
#####################################
# Using this custom reference database:
trnl.ref <-'/fast2/def-ilafores/envBarcodeMiner/trnL_CGAAATCGGTAGACGCTACG/hits.lineage.Genus.fa'

# check sequence length distribution
trnl.seq <- readDNAStringSet(trnl.ref, format = 'fasta')
seqlen <- width(trnl.seq) ; seqlen %>% length
ggplot(data.frame(length = seqlen), aes(x = length)) +
  geom_histogram(binwidth = 10, fill = "blue", color = "black") +
  labs(title = "Distribution of Sequence Lengths", x = "Sequence Length", y = "Frequency") +
  theme_minimal()

# filter within expected limits for trnL (???)
filtered_sequences <- trnl.seq[seqlen >= 200 & 
                                 width(trnl.seq) <= 800]
filtered_sequences %>% width %>% sort %>% hist
filt_trnL.path <- paste0(path.tax,'/filtered_trnl_ref.fa')
writeXStringSet(filtered_sequences, filt_trnL.path)

### ASSIGN TAXONOMY ###
taxa <- assignTaxonomy(seqtab.nochim, filt_trnL.path, 
                       multithread = min(ncores, 48), tryRC = TRUE, verbose = TRUE)
taxa[taxa == ""] <- NA # some are NA, others ""
taxa[is.na(taxa)] <- 'Unclassified' # otherwise Tax_glom flushes out the NAs .
write_rds(taxa, paste0(path.tax,'/taxonomy.RDS'))
taxa <- read_rds(paste0(path.tax,'/taxonomy.RDS'))
taxa %>% View


