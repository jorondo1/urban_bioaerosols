# This pipeline is optimised for using only the forward read, as trnL genes
# have extremely variable length and often cannot be merged using 300x300 reads.
# The most important step is deciding at what length to trim, where there is a 
# tradeoff between choosing a lower truncLen to keep more reads but ending up
# with shorter reads and thus lower resolution for taxonomic assignemnt.
# Visualisation is recommended (included in step 3).

#  ml StdEnv/2023 r/4.4.0 mugqic/cutadapt/2.10
library(dada2)
library(pacman)
p_load(dada2, tidyverse, Biostrings, ShortRead, parallel,
       update = FALSE)
source('scripts/myFunctions.R')

# CONFIG
barcode <- 'trnL_passive'
FWD <- "CGAAATCGGTAGACGCTACG"
REV <- "GGGGATAGAGGGACTTGAAC"

ncores <- 60
path_data <- paste0('data/',barcode)

# Raw data:
path_raw <- paste0(path_data, '/0_raw')
if(!dir.exists(path_raw)) dir.create(path_raw)

# directory for output figures:
path_out <- paste0(path_data, '/_out')
if(!dir.exists(path_out)) dir.create(path_out)

fnFs <- sort(list.files(path_raw, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path_raw, pattern="_R2_001.fastq", full.names = TRUE))

get.sample.name <- function(fname, pattern) {
  base <- tools::file_path_sans_ext(basename(fname))
  sub("_.*", "", base) 
}
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

### Filter and trim , forwards only from hereon out
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))

# Check quality
plotQualityProfile(cutFs[10:20]) %>% 
  ggsave(filename = file.path(path_out, 'qualPlot_raw.pdf'),
         bg = 'white', width = 2400, height = 2400, 
         units = 'px', dpi = 220, plot = .)

# Filter samples; define out files
filtFs <- file.path(path_data, 
                    "3_filtered", 
                    str_remove(basename(cutFs), '.gz'))

names(filtFs) <- sample.names

out <- filterAndTrim(cutFs, filtFs,
                     maxEE=4,
                     truncQ = 2,
                     minLen = 100, 
                     trimRight = 20,
                     rm.phix = TRUE, 
                     compress = FALSE, 
                     multithread = ncores)

plotQualityProfile(filtFs[10:20]) %>% 
  ggsave(filename = file.path(path_out, 'qualPlot_filtered.pdf'),
       bg = 'white', width = 2400, height = 2400, 
       units = 'px', dpi = 220, plot = .)

# Interestingly the quality profiles are a lot worse compared to the
# active traps run

##################################
# 3.2. LENGTH TRUNCATION ########
################################

# Number of sequences per ASV 
all_lengths <- integer(0)
for (file in filtFs) {
  sequences <- readDNAStringSet(file, format = 'fastq')  # Use AAStringSet for protein sequences
  seq_lengths <- width(sequences)
  all_lengths <- c(all_lengths, seq_lengths)
}

# Number of sequences per ASV length (summarise)
asv_len_count <- tibble(all_lengths) %>% 
  group_by(all_lengths) %>% 
  summarise(n = n()) 

# Visualize the cumulative number of sequences at various lengths.
# Remember that truncLen will not only truncate the reads at the
# specified length; it will also discards any shorter reads.

asv_len_plot<- asv_len_count %>% 
  arrange(all_lengths, n) %>% 
  mutate(cum_count = cumsum(n)) %>% 
  ggplot(aes(x = all_lengths, y = cum_count)) +
  geom_line() +
  labs(x = 'ASV length', y = 'Cumulative sum of sequences at that length')

ggsave(filename = file.path(path_out, 'asv_length_summary.pdf'),
       bg = 'white', width = 2400, height = 2400, 
       units = 'px', dpi = 220, plot = asv_len_plot)

# Play around with the proportions of reads kept depending on length
try_truncLen <- 258
asv_len_count %>% 
  mutate(group = case_when(all_lengths >= try_truncLen ~ 'Keep',
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
  paste0("3.2_filt_trunc_",try_truncLen),
  str_replace(basename(filtFs[file.exists(filtFs)]), 'fastq', 'fastq.gz')
) %>% 
  setNames(names(filtFs)[file.exists(filtFs)]) 

out_trunc <- filterAndTrim(filtFs, filtFs_survived_trunc, 
                     truncLen = try_truncLen,
                     compress = TRUE, 
                     multithread = ncores)  

####################################
# 4. ERROR MODEL AND CLEANUP ########
######################################

# Filtering with minlen may yield empty samples (e.g. neg. controls);
# list files that did survive filtering. We also want to let go of 
# samples with too few reads, which happened in this case (some had <100)

min_reads <- 1000  # arbitrary

filtFs_survived <- tibble(path = filtFs_survived_trunc) %>%
  filter(file.exists(path)) %>%                          # Keep existing files
  mutate(sample = sub("\\.gz$", "", basename(path))) %>%  # Remove .gz suffix
  #filter(sample %in% rownames(out_trunc)[out_trunc[, 2] >= min_reads]) %>%  # Filter reads
  filter(sample %in% rownames(out_trunc)) %>% 
  pull(path)                                             # Extract paths

# Learn errors from the data
errF <- learnErrors(filtFs_survived, multithread = ncores)

error_plots <- plotErrors(errF, nominalQ = TRUE)
ggsave(plot = error_plots,
       filename = file.path(path_out, 'errors.pdf'),
       bg = 'white', width = 2400, height = 2400, 
       units = 'px', dpi = 220)

# Infer sample composition
dadaFs <- dada(filtFs_survived, err = errF, pool = 'pseudo', 
               multithread = min(ncores, 36))

#####################
### SEQUENCE TABLE ###
#######################
path.tax <- file.path(path_data, "4_taxonomy")
if(!dir.exists(path.tax)) dir.create(path.tax)

seqtab <- makeSequenceTable(dadaFs)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(
  seqtab, method="consensus", multithread = ncores, verbose = TRUE
)

dim(seqtab); dim(seqtab.nochim) 
sum(seqtab); sum(seqtab.nochim)

# WRITE OUT
write_rds(seqtab.nochim, paste0(path.tax,'/seqtab.RDS'))

### TRACK PIPELINE READS
track_change <- track_dada(out.N = out.N, out = out,
                           dadaFs = dadaFs,
                           seqtab = seqtab,
                           seqtab.nochim = seqtab.nochim)

track_plot <- track_change %>% 
  plot_track_change() 

ggsave(plot = track_plot, filename = file.path(path_out, 'track.pdf'),
       bg = 'white', width = 2400, height = 2400, 
       units = 'px', dpi = 220)

###################################
# TRNL: CUSTOM TAXONOMY DATABASE ###
#####################################
# Using this custom reference database:
trnl.ref <- paste0(path_data,'/trnL_hits.lineage.filtered.Genus.fa')

# check sequence length distribution
trnl.seq <- readDNAStringSet(trnl.ref, format = 'fasta')
seqlen <- width(trnl.seq) ; seqlen %>% length
trnL_seqlen_plot<- ggplot(data.frame(length = seqlen), aes(x = length)) +
  geom_histogram(binwidth = 10, fill = "blue", color = "black") +
  labs(title = "Distribution of Sequence Lengths", x = "Sequence Length", y = "Frequency") +
  theme_minimal()

ggsave(plot = trnL_seqlen_plot, filename = file.path(path_out, 'trnL_seqlen.pdf'),
       bg = 'white', width = 2400, height = 2400, 
       units = 'px', dpi = 220)

# filter within expected limits for trnL (???)
filtered_sequences <- trnl.seq[seqlen >= 200 & 
                                 width(trnl.seq) <= 800]
filtered_sequences %>% width %>% sort %>% hist
filt_trnL.path <- paste0(path_data,'/filtered_trnl_ref.fa')
writeXStringSet(filtered_sequences, filt_trnL.path)

### ASSIGN TAXONOMY ###

# seqtab.nochim <- read_rds('data/trnL/4_taxonomy_E42_100/seqtab.RDS')
# Filter out sequences smaller than 200bp
keep <- nchar(colnames(seqtab.nochim)) >= 100
seqtab200 <- seqtab.nochim[, keep, drop = FALSE]
dim(seqtab200)

# Assign taxonomy
taxa <- assignTaxonomy(seqtab200, filt_trnL.path, 
                       multithread = min(ncores, 48), 
                       taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
                       tryRC = TRUE, verbose = TRUE)
taxa[taxa == ""] <- NA # some are NA, others ""
taxa[is.na(taxa)] <- 'Unclassified' # otherwise Tax_glom flushes out the NAs .
write_rds(taxa, paste0(path.tax,'/taxonomy.RDS'))
# taxa <- read_rds(paste0(path.tax,'/taxonomy.RDS'))
# taxa %>% View