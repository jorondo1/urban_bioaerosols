# Reproducing Sarah's script

library(pacman)
p_load(dada2, tidyverse, Biostrings, ShortRead, parallel)

# List raw files 
path_dbio <- '/Volumes/DBio_Rech_Data/LaforestLapointeI/PROJECTS/SARAH_POIRIER/DATA/trnL'
path_dbio_sarah <- paste0(path_dbio,'/trnl_SASS_sarah_replicate')
path_raw <- paste0(path_dbio,'/trnL_SASS_samples_2022/raw')

fnFs <- sort(list.files(path_raw, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path_raw, pattern="_R2_001.fastq", full.names = TRUE))

########################
# 1. N-FILTERING ########
##########################
# Inspect qc
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Identify primers
FWD <- "CGAAATCGGTAGACGCTACG"
REV <- "GGGGATAGAGGGACTTGAAC"

# Verify the presence and orientation of these primers in the data
allOrientsTrnL <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dnaTrnL <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orientsTrnL <- c(Forward = dnaTrnL, Complement = Biostrings::complement(dnaTrnL), Reverse = Biostrings::reverse(dnaTrnL), 
                   RevComp = Biostrings::reverseComplement(dnaTrnL))
  return(sapply(orientsTrnL, toString))  # Convert back to character vector
}
FWD.orientsTrnL <- allOrientsTrnL(FWD)
REV.orientsTrnL <- allOrientsTrnL(REV)

### PRE-FILTER (no qc, just remove Ns)
fnFs.filtN <- file.path(path_dbio_sarah, "1_filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path_dbio_sarah, "1_filtN", basename(fnRs))
outTrnL_N <- filterAndTrim(fnFs, fnFs.filtN,
                         fnRs, fnRs.filtN, 
                         maxN = 0, 
                         multithread = TRUE)
head(outTrnL_N)

###########################
# 2. PRIMER REMOVAL ########
#############################

# Count the number of times the primers appear in the forward and reverse read, 
# while considering all possible primer orientations
primerHitsTrnL <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orientsTrnL, primerHitsTrnL, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orientsTrnL, primerHitsTrnL, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orientsTrnL, primerHitsTrnL, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orientsTrnL, primerHitsTrnL, fn = fnRs.filtN[[1]])) 

# Primers are present so we need to remove them
FWD.lengthTrnL <- 20 #FWD has 20 nucleotides
REV.lengthTrnL <- 20 #REV has 20 nucleotides

# Place filtered files in filtered/ subdirectory
filtFsTrnL <- file.path(path_dbio_sarah, "filtered", paste0(sample.namesTrnL, "_F_filt.fastq.gz"))
filtRsTrnL <- file.path(path_dbio_sarah, "filtered", paste0(sample.namesTrnL, "_R_filt.fastq.gz"))
names(filtFsTrnL) <- sample.namesTrnL
names(filtRsTrnL) <- sample.namesTrnL

outTrnL <- filterAndTrim(fnFsTrnL, filtFsTrnL, fnRsTrnL, filtRsTrnL,
                         trimLeft = c(FWD.lengthTrnL, REV.lengthTrnL), #trimLeft to remove primers
                         truncLen=c(280, 160), # Trim parameter
                         maxN=0, maxEE=c(2,2), truncQ=2,
                         rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE)

####################################
# 4. ERROR MODEL AND CLEANUP ########
######################################
errFTrnL <- learnErrors(filtFsTrnL, multithread=TRUE)
errRTrnL <- learnErrors(filtRsTrnL, multithread=TRUE)

plotErrors(errFTrnL, nominalQ=TRUE)
plotErrors(errRTrnL, nominalQ=TRUE)

# Sample inference
dadaFsTrnL <- dada(filtFsTrnL, err=errFTrnL, multithread=TRUE)
dadaRsTrnL <- dada(filtRsTrnL, err=errRTrnL, multithread=TRUE)
dadaFsTrnL[[1]]

# Merge paired reads and create sequence table
concatTrnL <- mergePairs(dadaFsTrnL, filtFsTrnL, dadaRsTrnL, filtRsTrnL, 
                         returnRejects = TRUE, # NOPE
                         justConcatenate = TRUE, # NOPE
                         verbose=TRUE)
seqtabTrnL <- makeSequenceTable(concatTrnL)
dim(seqtabTrnL) # 220 24668
table(width(getSequences(seqtabTrnL))) %>% sort %>% plot # distrib of seq len

# Clean chimeras
seqtab.nochimTrnL <- removeBimeraDenovo(
  seqtabTrnL, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochimTrnL) # 220 7986
table(nchar(getSequences(seqtab.nochimTrnL)))

write_rds(seqtab.nochimTrnL, 'data/sarah/seqtab.nochimTrnL.rds')

### TRACK PIPELINE READS
getN <- function(x) sum(getUniques(x))
track_trnL <- cbind(outTrnL_N, outTrnL[,2], sapply(dadaFsTrnL, getN), sapply(dadaRsTrnL, getN), sapply(concatTrnL, getN),
               rowSums(seqtab.nochimTrnL))

colnames(track_trnL) <- c("input", "removeNs", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_trnL) <- sample.names

track_change_trnL <- track_trnL %>% data.frame %>% 
  rownames_to_column('Sample') %>% 
  tibble %>% 
  mutate(lost_Ns = (input-removeNs)/input,
         lost_filt = (removeNs-filtered)/removeNs,
         lost_noise = (filtered-denoisedR)/filtered,
         lost_merged = (denoisedR-merged)/denoisedR, # Proportion of reads lost to merging
         prop_chimera = (merged-nonchim)/merged) # proportion of chimeras

track_change_trnL %>%
  filter(lost_merged>0.10) %>% dim # Half of our samples lose more than 15% reads at merging; see https://github.com/benjjneb/dada2/issues/1651#issuecomment-1344934512

# Check chimera proportion
track_change_trnL %>% 
 # filter(prop_chimera>=0) %>% 
  ggplot(aes(x = NA, y = lost_filt ))+
  geom_jitter() + theme_minimal() # overall very low chimeric rate



rbind(track_change %>% mutate(truncLen = 'cutadapt+minlen50'), 
      track_change_trnL %>% mutate(truncLen = 'truncLen+manualPrimerRemoval')) %>% 
  ggplot(aes(x = truncLen, y = prop_chimera ))+ 
  geom_boxplot() + theme_minimal()

rbind(track_change %>% mutate(test = 'Jo'),
      track_change_trnL %>% mutate(test = 'Sarah')) %>% 
  pivot_longer(names_to = 'Step')











