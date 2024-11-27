##########################
### PROCESSING DADA2 ####
########################
#
### Verify the presence and orientation of these primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dnaTrnL <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orientsTrnL <- c(Forward = dnaTrnL, Complement = Biostrings::complement(dnaTrnL), Reverse = Biostrings::reverse(dnaTrnL), 
                   RevComp = Biostrings::reverseComplement(dnaTrnL))
  return(sapply(orientsTrnL, toString))  # Convert back to character vector
}

### Count the number of times the primers appear in the forward and reverse read, 
# while considering all possible primer orientations
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

### CUTADAPT
run_cutadapt <- function(i) {
  system2(
    cutadapt, args = c(
      R1.flags, R2.flags, 
      "-n", 2, 
      "-m", 21, '-M', 300, # see https://github.com/benjjneb/dada2/issues/2045#issuecomment-2449416862
      "-o", fnFs.cut[i], 
      "-p", fnRs.cut[i],
      fnFs.filtN[i], fnRs.filtN[i])
  )
}

### PLOTS 

topTaxa <- function(psmelt, taxLvl, topN) {
  psmelt %>% 
    group_by(!!sym(taxLvl)) %>% # group by tax level
    filter(relAb != 'NaN') %>% 
    summarise(relAb = mean(relAb)) %>% # find top abundant taxa
    arrange(desc(relAb)) %>% 
    mutate(aggTaxo = as.factor(case_when( # aggTaxo will become the plot legend
      row_number() <= topN ~ !!sym(taxLvl), #+++ We'll need to manually order the species!
      row_number() > topN ~ 'Others'))) # +1 to include the Others section!
}

