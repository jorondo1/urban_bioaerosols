library(pacman)
p_load(phyloseq, tidyverse, ape, ggtree, RColorBrewer, paletteer)
urbanbio.path <- '~/Desktop/ip34/urbanBio'
source(file.path(urbanbio.path, 'scripts/myFunctions.R'))
ps.ls <- read_rds(file.path(urbanbio.path, 'data/ps.ls.rds'))

barcode <- '16S'
barcode.path <- file.path(urbanbio.path,'data',barcode)
ps <- ps.ls[[barcode_mapping[barcode]]]

# If we masked ASVs to infer tree, there will potentially be duplicates.
# First, create a mapping table between ASVs and masked ASVs (mASV)
ASV.fa <- readDNAStringSet(file.path(barcode.path, "4_taxonomy/asv.fa"))
mASV.fa <- readDNAStringSet(file.path(barcode.path, "5_tree/asv_alignment_masked.fa"))

ASV_mapping <- data.frame(
  ASV_id = names(ASV.fa),
  ASV = as.character(ASV.fa),
  mASV = as.character(mASV.fa),
  stringsAsFactors = FALSE
)

# Verify if taxonomy reconciles with mASVs belonging to multiple ASVs
ASV_mapping_tax <- ps %>% tax_table %>% 
  data.frame %>% 
  rownames_to_column('ASV') %>% 
  right_join(ASV_mapping, by = 'ASV')

ASV_mapping_tax %>% 
  group_by(mASV) %>% 
  filter(n() > 1) %>%
  mutate(n = n()) %>% filter(n == 6) 

# Filter duplicate mASV by keeping the most conservative taxonomy

# Function to find the lowest common taxonomic rank when two or more ASVs share
# the same mASV
find_lowest_common_taxonomy <- function(group) {
  ranks <- taxRanks
  common_row <- group[1, ] # select a single row (doesn't matter which)
  
  for (rank in ranks) {
    unique_values <- unique(group[[rank]]) # extract unique values at taxrank
    unique_non_na <- unique_values[!is.na(unique_values)]
    if (length(unique_non_na) == 1) {
      # If there's only one unique non-NA value, it's the common rank
      common_row[[rank]] <- unique_non_na
    } else {
      # If values disagree, set all subsequent ranks to NA
      common_row[[rank]] <- NA
    }
  }
  
  return(common_row)
}

# Create a subset mASV taxonomy table
mASV_taxonomy <- ASV_mapping_tax %>%
  group_by(mASV) %>%
  group_modify(~ find_lowest_common_taxonomy(.x)) %>%
  ungroup()%>% 
  dplyr::select(-ASV_id, -ASV) 

# Map the conservative taxonomy to the complete ASV table via mASV
ASV_mapping_mTax <- mASV_taxonomy %>% 
  left_join(ASV_mapping, ., by = 'mASV')

# Update/subset/aggregate the phyloseq object with this new taxonomy
mASV_table <- otu_table(ps) %>% 
  data.frame %>% t %>% data.frame(check.names = FALSE) %>% 
  rownames_to_column('ASV') %>% 
  left_join(ASV_mapping_mTax, by = 'ASV') %>%
  dplyr::select(-ASV, -ASV_id) %>% 
  group_by(mASV) %>% 
  summarise(across(where(is.numeric), sum)) %>% # sum seqs of same mASVs
  column_to_rownames('mASV')

# To use ASVs as tip labels, load the fasta
relabel_tips <- function(ASV.fa, tree.path) {
  require(Biostrings)
  fasta <- ASV.fa
  tree <- read.tree(tree.path)
  mapping <- data.frame(
    ASV_label = names(fasta), # Sequence IDs
    ASV = as.character(fasta) # Sequences
  )
  tree$tip.label <- mapping$ASV[match(tree$tip.label, mapping$ASV_label)]
  return(tree)
}

# Import tree
tree <- relabel_tips(ASV.fa = mASV.fa, 
                     tree.path = file.path(barcode.path, "5_tree/tree_masked.treefile"))

# handle duplicate taxa names
duplicate_taxa <- taxa_names(tree)[duplicated(taxa_names(tree))]
unique_tree <- prune_taxa(setdiff(taxa_names(tree), duplicate_taxa), tree)

# Create updated phyloseq object
ps_mASV <- phyloseq(
  tax_table(mASV_taxonomy %>% column_to_rownames('mASV') %>% as.matrix),
  otu_table(mASV_table, taxa_are_rows = TRUE),
  sample_data(ps@sam_data),
  phy_tree(unique_tree)
)

saveRDS(ps_mASV, file.path(urbanbio.path, 'data/ps_mASV.rds'), compress = 'gzip', )

# Load the Newick tree
# tree_ITS <- read.tree(file.path(urbanbio.path, 'data/ITS/5_tree/final_tree.raxml.support'))
# plot(tree)
# ps.ls[['FUNG']] <- merge_phyloseq(ps.ls[['FUNG']], phy_tree(tree_ITS))
# 
# # Load the Newick tree
# tree_trnL <- read.tree(file.path(urbanbio.path, 'data/trnL/5_tree/final_tree.raxml.support'))
# plot(tree)
# ps.ls[['PLAN']] <- merge_phyloseq(ps.ls[['PLAN']], phy_tree(tree_trnL))

# Overwrite rds object
# saveRDS(ps.ls, file.path(urbanbio.path,'data/ps.ls.rds'))


### Tests
which_taxrank <- "Class"
ps_filt <- subset_taxa(ps_mASV, !is.na(tax_table(physeq)[, which_taxrank]))
tree <- phy_tree(ps_filt)
taxonomy <- tax_table(ps_filt)

# Extract
order_mapping <- data.frame(
  tip_label = rownames(taxonomy),
  Taxrank = as.character(taxonomy[, which_taxrank])
) 
#order_mapping$Taxrank[is.na(order_mapping$Taxrank)] <- "Unknown"

# Add taxonomic information to the tree in ggtree format
tree_gg <- ggtree(tree)
tree_gg <- tree_gg %<+% order_mapping

nTaxa <- order_mapping$Taxrank %>% unique %>% length
expanded_palette <- colorRampPalette(paletteer_d("ggthemes::Tableau_20"))(nTaxa) 

tree_plot <- tree_gg +
  geom_tree() + # Plot the tree structure
  geom_tippoint(aes(color = Taxrank), size = 1) + # Color tip nodes by Taxrank
  theme_tree2() + # Clean tree theme
  theme(legend.position = "right") + # Add legend for Taxrank
  scale_color_manual(values = expanded_palette) +
  geom_tiplab(aes(label = NA)) +# Suppress tip labels
  labs(color = which_taxrank)

tree_plot
ggsave(paste0('~/Desktop/ip34/urbanBio/out/tree_',barcode,'_',which_taxrank,'_masked.pdf'),
       bg = 'white', width = 2000, height = 2400, 
       units = 'px', dpi = 300)








