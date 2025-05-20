library(pacman)
p_load(dada2, tidyverse, magrittr, RColorBrewer, ggdist, tidyquant,
       phyloseq, ggh4x, Biostrings, ggridges, optparse, treemapify)

source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/tax_glom2.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/community_functions.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/rarefy_even_depth2.R')
source('~scripts/myFunctions.R')

tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
ps.ls <- readRDS('data/ps.ls.rds')

########################################
# Mean abundance of genera per family ###
##########################################

which_taxrank <- 'Family'
melted <- ps.ls$PLAN %>% 
  tax_glom2(taxrank = which_taxrank) %>% 
  psflashmelt() %>%
  filter(Abundance != 0) %>% 
  group_by(Sample) %>% 
  mutate(relAb = Abundance/sum(Abundance),
         time = factor(time, levels=c('Spring', 'Summer', 'Fall'))) %>% 
  ungroup

(top_taxa <- topTaxa(melted, which_taxrank, 10))

(top_families <- top_taxa[1:7,]$Family)

melt_PLAN <- psflashmelt(ps.ls$PLAN)
melt_family <- melt_PLAN %>% 
  filter(Abundance != 0 &
           Family %in% top_families)

# Mean relative abundance of genera within families across samples
family_genera <- melt_family %>% tibble %>% 
  group_by(Sample, Family, Genus) %>%
  reframe(Abundance = sum(Abundance)) %>%
  group_by(Sample, Family) %>%
  mutate(relAb = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  group_by(Family, Genus) %>%
  reframe(mean_relAb = mean(relAb)) %>%
  group_by(Family) %>%
  mutate(relAbNorm = mean_relAb / sum(mean_relAb)) %>%
  ungroup()


# Treemap
plots <- family_genera %>%
  split(.$Family) %>%
  imap(~ {
    # Check if "Unclassified" is present
    has_unclassified <- "Unclassified" %in% .x$Genus
    
    # Define a unique palette
    unique_colors <- scales::hue_pal()(length(unique(.x$Genus)) - has_unclassified)
    color_mapping <- setNames(
      c(if (has_unclassified) 'grey' else NULL, unique_colors),
      c(if (has_unclassified) "Unclassified" else NULL, setdiff(unique(.x$Genus), "Unclassified"))
    )
    
    # Create the treemap for this family
    ggplot(.x, aes(area = relAbNorm, fill = Genus, label = Genus)) +
      geom_treemap() +
      geom_treemap_text(
        aes(label = ifelse(relAbNorm > 0.01, Genus, "")), # Label only larger proportions
        size = 8,
        place = "centre",
        grow = TRUE
      ) +
      scale_fill_manual(values = color_mapping) +
      labs(title = .y) +
      theme_minimal() +
      theme(legend.position = "none")
  })

wrap_plots(plots, ncol = 1) +
  plot_annotation(title = "Genus Proportions by Family")

ggsave('out/composition/family_genus_prop.pdf', bg = 'white', 
       width = 1200, height = 1600, 
       units = 'px', dpi = 180)


######################################
### ASV length vs. classifiability ####
########################################

tax_tibble <- ps.ls$PLAN@tax_table %>%
  data.frame() %>% 
  rownames_to_column('OTU') %>% tibble

length_classification <- tax_tibble %>% 
  rowwise() %>% 
  mutate(classification_resolution = {
    # Get the vector of column names for non-"unclassified" values
    valid_cols <- tax_ranks[unlist(across(all_of(tax_ranks), ~ .x != "Unclassified"))]
    # Pick the last column name or NA if none exists
    if (length(valid_cols) > 0) tail(valid_cols, 1) else 'No_classification'
  }) %>%
  ungroup() %>% 
  mutate(asv_len = nchar(OTU),
         classification_resolution = factor(classification_resolution,
                                  levels = c('No_classification',tax_ranks))) %>% 
  select(classification_resolution, asv_len)

# classification resolution of our ASVs along their length.
# This shows the length distribution at the highest classification 
# (so each colour represents a distinct set of ASVs)

# First, count ASVs per highest classification
count_labels <- length_classification %>% 
  dplyr::count(classification_resolution) %>% 
  mutate(label = paste0(classification_resolution, " (n = ", n, ")")) %>%
  select(-n) %>% 
  deframe()

length_classification %>% 
  ggplot(aes(x = asv_len, y = classification_resolution)) +
  geom_boxplot(position = position_nudge(y = 0.2), width = 0.4) +
  geom_density_ridges2(
    aes(fill = classification_resolution ),
    alpha = 0.6, show.legend = FALSE,
    scale = 3
    ) +
  scale_y_discrete(expand = c(0.01, 0),
                   labels = count_labels) +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(x = 'ASV length', y = 'Highest classification resolution')+
  theme_classic() 

ggsave(paste0('out/asv_processing/highest_classif_res.pdf'),
       bg = 'white', width = 3600, height = 2000, 
       units = 'px', dpi = 220)

#####################################################
### proportion of unclassified sequences by taxlev ###
#######################################################

melted_all <- melt_PLAN %>% 
  select(all_of(tax_ranks), Abundance, Sample, time, city) %>% 
      mutate(across(all_of(tax_ranks),~ replace_na(., "Unclassified")))
  
# # Extract proportion of unclassified ASV for each rank per sample/approach
test <- map_dfr(tax_ranks, function(taxRank) {
  melted_all %>%
    group_by(Sample, !!sym(taxRank), time, city) %>%
    summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
    group_by(Sample) %>%
    mutate(relAb = Abundance/sum(Abundance)) %>%
    filter(!!sym(taxRank) == "Unclassified") %>%
    mutate(Rank = factor(taxRank, levels = tax_ranks))
})

# Boxplot
test %>% 
  filter(#marker %in% tests_list &
           time != 'NA') %>% 
  ggplot(aes(x = time, y = relAb, fill = time)) +
  geom_violin(linewidth = 0.2, outliers = FALSE) +
  labs(y = 'Proportion of unclassified ASVs') +
  facet_grid(Rank~city, scales = 'free_x') +
  theme_light() +
  labs(x = '') + 
  theme(
    axis.text.x = element_blank()
  ) +
  scale_fill_brewer(palette = "Set2", na.value = "beige")

ggsave('out/asv_processing/unclassified_rate.pdf', bg = 'white', 
       width = 1600, height = 1600, 
       units = 'px', dpi = 180)

# or ridge plot 
test %>%
  filter(Rank != 'Kingdom' &
     #      marker %in% tests_list &
           time != 'NA') %>% 
  ggplot(aes(x = relAb, y = Rank, fill = time)) +
  geom_density_ridges(
    linewidth = 0.2,
    scale = 0.9, 
    alpha = 0.7
  ) +
  facet_grid(~ city, scales = 'free_x') +
 theme_light() +
  labs(
    x = 'Proportion of unclassified ASVs',
    y = 'Taxonomic Rank'
  ) +
  theme(
    axis.title.y = element_text(vjust = 1.5),
    axis.title.x = element_text(vjust = -1)
    ) +
  scale_fill_brewer(palette = "Set2", na.value = "beige")

ggsave('out/asv_processing/unclassified_ridge.pdf', bg = 'white', width = 1600, height = 1200, 
       units = 'px', dpi = 180)


# Look at sequence length in unclassified data :
melted_length <- melt_PLAN %>% 
  filter(Abundance != 0) %>% 
  mutate(seqlen = str_length(OTU))

classified_bool <- map_dfr(tax_ranks, function(taxRank) {
  melted_length %>%
    mutate(classified = case_when(!!sym(taxRank) == 'Unclassified' ~ FALSE,
                                  TRUE~ TRUE),
           Rank = factor(taxRank, levels = tax_ranks))
  })
  
  
classified_bool %>% 
  mutate(classified = factor(classified, levels = c(TRUE, FALSE))) %>% 
  ggplot(aes(x = seqlen, y = Rank, fill = classified)) +
    geom_density_ridges(
      aes(height = sqrt(after_stat(count))),
      linewidth = 0.1,
      stat = 'binline',
      binwidth = 2,
      scale=1, alpha = 0.8
      ) + theme_minimal() +
  scale_fill_manual(values = c('skyblue', 'indianred')) +
  labs(y = 'square root of read count')

ggsave('out/asv_processing/classification_seqlen.pdf', 
         bg = 'white', width = 1200, height = 1200, 
         units = 'px', dpi = 180)
  






















