library(pacman)
p_load(tidyverse, phyloseq, ggh4x)

# Diversity functions:
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/community_functions.R")
source('scripts/0_config.R') # Variable naming and such

betadiv_full_ls <- read_rds('data/diversity/beta_diversity_full.ls.rds')

aitch_test <- betadiv_full_ls$BACT$robust.aitchison

# Create sample name groups by metadata variable combinations
group_sample_names <- function(df, col_names){
  compact(
    split(
      rownames(df), 
      as.list(df[col_names])
    )
  )
}

# Which variables to group by
group_by_vars <- c('city', 'time')

# Loop over barcodes
dist_long <- 
  imap(kingdoms, function(barcode_name, barcode) {
    
    # Extract diss/dist lists
    barcode.ls <- betadiv_full_ls[[barcode]]
    
    # map over all dist/dissimilarities elements
    imap(barcode.ls, function(dist.ls, dist_metric) {
      
      # Group sample names by metadata
      samples_grouped <- group_sample_names(dist.ls$metadata, group_by_vars)
      
      # For each group, compile long dataset
      imap(samples_grouped, function(samples, sample_group){
        
        # Compile distance pairs 
        compile_dist_pairs(
          dist.mx = dist.ls$dist.mx,
          sample_subset = samples
        ) %>% 
          tibble() %>% 
          mutate(group = sample_group) %>% 
          separate(group, into = group_by_vars) # recreate grouping vars
        # Compile lists :
      }) %>% 
        list_rbind() %>% 
        mutate(Dist_metric = dist_metric)
    }) %>% 
      list_rbind() %>% 
      mutate(Barcode = barcode_name)
  }) %>% 
  list_rbind() %>% 
  mutate(time = factor(time, levels = periods))

# PLot !

dist_long %>% 
  filter(Dist_metric == 'bray') %>% 
  ggplot(aes(x = time, y = Distance, fill = Barcode)) +
  geom_violin(linewidth = 0.2, 
              draw_quantiles = 0.5) +
  facet_nested(city~time, scales = 'free') +
  scale_fill_manual(values = c('red3', 'palegoldenrod', 'forestgreen')) +
  theme_light() +
  theme(
    strip.background = element_rect(fill = 'grey50'),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.5,0.5),
    legend.title = element_blank(),
    legend.background = element_rect(
      fill = "white",    
      color = "black",   
      linewidth = 0.5    
    )
  ) +
  labs(y = 'Bray-Curtis dissimilarity') +
  ylim(0,1)

ggsave('out/_SUPP/diss_within.pdf', bg = 'white', width = 2200, height = 1200, 
       units = 'px', dpi = 220)


## GLM logit link with RE ?
lme4::
  
  
  
  ### GRAVEYARD
  
  ## Are there differences between seasons within city?
  library(rstatix)

# Mann-whitney test despite samples not being independent (sites sampled multiple times)
dist_long %>% 
  filter(Dist_metric == 'robust.aitchison') %>% 
  group_by(city, Barcode) %>% 
  wilcox_test(Distance ~ time)

# Paired test by subsetting

test_U <- dist_long %>% 
  # create site_id pair
  mutate(siteID_pair = paste0(
    str_extract(Sample1, "[^-]+$"),
    "_",
    str_extract(Sample2, "[^-]+$")
  )) 



