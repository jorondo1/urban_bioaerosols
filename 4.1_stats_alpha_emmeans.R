library(pacman)
p_load(mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
       tidyverse, magrittr, kableExtra, rstatix, parallel,
       lmerTest, sandwich, lmtest, emmeans,
       ggpubr)
source('scripts/0_config.R')

#################################################################
# 1. Alpha diversity across seasons - Clustered standard errors ###
#####################################################################

alphadiv.df <- readRDS('data/diversity/alpha_diversity.rds')

emm_fit <- map(kingdoms, function(barcode) {
  map(cities, function(city) {
    
    lmfit <- alphadiv.df %>% 
      filter(Barcode == barcode & City == city) %>% 
      lm(Shannon ~ time, data = .)
    
    # All contrasts
    emm <- emmeans(lmfit, ~ time, vcov = vcovCL(lmfit, cluster = ~site_id))
    
    # Process output
    data.frame(pairs(emm)) %>% 
      separate(contrast, into = c('group1', 'group2'), sep = " - ") %>% 
      mutate(Barcode = barcode, City = city) %>% 
      as_tibble()  
    
  }) %>% list_rbind()
}) %>% list_rbind()

# Adjust p-values 
emm_fit %<>% 
  mutate(p.adjusted = p.adjust(p.value, method = 'holm'),
         # Convert to y-position (adjust based on your data range)
         # Create comparison labels for geom_signif
         p_label = case_when(
           p.adjusted <0.0001 ~"****",
           p.adjusted < 0.001 ~ "***",
           p.adjusted < 0.01 ~ "**",
           p.adjusted < 0.05 ~ "*",
           TRUE ~ "ns"
         )
  )

write_rds(emm_fit, 'out/stats/alpha_emm_fit.RDS', compress = 'gz')

city <- 'Quebec'

adiv_data <- alphadiv.df %>% 
  filter(City==city) %>% 
  select(Shannon, time, Barcode)

max_div <- adiv_data %>% 
  group_by(Barcode) %>% 
  summarise(y_pos = max(Shannon) + 0.1)

period_levels <- intersect(periods, unique(adiv_data$time))

pval_data <- emm_fit %>% 
  filter(City==city & p.adjusted < 0.05) %>% 
  left_join(max_div, by = 'Barcode') %>%
  group_by(Barcode) %>%
  # Rank comparisons within each Barcode and offset accordingly
  mutate(
    rank = row_number(),
    y_pos = case_when(
      rank == 1 ~ y_pos + 0.1,
      TRUE ~ y_pos + (rank)*0.3
      ),
    x_min = as.numeric(factor(group1, levels = period_levels)),
    x_max = as.numeric(factor(group2, period_levels))
  ) %>%
  ungroup() 

adiv_data %>% 
  ggplot(aes(x = time, y = Shannon, fill = time)) + 
  scale_fill_manual(values = period_colours) +
  geom_violin(width = 1, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_bracket(
    data = pval_data,
    aes(xmin = x_min, xmax = x_max, y.position = y_pos, label = p_label),
    inherit.aes = FALSE,
    step.increase = 0,
    tip.length = 0.01
  ) +
  facet_grid(. ~ Barcode) +
  labs(fill = 'Sampling\nperiod', tag = "D") +
  theme(
    axis.title.x = element_blank()
  )

