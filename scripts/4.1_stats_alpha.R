#############
# X. SETUP ###
###############

library(pacman)
p_load(tidyverse, magrittr, vegan, kableExtra, rstatix, parallel)
source('scripts/0_config.R') # Variable naming and such
source('scripts/myFunctions.R')

######################################
# 1. Alpha diversity across seasons ###
########################################

alphadiv.df <- read_rds('data/diversity/alpha_diversity.rds')

wilcox_e_p_combo <- function(df, grouping_vars, # a vector
                             model) {
  model <- as.formula(model)
  wilcox_p <- df %>% 
    group_by(across(all_of(grouping_vars))) %>% 
    wilcox_test(formula = model) %>% 
    select(-n1, -n2, -.y.) 
  
  wilcox_e <- df %>% 
    group_by(across(all_of(grouping_vars))) %>% 
    wilcox_effsize(formula = model, 
                   ci = TRUE, 
                   ci.type = 'norm')
  
  full_join(wilcox_p, wilcox_e, 
            by = c(grouping_vars, 'group1', 'group2'))
  
}

# Run on all indices
indices <- c('Richness', 'Shannon', 'Simpson', 'Tail')

wilcox.ls <- mclapply(indices, function(var) {
  wilcox_e_p_combo(
    df = alphadiv.df,
    grouping_vars = c('City', 'Barcode'),
    model = as.formula(paste(var, "~ time")))  # Dynamic formula
}, mc.cores = parallel::detectCores()) 

names(wilcox.ls) <- indices

map(indices, function(idx) {
  
  # Effsize/pval Plot 
  wilcox.ls[[idx]] %>% 
    mutate(group_pair = paste0(group1,'_',group2)) %>% 
    ggplot(aes(x = log10(p.adj), y = effsize, 
               colour = Barcode, shape = City)) +
    geom_vline(aes(xintercept = log10(0.05), linetype = "p = 0.05"), 
               color = "red", linewidth = 0.3, alpha = 0.5) +
    geom_vline(aes(xintercept = log10(0.01), linetype = "p = 0.01"), 
               color = "blue", linewidth = 0.3, alpha = 0.5) +
    geom_errorbar(aes(y = effsize,
                      x = log10(p.adj),
                      ymin = conf.low, 
                      ymax = conf.high), 
                  width=.4, linewidth = 0.2,
                  inherit.aes = FALSE) +
    geom_point(size = 3) +
    facet_grid(group_pair~.) +
    xlim(NA, 0) +
    scale_linetype_manual(
      values = c('dashed', 'dashed'),
      name = 'Thresholds'
    ) +
    labs(x = 'Log(adjusted p-value)',
         y = 'Effect size (r)', 
         colour = 'Barcode',
         caption = 'The r value varies from 0 to close to 1. The interpretation values for r commonly in published litterature are:
       0.10 - < 0.3 (small effect), 0.30 - < 0.5 (moderate effect) and >= 0.5 (large effect).')
  
  ggsave(paste0('out/stats/wilcox_',idx,'.pdf'),
         bg = 'white', width = 2400, height = 2000, 
         units = 'px', dpi = 300)
  
  # Output table
  wilcox.ls[[idx]] %>% 
    select(-p, -statistic, -.y.) %>%
    mutate(across(where(is.numeric), ~round(.,3))) %>% 
    kable("html",
          caption = '') %>%
    kable_styling(full_width = FALSE) %>% 
    save_kable(file = paste0('out/stats/wilcox_',idx,'.html'))
  
})

######################################
# 1. Alpha diversity across seasons and NVDI ###
########################################
