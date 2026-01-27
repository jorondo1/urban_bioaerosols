#############
# X. SETUP ###
###############

library(pacman)
p_load(mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
       tidyverse, magrittr, vegan, kableExtra, rstatix, parallel)

######################################
# 1. Alpha diversity across seasons ###
########################################

alphadiv.df <- readRDS('data/diversity/alpha_diversity.rds')

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
# 2. Alpha diversity model within plot variation ###
########################################

p_load(lmerTest, sandwich, lmtest, emmeans)

# Random effect
lmerfit <- alphadiv.df %>% 
  filter(Barcode=="Bacteria") %>% 
  lmer(Shannon ~  (1|site_id),
     data = .) # singular fit
summary(lmerfit)
# null variance for site_id (RE) suggests very similar diverisity within cities

# Remember that we have 1-3 points per site only
alphadiv.df %>% 
  group_by(Barcode, site_id) %>% 
  summarise(N = n()) %>% 
  select(-site_id) %>% 
  table()

alphadiv.df %>% 
  ggplot(aes(x = site_id, y = Shannon, colour = City)) +
  geom_boxplot() +
  facet_grid(Barcode~., scales = 'free')

# What if we only keep sites with 3 samples
site_3 <- alphadiv.df %>% 
  group_by(Barcode, site_id) %>% 
  summarise(N = n()) %>% filter(N==3)

# Regular linear model with time as an explanatory
lmfit <- alphadiv.df %>% 
  filter(Barcode=="Bacteria") %>% 
  lm(Shannon ~ time + City,
       data = .)
summary(lmfit)

# Cluster-robust (by site) standard errors
model <- alphadiv.df %>% 
  filter(Barcode=="Plant particles" & City == 'Sherbrooke') %>%
  lm(Shannon ~ time, data = .)
coeftest(model, vcov = vcovCL, cluster = ~site_id)
coefci(model, vcov = vcovCL, cluster = ~site_id)

# All contrasts
emm <- emmeans(model, ~ time, vcov = vcovCL(model, cluster = ~site_id))
pairs(emm)

# That works, see next script
