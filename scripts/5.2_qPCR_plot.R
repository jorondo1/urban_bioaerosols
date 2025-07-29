# Boxplot of qPCR bacterial loads

#############
# X. SETUP ###
###############
library(pacman)
p_load(tidyverse, phyloseq, readxl, ggpubr, rstatix, scales)

source('scripts/0_config.R') # Variable naming and such
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/community_functions.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/psflashmelt.R')

ps <- read_rds('data/ps.ls.rds')[['BACT']]

####################
# qPCR quantification
##########################

# Add qpcr data
parse_bacterial_load <- function(input_tibble) {
  require(magrittr, dplyr)
  input_tibble %>% 
    filter(is.na(note)) %>% 
    mutate(copy_number = as.numeric(copy_number))
}

bact_load <- read_xlsx('data/metadata/load_bacteria_16S.xlsx') %>% 
  set_names(c('Sample', 'replicate', 'copy_number', 'note')) %>% 
  parse_bacterial_load() %>% 
  select(-note)

# filter only samples from dataset
bact_load_samdat <- samdat_as_tibble(ps) %>% 
  select(Sample, time, city) %>% 
  left_join(bact_load,
            by = 'Sample') %>% 
  filter(!is.na(copy_number))

# remove insanely extreme values where 
bact_load_noExtremes <- bact_load_samdat %>% 
  group_by(Sample, city, time) %>%
  summarise(mean_copy_number = mean(copy_number),
         sd_copy_number = sd(copy_number),
         .groups = 'drop') %>% 
  filter(!mean_copy_number > 5e+06) %>% 
  mutate(time = factor(time, levels=c('Spring', 'Summer', 'Fall')),
         mean_copy_number = log10(mean_copy_number))

# By City
wilcox_byCity <- bact_load_noExtremes %>%
  group_by(city) %>%
  wilcox_test(mean_copy_number ~ time) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "time") # Automatically finds positions for the bars

P_city <- bact_load_noExtremes %>%
  ggplot(aes(x = time, y = mean_copy_number, fill = time)) +
  geom_violin(width = 1, alpha = 0.5, linewidth = 0.3) +
  geom_boxplot(width = 0.1, outlier.shape = NA, linewidth = 0.3)+  # Add the pre-computed p-values
  stat_pvalue_manual(
    wilcox_byCity,
    label = "p.adj.signif",
    inherit.aes = FALSE,
    hide.ns = FALSE # Optional: hides non-significant comparisons
  ) +
  facet_grid(.~city) +
  #scale_y_continuous(labels = label_scientific()) +
  scale_fill_manual(values = period_colours) +
  labs(fill = 'Sampling period')

# By City
wilcox_byTime <- bact_load_noExtremes %>%
  group_by(time) %>%
  wilcox_test(mean_copy_number ~ city) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "city") # Automatically finds positions for the bars

P_time<- bact_load_noExtremes %>%
  ggplot(aes(x = city, y = mean_copy_number, fill = city)) +
  geom_violin(width = 1, alpha = 0.5, linewidth = 0.3) +
  geom_boxplot(width = 0.1, outlier.shape = NA, linewidth = 0.3)+  # Add the pre-computed p-values
  stat_pvalue_manual(
    wilcox_byTime,
    label = "p.adj.signif",
    inherit.aes = FALSE,
    hide.ns = FALSE # Optional: hides non-significant comparisons
  ) +
  scale_fill_manual(values = city_colours) +
  facet_grid(.~time) +
  labs(fill = 'City')

P_city / P_time &
  theme(strip.text = element_text(color = "black",size = 14,face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = c(0.01, 0.13),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(face = "bold"),
        legend.justification = c(0, 0.5)) &
  labs(y = 'Mean sample copy number (log10 scale)')

ggsave(paste0('out/_MAIN/7_qPCR_load_test.pdf'),
       bg = 'white', width = 2000, height =1800, 
       units = 'px', dpi = 160)

# fung_load <- read_xlsx('data/metadata/load_fungi_ITS.xlsx') %>% 
#   set_names(c('Sample', 'replicate', 'copy_number_16S','copy_number', 'note')) %>%
#   parse_bacterial_load() %>% 
#   group_by(Sample) %>% 
#   summarise(mean_copy_number = mean(copy_number),
#             sd_copy_number = sd(copy_number),
#             n = n())
# 
# fung_load%>%  filter(n == 1) # many, not ideal
# 