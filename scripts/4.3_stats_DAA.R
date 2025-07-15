#############
# X. SETUP ###
###############

library(pacman)
p_load(tidyverse, magrittr, vegan, kableExtra, ANCOMBC, parallel, MetBrewer)
source('scripts/0_config.R') # Variable naming and such
source('scripts/myFunctions.R')
source(url('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/tax_glom2.R'))

ps.ls <- read_rds('data/ps.ls.rds')
taxLvl <- 'Genus'
ps_glom.ls <- lapply(ps.ls, tax_glom2, taxrank = taxLvl)

run_ancom <- function(ps){
  require(ANCOMBC)
  
  # Reorder levels for Dunnett test
  sample_data(ps)$time <- factor(
    sample_data(ps)$time, 
    levels = c('Summer', 'Spring', 'Fall') # make Summer the reference
    )
  
  ancom_out <- ancombc2(
    data = ps, 
    #tax_level= "taxLvl", # we do it ourselves above to retain the taxLvl level ids
    prv_cut = 0.20, 
    fix_formula="city + time + vegetation_index_NDVI_landsat", 
    group = "time", # specify group if >=3 groups exist, allows structural zero detection 
    struc_zero = TRUE,
    dunnet = TRUE,
    alpha = 0.01,
    verbose = TRUE,
    n_cl = 8 # cores for parallel computing
  )
  
  return(ancom_out)
}

ancom_out_dunnett.ls <- lapply(ps_glom.ls, run_ancom)

write_rds(ancom_out_dunnett.ls, paste0('data/ancom_out_dunnett_',taxLvl,'.RDS'))

# write_rds(ancom_out, 'data/ancom_out_ASV_BACT.RDS') # run by mistake lol


# processing ancombc out 
# ##########
signif_threshold <- 0.01
domain <- 'BACT'

ancom_out_filtered <- tibble(ancom_out_dunnett.ls[[domain]]$res_dunn) %>%
  dplyr::select(-starts_with('W_'), -starts_with('p_'), -starts_with('diff_')) %>% 
  # Remove taxon for which no comparison passed the ss 
  filter(!if_all(.cols = contains('passed_ss_'), 
                 .fns = ~ .x == FALSE)) %>%
  # Then only keep taxa where at least one has q < 0.01
  filter(!if_all(.cols = contains('q_'), 
                 .fns = ~ .x > signif_threshold)) 

# prep data for plotting
ancom_out_long <- ancom_out_filtered %>% 
  # Keep only if one LFC is > 1.3
  filter(!if_all(.cols = contains('lfc_'), 
                 .fns = ~ abs(.x) <1.5)) %>% 
  # long format 
  pivot_longer(cols = -taxon, 
               names_to = c(".value", "Group"), 
               names_pattern = "(lfc|se|q|passed_ss)_(.+)", 
               values_drop_na = TRUE) %>% 
  # lfc become 0 when q > threshold for plotting purposes
  mutate(
    across(c(lfc,se), ~ case_when(q > signif_threshold ~ 0, TRUE ~ .x)),
    textcolour = case_when(lfc==0 ~ "white", TRUE ~ "black"),
    Group = factor(
      case_when(Group == 'timeSpring' ~ 'Spring',
                Group == 'timeFall' ~ 'Fall')
    ), 
    q=q, 
    .keep = 'unused'
  ) %>% 
  left_join(ps_glom.ls[[domain]] %>% # identifier \ species association table
               tax_table() %>% data.frame() %>% 
               select(taxLvl, Class) %>% tibble(),
             join_by(taxon == !!sym(taxLvl)))

### --- PLOT
taxLvls <- ancom_out_long %>% 
  arrange(desc(Class), desc(taxon)) %$% taxon %>% unique

ancom_out_long %<>% 
  # reorder taxa by taxLvl
  mutate(taxon = factor(taxon, levels = taxLvls))

p_main <- ancom_out_long %>% 
  ggplot(aes(x = Group, y = taxon, fill = lfc)) +
  geom_tile() +
  scale_fill_gradient2(low = met.brewer("Cassatt1")[1], 
                       mid = "white", 
                       high = met.brewer("Cassatt1")[8], 
                       midpoint = 0) +
  geom_text(aes(Group, taxon, label = round(lfc, 2), color=textcolour)) +
  scale_color_identity(guide = FALSE) + 
  theme(axis.text.x = element_text(margin = margin(t = 5, r = 5, b = 5, l = 5)),
        axis.text.y = element_blank()
        #legend.position = 'bottom'
        ) +
  labs(x = '', y = '', fill = "Log2 fold-change\nrelative to Summer")


# taxLvl -coloured tile:
p_tile <- ancom_out_long %>% 
  ggplot(aes(x = '1', y = taxon, fill = Class)) +
  geom_tile() + #theme_void() + 
  theme(axis.text.x = element_blank(), 
        axis.title = element_blank(),
        axis.text.y = element_text(hjust = 1)) +
  scale_fill_met_d(name = 'Signac', direction = -1)


p_tile + p_main  + plot_layout(
      guides = "collect",
      design = "ABBBBBBBBBBBBB")
    
    
    
    
    
    
    
    
    
    