#############
# X. SETUP ###
###############

library(pacman)
p_load(mgx.tools, # devtools::install_github("jorondo1/mgx.tools")
       tidyverse, magrittr, vegan, kableExtra, 
       patchwork, gridExtra, cowplot, gtable,
       ANCOMBC, parallel, MetBrewer, phyloseq,
       update = FALSE)
source('scripts/0_config.R') # Variable naming and such

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

# ancom_out_dunnett.ls <- lapply(ps_glom.ls, run_ancom)
# write_rds(ancom_out_dunnett.ls, paste0('data/ancom_out_dunnett_',taxLvl,'.RDS'))
# write_rds(ancom_out, 'data/ancom_out_ASV_BACT.RDS') # run by mistake lol

ancom_out_dunnett.ls <- read_rds('data/ancom_out_dunnett_Genus.RDS')

############################
# processing ancombc out 
# ##########################
#ancom_out_dunnett.ls <- readRDS('data/ancom_out_dunnett_Genus.RDS')

signif_threshold <- 0.01

# Manual p-value corrections, because we are testing multiple datasets and
# ancom will only apply correction per dataset (since one ancom function per
# dataset is executed )

# Number of tests performed, multiplied by two because we have 3 groups (two comparisons to the ref group)
num_tests <- 2*sum(sapply(ancom_out_dunnett.ls, 
                          function(x) nrow(x$res_dunn)))

ancom_out_filt.ls <- map(ancom_out_dunnett.ls, function(ancom_out) {
  
  tibble(ancom_out$res_dunn) %>%
    dplyr::select(-starts_with('W_'), -starts_with('q_'), -starts_with('diff_')) %>% 
    # Manual apply holm test with the TOTAL number of statistical tests
    mutate(across(
      starts_with("p_"),
      ~ p.adjust(.x, method = "holm", n = num_tests),
      .names = "{.col}"
    )) %>%
    rename_with(~ sub("^p_", "q_", .), starts_with("p_")) %>% 
    select(-starts_with('p_')) %>% 
    # Remove taxon for which no comparison passed the ss 
    filter(!if_all(.cols = contains('passed_ss_'), 
                   .fns = ~ .x == FALSE)) %>%
    # Then only keep taxa where at least one has q < 0.01
    filter(!if_all(.cols = contains('q_'), 
                   .fns = ~ .x > signif_threshold)) %>% 
    rowwise() %>% 
    # Remove taxon where no comparison is both <= signif_threshold & passed_ss TRUE
    filter(any(
      c_across(contains('q_')) <= signif_threshold
      & c_across(contains('passed_ss_')) == TRUE
    )) %>% ungroup()
})

# prep data for plotting
ancom_out_long.ls <- imap(ancom_out_filt.ls, function(ancom_filt, domain) {
  ancom_filt %>% 
    # Keep only if one LFC is > 1.3
    filter(!if_all(.cols = contains('lfc_'), 
                   .fns = ~ abs(.x) <1)) %>% 
    # long format 
    pivot_longer(cols = -taxon, 
                 names_to = c(".value", "Group"), 
                 names_pattern = "(lfc|se|q|passed_ss)_(.+)", 
                 values_drop_na = TRUE) %>% 
    # lfc become 0 when q > threshold for plotting purposes
    mutate(
      across(c(lfc,se), ~ case_when(q > signif_threshold ~ 0, TRUE ~ .x)),
      textcolour = case_when(lfc==0 | passed_ss == FALSE ~ "white", TRUE ~ "black"),
      Group = factor(
        case_when(Group == 'timeSpring' ~ 'Spring',
                  Group == 'timeFall' ~ 'Fall')
      ), 
      q=q, 
      .keep = 'unused'
    ) %>% 
    left_join(ps_glom.ls[[domain]] %>% # identifier \ species association table
                tax_table() %>% data.frame() %>% 
                select(all_of(taxLvl), Class) %>% tibble(),
              join_by(taxon == !!sym(taxLvl))) %>% 
    mutate(taxon = str_replace(taxon, '_gen_Incertae_sedis', '*')) %>% 
    filter(lfc!=0) 
})

write_rds(ancom_out_long.ls, 'data/DAA_out_long.RDS')
