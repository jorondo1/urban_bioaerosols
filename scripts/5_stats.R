# 1. Alpha diversity tests
# 2. Beta diversity tests
# 3. Differential abundance tests

#############
# X. SETUP ###
###############

library(pacman)
p_load(tidyverse, magrittr, vegan, kableExtra)
source('scripts/0_config.R') # Variable naming and such
source('scripts/myFunctions.R')

#######################
# 1. Alpha diversity ###
#########################

alphadiv.df <- read_rds('data/diversity/alpha_diversity.rds')


######################
# 2. Beta diversity ###
########################

# Data from scripts/3_metrics.R
betadiv.ls <- read_rds('data/diversity/beta_diversity.ls.rds')

# Permanova-generarting function
# One permanova by city*barcode
iterate_permanova <- function(pcoa.ls, 
                              metric = 'bray',
                              vars, 
                              partType) {
  partType_choices = c(NULL, "terms", "margin", "onedf")
  if (!partType %in% partType_choices){
    stop(paste0("partType must be one of the following: ", paste0(partType_choices, collapse = ', ')))
  }
  
  # Loop by city
  map(cities, function(city) {
    
    # Loop by barcode (target all sublists)
    imap(pcoa.ls[[city]], function(pcoa_barcode.ls, barcode) {
      
      pcoa_barcode.ls <- pcoa_barcode.ls[[metric]]
      dist.mx <- pcoa_barcode.ls$dist.mx
      samData <- pcoa_barcode.ls$metadata
      vars <- intersect(vars, colnames(samData))
      
      formula <- as.formula(paste("dist.mx ~", paste(vars, collapse = " + ")))
      
      res <- adonis2(formula = formula,
                     permutations = 1000,
                     data = samData,
                     by = partType,
                     na.action = na.exclude,
                     parallel = 8)
      
      vars_res <- c(vars, 'Residual', 'Shared') # add residual to pickup value in summary table
      
      shared_var = 1 - sum(res$R2[1:length(vars)]) - res$R2[length(vars)+1] # variance explained by variable interaction
      
      # Parse variables from summary into a tibble
      lapply(seq_along(vars_res), function(i) {
        # each iteration makes a one-row tibble
        tibble(
          City = city,
          Barcode = barcode,
          variable = vars_res[i],
          # shared is "absent" from summary, and it's the last value in vars_res
          R2 = ifelse( # if it's the last iteration on vars_res, 
            i == length(vars_res), shared_var,
            res$R2[i]
          ),
          p = res$`Pr(>F)`[i]) 
        
      }) %>% list_rbind # bind all rows into one tibble
    }) %>% list_rbind # bind all barcode tibbles
  }) %>% list_rbind %>% # bind all city tibbles
    # Reorder factors 
    mutate(variable = factor(variable, levels = c(model_vars, "Shared",'Residual')),
           p_sig = as.factor(case_when(p < 0.05 ~ 'yes', TRUE~ 'no')))
}

model_vars <- c('concDNA' ,'median_income', 'population_density',
                'mean_temperature', 'vegetation_index_NDVI_landsat', 
                'mean_wind_speed', 
                'precip','mean_relative_humidity' , 'date'
                )

perm_out <- list(
  margin = iterate_permanova(betadiv.ls,
                              metric = 'bray',
                              model_vars, 
                              partType = 'margin'),
  terms = iterate_permanova(betadiv.ls,
                            metric = 'bray',
                            model_vars, 
                            partType = 'terms'))

# Kable the output tibbles
map(c('terms', 'margin'), function(effect_type){
  map(cities, function(city) {
    
    # Pivot wide to have one set of columns by barcode
    perm_out[[effect_type]] %>% 
      filter(City == city) %>% 
      select(-City, -p_sig) %>% 
      pivot_wider(values_from = c('R2', 'p'),
                  names_from = Barcode,
                  id_cols = c('variable'),
                  names_glue = "{.value}_{Barcode}",
                  names_vary = "slowest" # control column order, so ordered by names_from and not by values_from
                  ) %>% 
      mutate(across(where(is.numeric), ~ round(., 3))) %>% 
      
      # Build table
      kable("html",
            caption = paste('perMANOVA by', effect_type)) %>%
      kable_styling(full_width = FALSE) %>% 
      add_header_above( # Secondary header
        header = c(
        "Variables" = 1,
        rep(c("R<sup>2</sup>" = 1,
              "<i>p</i>" = 1), 3)
        ), 
        bold = FALSE,
        align = c('l', rep('c', 6)), # align headers left
        escape = FALSE # allows html in variable names
      ) %>%
      add_header_above( # Top header
        header = c(
        setNames(1, city),
        "Bacteria" = 2,
        "Fungi" = 2,
        "Pollen" = 2
      ), align = c('l', rep('c',3))) %>% 
      row_spec(0, extra_css = "display: none;") %>% # Remove original header 
  
      # Save !
      save_kable(file = paste0('out/stats/perm_',city,'_',effect_type,'.html'))
  })
})


# Let's make a cool plot
perm_out$terms %>% 
  ggplot(aes(x = City, y = R2, #alpha = p_sig,
             fill = variable)) +
  geom_col() +
  facet_grid(~Barcode) +
  scale_fill_brewer(palette = 'Set3') +
  #  scale_alpha_manual(values = c("no" = 0.4, "yes" = 1)) +  # Map alpha values explicitly
  theme_light() +
  guides(alpha = 'none')

ggsave('out/perMANOVA_margin.pdf',
       bg = 'white', width = 2400, height = 2000, 
       units = 'px', dpi = 300)








# div ~ city + veg_dens*med_income
# same with distance

# Random effects? MiRKAT (not db) ou brms bayesian?!

