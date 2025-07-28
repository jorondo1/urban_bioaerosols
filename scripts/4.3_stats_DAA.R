#############
# X. SETUP ###
###############

library(pacman)
p_load(tidyverse, magrittr, vegan, kableExtra, 
       patchwork, gridExtra, cowplot, gtable,
       ANCOMBC, parallel, MetBrewer, phyloseq)
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
    prv_cut = 0.30, 
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

ancom_out_dunnett.ls <- read_rds('data/ancom_out_dunnett_Genus20.RDS')

############################
# processing ancombc out 
# ##########################

signif_threshold <- 0.01

ancom_out_filt.ls <- map(ancom_out_dunnett.ls, function(ancom_out) {
  
  tibble(ancom_out$res_dunn) %>%
    dplyr::select(-starts_with('W_'), -starts_with('p_'), -starts_with('diff_')) %>% 
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

###############################
#### ---- WATERFALL PLOT
################################
legend_labels <- c(#'' = 'white',
  'Higher in Fall' = met.brewer("VanGogh2")[5],
  'Higher in Spring' = met.brewer("VanGogh2")[7],
  'Lower in Fall' = met.brewer("VanGogh2")[2],
  'Lower in Spring' = met.brewer("VanGogh2")[4])

# --- 1. Dataframes for plots
waterfall_plot_df.ls <- imap(ancom_out_long.ls, function(ancom_long, domain) {
  
  # Filter fungi because there are too many abundant genera
  if(domain == 'FUNG') {
    remove_taxa <- ancom_long %>% 
      group_by(taxon) %>% 
      summarise(max_lfc = max(abs(lfc))) %>% 
      filter(max_lfc < 1.5) %>% 
      pull(taxon)
  } else {
    remove_taxa <- c('')
  }
  
  ancom_long %<>% 
    filter(!taxon %in% remove_taxa) 
  
  # Order them considering the mean of non-null lfcs
  taxLvls2 <- ancom_long %>%
    group_by(taxon) %>% 
    summarise(mean_lfc = mean(lfc)) %>% 
    arrange(desc(mean_lfc)) %$% taxon %>% unique()
  
  ancom_long %>% 
    # reorder taxa by taxLvl
    mutate(taxon = factor(taxon, levels = taxLvls2),
           Group = factor(Group, levels = c('Spring', 'Fall')),
           lfc_cat = factor(
             case_when(
               sign(lfc) == 1 & Group == 'Spring' ~ 'Higher in Spring',
               sign(lfc) == -1 & Group == 'Spring' ~ 'Lower in Spring',
               sign(lfc) == 1 & Group == 'Fall' ~ 'Higher in Fall',
               sign(lfc) == -1 & Group == 'Fall' ~ 'Lower in Fall',
               TRUE ~ '' # just to keep the 0s, otherwise they need their own fill aes 
             ), levels = names(legend_labels)
           )
    )
})

# --- 2. Waterfall plots
wf_plots.ls <- imap(waterfall_plot_df.ls, function(plot_df, domain) {
  
  # Background tile
  bg_waterfall_data <- plot_df %>%
    distinct(taxon) %>% # Get unique taxa
    mutate(
      taxon_index = as.numeric(as.factor(taxon)),
      bg_color_id = (taxon_index %% 2) # 0 for pale grey, 1 for white
    )
  bg_colors <- c("0" = "grey80", "1" = "white")
  
  taxon_levels <- levels(bg_waterfall_data$taxon)
  num_taxa <- length(taxon_levels)
  
  plot_df %>% 
    ggplot(aes(x = lfc, y = taxon,
               fill = lfc_cat)) +
    geom_rect(data = bg_waterfall_data, 
              aes(xmin = -Inf, xmax = Inf, 
                  ymin = as.numeric(taxon) - 0.5, 
                  ymax = as.numeric(taxon) + 0.5), 
              fill = bg_colors[factor(bg_waterfall_data$bg_color_id)], # DIRECTLY SET FILL
              alpha = 0.5, 
              inherit.aes = FALSE) +    # Crucial: Don't inherit fill/group from main aes
    geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
    geom_errorbar(aes(x = lfc,
                      xmin = lfc - se, 
                      xmax = lfc + se), 
                  width = 0.2,
                  linewidth = 0.3,
                  position = position_dodge(width = 0.9), 
                  color = "grey20") +
    scale_fill_manual(values = legend_labels,
                      limits = names(legend_labels),
                      drop = FALSE) +
    theme_minimal() +
    theme(legend.position = 'none',
          axis.title = element_blank(),
          axis.text.y = element_blank()
    ) +
    guides(fill = guide_legend(override.aes = list(shape = 21)))  # Ensure shapes show in legend
  
})

# --- 3. Class tiles for left side of the plot
tile_palettes <- list(
  'BACT' = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
             "#FF7F00", "#FFFF33", "#A65628", "#F781BF"),
  'FUNG' = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
             "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"),
  'PLAN' = c("#FFFFB3", "#BEBADA", "#FB8072", "#FDB462")
)


tile_plots.ls <- imap(waterfall_plot_df.ls, function(plot_df, barcode) {
  plot_df %>% 
    ggplot(aes(x = '1', y = taxon, fill = Class)) +
    geom_tile(width = 2) + theme_minimal()+
    theme(axis.text.x = element_blank(), 
          axis.title = element_blank(),
          axis.text.y = element_text(hjust = 1),
          legend.position = 'none') +
    scale_fill_manual(values = tile_palettes[[barcode]])
})

# --- 4. Extract legends
# Common LFC legend :
common_legend_plot <- waterfall_plot_df.ls$BACT %>%
  ggplot(aes(x = 1, y = 1, fill = lfc_cat)) +
  geom_tile() +
  labs(fill = "Log fold-changes in absolute abundances relative to summer") +
  scale_fill_manual(values = legend_labels,
                    limits = names(legend_labels),
                    drop = FALSE) +
  theme(legend.position = "bottom",
        legend.title.position = 'top',
        # Add some margin to the bottom to push the title up if needed
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0))

common_legend <- get_plot_component(common_legend_plot, "guide-box-bottom")

# Distinct Class legends: 
class_legends_grobs <- imap(waterfall_plot_df.ls, function(waterfall_plot_df, barcode) {
  dummy_plot <- waterfall_plot_df %>%
    ggplot(aes(x = lfc, y = taxon, fill = Class)) +
    geom_bar(stat = "identity") +
    labs(fill = paste(kingdoms[barcode], 'classes')) +
    scale_fill_manual(values = tile_palettes[[barcode]]) +
    theme_void() + # Important for extracting just the legend
    theme(
      legend.position = "right", 
      legend.justification = "left",
      #  legend.direction = "vertical",
      legend.box.just = "left", # Align legend box to the left
      legend.box.margin = margin(0,0,0,0) # Remove any default margins
    )
  
  get_plot_component(dummy_plot, "guide-box-right") # get the legend
})

# --- 5. Manual alignment of legend labels 
lg_bact <- class_legends_grobs$BACT
lg_fung <- class_legends_grobs$FUNG
lg_plan <- class_legends_grobs$PLAN
# 
# # Must be at least as 
# max_lg_width <- max(lg_bact$widths, lg_fung$widths, lg_plan$widths)
# lg_bact <- patchwork::wrap_elements(lg_bact) & plot_layout(widths = max_lg_width)
# lg_fung <- patchwork::wrap_elements(lg_fung) & plot_layout(widths = max_lg_width)
# lg_plan <- patchwork::wrap_elements(lg_plan) & plot_layout(widths = max_lg_width)

stacked_distinct_legends <- plot_grid(
  lg_bact,
  lg_fung,
  lg_plan,
  ncol = 1,
  # Use NULL for rel_heights to let cowplot auto-distribute based on content,
  # or provide values like c(4,4,8) if you want specific proportions.
  # For legends, equal spacing is often fine.
  rel_heights = c(1, 1, 1) # Equal heights for simplicity within this stack
)

# --- 6. Define the plots without legends 
P1 <- tile_plots.ls$BACT + labs(title = 'A')
P2 <- tile_plots.ls$FUNG + labs(title = 'B') 
P3 <- tile_plots.ls$PLAN + labs(title = 'C')

P1_combined <- (P1 + wf_plots.ls$BACT + plot_layout(widths = c(1, 20))) & theme(legend.position = "none")
P2_combined <- (P2 + wf_plots.ls$FUNG + plot_layout(widths = c(1, 20))) & theme(legend.position = "none")
P3_combined <- (P3 + wf_plots.ls$PLAN + plot_layout(widths = c(1, 20))) & theme(legend.position = "none")

# --- 7. Grid layout for three plots

# Proportions :
design_layout <- "
acd
acd
acd
bcd
"

# Now, define the full top panel using wrap_plots with this design
# Pass the elements in the order corresponding to the 'a', 'b', 'c', 'd' in the design string
top_panel <- wrap_plots(
  a = P1_combined,
  b = P3_combined,
  c = P2_combined,
  d = stacked_distinct_legends,
  design = design_layout,
  widths = c(21, 21, 12), # Col 1 (for a/b), Col 2 (for c), Col 3 (for d)
  heights = c(4, 4)      # Row 1 (for a/c/d top half), Row 2 (for b/c/d bottom half)
)

# --- 8. Assemble the complete plot : 
top_panel / common_legend +
  plot_layout(
    heights = c(14, 1) # Total height of top_panel and common_legend
  ) & # Adjust the plot letter identifier placement :
  theme(
    plot.title = element_text(
      hjust = 0,
      margin = margin(l = -85, b = -20, unit = "pt")),
    panel.grid.major.y = element_blank(), # Remove major Y grid lines
    panel.grid.minor = element_blank(), # Remove minor Y grid lines
    
  )

# SAVE !
ggsave(paste0('out/DAA/waterfall_composite.pdf'), 
       bg = 'white', width = 2000, height = 2400, 
       units = 'px', dpi = 200)


theme(legend.position = 'bottom',
      ,
      legend.title.position = 'top',
      axis.title = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
)


#########################
### --- HEATMAP PLOT
###########################
taxLvls <- ancom_out_long.ls$BACT %>% 
  arrange(desc(Class), desc(taxon)) %$% taxon %>% unique

ancom_out_long.heatmap <- ancom_out_long.ls$BACT  %>% 
  # reorder taxa by taxLvl
  mutate(taxon = factor(taxon, levels = taxLvls))

p_main <- ancom_out_long.heatmap %>% 
  ggplot(aes(x = Group, y = taxon, fill = lfc)) +
  geom_tile() +
  scale_fill_gradient2(low = met.brewer("Homer2")[1], 
                       mid = "white", 
                       high = met.brewer("Homer2")[6], 
                       midpoint = 0) +
  geom_text(aes(Group, taxon, label = round(lfc, 2), color=textcolour)) +
  scale_color_identity(guide = 'none') + 
  theme_minimal()+
  theme(axis.text.x = element_text(margin = margin(t = 5, r = 5, b = 5, l = 5)),
        axis.text.y = element_blank()
        #legend.position = 'bottom'
  ) +
  labs(x = '', y = '', fill = "Log2 fold-change\nrelative to Summer")

# taxLvl -coloured tile:
p_tile <- ancom_out_long.heatmap %>% 
  ggplot(aes(x = '1', y = taxon, fill = Class)) +
  geom_tile() + theme_minimal()+
  theme(axis.text.x = element_blank(), 
        axis.title = element_blank(),
        axis.text.y = element_text(hjust = 1)) +
  scale_fill_met_d(name = 'Signac', direction = -1)

p_tile + p_main  + plot_layout(
  guides = "collect",
  design = "ABBBBBBBBBBBBB") 

ggsave(paste0('out/DAA/heatmap_',domain,'.pdf'), 
       bg = 'white', width = 1400, height = 2000, 
       units = 'px', dpi = 220)
