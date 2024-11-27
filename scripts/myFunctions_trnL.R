
topTaxa <- function(psmelt, taxLvl, topN) {
  psmelt %>% 
    group_by(!!sym(taxLvl)) %>% # group by tax level
    filter(relAb != 'NaN') %>% 
    summarise(relAb = mean(relAb)) %>% # find top abundant taxa
    arrange(desc(relAb)) %>% 
    mutate(aggTaxo = as.factor(case_when( # aggTaxo will become the plot legend
      row_number() <= topN ~ !!sym(taxLvl), #+++ We'll need to manually order the species!
      row_number() > topN ~ 'Others'))) # +1 to include the Others section!
}
