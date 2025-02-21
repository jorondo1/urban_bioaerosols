library(pacman)
p_load(tidyverse, magrittr)

urbanbio.path <- '~/Desktop/ip34/urbanBio'
source(file.path(urbanbio.path, 'scripts/myFunctions.R'))

ps.ls <- read_rds(file.path(urbanbio.path, 'data/ps.ls.rds'))
ps_rare.ls <- read_rds(file.path(urbanbio.path, 'data/ps_rare.ls.rds'))




# median vs vegetation index 
ps.ls$BACT@sam_data %>% as_tibble %>% filter(city != 'None') %>%
  ggplot(aes(x = median_income, y = vegetation_index_NDVI_landsat, colour = city)) + 
  geom_point() + geom_smooth(method = 'lm') +
  theme_light()


