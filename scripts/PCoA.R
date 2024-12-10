library(pacman)
p_load(tidyverse)
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/community_functions.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/tax_glom2.R")
source("https://github.com/jorondo1/misc_scripts/raw/refs/heads/main/rarefy_even_depth2.R")

urbanbio.path <- '~/Desktop/ip34/urbanBio'
ps.ls <- read_rds(file.path(urbanbio.path, 'data/ps.list.rds'))

ps_rare.ls <- lapply(ps.ls, function(ps) {
  set.seed(1234)
  rarefy_even_depth2(ps, ncores = 7)
})

pcoa.ls <- lapply(ps.ls, function(ps) {
  compute_pcoa(ps, dist = 'bray')
})

pcoa_rare.ls <- lapply(ps_rare.ls, function(ps) {
  compute_pcoa(ps, dist = 'bray')
})
