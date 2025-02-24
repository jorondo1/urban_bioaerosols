# Migrate here from ps.ls creation
# Before rarefaction, check contaminants (decontam) and impossible trnLs 
# Produce tables 


# Rarefy phyloseq tables
ps_rare.ls <- lapply(ps.ls, function(ps) {
  set.seed(1234)
  prune_samples(sample_sums(ps) >= 2000, ps) %>% 
    rarefy_even_depth2(ncores = 7) 
})
write_rds(ps_rare.ls, file.path(urbanbio.path, 'data/ps_rare.ls.rds'),
          compress = 'gz')

