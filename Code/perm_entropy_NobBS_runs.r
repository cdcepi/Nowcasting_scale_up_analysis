library(tidyverse)
library(statcomp)
library(here)

load(here("Data/nobbs_dyanmic.Rdata"))
load(here("Data/nobbs_default.Rdata"))
load(here("Data/nobbs_fixed.Rdata"))

dat <-nobbs_fixed %>% #nobbs_default %>% #nobbs_dyanmic %>% #
  #filter(win!="NULL") %>%
  dplyr::select(state, model_end_date, win) %>%
  distinct()
state_for_run <-sort(unique(dat$state))
unique_dat <-dat %>%
  filter(state==state_for_run[7]) %>%
  mutate(win=as.numeric(ifelse(win=="NULL", NA, win)))

# ## Rd in data
load("./Data/dat_COVID.Rdata")
dat_wkly <-dat_COVID %>% #Subset to state
  filter(delay_wk >= 0,
         state==state_for_run[7]) %>% #re-run FL when get a chance
  arrange(report_wk)
#rm(dat_COVID)

# load(here("Data/dat_denv.Rdata")) 

##############################################################################
##############################################################################

## Entropy
perm_entro_full <-list()
for(i in 1:nrow(unique_dat)){
  dat_sub <-dat_wkly %>% #Creates subsets of datasets for running nowcasts
    filter(report_wk <= unique_dat$model_end_date[i]) %>%
    group_by(state, onset_wk) %>%
    count() 
  perm_entro_full[[i]] <-permutation_entropy(ordinal_pattern_distribution(dat_sub$n,4))
} 

perm_entro_part <-list()
for(i in 1:nrow(unique_dat)){
  dat_sub <-dat_wkly %>% #Creates subsets of datasets for running nowcasts
      filter(report_wk <= unique_dat$model_end_date[i]) %>%
      group_by(state, onset_wk) %>%
      count() %>%
      tail(.,  12) 
  perm_entro_part[[i]] <-permutation_entropy(ordinal_pattern_distribution(dat_sub$n,4))
}

parms_pe <-bind_cols(unique_dat, do.call(rbind.data.frame, perm_entro_full), do.call(rbind.data.frame, perm_entro_part)) %>%
  rename(pe_full=4, pe_part=5)
save(parms_pe, file=here(paste0("Data/nowcast_fixed_parameters_with_entropy_", state_for_run[7],".rds")))
#
#"Data/nowcast_dynamic_parameters_with_entropy_", state_for_run[7],".rds")))
#"Data/nowcast_default_parameters_with_entropy_", state_for_run[7],".rds")))
