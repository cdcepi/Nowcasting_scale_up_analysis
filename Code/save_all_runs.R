library(tidyverse)
library(here)

## Rd in data
load(here("Data/dat_COVID.Rdata")) ### COVID data
dat_COVID <-dat_COVID %>%
  filter(delay_wk >= 0)
load(here("Data/dat_denv.Rdata")) ### Dengue data
# load("Data/dat_eventual_nobbs.Rdata")

### All data
dat_all <-bind_rows(dat_denv, dat_COVID)

## Reported cases for NobBS
dat_eventual_nobbs <-dat_all %>% 
  group_by(state, onset_wk) %>%
  count() %>%
  rename(eventual=n)
save(dat_eventual_nobbs, file='./Data/dat_eventual_nobbs.Rdata')

##############################################################################
##############################################################################

## Functions for coverage and log scores
covs <-function(dat){
  dat <-dat %>%
    mutate(cov95=ifelse(eventual >= q2.5 & eventual <= q97.5, 1, 0), #Create coverage metrics
           cov50=ifelse(eventual >= q25 & eventual <= q75, 1, 0)) 
}

logs <-function(dat){
  dat <-dat %>%
    mutate(log_score=dnbinom(eventual, mu=mean, size=dispersion, log=T)) #Create log score
}

##############################################################################
##############################################################################

datapath <-here("Data")
state_for_run <-sort(unique(dat_eventual_nobbs$state))

## NobBS with default parameters
filenames <-list.files(path=datapath, pattern="nowcast_sum_dat.+default_parms")
mods_default <-list()
for (i in c(1:length(filenames))) {
  load(paste0(datapath, "/", filenames[i]))
  mods_default[[i]] <-casts_dat  
  mods_default[[i]] <-do.call(rbind.data.frame, mods_default[[i]]) %>% #makes one DF per jurisdiction
    mutate(state=state_for_run[[i]],
           delay_mult="default selection", 
           wind_mult="default selection",
           win=as.character(win),
           method="NobBS",
           model_parms="default selection",
           multi="default selection",
           id=paste0(model_end_date, model_parms, state)) 
}
nobbs_default <-do.call(rbind.data.frame, mods_default) %>%
  left_join(., dat_eventual_nobbs, by=c("onset_wk", "state")) %>%
  mutate(eventual=ifelse(is.na(eventual), 0, eventual)) %>%
  rename(q2.5=lower_95,
         q97.5=upper_95,
         q25=lower_50,
         q75=upper_50) %>%
  covs(.) %>%
  logs(.)

save(nobbs_default, file='./Data/nobbs_default.Rdata')

## NobBS with dynamic daramters
filenames <-list.files(path=datapath, pattern="nowcast_sum_dat")
filenames <-filenames[!grepl("nowcaster", filenames)] 
filenames <-filenames[!grepl("default", filenames)] 
filenames <-filenames[!grepl("fixed", filenames)] 

mods <-list()
for (i in c(1:length(filenames))) {
  load(paste0(datapath, "/", filenames[i]))
  mods[[i]] <-casts_dat  
  mods[[i]] <-do.call(c, mods[[i]]) #groups quantile runs
  mods[[i]] <-do.call(c, mods[[i]]) #groups delay multiplier runs
  mods[[i]] <-do.call(c, mods[[i]]) #groups window multiplier runs
  mods[[i]] <-do.call(rbind.data.frame, mods[[i]]) %>% #makes one DF per jurisdiction
    mutate(state=state_for_run[[i]],
           method="NobBS",
           model_parms=paste0("quantile=", quant, " delay multipler=", delay_mult,", " , "window multipler=", wind_mult),
           multi=paste0("delay multipler=", delay_mult,", " , "window multipler=", wind_mult),
           model_parms=ifelse(win=="NULL", paste0(model_parms, "*"), model_parms),
           id=paste0(model_end_date, model_parms, state),
           wind_mult=as.character(wind_mult)
  )
}
nobbs_dyanmic_og <-do.call(rbind.data.frame, mods) %>%
  left_join(., dat_eventual_nobbs, by=c("onset_wk", "state")) %>%
  mutate(eventual=ifelse(is.na(eventual), 0, eventual)) %>%
  covs(.) %>%
  logs(.)


## Nowcaster 
state_sub <-state_for_run[2:3]
filenames <-list.files(path=datapath, pattern="sum_dat_.+nowcaster") 
filenames <-filenames[2:3]
mods <-list()
mods_nowcaster <-list()
for (i in c(1:length(filenames))) {
  load(paste0(datapath, "/", filenames[i]))
  mods_nowcaster[[i]] <-casts_dat  
  mods_nowcaster[[i]] <-do.call(rbind.data.frame, mods_nowcaster[[i]]) %>% #makes one DF per jurisdiction
    mutate(state=state_sub[[i]],
           onset_wk=onset_wk+7,
           model_end_date=model_end_date+7,
           model_parms=paste0("quantile=", quant, " delay multipler=", delay_mult,", " , "window multipler=", wind_mult),
           id=paste0(model_end_date, model_parms, state)) %>%
    filter(win!="NULL") %>%
    rename(mean=Mean,
           median=Median,
           q1=q01,
           q2.5=LI,
           q97.5=LS,
           q25=LIb,
           q75=LSb)
}
nowcaster <-do.call(rbind.data.frame, mods_nowcaster)
save(nowcaster, file='./Data/nowcaster.Rdata')

## NobBS for comparison
filenames <-list.files(path=datapath, pattern="NobBS_for_.+nowcast_sum_dat")
state_sub <-state_for_run[2:3]
mods_nobbs <-list()
for (i in c(1:length(filenames))) {
  load(paste0(datapath, "/", filenames[i]))
  mods_nobbs[[i]] <-casts_dat
  mods_nobbs[[i]] <-do.call(rbind.data.frame, mods_nobbs[[i]]) %>% #makes one DF per jurisdiction
    mutate(state=state_sub[[i]],
           method="NobBS",
           model_parms=paste0("quantile=", quant, " delay multipler=", delay_mult,", " , "window multipler=", wind_mult),
           id=paste0(model_end_date, model_parms, state)) %>%
    filter(win!="NULL")
}
nobbs_for_compare <-do.call(rbind.data.frame, mods_nobbs) %>%
  left_join(., dat_eventual_nobbs, by=c("onset_wk", "state")) %>%
  mutate(eventual=ifelse(is.na(eventual), 0, eventual)) %>%
  covs(.) %>%
  logs(.)
save(nobbs_for_compare, file='./Data/nobbs_for_compare.Rdata')

## Epinowcast
datapath <-here("Data/epinowcasts")
files_MI <-list.files(path=datapath, pattern="_MI_")
files_ID <-list.files(path=datapath, pattern="_ID_")

mods <-list()
for (h in c(1:length(files_MI))) {
  load(paste0(datapath, "/", files_MI[h]))
  mods[[h]] <-casts }

samples_MI <-list()
for (h in c(1:length(files_MI))) {
  samples_MI[[h]] <-summary(mods[[h]], type="nowcast_samples")} #
for (h in c(1:length(files_MI))) {
  samples_MI[[h]] <-samples_MI[[h]] %>%
    mutate(reference_date=as.Date(reference_date), 
       onset_wk=as.Date(cut(reference_date, "week", start.on.monday = FALSE)) + 6) %>%
  group_by(.iteration, .chain, onset_wk) %>%
  summarize(sum_wk=sum(sample)) %>% # Finds the sum of predictions over each week
  group_by(onset_wk) %>%
  summarize(mean = mean(sum_wk), # Estimates the quantiles of those samples
            median = median(sum_wk),
            q1 = quantile(sum_wk, probs = 0.01),
            q99 = quantile(sum_wk, probs = 0.99),
            q2.5 = quantile(sum_wk, probs = 0.025),
            q97.5 = quantile(sum_wk, probs = 0.975),
            q25 = quantile(sum_wk, probs = 0.25),
            q75 = quantile(sum_wk, probs = 0.75)) %>%
  mutate(model_end_date=max(onset_wk),
         method="Epinowcast", 
         state="MI")
}
casts_MI <-do.call(rbind.data.frame, samples_MI)
save(casts_MI, file="./Data/epinowcast_preds_MI.Rdata")

samples_ID <-list()
for (h in c(1:length(files_ID))) {
  load(paste0(datapath, "/", files_ID[h]))
  mods[[h]] <-casts }
for (h in c(1:length(files_ID))) {
  samples_ID[[h]] <-summary(mods[[h]], type="nowcast_samples") %>%
    mutate(reference_date=as.Date(reference_date), 
           onset_wk=as.Date(cut(reference_date, "week", start.on.monday = FALSE)) + 6) %>%
    group_by(.iteration, .chain, onset_wk) %>%
    summarize(sum_wk=sum(sample, na.rm=T)) %>% # Finds the sum of predictions over each week
    group_by(onset_wk) %>%
    summarize(mean = mean(sum_wk), # Estimates the quantiles of those samples
              median = median(sum_wk),
              q1 = quantile(sum_wk, probs = 0.01),
              q99 = quantile(sum_wk, probs = 0.99),
              q2.5 = quantile(sum_wk, probs = 0.025),
              q97.5 = quantile(sum_wk, probs = 0.975),
              q25 = quantile(sum_wk, probs = 0.25),
              q75 = quantile(sum_wk, probs = 0.75)) %>%
    mutate(model_end_date=max(onset_wk),
           method="Epinowcast", 
           state="ID")
}
casts_ID <-do.call(rbind.data.frame, samples_ID)
save(casts_ID, file="./Data/epinowcast_preds_ID.Rdata")

load(here("Data/dat_eventual_nobbs.Rdata"))
epinow <-bind_rows(casts_MI, casts_ID) %>%
  left_join(., dat_eventual_nobbs, by=c("onset_wk", "state")) %>%
  mutate(eventual=ifelse(is.na(eventual), 0, eventual),
         cov95=ifelse(eventual >= q2.5 & eventual <= q97.5, 1, 0), #Create coverage metrics
         cov50=ifelse(eventual >= q25 & eventual <= q75, 1, 0)) 
save(epinow, file="./Data/epinowcast_all_preds.Rdata")

##############################################################################
##############################################################################

## Runs with default delay at 12 weeks
filenames <-list.files(path=datapath, pattern="fixed.+nowcast_sum_dat")
filenames <-filenames[!grepl("24", filenames)] 
mods <-list()
for (i in c(1:length(filenames))) {
  load(paste0(datapath, "/", filenames[i]))
  mods[[i]] <-casts_dat  
  mods[[i]] <-do.call(c, mods[[i]]) #groups quantile runs
  mods[[i]] <-do.call(rbind.data.frame, mods[[i]]) %>% #makes one DF per jurisdiction
    mutate(state=state_for_run[[i]],
           win=as.character(win),
           quant="Fixed delay at 12 wks",
           model_parms=paste0(quant,
             ", " , "window multipler=", wind_mult),
           multi=paste0("fixed delay, window multipler=", wind_mult),
           model_parms=ifelse(win=="NULL", paste0(model_parms, "*"), model_parms),
           wind_mult=as.character(wind_mult),
           method="NobBS",
           id=paste0(model_end_date, model_parms, state)) 
}
nobbs_fixed_12wks <-do.call(rbind.data.frame, mods) %>%
  left_join(., dat_eventual_nobbs, by=c("onset_wk", "state")) %>%
  mutate(eventual=ifelse(is.na(eventual), 0, eventual)) %>%
  covs(.) %>%
  logs(.)

## Runs at 24 weeks
filenames <-list.files(path=datapath, pattern="fixed 24.+nowcast_sum_dat")
mods <-list()
for (i in c(1:length(filenames))) {
  load(paste0(datapath, "/", filenames[i]))
  mods[[i]] <-casts_dat  
  mods[[i]] <-do.call(c, mods[[i]]) #groups quantile runs
  mods[[i]] <-do.call(rbind.data.frame, mods[[i]]) %>% #makes one DF per jurisdiction
    mutate(state=state_for_run[[i]],
           win=as.character(win),
           quant="Fixed delay at 24 wks",
           model_parms=paste0(quant,
             ", " , "window multipler=", wind_mult),
           multi=paste0("fixed delay, window multipler=", wind_mult),
           model_parms=ifelse(win=="NULL", paste0(model_parms, "*"), model_parms),
           wind_mult=as.character(wind_mult),
           method="NobBS",
           id=paste0(model_end_date, model_parms, state)) 
}
nobbs_fixed_24wks <-do.call(rbind.data.frame, mods) %>%
  left_join(., dat_eventual_nobbs, by=c("onset_wk", "state")) %>%
  mutate(eventual=ifelse(is.na(eventual), 0, eventual), 
         delay=ifelse(delay=="NULL", NA, 24)) %>%
  covs(.) %>%
  logs(.)

nobbs_fixed <-bind_rows(nobbs_fixed_12wks, nobbs_fixed_24wks)
save(nobbs_fixed, file='./Data/nobbs_fixed.Rdata')

##############################################################################
##############################################################################

# Permutation entropy
filenames <-list.files(path=datapath, pattern="fixed.+entropy")
mods <-list()
for (i in c(1:length(filenames))) {
  load(paste0(datapath, "/", filenames[i]))
  mods[[i]] <-as.data.frame(parms_pe) }
nobbs_entro_fixed <-do.call(rbind.data.frame, mods) %>%
  mutate(delay_mod="fixed")
save(nobbs_entro_fixed, file='./Data/nobbs_entro_fixed.Rdata')

filenames <-list.files(path=datapath, pattern="default.+entropy")
mods <-list()
for (i in c(1:length(filenames))) {
  load(paste0(datapath, "/", filenames[i]))
  mods[[i]] <-as.data.frame(parms_pe) }
nobbs_entro_default <-do.call(rbind.data.frame, mods) %>%
  mutate(delay_mod="default")
save(nobbs_entro_default, file='./Data/nobbs_entro_default.Rdata')

filenames <-list.files(path=datapath, pattern="dynamic.+entropy")
mods <-list()
for (i in c(1:length(filenames))) {
  load(paste0(datapath, "/", filenames[i]))
  mods[[i]] <-as.data.frame(parms_pe) }
nobbs_entro_dynamic <-do.call(rbind.data.frame, mods) %>%
  mutate(delay_mod="dynamic")
save(nobbs_entro_dynamic, file='./Data/nobbs_entro_dynamic.Rdata')
