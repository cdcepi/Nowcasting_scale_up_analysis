library(tidyverse)
library(INLA)
library(here)

gc()
source(here("Code/updated_nowcaster_source_code.R"))

## Rd in data
load(here("Data/dat_COVID.Rdata"))
state_for_run <-sort(unique(dat_COVID$state))
state_for_run <-state_for_run[2:3]

# Subset to state
dat_wkly <-dat_COVID %>%
  filter(delay_wk >= 0,
         state==state_for_run[2]) %>%
  arrange(report_wk)

##############################################################################
##############################################################################

## Parameters
delay_min <-4 #Use if selected delay < value
delay_max <-12 #Use if selected delay > value
cutpoints <-c(0.95) #quantile threshold for selecting max delay
delay_mult <-c(1) #delay multiplier
wind_mult <-c(1) #window multiplier

# Look at nowcast performance at: from 3 months after start to 8 months before end date, in weekly increments 
# Note that for the last week or reported data, delays were at 8 months at 90th quantile for all states included
## Create datasets for nowcasts and parameters for each dataset
end_date_choice <-seq(min(dat_wkly$report_wk)+4*3*7, 
                      max(dat_wkly$report_wk)-4*8*7,
                      by="week")

dat_sub <-list()
wind_thresh <-list()
parms <-list()
casts <-list()
casts_dat <-list()

for(i in 1:length(end_date_choice)){
  
  start_time <-Sys.time()
  
  dat_sub[[i]] <-dat_wkly %>% #Creates subsets of datasets for running nowcasts
    filter(report_wk <= end_date_choice[[i]])
  
  wind_thresh[[i]] <-length(seq(min(as.Date(dat_sub[[i]]$onset_wk)), max(as.Date(as.Date(dat_sub[[i]]$report_wk))), 
                                by="week"))
  
  #These parameter values should be run for each nowcast
  parms[[i]] <-dat_sub[[i]] %>% #Finds delays within those data subsets
    filter(report_wk == max(report_wk)) %>%
    summarise(quants = cutpoints, 
              delays = round(quantile(delay_wk, cutpoints))) %>% #Manipulates delays by multiplier
    mutate(delays=ifelse(delays<delay_min, delay_min,
                             ifelse(delays>delay_max, delay_max, delays)),
               delays=delays*delay_mult,
               delays=ifelse(delays >= wind_thresh[[i]], wind_thresh[[i]]-1, delays)) %>% 
    mutate(multi=wind_mult, #Manipulates window by multiplier
           winds=round(ifelse(multi==1, delays*multi + 1,
                              delays*multi)), 
                 winds=ifelse(winds > wind_thresh[[i]], wind_thresh[[i]], winds),
                 delays=unname(delays)) %>%
    dplyr::select(-multi)
  
  ## Run nowcasts with parameter values and historic data
  casts[[i]] <-nowcasting_inla(dataset = dat_sub[[i]], 
                               date_onset = onset_wk, #onset, 
                               date_report = report_wk, #report, 
                               data.by.week = T, 
                               Dmax=parms[[i]]$delays, 
                               Wdw=parms[[i]]$winds)
  
  casts_dat[[i]]  <-casts[[i]]$total %>%        
    rename(onset_wk=dt_event) %>%
    mutate(delay=parms[[i]]$delays,
           win=ifelse(parms[[i]]$winds > wind_thresh[[i]], "NULL", parms[[i]]$winds),
           quant=cutpoints,
           delay_mult=delay_mult,
           wind_mult=wind_mult,
           model_end_date=max(onset_wk),
           model_parms=paste0("delay=",delay,", " , "window=", win),
           model_run=paste0("delay=",delay,", " , "window=", win, " for nowcast on ", model_end_date),
           method="Nowcaster", 
           run_time=Sys.time()-start_time
           )
} 
save(casts, file=here(paste0("Data/nowcast_samples_", state_for_run[2],"_nowcaster.rds")))
save(casts_dat, file=here(paste0("Data/nowcast_sum_dat_", state_for_run[2],"_nowcaster.rds")))
