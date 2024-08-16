library(tidyverse)
library(rjags)
library(here)

source(here("Code/NobBS_dynamic delays_source code.R"))

# ## Rd in data
load("./Data/dat_COVID.Rdata")
state_for_run <-sort(unique(dat_COVID$state))

# Subset to state
dat_wkly <-dat_COVID %>%
  filter(delay_wk >= 0,
         state==state_for_run[6]) %>%
  arrange(report_wk)

# load("./Data/dat_denv.Rdata")
# dat_wkly <-dat_denv %>%
#   filter(delay_wk >= 0)
# state_for_run <-sort(unique(dat_wkly$state))

##############################################################################
##############################################################################

## Static parameters
delay_min <-4 #Use if selected delay < value
delay_max <-12 #Use if selected delay > value

## Dynamic parameters 
cutpoints <-c(0.90, 0.95, 0.99) #quantile threshold for selecting max delay
delay_mult <-3#c(1, 2) #delay multiplier
wind_mult <-3#c(1, 1.5, 2) #window multiplier

# Look at nowcast performance at: from 3 months after start to 8 months before end date, in weekly increments 
# Note that for the last week or reported data, delays were at 8 months at 90th quantile for all states included
## Create datasets for nowcasts and parameters for each dataset
end_date_choice <-seq(min(dat_wkly$report_wk)+4*3*7, 
                      max(dat_wkly$report_wk)-4*8*7,
                      by="week")
# end_date_choice <-end_date_choice[seq(1, length(end_date_choice), 4)] #for PR only

dat_sub <-list()
delay_dat <-list()
wind_thresh <-list()

report_delays <-list()
parms <-list()
casts <-list()
casts_rhat <-list()
casts_dat <-list()

for(g in 1:length(end_date_choice)){
  dat_sub[[g]] <-dat_wkly %>% #Creates subsets of datasets for running nowcasts
    filter(report_wk <= end_date_choice[[g]])
  
  wind_thresh[[g]] <-length(seq(min(as.Date(dat_sub[[g]]$onset_wk)), max(as.Date(as.Date(dat_sub[[g]]$report_wk))), 
                                by="week")) 
  
  delay_dat[[g]] <-list()
  report_delays[[g]] <-list()
  parms[[g]] <-list()
  casts[[g]] <-list()
  casts_rhat[[g]] <-list()
  casts_dat[[g]] <-list()
  
 for(h in 1:length(cutpoints)){
   delay_dat[[g]][[h]] <-dat_sub[[g]] %>% #Finds delays within those data subsets
    filter(report_wk == max(report_wk)) %>%
    summarise(quants = cutpoints[[h]], 
              delays = round(quantile(delay_wk, cutpoints[[h]])))  
  
  report_delays[[g]][[h]] <-list()
  parms[[g]][[h]] <-list()
  casts[[g]][[h]] <-list()
  casts_rhat[[g]][[h]] <-list()
  casts_dat[[g]][[h]] <-list()
  
  for(i in 1:length(delay_mult)){
    report_delays[[g]][[h]][[i]] <-delay_dat[[g]][[h]] %>% #Manipulates delays by multiplier
      mutate(delays=ifelse(delays<delay_min, delay_min,
                           ifelse(delays>delay_max, delay_max, delays)),
             delays=delays*delay_mult[[i]],
             delays=ifelse(delays >= wind_thresh[[g]], wind_thresh[[g]]-1, delays))

    parms[[g]][[h]][[i]] <-list()
    casts[[g]][[h]][[i]] <-list()
    casts_rhat[[g]][[h]][[i]] <-list()
    casts_dat[[g]][[h]][[i]] <-list()

    for(j in 1:length(wind_mult)){ #These parameter values should be run for each nowcast
      parms[[g]][[h]][[i]][[j]] <-report_delays[[g]][[h]][[i]] %>% #Manipulates window by multiplier
        mutate(multi=wind_mult[[j]],
               winds=round(ifelse(multi==1, delays*multi + 1,
                            delays*multi)), 
               winds=ifelse(winds > wind_thresh[[g]], NA, winds),
               delays=unname(delays)) %>%
        dplyr::select(-multi)

      ## Run nowcasts with parameter values and historic data
      casts[[g]][[h]][[i]][[j]] <-now_cast(df=dat_sub[[g]],
                                      delay=parms[[g]][[h]][[i]][[j]]$delays,
                                      win=parms[[g]][[h]][[i]][[j]]$winds) #%>%
      
      casts_rhat[[g]][[h]][[i]][[j]] <-mcmcr::rhat(casts[[g]][[h]][[i]][[j]], by="term")$sum.n #Rhat values
      
      casts_dat[[g]][[h]][[i]][[j]] <-extract_dat(df=dat_sub[[g]], dat=as.data.frame(as.matrix(casts[[g]][[h]][[i]][[j]]))) %>%
        mutate(delay=parms[[g]][[h]][[i]][[j]]$delays,
               win=ifelse(is.na(parms[[g]][[h]][[i]][[j]]$winds), "NULL", parms[[g]][[h]][[i]][[j]]$winds),
               quant=cutpoints[[h]],
               delay_mult=delay_mult[[i]],
               wind_mult=wind_mult[[j]],
               model_end_date=max(onset_wk),
               model_parms=paste0("delay=",delay,", " , "window=", win),
               model_run=paste0("delay=",delay,", " , "window=", win, " for nowcast on ", model_end_date)) %>%
        bind_cols(., casts_rhat[[g]][[h]][[i]][[j]]) %>%
        rename(rhat=19)
    } #Closes window and nowcast run loop
  } #Closes delay loop
 } #Closes cutpoint loop
} #Closes data subset loop

save(casts, file=here(paste0("Data/nowcast_object_", state_for_run[6],".rds")))
save(casts_dat, file=here(paste0("Data/nowcast_sum_dat_", state_for_run[6],".rds")))
