library(tidyverse)
library(epinowcast)
library(cmdstanr)
library(here)

# Rd in data
load(here("Data/dat_COVID.Rdata"))
state_for_run <-sort(unique(dat_COVID$state))
state_for_run <-state_for_run[2:3]

# Subset to state: Note that only MI and ID are included here
dat_wkly <-dat_COVID %>%
  mutate(delay_day=difftime(report_dt, onset_dt, units="day"), #Estimates difference between report and delay in day
         delay_day=as.numeric(gsub("days", "", delay_day))
  ) %>%
  filter(delay_day >= 0,
         state==state_for_run[2]) %>%
  arrange(report_wk) %>%
  filter(delay_day <= 90) 
rm(dat_COVID)

# Look at nowcast performance at: from 3 months after start to 8 months before end date, in weekly increments 
# Note that for the last week or reported data, delays were at 8 months at 90th quantile for all states included
## Create datasets for nowcasts and parameters for each dataset
end_date_choice <-seq(min(dat_wkly$report_wk)+4*3*7, #min(dat_wkly$report_wk)+4*6*7, # 
                      max(dat_wkly$report_wk)-4*8*7,
                      by="week")

##############################################################################
##############################################################################

## Static parameters
delay_min <-4*7 #Use if selected delay < value
delay_max <-12*7 #Use if selected delay > value
cutpoints <-c(0.95) #quantile threshold for selecting max delay
delay_mult <-c(1) #delay multiplier
#wind_mult <-c(1) #window multiplier; not used here since entire window period is used 

model <-enw_model()#threads = TRUE) # Precompiling the model
options(mc.cores = parallel::detectCores())

##############################################################################
##############################################################################

## Set up needed data
dat_proc2 <-list()
expectation_module <-list()

for(i in 1:length(end_date_choice)){
  dat_sub <-dat_wkly %>% #Creates subsets of datasets for running nowcasts
    filter(report_wk <= end_date_choice[[i]]) %>%
    rename(reference_date=onset_dt, 
           report_date=report_dt)
  
  wind_thresh <-length(seq(min(as.Date(dat_sub$reference_date)), max(as.Date(dat_sub$report_date)), 
                                by="day"))
  
  delay_dat <-dat_sub %>% #Finds delays within those data subsets
    filter(report_wk == max(report_wk)) %>%
    summarise(quants = cutpoints,
              delays = round(quantile(delay_day, cutpoints))) %>% #Manipulates delays by multiplier
    mutate(delays=ifelse(delays<delay_min, delay_min,
                         ifelse(delays>delay_max, delay_max, delays)),
           delays=delays*delay_mult,
           delays=ifelse(delays >= wind_thresh, wind_thresh-1, delays))
  
  dat_sub <-dat_sub %>%
    filter(report_date >= max(report_date)-(delay_dat$delays)) #Truncation to help with run time

  # Set up data for processing with epinowcasts
  dat_proc2[[i]] <-dat_sub %>%
    enw_linelist_to_incidence() %>%
    enw_add_cumulative() %>%
    enw_preprocess_data(., max_delay=delay_dat$delays) 

    # # Expectation model [this basically models the epidemic curve]
  expectation_module[[i]] <-enw_expectation(~ 0 + (1 | day),
                                            latent_reporting_delay = 0.01,
                                            data = dat_proc2[[i]])
  }

rm(dat_wkly, dat_sub, delay_dat)

## Run nowcast

start_all <-Sys.time()
for(i in 38:length(end_date_choice)){
  casts <-epinowcast(data=dat_proc2[[i]],
                          expectation = expectation_module[[i]],
                          fit = enw_fit_opts(
                            save_warmup = FALSE, 
                            chains = 8, 
                            iter_sampling =1000, iter_warmup = 250, 
                            show_messages =FALSE),
                          model = model)
  save(casts, file=here(paste0("Data/epinowcasts/daily_samples_", state_for_run[2],"_epinowcast_run", i,".RDS")))
}
total_time <-(start_all-Sys.time())
