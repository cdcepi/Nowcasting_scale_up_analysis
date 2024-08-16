library(tidyverse)
library(here)
library(scoringutils)

#### Nowcaster
load('./Data/nowcaster.Rdata')
load(here('Data/reports_nowcaster_dates.Rdata'))
reports_nowcaster_dates <-reports_nowcaster_dates %>%
  rename(eventual=n)

#### NobBS
load('./Data/nobbs_for_compare.Rdata')
nobbs_for_compare <-nobbs_for_compare %>%
  group_by(state, model_end_date) %>%
  do(tail(., n=3)) 

#### Epinowcast
load('./Data/epinowcast_all_preds.Rdata')
load('./Data/epinowcast_model_fit.Rdata')
 
all_mods <-all_mods %>%
  mutate(id=paste0(state,".", model_end_date))

times <-all_mods %>%
  dplyr::select(id, run_time)
epinow <-epinow %>%
   mutate(id=paste0(state,".", model_end_date)) %>%
   group_by(state, model_end_date) %>%
   do(tail(., n=3)) %>%
  left_join(., times, "id")
 
## Score
covs <-function(dat){
  dat <-dat %>%
    mutate(cov95=ifelse(eventual >= q2.5 & eventual <= q97.5, 1, 0), #Create coverage metrics
           cov50=ifelse(eventual >= q25 & eventual <= q75, 1, 0)) 
}

all_mods_compare_no_filter <-nowcaster %>%
  group_by(state, model_end_date) %>%
  do(tail(., n=3)) %>%
  left_join(., reports_nowcaster_dates, by=c("onset_wk"="in_wk", "state")) %>%
  ungroup() %>%
  covs(.) %>%
  bind_rows(., nobbs_for_compare) %>%
  mutate(run_time=as.numeric(gsub("secs","", run_time))) %>%
  bind_rows(., epinow) %>%
  dplyr::select(-id, -dispersion, -rhat, -max_rhat,-rhat_flag, -Time) 
save(all_mods_compare_no_filter, file="./Data/all_nowcasts_mods_for_comparison_no_filter.Rdata")

vars_select <-c("median","q1", "q99", "q2.5", "q97.5", "q25", "q75",
                "onset_wk", "model_end_date", #"model_run",
                "state", "method", "onset_wk", "true_value")

interval_scores_no_filter <-all_mods_compare_no_filter %>%
  rename(true_value=eventual) %>%
  group_by(state, model_end_date, method) %>%
  ungroup() %>%
  select(all_of(vars_select)) %>%
  rename(q50=median) %>%
  pivot_longer(q50:q75) %>%
  mutate(name=as.numeric(gsub("q", "", name)),
         name=name/100) %>%
  rename(quantile=name,
         prediction=value,
         model=method) %>%
  score() %>%
  mutate(form="linear_scale")

save(interval_scores_no_filter, file="./Data/interval_scores_no_filter.Rdata")
