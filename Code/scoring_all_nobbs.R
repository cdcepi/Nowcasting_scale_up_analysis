library(tidyverse)

## Threshold for log scores
log_thresh <-log(1/10E3)

## PE data
load(here("Data/nobbs_entro_dynamic.Rdata"))
load(here("Data/nobbs_entro_default.Rdata"))
load(here("Data/nobbs_entro_fixed.Rdata"))

### Nowcast data
load(here("Data/nobbs_dyanmic.Rdata"))
nobbs_dyanmic <-nobbs_dyanmic %>%
  mutate(delay_mult=paste0("selected delay x ", delay_mult)) %>%
  mutate(quant=format(quant, nsmall = 2),
         quant=paste0(quant, " quantile")
  ) %>%
  mutate(win=as.numeric(ifelse(win=="NULL", NA, win))) %>%
  left_join(., nobbs_entro_dynamic, by=c("state", "model_end_date", "win"))

load(here("Data/nobbs_default.Rdata"))
nobbs_default <-nobbs_default %>% 
  mutate(win=as.numeric(win)) %>%
  left_join(., nobbs_entro_default, by=c("state", "model_end_date", "win")) %>%
  mutate(quant=" ")

load(here("Data/nobbs_fixed.Rdata"))
nobbs_fixed <-nobbs_fixed %>% 
  mutate(win=as.numeric(ifelse(win=="NULL", NA, win))) %>%
  left_join(., nobbs_entro_fixed, by=c("state", "model_end_date", "win"))


##############################################################################################
##############################################################################################

all_nobbs <-bind_rows(nobbs_dyanmic, nobbs_default) %>% #set factor order for multipliers
  bind_rows(., nobbs_fixed) %>%
  group_by(state, multi, quant, model_parms, model_end_date, model_run, method) %>%
  do(tail(., n=3)) %>%  #Limit eval metrics to last 3 nowcasts made
  mutate(id=paste0(model_end_date, "_", model_parms, "_", state),
         quant=as.factor(quant), 
         quant=factor(quant, levels=c(" ", "Fixed delay at 12 wks", "Fixed delay at 24 wks",  
                                      "0.90 quantile", "0.95 quantile", "0.99 quantile"))) %>%
  mutate(log_score_clean=ifelse(log_score <= log_thresh, log_thresh, log_score), 
         failed=ifelse(log_score <= log_thresh, 1, 0)) %>%
         #failed=ifelse(mean(log_score) <= log_thresh, 1, 0)) %>%
  ungroup() %>%
  mutate(win=as.numeric(win), 
         delay_mod=ifelse(grepl("Fixed", quant), "fixed", 
                          ifelse(grepl("quantile", quant), "dynamic", "default")))
save(all_nobbs, file="Data/all_nobbs.Rdata")

per_rescore <-all_nobbs %>%
  dplyr::select(failed, state, multi, quant, model_parms, model_end_date, model_run, method, onset_wk) %>%
  ungroup()

ids_to_score <-per_rescore  %>%
  group_by(state, model_end_date, onset_wk) %>%  
  mutate(total=n(), 
         recodes=sum(failed),
         per=recodes/total) %>%
  filter(recodes==0) %>%
  mutate(score_id=paste0(state, ".", onset_wk, ".", model_end_date))

all_nobbs_scored <-all_nobbs %>%
  mutate(score_id=paste0(state, ".", onset_wk, ".", model_end_date)) %>%
  filter(failed==0, 
         score_id %in% ids_to_score$score_id) 

1-nrow(all_nobbs_scored)/nrow(all_nobbs) # percent failed
save(all_nobbs_scored, file="Data/all_nobbs_scored.Rdata")

##############################################################################################
##############################################################################################

n_for_cast <-nobbs_default %>% #Counts at time of nowcasts
  mutate(id=paste0(model_end_date, "_", model_parms, "_", state),
         id2=paste0(state, "_", model_end_date, "_", onset_wk),
         n.reported=ifelse(is.na(n.reported), 0, n.reported)) %>%
  rename(n_at_cast=n.reported) %>% 
  distinct() %>% 
  group_by(id) %>%
  do(tail(., n=3)) %>%
  ungroup() %>%
  select(id2, n_at_cast)

all_nobbs <-all_nobbs %>%
  mutate(id2=paste0(state, "_", model_end_date, "_", onset_wk),
         dynamic=ifelse(model_parms=="default selection", "default parameters", "dynamic parameters")) %>%
  left_join(., n_for_cast, by=c("id2")) %>%
  mutate(id3=sub("_[^_]+$", "", id2))

all_nobbs_mean_withF <-all_nobbs %>%
  mutate(range_50=q75-q25, 
         range_95=q97.5-q2.5) %>%
  group_by(state, model_parms, multi, quant, delay, model_end_date) %>%
  summarize(cov_95_mean=mean(cov95), 
            pe_full=mean(pe_full), 
            pe_part=mean(pe_part), 
            mean_ls=mean(log_score),
            mean_count=mean(n_at_cast),
            range_50=mean(range_50), 
            range_95=mean(range_95),
            range_50_cases=(range_50/(mean_count)) ,
            range_95_cases=(range_95/(mean_count)) ) %>%
  group_by(state) %>%
  mutate(failed_log_score=ifelse(mean_ls <= log_thresh, "Log scores <= -9.2", "Log scores > -9.2"),
         failed_cov=ifelse(cov_95_mean <= 0.95, "Failed", "Not failed"))
save(all_nobbs_mean_withFL, file="Data/all_nobbs_mean_withFL.Rdata")

all_nobbs_mean <-all_nobbs_mean_withFL %>%
  filter(state!="FL")
save(all_nobbs_mean, file="Data/all_nobbs_mean.Rdata")
