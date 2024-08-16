library(tidyverse)
library(here)

load(here("Data/nowcaster.Rdata"))
nowcaster <-nowcaster %>%
  filter(win!="NULL") %>%
  select(state, onset_wk) %>%
  distinct() 

load(here("Data/dat_COVID.Rdata"))
dat_COVID <-dat_COVID %>%
  filter(delay_wk >= 0) %>%
  arrange(report_dt) 
states <-sort(unique(dat_COVID$state))
states <-states[2:3]

state_counts <-list()
pred_weeks <-list()

end_dates <-list()
seq_dates <-list()
for(h in 1:length(states)){
  state_counts[[h]] <-dat_COVID %>% 
    filter(state==states[h]) %>%
    mutate(in_wk=NA)
  
  pred_weeks[[h]] <-nowcaster %>%
    filter(state==states[h])
  
  end_dates[[h]] <-sort(unique(nowcaster$onset_wk))
  end_dates[[h]] <-append(min(nowcaster$onset_wk) - 7, end_dates[[h]])
  end_dates[[h]] <-append(min(nowcaster$onset_wk) - 7*2, end_dates[[h]])
  
  state_counts[[h]] <-state_counts[[h]]%>%
    filter(onset_dt <=max(end_dates[[h]]))
  
  seq_dates[[h]] <-list()

  for(i in 1:length(end_dates[[h]])) {
    seq_dates[[h]][[i]] <-seq(as.Date(end_dates[[h]][[i]]-6), as.Date(end_dates[[h]][[i]]), by="day")
  }
}

for(h in 1:length(states)){ 
  for(i in 1:length(end_dates[[h]])) {
    for(j in 1:nrow(state_counts[[h]])){
      if(state_counts[[h]][j, "onset_dt"]  %in% seq_dates[[h]][[i]]) 
        state_counts[[h]][j, "in_wk"] <-end_dates[[h]][[i]]
      }
  }
  }
   
for(h in 1:length(states)){ 
  state_counts[[h]]$in_wk <-as.Date(state_counts[[h]]$in_wk, origin='1970-01-01')
  
  state_counts[[h]] <-state_counts[[h]] %>%
    group_by(in_wk) %>%
    count()
}

reports_nowcaster_dates <-do.call(rbind.data.frame, state_counts) 
save(reports_nowcaster_dates, file="./Data/reports_nowcaster_dates.Rdata")