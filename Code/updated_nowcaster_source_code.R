### Function to set up dataframe with no age variable
data.w_no_age<-function(dataset,
                        trim.data,
                        date_onset,
                        date_report,
                        K=0,
                        silent=F){
  if(!silent){
    ## Last digitation date considered
    if(missing(trim.data)){
      trim.data <-  0
      warning("Using default, no trimming out of the data")
    } else {
      warning("Using default, trimming out ",
              trim.data ,
              " days of data",
              call. = T)
    }
  }else{
    ## Trim.data
    trim.data <-  0
  }
  
  ## Transforming trim.data into weeks
  trim.data.w <- 7*trim.data
  
  ## K parameter of forecasting
  K.w<-7*K
  
  ## Maximum date to be considered on the estimation
  DT_max <- max(dataset |>
                  dplyr::pull(var = {{date_report}}),
                na.rm = T) - trim.data.w + K.w
  
  ## Last day of the week for the digitation date calculation
  DT_max_diadasemana <- as.integer(format(DT_max, "%w"))
  
  
  ## Accounting for the maximum of days on the last week to be used
  data_w <- dataset |>
    dplyr::rename(date_report = {{date_report}},
                  date_onset = {{date_onset}}) |>
    dplyr::filter(date_report <= DT_max ) |>
    dplyr::mutate(
      ## Altering the date for the first day of the week
      date_onset = date_onset -
        as.integer(format(date_onset, "%w")) -
        (6-DT_max_diadasemana),
      date_report = date_report -
        as.integer(format(date_report, "%w")) -
        (6-DT_max_diadasemana),
      Delay = as.numeric(date_report - date_onset) / 7) |>
    dplyr::filter(Delay >= 0)
  
  # Returning the data
  return(data_w)
}


### Model function
nowcasting_no_age <- function(dataset,
                              zero_inflated=FALSE){
  
  ## Check for the zero-inflated
  if (zero_inflated){
    family <- "zeroinflatednbinomial1"
    control.family <- list(
      hyper = list("theta1" = list(prior = "loggamma",
                                   param = c(0.01, 0.01)),
                   "theta2" = list(prior = "gaussian",
                                   param = c(0, 0.4)))
    )
  } else {
    family <- 'nbinomial'
    control.family <- list(
      hyper = list("theta" = list(prior = "loggamma",
                                  param = c(0.001, 0.001)))
    )
  }
  
  index.missing <- which(is.na(dataset$Y))
  
  ## Model equation: intercept + f(time random effect) + f(Delay random effect)
  ## Y(t) ~ 1 + rw2(t) + rw1(delay),
  ## prec(rw2) ~ logGamma(10e-3, 10e-3), prec(rw1) ~ logGamma(10e-3, 10e-3)
  model <- Y ~ 1 +
    f(Time,
      model = "rw2",
      hyper = list("prec" = list(prior = "loggamma",
                                 param = c(0.001, 0.001))
      )) +
    f(delay, model = "rw1",
      hyper = list("prec" = list(prior = "loggamma",
                                 param = c(0.001, 0.001)))
    )
  
  ## Running the Negative Binomial model in INLA
  output0 <- INLA::inla(model,
                        family = family,
                        data = dataset,
                        control.predictor = list(link = 1, compute = T),
                        control.compute = list( config = T, waic=F, dic=F),
                        control.family = control.family
  )
  
  # }
  
  ## Algorithm to get samples for the predictive distribution for the number of cases
  
  ## Step 1: Sampling from the approximate posterior distribution using INLA
  srag.samples0.list <- INLA::inla.posterior.sample(n = 1000, output0)
  
  ## Give a parameter to trajectories, 
  
  ## Step 2: Sampling the missing triangle from the likelihood using INLA estimates
  vector.samples0 <- lapply(X = srag.samples0.list,
                            FUN = function(x, idx = index.missing){
                              if(zero_inflated){
                                unif.log <- as.numeric(runif(idx,0,1) < x$hyperpar[2])
                              }else{
                                unif.log = 1
                              }
                              stats::rnbinom(n = idx,
                                             mu = exp(x$latent[idx]),
                                             size = x$hyperpar[1]
                              ) * unif.log
                            } )
  
  ## Step 3: Calculate N_{a,t} for each triangle sample {N_{t,a} : t=Tactual-Dmax+1,...Tactual}
  
  gg.age <- function(x, dados.gg, idx){
    data.aux <- dados.gg
    Tmin <- min(dados.gg$Time[idx])
    data.aux$Y[idx] <- x
    data.aggregated <- data.aux |>
      ## Selecionando apenas os dias faltantes a partir
      ## do domingo da respectiva ultima epiweek
      ## com dados faltantes
      dplyr::filter(Time >= Tmin  ) |>
      dplyr::group_by(Time, dt_event) |>
      dplyr::summarise(
        Y = sum(Y), .groups = "keep"
      )
    data.aggregated
  }
  
  ## Step 4: Applying the age aggregation on each posterior
  tibble.samples.0 <- lapply( X = vector.samples0,
                              FUN = gg.age,
                              dados = dataset,
                              idx = index.missing)
  
  srag.pred.0 <- dplyr::bind_rows(tibble.samples.0, .id = "sample")
  
  return(srag.pred.0)
  
}

### Test function for extracting dispersion
disper_test <- function(dataset,
                              zero_inflated=FALSE){
  
  ## Check for the zero-inflated
  if (zero_inflated){
    family <- "zeroinflatednbinomial1"
    control.family <- list(
      hyper = list("theta1" = list(prior = "loggamma",
                                   param = c(0.01, 0.01)),
                   "theta2" = list(prior = "gaussian",
                                   param = c(0, 0.4)))
    )
  } else {
    family <- 'nbinomial'
    control.family <- list(
      hyper = list("theta" = list(prior = "loggamma",
                                  param = c(0.001, 0.001)))
    )
  }
  
  index.missing <- which(is.na(dataset$Y))
  
  ## Model equation: intercept + f(time random effect) + f(Delay random effect)
  ## Y(t) ~ 1 + rw2(t) + rw1(delay),
  ## prec(rw2) ~ logGamma(10e-3, 10e-3), prec(rw1) ~ logGamma(10e-3, 10e-3)
  model <- Y ~ 1 +
    f(Time,
      model = "rw2",
      hyper = list("prec" = list(prior = "loggamma",
                                 param = c(0.001, 0.001))
      )) +
    f(delay, model = "rw1",
      hyper = list("prec" = list(prior = "loggamma",
                                 param = c(0.001, 0.001)))
    )
  
  ## Running the Negative Binomial model in INLA
  output0 <- INLA::inla(model,
                        family = family,
                        data = dataset,
                        control.predictor = list(link = 1, compute = T),
                        control.compute = list( config = T, waic=F, dic=F),
                        control.family = control.family
  )
  
  converge <-output0$mode$mode.status
  
  hypers <-output0$summary.hyperpar %>%
    rownames_to_column("parameter") %>%
    bind_cols(., converge) %>%
    rename(mode.status=8)
  
  return(hypers)
  
}

### Model output function
nowcasting.summary <- function(trajetory, age = F){
  
  total.summy <- trajetory |>
    dplyr::group_by(Time, dt_event, sample) |>
    dplyr::summarise(Y = sum(Y, na.rm = T)) |>
    dplyr::group_by(Time, dt_event) |>
    dplyr::summarise(Mean = mean(Y, na.rm = T), #added
                     Median = stats::median(Y, na.rm = T),
                     LI = stats::quantile(Y, probs = 0.025, na.rm = T),
                     LS = stats::quantile(Y, probs = 0.975, na.rm = T),
                     LIb = stats::quantile(Y, probs = 0.25, na.rm = T),
                     LSb = stats::quantile(Y, probs = 0.75, na.rm = T),
                     q01 = stats::quantile(Y, probs = 0.01, na.rm = T),
                     q99 = stats::quantile(Y, probs = 0.99, na.rm = T),
                     .groups = "drop")

  
  if(age){
    age.summy <- trajetory |>
      dplyr::group_by(Time, dt_event, fx_etaria, fx_etaria.num) |>
      dplyr::summarise(Mean = mean(Y, na.rm = T), #added,
                       Median = stats::median(Y, na.rm = T),
                       LI = stats::quantile(Y, probs = 0.025, na.rm = T),
                       LS = stats::quantile(Y, probs = 0.975, na.rm = T),
                       LIb = stats::quantile(Y, probs = 0.25, na.rm = T),
                       LSb = stats::quantile(Y, probs = 0.75, na.rm = T),
                       q01 = stats::quantile(Y, probs = 0.01, na.rm = T),
                       q99 = stats::quantile(Y, probs = 0.99, na.rm = T),
                       .groups = "drop")
    
    output <- list()
    output$total <- total.summy
    output$age <- age.summy
    
  }else{
    output<- list()
    output$total <- total.summy
    
  }
  
  
  return(output)
}

### Nowcasting call
nowcasting_inla <- function(dataset,
                            bins_age="SI-PNI",
                            trim.data=0,
                            Dmax = 15,
                            wdw = 30,
                            age_col,
                            date_onset,
                            date_report,
                            data.by.week = FALSE,
                            # return.age = NULL,
                            silent = F,
                            K = 0,
                            trajectories = F,
                            zero_inflated = F,
                            ...){
  
  dots<-list(...)
  
  ## Safe tests
  if(missing(dataset)){
    stop("Dataset is missing!")
  }
  if(ncol(dataset) < 2){
    if(!missing(age_col) & ncol(dataset) < 3){
      stop("Missing 'age_col'! \n
           Dataset does not have 3 columns!")
    }else{
      stop("Dataset does not have 2 columns!")
    }
    
  }
  if(missing(date_onset) | missing(date_report)){
    stop("date_onset or date_report missing! Please give a column name for each of this parameters")
  }
  if(K < 0 ){
    stop("K less than 0, we cannot produce backcasting! \n
         Please set the K to anything greater than 0 to Forecasting")
  }
  
  ## Warnings
  if(missing(silent) | silent == FALSE){
    
    ## K parameter and trim.data warnings
    if (K > 0 & trim.data != 0){
      warning(paste0("Using K = ", K, " and trim.data = ", trim.data, ", is that right?"))
    }else{
      if(K > 0){
        message(paste0("Forecasting with K ", K, " ahead"))
      }else{
        message("Nowcasting only")
      }
    }
    ## Missing age column warning
    if(missing(bins_age)){
      bins_age <- "SI-PNI"
      warning("Using 'SI-PNI' age bins!")
    }else{
      bins_age<-bins_age
      message("Using age bins inputed")
    }
    ## Missing trim.data warning
    if(missing(trim.data)){
      trim.data <- 0
      warning("Using default to trim dates, trim.data = 0")
    }else{
      trim.data<-trim.data
      message("Using trim.data inputed")
    }
    ## Missing Dmax warning
    if(missing(Dmax)){
      Dmax <- 15
      warning("Using default to maximum delay, Dmax = 15")
    }else{
      Dmax<-Dmax
      message("Using Dmax inputed")
    }
    ## Missing wdw warning
    if(missing(wdw)){
      wdw <- 30
      warning("Using default to window of action, wdw = 30")
    }else{
      wdw<-wdw
      message("Using wdw inputed")
    }
    ## Missing data.by.week warning
    if(missing(data.by.week)){
      data.by.week <- FALSE
      warning("Using default to returning option for the data, data.by.week = FALSE")
    }else{
      data.by.week<-data.by.week
      message("Returning data.by.week")
    }
    # ## Missing return.age warning
    # if(missing(return.age)){
    #   return.age <- TRUE
    #   warning("Using default to returning estimate by age, return.age = TRUE")
    # }
    ## Missing age_col warning
    if(missing(age_col)){
      warning("Age_col missing, nowcasting with unstructured model")
    }else{
      message("Age col inputed, nowcasting with structured model")
    }
    
    if(missing(trajectories) | trajectories == FALSE){
      warning("Not returning trajectories")
    }else{
      trajectories = TRUE
      message("Trajectories returned")
      
      if(!missing(age_col) & !missing(zero_inflated)){
        zero_inflated<-FALSE
        warning("age_col parsed, zero_inflated ignored!")
      }
    }
  }
  # else{
  #     bins_age<-bins_age;
  #     trim.data<-trim.data;
  #     Dmax<-Dmax;
  #     wdw<-wdw;
  #     data.by.week<-data.by.week;
  #     zero_inflated<-zero_inflated
  #   }
  
  
  
  ## Objects for keep the nowcasting
  ## Filtering out cases without report date
  if(missing(age_col)){
    data<-dataset |>
      dplyr::select({{date_report}}, {{date_onset}})  |>
      tidyr::drop_na({{date_report}})
  } else {
    data <- dataset  |>
      dplyr::select({{date_report}}, {{date_onset}}, {{age_col}})  |>
      tidyr::drop_na({{date_report}})
  }
  
  ## Filtering data to the parameters setted above
  if(missing(age_col)){
    data_w<-data.w_no_age(dataset = data,
                          trim.data = trim.data,
                          date_onset = {{date_onset}},
                          date_report = {{date_report}},
                          K = K,
                          silent = silent)
  }else {
    data_w <- data.w(dataset = data,
                     bins_age = bins_age,
                     trim.data = trim.data,
                     age_col = {{age_col}},
                     date_onset = {{date_onset}},
                     date_report = {{date_report}},
                     K = K,
                     silent = silent)
  }
  
  ## Parameters of Nowcasting estimate
  Tmax <- max(data_w |>
                dplyr::pull(var = date_onset))
  
  ## Data to be entered in Nowcasting function
  ##
  if(missing(age_col)){
    data.inla <- data_w  |>
      ## Filter for dates
      dplyr::filter(date_onset >= Tmax - 7 * wdw,
                    Delay <= Dmax)  |>
      ## Group by on Onset dates, Amounts of delays and Stratum
      dplyr::group_by(date_onset, delay = Delay)  |>
      ## Counting
      dplyr::tally(name = "Y")  |>
      dplyr::ungroup()
  } else {
    data.inla <- data_w  |>
      ## Filter for dates
      dplyr::filter(date_onset >= Tmax - 7 * wdw,
                    Delay <= Dmax)  |>
      ## Group by on Onset dates, Amounts of delays and Stratum
      dplyr::group_by(date_onset, delay = Delay, fx_etaria)  |>
      ## Counting
      dplyr::tally(name = "Y")  |>
      dplyr::ungroup()
  }
  
  
  ## Auxiliary date table
  if(K==0){
    dates<-unique(data.inla |>
                    dplyr::pull(var = date_onset - 7*trim.data))
  } else {
    ## This is done to explicitly say for the forecast part that its date of onset is the present date
    date_k<-(max(data.inla$date_onset + 7*K - 7*trim.data))
    dates<-c(unique(data.inla$date_onset), date_k)
  }
  
  ## To make an auxiliary date table with each date plus an amount of dates  to forecast
  tbl.date.aux <- tibble::tibble(
    date_onset = dates
  )  |>
    tibble::rowid_to_column(var = "Time")
  
  ## Joining auxiliary date tables
  data.inla <- data.inla  |>
    dplyr::left_join(tbl.date.aux)
  
  ## Time maximum to be considered
  Tmax.id <- max(data.inla$Time)
  
  # Auxiliary date table on each stratum, By age
  if(missing(age_col)){
    tbl.NA <-
      expand.grid(Time = 1:(Tmax.id+K),
                  delay = 0:Dmax)  |>
      dplyr::left_join(tbl.date.aux, by = "Time")
  } else{
    tbl.NA <-
      expand.grid(Time = 1:(Tmax.id+K),
                  delay = 0:Dmax,
                  fx_etaria = unique(data.inla$fx_etaria)
      ) |>
      dplyr::left_join(tbl.date.aux, by = "Time")
  }
  
  ## Joining the auxiliary date table by Stratum
  if(missing(age_col)){
    data.inla <- data.inla  |>
      dplyr::full_join(tbl.NA) |>  #View()
      dplyr::mutate(
        Y = ifelse(Time + delay > Tmax.id, as.numeric(NA), Y),
        ## If Time + Delay is greater than Tmax, fill with NA
        Y = ifelse(is.na(Y) & Time + delay <= Tmax.id, 0, Y ),
        ## If Time + Delay is smaller than Tmax AND Y is NA, fill 0
      )  |>
      dplyr::arrange(Time, delay) |>
      dplyr::rename(dt_event = date_onset) |>
      tidyr::drop_na(delay)
  }else {
    data.inla <- data.inla  |>
      dplyr::full_join(tbl.NA)  |>   #View()
      dplyr::mutate(
        Y = ifelse(Time + delay > Tmax.id, as.numeric(NA), Y),
        ## If Time + Delay is greater than Tmax, fill with NA
        Y = ifelse(is.na(Y) & Time + delay <= Tmax.id, 0, Y ),
        ## If Time + Delay is smaller than Tmax AND Y is NA, fill 0
      )  |>
      dplyr::arrange(Time, delay, fx_etaria) |>
      dplyr::rename(dt_event = date_onset) |>
      tidyr::drop_na(delay)
  }
  ## Precisamos transformar essa datas de volta no valor que é correspondente delas,
  ## a ultima data de primeiro sintomas foi jogada pra até uma semana atrás
  
  
  if(missing(age_col)){
    
    if(zero_inflated){
      ## Nowcasting estimate
      sample.now <- nowcasting_no_age(dataset = data.inla,
                                      zero_inflated = T)
      hyper_parms <-list(disper_test(dataset = data.inla,
                                 zero_inflated = T))
      names(hyper_parms) <-"hyper_parms"

    }else{
      ## Nowcasting estimate
      sample.now <- nowcasting_no_age(dataset = data.inla,
                                      zero_inflated = F)
      hyper_parms <-list(disper_test(dataset = data.inla,
                                      zero_inflated = F))
      names(hyper_parms) <-"hyper_parms"
    }
    
    ## Summary on the posteriors of nowcasting
    now_summary<-nowcasting.summary(trajetory = sample.now,
                                    age = F)
    
    l<-1
  } else {
    
    if(zero_inflated){
      ## Nowcasting estimate
      sample.now <- nowcasting_age(dataset = data.inla,
                                   zero_inflated = T)
    }else{
      ## Nowcasting estimate
      sample.now <- nowcasting_age(dataset = data.inla,
                                   zero_inflated = F)
    }
    
    ## Summary on the posteriors of nowcasting
    now_summary<-nowcasting.summary(trajetory = sample.now,
                                    age = T)
    
    l<-0
  }
  
  ## Objects to be returned
  
  if(data.by.week){
    
    now_summary[[3-l]]<-data_w |>
      dplyr::group_by(date_onset) |>
      dplyr::summarise(observed = dplyr::n(),
                       Delay = Delay)
    
    names(now_summary)[3-l]<-"data"
    
    if(trajectories){
      now_summary[[4-l]]<-sample.now
      names(now_summary)[4-l]<-"trajectories"
      
    }
  } else {
    if(trajectories){
      now_summary[[3-l]]<-sample.now
      names(now_summary)[3-l]<-"trajectories"
      
      
    }
  }
  
  ## Final object returned
  
  now_summary <-append(now_summary, hyper_parms) 
  
  return(now_summary)
  
}
