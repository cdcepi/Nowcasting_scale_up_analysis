### Find delays
find_delay_params <-function(df){
      
     parms <-df %>%
            filter(report_wk == max(report_wk)) %>%
            filter(delay_wk >= 0) %>%
            summarise(quants = cutpoints, 
                      delays = round(quantile(delay_wk, cutpoints))) %>% 
            mutate(delays=ifelse(delays<delay_min, delay_min, 
                                 ifelse(delays<delay_max, delay_max, delays)), 
                   delays=delays*delay_mult,
                   winds=(delays*wind_mult) + offset)
     return(parms)
}

### Function for recoding window periods to NULL
recode_null_function <-function(x) {if (!is.na(x)) return(x) else return(NULL)}

### Function for nowcast with dynamic params
now_cast <-function(df, delay, win) {
  
  win <-recode_null_function(win)
  delay <-recode_null_function(delay) 
  
  data=df
  now = max(df$report_wk) 
  units = "1 week"
  onset_date = "onset_wk"
  max_D=delay
  report_date = "report_wk"
  moving_window=win
  cutoff_D=T
  proportion_reported=1
  
  quiet=TRUE
  specs=list(
    dist=c("NB"),
    alpha1.mean.prior=0,
    alpha1.prec.prior=0.001,
    alphat.shape.prior=0.001,
    alphat.rate.prior=0.001,
    beta.priors=NULL,
    param_names=NULL,
    conf=0.95,
    dispersion.prior=NULL,
    nAdapt=3500, #3500,
    nChains=1,
    nBurnin=1000,
    nThin=1,
    nSamp=10000)
  
  # Check that "now" is entered as a Date
  if(inherits(now, "Date")==FALSE){
    stop("'Now' argument must be of datatype Date (as.Date)")
  }
  
  # Check that "now" is possible in the sequence of reporting data
  if(dplyr::last(seq(unique(data[,onset_date])[1],now,by=units))!=now){
    stop("The date `now` is not possible to estimate: the possible nowcast dates are seq(unique(data[,onset_date])[1],now,by=units).")
  }
  
  # Print date
  message(paste("Computing a nowcast for ",now))
  # Define "T", the length of dates between the first date of data and "now", making sure that "T" is unaffected by skipped-over dates in the time series
  # If the moving window is specified, "T" considers only the dates within the moving window; otherwise considers all historical data
  now.T <- ifelse(is.null(moving_window),length(seq(min(data[,onset_date]),as.Date(now),by=units)),
                  moving_window)
  
  # Check the default arguments
  if (is.null(moving_window)) {
    moving_window <- now.T
  }
  if (is.null(max_D)) {
    max_D <- now.T-1 # ifelse(is.null(moving_window),now.T-1,moving_window-1)
  }
  if (is.null(cutoff_D)) {
    cutoff_D <- TRUE
  }
  if(quiet==TRUE){
    progress.bar <- "none"
  }
  if(quiet==FALSE){
    progress.bar <- "text"
  }
  
  # Check that proportion_reported is between 0,1
  if (proportion_reported > 1 | proportion_reported<=0){
    stop("The proportion_reported must be a number between (0,1].")
  }
  
  # Manipulate the control arguments
  if ("Poisson"%in%(specs[["dist",exact=TRUE]])) { # if no distribution specified, take Poisson as default
    specs$dist <- "Poisson"
  }
  if (is.null(specs[["dist",exact=TRUE]])) {
    specs$dist <- "Poisson"
  }
  if (is.null(specs[["alpha1.mean.prior",exact=TRUE]])) {
    specs$alpha1.mean.prior <- 0
  }
  if (is.null(specs[["alpha1.prec.prior",exact=TRUE]])) {
    specs$alpha1.prec.prior <- 0.001
  }
  if (is.null(specs[["alphat.shape.prior",exact=TRUE]])) {
    specs$alphat.shape.prior <- 0.001
  }
  if (is.null(specs[["alphat.rate.prior",exact=TRUE]])) {
    specs$alphat.rate.prior <- 0.001
  }
  if (is.null(specs[["beta.priors",exact=TRUE]])) {
    specs$beta.priors <- rep(0.1, times=(max_D)+1)
  }
  if (is.null(specs[["param_names",exact=TRUE]])&(specs[["dist"]]=="Poisson")) {
    specs$param_names <- c( "lambda","alpha","beta.logged","tau2.alpha","sum.n")
  }
  if (is.null(specs[["param_names",exact=TRUE]])&(specs[["dist"]]=="NB")) {
    specs$param_names <- c( "lambda","alpha","beta.logged","tau2.alpha","sum.n","r")
  }
  if (is.null(specs[["conf",exact=TRUE]])) {
    specs$conf <- 0.95
  }
  if (is.null(specs[["dispersion.prior",exact=TRUE]])&(specs[["dist"]]=="NB")) {
    specs$dispersion.prior <- c(0.001,0.001)
  }
  if (is.null(specs[["nAdapt",exact=TRUE]])) {
    specs$nAdapt <- 1000
  }
  if (is.null(specs[["nChains",exact=TRUE]])) {
    specs$nChains <- 1
  }
  if (is.null(specs[["nBurnin",exact=TRUE]])) {
    specs$nBurnin <- 1000
  }
  if (is.null(specs[["nThin",exact=TRUE]])) {
    specs$nThin <- 1
  }
  if (is.null(specs[["nSamp",exact=TRUE]])) {
    specs$nSamp <- 10000
  }
  
  # Warnings
  if(max_D>(moving_window-1)){
    stop("Maximum delay cannot be greater than the length of the moving window minus 1 time unit")
  }
  
  # Prep the data: filter only to observable cases reported at or before "now"
  unit.num <- switch(units, "1 day"=1,"1 week"=7)
  w.days <- max((moving_window-1)*unit.num,(now.T-1)*unit.num) # moving window converted to days
  
  realtime.data <- subset(data,(data[,onset_date]<=now) & (data[,onset_date]>=now-w.days) & (data[,report_date]<=now) & (data[,report_date]>=now-w.days))
  realtime.data$week.t <- (as.numeric(realtime.data[,onset_date]-min(realtime.data[,onset_date]))/unit.num)+1
  realtime.data$delay <- as.numeric(realtime.data[,report_date]-realtime.data[,onset_date])/unit.num
  
  if(cutoff_D==FALSE){
    realtime.data$delay <- ifelse(realtime.data$delay>=max_D,max_D,realtime.data$delay)
  }
  
  if(length(unique(realtime.data$week.t))!=now.T){
    warning("Warning! The line list has zero case reports for one or more possible onset dates at one or more delays. Proceeding under the assumption that the true number of cases at the associated delay(s) and week(s) is zero.")
  }
  
  # Build the reporting triangle, fill with NAs where unobservable
  reporting.triangle <- matrix(NA, nrow=now.T,ncol=(max_D+1))
  
  for(t in 1:now.T){
    for(d in 0:max_D){
      reporting.triangle[t,(d+1)] <- nrow(realtime.data[which(realtime.data$week.t==t & realtime.data$delay==d),])
      if(now.T < (t+d)){
        reporting.triangle[t,(d+1)] <- NA
      }
    }
  }
  
  # Run the JAGS model
  
  if(specs[["dist"]]=="Poisson"){
    params=c( "lambda","alpha","beta.logged","tau2.alpha","n","sum.n","sum.lambda")
  }
  if(specs[["dist"]]=="NB"){
    params=c( "lambda","alpha","beta.logged","tau2.alpha","n","sum.n","sum.lambda","r")
  }
  nAdapt = specs[["nAdapt"]] #default = 1000
  nChains = specs[["nChains"]] # default=1
  nBurnin = specs[["nBurnin"]] # default=1000
  nThin = specs[["nThin"]] # default=1
  nKeep = specs[["nSamp"]] # default=10,000
  nIter = nKeep * nThin
  
  if(specs[["dist"]]=="Poisson"){
    dataList = list(Today = now.T,
                    D = max_D,
                    n = reporting.triangle,
                    alpha1.mean.prior=specs$alpha1.mean.prior,
                    alpha1.prec.prior=specs$alpha1.prec.prior,
                    alphat.rate.prior=specs$alphat.rate.prior,
                    alphat.shape.prior=specs$alphat.shape.prior,
                    beta.priors=specs$beta.priors)
  }
  
  if(specs[["dist"]]=="NB"){
    dataList = list(Today = now.T,
                    D = max_D,
                    n = reporting.triangle,
                    alpha1.mean.prior=specs$alpha1.mean.prior,
                    alpha1.prec.prior=specs$alpha1.prec.prior,
                    alphat.rate.prior=specs$alphat.rate.prior,
                    alphat.shape.prior=specs$alphat.shape.prior,
                    beta.priors=specs$beta.priors,
                    dispersion.prior.shape=specs$dispersion.prior[1],
                    dispersion.prior.rate=specs$dispersion.prior[2])
  }
  
  JAGSmodPois <-(("Code/nowcastPois.txt")) #system.file("JAGS", "nowcastPois.txt", package="NobBS")) # file.path(path.package('NobBS'),"nowcastPois.txt")
  JAGSmodNB <- (("Code/nowcastNB.txt")) #system.file("JAGS", "nowcastNB.txt", package="NobBS")) #file.path(path.package('NobBS'),"nowcastNB.txt")
  
  nowcastmodel = jags.model(
    file = ifelse(specs[["dist"]]=="Poisson",JAGSmodPois,JAGSmodNB),
    data = dataList,
    n.chains = nChains,
    n.adapt = nAdapt,
    inits=list(.RNG.seed=1,.RNG.name="base::Super-Duper"),
    quiet=quiet)
  
  update( object = nowcastmodel, n.iter = nBurnin , progress.bar = progress.bar)
  
  lambda.output = coda.samples(
    model = nowcastmodel,
    variable.names =  if("sum.n"%in%specs$param_names) c(specs$param_names) else c(specs$param_names,"sum.n"),
    n.iter = nIter,
    thin = nThin,
    quiet=quiet,
    progress.bar=progress.bar)
  
  mymod.mcmc <- as.mcmc(lambda.output)
  #mymod.dat <- as.data.frame(as.matrix(mymod.mcmc))
  
  #return(mymod.dat)#$estimates)
  #return(nowcastmodel)
  return(mymod.mcmc)
}

### Function for extracting last two weeks of each nowcast
extract_dat <-function(dat, df){
  
  r <-dat %>% 
    dplyr::select(starts_with("r")) %>%
    colMeans() 
  
  dat_preds <-dat %>% 
        dplyr::select(starts_with("sum")) 
  
  onset_wks <-seq(max(df$report_wk)-(ncol(dat_preds)-1)*7, max(df$report_wk), by="week")
 
  dat_preds %>%
    gather("variable", "value") %>% 
    mutate(variable = gsub('[^[:alnum:] ]', "", variable),
           variable = as.numeric(gsub('sumn', "", variable))) %>%
    group_by(variable) %>% 
    arrange(desc(variable)) %>%
    summarize(mean = mean(value),
              median = median(value),
              q1 = quantile(value, probs = 0.01),
              q99 = quantile(value, probs = 0.99),
              q2.5 = quantile(value, probs = 0.025),
              q97.5 = quantile(value, probs = 0.975),
              q25 = quantile(value, probs = 0.25),
              q75 = quantile(value, probs = 0.75)) %>%
    dplyr::select(-variable) %>%
    bind_cols(., onset_wks, r) %>%
    rename(onset_wk=9,
           dispersion=10)
  
}
