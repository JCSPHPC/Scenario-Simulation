#
# Successful control scenarios as 95% of simulation to control an epidemic.
#
# Edward McNeil
# November 2022
#

# Load required libraries
library(knitr)
library(kableExtra)

# Set the global paths for data and output
Path.data <- "F:/Edward/Kwok/Simulation/data"
Path.out <- "F:/Edward/Kwok/Simulation/output"

folders <- dir(Path.data)

# ff is the index vector for the folders since SARS has a different optimal CTR value of 0.7 (scenario 8))
# The others all use 100% (scenario 11)
f_table <- function(ff) {
  
  Cases <- list()
  
  S1 <- list()
  S2 <- list()
  S3 <- list()
  S4 <- list()
  S5 <- list()
  S6 <- list()
  CFR <- list()
  IFR <- list()
  S8 <- list()
  S9 <- list()
  S10 <- list()
  
  i <- 0
  
  # Loop through the data folders
  for(f in folders[ff]) {
    i <- i + 1  
    
    # Change the global paths
    path.data <- paste(Path.data, f, sep="/")
    path.out <- paste(Path.out, f, sep="/")
    
    # Import the data
    setwd(path.data)
    cat("\nImporting data for", f, "...\n")
    
    source("F:/Edward/Kwok/Simulation/import.R")
    
    cat("\nEstimating results for", f, "...\n")
    source("F:/Edward/Kwok/Simulation/sim.plot.R")
    
    # Call the function (one with 95% CI and one without)
    if(f=="sars")
      times <- c(5, 10, 15, 20)
    else
      times <- c(5, 10, 15, 20, 40, 60, 80)
    
    Cases[[i]] <- Sim.plot(data, param='total_case', ci=TRUE, times=times, plot=FALSE)
    
    
    # Parameter estimates
    cat("\nEstimating parameters for", f, "...\n")
    
    # Load the required functions  
    source("F:/Edward/Kwok/Simulation/estimate.R")
    
    
    # Calculate the estimates
    
    # 1. Median (95% CI) DAILY number of individuals traced in each contact tracing strategy. I would like to check the capacity for the community quarantine centre. (this one is correct; Mean of NON-ZERO values in column of n_quarantine_events in each excel file then choose 50th, 2.5th and 97.5th from 100 excel files estimates) - average daily resources spent 
    
    quarantine_mean <- Mean(data, 'n_quarantine_events')
    S1[[i]] <- Out(data, quarantine_mean)
    
    
    # 2. Mean, maximum and minimum DAILY number of individuals traced in each contact tracing strategy (taking the mean within each simulation/excel file). I would like to check the capacity for the community quarantine centre.
    quarantine_sum <- Sum(data, 'n_quarantine_events')
    S2[[i]] <- Out(data, quarantine_sum, digits=0)
    
    
    # 3. The actual number of contact to be traced per case (formula = sum the column of n_quarantine_events OVER the last value of total_case variable) in given contact tracing strategies in the early phase of the future outbreak of novel pathogens with 95%CI by ranking 100 simulations and values of  97.5th and 2.5th
    
    total_case <- Max(data, 'total_case')
    
    D <- tibble(scenario=rep(1:nscen, times=unlist(lapply(data, length))),
                total_case=unlist(total_case),
                quarantine_sum=unlist(quarantine_sum),
                contacts.traced.per.case=quarantine_sum/total_case)
    
    fdig <- paste0("%.", 1, "f")
    
    S3[[i]] <- group_by(D, scenario) %>%
      summarise(min = sprintf(fdig, min(contacts.traced.per.case, na.rm=TRUE)),
                q.025 = sprintf(fdig, quantile(contacts.traced.per.case, probs=0.025, na.rm=TRUE)),
                q.50 = sprintf(fdig, quantile(contacts.traced.per.case, probs=0.5, na.rm=TRUE)),
                q.975 = sprintf(fdig, quantile(contacts.traced.per.case, probs=0.975, na.rm=TRUE)), 
                max = sprintf(fdig, max(contacts.traced.per.case, na.rm=TRUE)), .groups='drop')
    
    # 4. Median number of deaths in given contact tracing strategies in the early phase of the future outbreak of novel pathogens with 95%CI by ranking 100 simulations and values of  50th,  97.5th and 2.5th (using MAX of n_death in each excel file)
    
    total_death <- Max(data, 'total_death')
    S4[[i]] <- Out(data, total_death, digits=1)
    
    
    # 5. Maximum, minimum and mean number of DAILY hospital bed occupied in given contact tracing strategies in the early phase of the future outbreak of novel pathogens with 95%CI by ranking 100 simulations and values of 50th, 97.5th and 2.5th (using MAX, MIN and mean value of NON-ZERO value in n_hospital and obtain the 50th, 97.5th and 2.5th from 100 realization)
    
    n_hospital <- Mean(data, 'n_hospital')
    S5[[i]] <- Out(data, n_hospital, digits=1)
    
    
    # 6. Maximum, minimum and mean number of DAILY ICU admission in given contact tracing strategies in the early phase of the future outbreak of novel pathogens with 95%CI by ranking 100 simulations and values of  97.5th and 2.5th  (using MAX, MIN and mean value of NON-ZERO values in n_critical and obtain the 50th, 97.5th and 2.5th from 100 realization
    
    n_critical <- Mean(data, 'n_critical')
    S6[[i]] <-  Out(data, n_critical, digits=1)
    
    
    # 7.Median CFR/IFR given contact tracing strategies using i) total no. of death/ total no. of CASES ii) total no. of death/ total no. of infection Using the formula i)  MAX of n_death/MAX of total case ii)  MAX of n_death/MAX of total infected
    
    total_infected <- Max(data, 'total_infected')
    
    D <- tibble(scenario=rep(1:nscen, times=unlist(lapply(data, length))),
                total_case=unlist(total_case),
                total_death=unlist(total_death),
                total_infected=unlist(total_infected),
                CFR=total_death/total_case,
                IFR=total_death/total_infected)
    
    CFR[[i]] <- group_by(D, scenario) %>%
      summarise(min = sprintf(fdig, 100*min(CFR, na.rm=TRUE)),
                q.025 = sprintf(fdig, 100*quantile(CFR, probs=0.025, na.rm=TRUE)),
                q.50 = sprintf(fdig, 100*quantile(CFR, probs=0.5, na.rm=TRUE)),
                q.975 = sprintf(fdig, 100*quantile(CFR, probs=0.975, na.rm=TRUE)), 
                max = sprintf(fdig, 100*max(CFR, na.rm=TRUE)), .groups='drop')
    
    IFR[[i]] <- group_by(D, scenario) %>%
      summarise(min = sprintf(fdig, 100*min(IFR, na.rm=TRUE)),
                q.025 = sprintf(fdig, 100*quantile(IFR, probs=0.025, na.rm=TRUE)),
                q.50 = sprintf(fdig, 100*quantile(IFR, probs=0.5, na.rm=TRUE)),
                q.975 = sprintf(fdig, 100*quantile(IFR, probs=0.975, na.rm=TRUE)), 
                max = sprintf(fdig, 100*max(IFR, na.rm=TRUE)), .groups='drop')
    
    
    #8. Number of cases
    S8[[i]] <- Out(data, total_case, digits=0)
    
    
    #9. Number of infections
    S9[[i]] <- Out(data, total_infected, digits=0)
    
    #10. Median duration of containment
    options(warn = -1)
    
    durations <- lapply(data, function(x) 
      lapply(x, function(y) 
        y$time[1 + try(max(which(diff(y$total_case) !=0)), silent=TRUE)]))
    options(warn = 0)
    
    D <- tibble(scenario=rep(1:nscen, times=unlist(lapply(data, length))),
                duration=unlist(durations)) 
    
    fdig <- paste0("%.", 0, "f")
    
    S10[[i]] <- group_by(D, scenario) %>%
      summarise(min = sprintf(fdig, round(min(duration, na.rm=TRUE))),
                q.025 = sprintf(fdig, round(quantile(duration, probs=0.025, na.rm=TRUE))),
                q.50 = sprintf(fdig, round(quantile(duration, probs=0.5, na.rm=TRUE))),
                q.975 = sprintf(fdig, round(quantile(duration, probs=0.975, na.rm=TRUE))), 
                max = sprintf(fdig, round(max(duration, na.rm=TRUE))), .groups='drop')
    
    cat("\n")
  }
  
  
  foo <- function(x){
    X <- x[x$Times=="10 X",]
    i <- min(which(X$mean>=0.95))
    x$ctr[i]
  }
  
  sapply(Cases, foo)
  
  
  lodc.A <- sapply(Cases, function(x) 100 * x$mean[x$Times=="10 X" & x$ctr==1.1]); lodc.A
  lodc.B <- sapply(Cases, function(x) 100 * x$mean[x$Times=="10 X" & x$ctr==0]); lodc.B
  
  octr <- sapply(Cases, foo); octr
  
  if(length(ff)==1)
    sc <- 8 # ctr=0.7
  else
    sc <- 11 # ctr=1.0
  
  ncases <- sapply(S8, function(x) x$q.50[x$scenario==sc]); ncases
  ncases.lower <- sapply(S8, function(x) x$q.025[x$scenario==sc]); ncases.lower
  ncases.upper <- sapply(S8, function(x) x$q.975[x$scenario==sc]); ncases.upper
  
  ndeaths <- sapply(S4, function(x) x$q.50[x$scenario==sc]); ndeaths
  ndeaths.lower <- sapply(S4, function(x) x$q.025[x$scenario==sc]); ndeaths.lower
  ndeaths.upper <- sapply(S4, function(x) x$q.975[x$scenario==sc]); ndeaths.upper
  
  ncontacts <- sapply(S2, function(x) x$q.50[x$scenario==sc]); ncontacts
  ncontacts.lower <- sapply(S2, function(x) x$q.025[x$scenario==sc]); ncontacts.lower
  ncontacts.upper <- sapply(S2, function(x) x$q.975[x$scenario==sc]); ncontacts.upper
  
  ncpc <- sapply(S3, function(x) x$q.50[x$scenario==sc]); ncpc
  ncpc.lower <- sapply(S3, function(x) x$q.025[x$scenario==sc]); ncpc.lower
  ncpc.upper <- sapply(S3, function(x) x$q.975[x$scenario==sc]); ncpc.upper
  
  nhospital <- sapply(S5, function(x) x$q.50[x$scenario==sc]); nhospital
  nhospital.lower <- sapply(S5, function(x) x$q.025[x$scenario==sc]); nhospital.lower
  nhospital.upper <- sapply(S5, function(x) x$q.975[x$scenario==sc]); nhospital.upper
  
  nicu <- sapply(S6, function(x) x$q.50[x$scenario==sc]); nicu
  nicu.lower <- sapply(S6, function(x) x$q.025[x$scenario==sc]); nicu.lower
  nicu.upper <- sapply(S6, function(x) x$q.975[x$scenario==sc]); nicu.upper
  
  nduration.0 <- sapply(S10, function(x) x$q.50[x$scenario==1]); nduration.0
  nduration.0.lower <- sapply(S10, function(x) x$q.025[x$scenario==1]); nduration.0.lower
  nduration.0.upper <- sapply(S10, function(x) x$q.975[x$scenario==1]); nduration.0.upper
  
  nduration.1 <- sapply(S10, function(x) x$q.50[x$scenario==sc]); nduration.1
  nduration.1.lower <- sapply(S10, function(x) x$q.025[x$scenario==sc]); nduration.1.lower
  nduration.1.upper <- sapply(S10, function(x) x$q.975[x$scenario==sc]); nduration.1.upper
  
  cfr <- sapply(CFR, function(x) x$q.50[x$scenario==sc]); cfr
  cfr.lower <- sapply(CFR, function(x) x$q.025[x$scenario==sc]); cfr.lower
  cfr.upper <- sapply(CFR, function(x) x$q.975[x$scenario==sc]); cfr.upper
  
  
  #
  bar <- function(x, low, up) {
    paste0(x, " (", low, " - ", up, ")")
  }
  
  bar(ncases, ncases.lower, ncases.upper)
  
  Table <- data.frame(lodc.A,
                      lodc.B,
                      octr,
                      ncases=bar(ncases, ncases.lower, ncases.upper),
                      ndeaths=bar(ndeaths, ndeaths.lower, ndeaths.upper),
                      ncontacts=bar(ncontacts, ncontacts.lower, ncontacts.upper),
                      ncpc=bar(ncpc, ncpc.lower, ncpc.upper),
                      nhospital=bar(nhospital, nhospital.lower, nhospital.upper),
                      nicu=bar(nicu, nicu.lower, nicu.upper),
                      nduration.0=bar(nduration.0, nduration.0.lower, nduration.0.upper),
                      nduration.1=bar(nduration.1, nduration.1.lower, nduration.1.upper),
                      CFR=bar(cfr, cfr.lower, cfr.upper))
  
  Table1 <- t(Table)
  Table1
}

Table1.flu <- f_table(1:4)
Table1.sars <- f_table(5)

Table1 <- cbind(Table1.flu, Table1.sars)

colnames(Table1) <- c("Influenza (2009)","Influenza (5x)","Influenza (1918)","Influenza (1968)","Sars (2002-2004)")

rownames(Table1) <- c("Likelihood of disease containment \n(social distancing measures alone)(%)",
                      "Likelihood of disease containment \n(+contact tracing at home/work)(%)",
                      "Optimal contact tracing ratio",
                      "Number of cases",
                      "Number of deaths",
                      "Number of contacts traced",
                      "Number of contacts traced per case",
                      "Number of hospital admissions",
                      "Number of ICU admissions",
                      "Duration of disease containment (days)",
                      " - ",
                      "Case fatality rate"
)

kable(Table1) %>%
  kable_styling()
