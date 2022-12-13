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

# Loop through the data folders
for(f in folders[5]) {

  # Change the global paths
  path.data <- paste(Path.data, f, sep="/")
  path.out <- paste(Path.out, f, sep="/")

  # Import the data
  setwd(path.data)
  cat("\nImporting data for", f, "...\n")

  # Load the plotting function
  source("F:/Edward/Kwok/Simulation/import.R")
  
  # Graphs
  cat("\nPlotting data for", f, "...\n")
  source("F:/Edward/Kwok/Simulation/sim.plot.R")
  
  # Call the function (one with 95% CI and one without)
  if(f=="sars")
    times <- c(5, 10, 15, 20)
  else
    times <- c(5, 10, 15, 20, 40, 60, 80)
  
  cases <- Sim.plot(data, param='total_case', ci=TRUE, times=times)
  ggsave(paste(path.out, "Figure 1.png", sep="/"), width=7, height=7, units="in")

  cases <- Sim.plot(data, param='total_case', ci=FALSE, times=times)
  ggsave(paste(path.out, "Figure 2.png", sep="/"), width=7, height=7, units="in")
  
  select(cases, -method) %>%
    mutate(lower=round(lower, 2),
           upper=round(upper, 2)) %>%
    write.csv(file=paste(path.out, "cases.csv", sep="/"), row.names=FALSE)

    # Parameter estimates
  cat("\nEstimating parameters for", f, "...\n")

  # Load the required functions  
  source("F:/Edward/Kwok/Simulation/estimate.R")
  
  # Calculate the estimates and send to an output file
  
  # 1. Median (95% CI) DAILY number of individuals traced in each contact tracing strategy. I would like to check the capacity for the community quarantine centre. (this one is correct; Mean of NON-ZERO values in column of n_quarantine_events in each excel file then choose 50th, 2.5th and 97.5th from 100 excel files estimates) - average daily resources spent 
  
  quarantine_mean <- Mean(data, 'n_quarantine_events')
  S1 <- Out(data, quarantine_mean)
  path <- paste(path.out, 'quarantine_mean', sep="/")
  write.csv(S1, paste0(path, ".csv"), row.names=FALSE)
  
  
  # 2. Mean, maximum and minimum DAILY number of individuals traced in each contact tracing strategy (taking the mean within each simulation/excel file). I would like to check the capacity for the community quarantine centre.
  quarantine_sum <- Sum(data, 'n_quarantine_events')
  S2 <- Out(data, quarantine_sum, digits=0)
  path <- paste(path.out, 'quarantine_sum', sep="/")
  write.csv(S2, paste0(path, ".csv"), row.names=FALSE)
  
  
  # 3. The actual number of contact to be traced per case (formula = sum the column of n_quarantine_events OVER the last value of total_case variable) in given contact tracing strategies in the early phase of the future outbreak of novel pathogens with 95%CI by ranking 100 simulations and values of  97.5th and 2.5th
  
  total_case <- Max(data, 'total_case')
  
  D <- tibble(scenario=rep(1:nscen, times=unlist(lapply(data, length))),
              total_case=unlist(total_case),
              quarantine_sum=unlist(quarantine_sum),
              contacts.traced.per.case=quarantine_sum/total_case)
  
  fdig <- paste0("%.", 1, "f")
  
  S3 <- group_by(D, scenario) %>%
    summarise(min = sprintf(fdig, min(contacts.traced.per.case, na.rm=TRUE)),
              q.025 = sprintf(fdig, quantile(contacts.traced.per.case, probs=0.025, na.rm=TRUE)),
              q.50 = sprintf(fdig, quantile(contacts.traced.per.case, probs=0.5, na.rm=TRUE)),
              q.975 = sprintf(fdig, quantile(contacts.traced.per.case, probs=0.975, na.rm=TRUE)), 
              max = sprintf(fdig, max(contacts.traced.per.case, na.rm=TRUE)), .groups='drop')
  
  write.csv(S3, paste(path.out, "contacts_traced_per_case.csv", sep="/"), row.names=FALSE)
  
  
  # 4. Median number of deaths in given contact tracing strategies in the early phase of the future outbreak of novel pathogens with 95%CI by ranking 100 simulations and values of  50th,  97.5th and 2.5th (using MAX of n_death in each excel file)
  
  total_death <- Max(data, 'total_death')
  S4 <- Out(data, total_death, digits=0)
  path <- paste(path.out, 'total_death', sep="/")
  write.csv(S4, paste0(path, ".csv"), row.names=FALSE)
  
  
  # 5. Maximum, minimum and mean number of DAILY hospital bed occupied in given contact tracing strategies in the early phase of the future outbreak of novel pathogens with 95%CI by ranking 100 simulations and values of 50th, 97.5th and 2.5th (using MAX, MIN and mean value of NON-ZERO value in n_hospital and obtain the 50th, 97.5th and 2.5th from 100 realization)
  
  n_hospital <- Mean(data, 'n_hospital')
  S5 <- Out(data, n_hospital, digits=0)
  path <- paste(path.out, 'n_hospital', sep="/")
  write.csv(S5, paste0(path, ".csv"), row.names=FALSE)
  
  
  # 6. Maximum, minimum and mean number of DAILY ICU admission in given contact tracing strategies in the early phase of the future outbreak of novel pathogens with 95%CI by ranking 100 simulations and values of  97.5th and 2.5th  (using MAX, MIN and mean value of NON-ZERO values in n_critical and obtain the 50th, 97.5th and 2.5th from 100 realization
  
  n_critical <- Mean(data, 'n_critical')
  S6 <-  Out(data, n_critical, digits=0)
  path <- paste(path.out, 'n_critical', sep="/")
  write.csv(S6, paste0(path, ".csv"), row.names=FALSE)
 
  
  # 7.Median CFR/IFR given contact tracing strategies using i) total no. of death/ total no. of CASES ii) total no. of death/ total no. of infection Using the formula i)  MAX of n_death/MAX of total case ii)  MAX of n_death/MAX of total infected
  
  total_infected <- Max(data, 'total_infected')
  
  D <- tibble(scenario=rep(1:nscen, times=unlist(lapply(data, length))),
              total_case=unlist(total_case),
              total_death=unlist(total_death),
              total_infected=unlist(total_infected),
              CFR=total_death/total_case,
              IFR=total_death/total_infected)
  
  CFR <- group_by(D, scenario) %>%
    summarise(min = sprintf(fdig, 100*min(CFR, na.rm=TRUE)),
              q.025 = sprintf(fdig, 100*quantile(CFR, probs=0.025, na.rm=TRUE)),
              q.50 = sprintf(fdig, 100*quantile(CFR, probs=0.5, na.rm=TRUE)),
              q.975 = sprintf(fdig, 100*quantile(CFR, probs=0.975, na.rm=TRUE)), 
              max = sprintf(fdig, 100*max(CFR, na.rm=TRUE)), .groups='drop')
  
  write.csv(CFR, paste(path.out, "CFR.csv", sep="/"), row.names=FALSE)
  
  IFR <- group_by(D, scenario) %>%
    summarise(min = sprintf(fdig, 100*min(IFR, na.rm=TRUE)),
              q.025 = sprintf(fdig, 100*quantile(IFR, probs=0.025, na.rm=TRUE)),
              q.50 = sprintf(fdig, 100*quantile(IFR, probs=0.5, na.rm=TRUE)),
              q.975 = sprintf(fdig, 100*quantile(IFR, probs=0.975, na.rm=TRUE)), 
              max = sprintf(fdig, 100*max(IFR, na.rm=TRUE)), .groups='drop')
  
  write.csv(IFR, paste(path.out, "IFR.csv", sep="/"), row.names=FALSE)
  
  
  #8. Number of cases
  S8 <- Out(data, total_case, digits=0)
  path <- paste(path.out, 'total_case', sep="/")
  write.csv(S8, paste0(path, ".csv"), row.names=FALSE)
  
  
  #9. Number of infections
  S9 <- Out(data, total_infected, digits=0)
  path <- paste(path.out, 'total_infected', sep="/")
  write.csv(S9, paste0(path, ".csv"), row.names=FALSE)
  
  
  #10. Median duration of containment
  options(warn = -1)
  
  durations <- lapply(data, function(x) 
    lapply(x, function(y) 
      y$time[1 + try(max(which(diff(y$total_case) !=0)), silent=TRUE)]))
  options(warn = 0)
  
  D <- tibble(scenario=rep(1:nscen, times=unlist(lapply(data, length))),
              duration=unlist(durations)) 
  
  fdig <- paste0("%.", 0, "f")
  
  S10 <- group_by(D, scenario) %>%
    summarise(min = sprintf(fdig, round(min(duration, na.rm=TRUE))),
              q.025 = sprintf(fdig, round(quantile(duration, probs=0.025, na.rm=TRUE))),
              q.50 = sprintf(fdig, round(quantile(duration, probs=0.5, na.rm=TRUE))),
              q.975 = sprintf(fdig, round(quantile(duration, probs=0.975, na.rm=TRUE))), 
              max = sprintf(fdig, round(max(duration, na.rm=TRUE))), .groups='drop')
  
  write.csv(S10, paste(path.out, "duration.csv", sep="/"), row.names=FALSE)
  
  cat("\n")
}


