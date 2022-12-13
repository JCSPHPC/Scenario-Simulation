###
library(dplyr)

# Estimation of various parameters

### Function to return the maximum value of a parameter from many realizations within a list
Max <- function(data, param){
  get_max <- function(data) {
    sapply(data, function(x) max(x[, param]))
  }
  lapply(data, get_max)
}


### Function to return the mean value of a non-zero parameter from many realizations within a list
Mean <- function(data, param){
  get_mean <- function(data) {
    sapply(data, function(x) mean(x[, param][x[, param]>0]))
  }
  lapply(data, get_mean)
}

### Function to return the sum of a parameter from many realizations within a list
Sum <- function(data, param){
  get_sum <- function(data) {
    sapply(data, function(x) sum(x[, param]))
  }
  lapply(data, get_sum)
}


# Function to write the output to a file.
Out <- function(data, x, digits=1) {
  nscen <- length(x)
  
  D <- data.frame(scenario=rep(1:nscen, times=unlist(lapply(data, length))),
                  x=unlist(x))
  
  fdig <- paste0("%.", digits, "f")
  
  group_by(D, scenario) %>%
    summarise(min = sprintf(fdig, min(x, na.rm=TRUE)),
              q.025 = sprintf(fdig, quantile(x, probs=0.025, na.rm=TRUE)),
              q.50 = sprintf(fdig, median(x, na.rm=TRUE)),
              q.975 = sprintf(fdig, quantile(x, probs=0.975, na.rm=TRUE)),
              max = sprintf(fdig, max(x, na.rm=TRUE)), .groups='drop') 
}

# End