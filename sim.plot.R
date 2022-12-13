library(dplyr, warn=FALSE)
library(tidyr, warn=FALSE)
library(binom, warn=FALSE)
library(ggplot2, warn=FALSE)
library(scales, warn=FALSE)
library(gghighlight, warn=FALSE)

# Begin graphing procedure

# User-defined function to create a line graph of Likelihood of control (y-axis) vs Contact tracing ratio (x-axis) for certain parameters.
# The x-axis can be specified by the user in future versions
Sim.plot <- function(data, seed=10, times=NULL, param='total_case', x=ctr, title=NULL, ci=TRUE, plot=TRUE) {
  
  # seed is the initial seed
  # times is the number of times to multiply the seed to obtain the threshold
  if(is.null(times)) times <- c(5, 10, 15, 20)
  
  # Define the thresholds (number of cumulative cases or other outcome measure. See param argument.)
  thresholds <- seed * times
  
  # Function to calculate the number of times the control measure was successful (param < threshold)
  control <- function(data) {
    rowSums(
      sapply(data, function(i) max(i[, param]) < thresholds))
  }
  
  # Contact tracing ratio
  ctr <- seq(0, by=0.1, length.out=length(data))
  
  # Apply the function to the scenario data
  output <- setNames(lapply(data, control), ctr)
  
  # Combine output and add the number of times
  result <- do.call(bind_rows, output)
  times <- paste(times, 'X')
  result <- mutate(result, Times=factor(times, levels=times))
  
  # Number of simulations in each scenario
  nsims <- unlist(lapply(data, length))
  
  # Pivot to long form and add 95% binomial confidence intervals
  results <- pivot_longer(result, -Times, names_to="ctr", values_to="x") %>%
    mutate(nsims=rep(nsims, times=length(times)), 
           p=binom.confint(x, nsims, method="wilson"),
           ctr=as.numeric(ctr)) %>%
    select(-c(x, nsims)) %>% unnest(p) %>%
    mutate(mean=round(mean, 2),
           ctrx=ctr+0.1,
           ctrx=if_else(ctrx>1.1, 0, ctrx))
  
  # Plot
  if(plot) {
    theme_set(theme_bw())
    
    theme_update(panel.grid.minor=element_blank())
    
    # Scaling for y-axis
    exp2 <- trans_new(
      name="exp2",
      transform = function(x) exp(exp(x)),
      inverse= function(x) log(log(x))
    )
    
    gg <- ggplot(results, aes(x=ctrx, y=mean, col=Times)) +
      geom_line(aes(group=Times), 
                size=1,
                position=position_dodge(width=0.05)) +
      scale_x_continuous(limits=c(-0.03, 1.21),
                         breaks=seq(from=0, to=length(ctr)/10, by=0.1),
                         labels=c('0*', seq(from=0, to=length(ctr)/10 - 0.1, by=0.1))) +
      labs(y="Probability of containment", 
           x="Contact tracing ratio", 
           title="") +
      theme(legend.position="none")
    
    if(f=="sars") {
      gg <- gg +
      scale_y_continuous(trans=exp2, 
                         limits = c(0.05, 1.02), 
                         breaks=c(0,0.5,0.75,0.85,0.9,0.95,1))
    }
    
    if(ci) {
      gg <- gg + geom_errorbar(aes(ymin=lower, ymax=upper),
                               width=0.05,
                               position=position_dodge(width=0.05))
    }
    
    ggh <- gg + gghighlight(use_direct_label = TRUE,
                            label_params = list(nudge_x = 0.1))
    print(ggh)
  }
  
  return(results)
}
