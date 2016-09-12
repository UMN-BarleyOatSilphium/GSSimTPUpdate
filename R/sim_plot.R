#' Plot a list of results
#' 
#' @description 
#' Creates a plot of the statistic of interest over the cycles in the simulation.
#' The input should be a list derived from parsing the results.
#' 
#' @param data.list A \code{list} of parsed data for a particular statistic.
#' @param legend.pos The position of the legend. Accepts character inputs such as
#' \code{"topleft"}.
#' @param just.data Logical whether to return the estimates of the mean and 
#' confidence interval
#' @param ... Additional arguments to pass to \code{plot}, including \code{xlim},
#' \code{xlab}, \code{ylim}, \code{ylab}, and \code{main}.
#' 
#' @import dplyr
#' 
#' @export
#' 
sim.plot <- function(data.list, legend.pos, just.data = FALSE, ... ){
  
  # Handle additional arguments
  other.args <- list(...)
  
  # Parse the other arguments
  xlim <- other.args$xlim
  ylim <- other.args$ylim
  xlab <- other.args$xlab
  ylab <- other.args$ylab
  main <- other.args$main
  
  # Extract the names from the data.list for use as a factor
  list.names.factor <- names(data.list) %>%
    as.factor()
  list.names.numeric <- as.numeric(list.names.factor)
  
  # Apply a function over the data.list to calculate the mean and a confidence
  ## interval at every cycle
  data.parameters <- lapply(X = data.list, FUN = function(data) {
    
    # Find the mean across iterations
    mu <- apply(X = data, MARGIN = 1, FUN = mean, na.rm = T)
    
    # Calculate a confidence interval based on a t-distribution
    CI <- apply(X = data, MARGIN = 1, FUN = function(cycle) {
      t.per <- qt(p = (1 - (0.05 / 2)), df = length(cycle) - 1)
      t.per * ( sd(cycle, na.rm = T) / sqrt(length(cycle)) ) }) 
    
    list(mu = mu, CI = CI) })
  
  # If just the data is requested, return it
  if (just.data) {
    return(data.parameters)
  }
  
  # Determine the number of cycles
  n.cycles <- data.parameters[[1]][[1]] %>%
    length()
  
  # Options for null other arguments
  if (is.null(xlim)) {
    xlim <- c(1, n.cycles)
  }
  
  if (is.null(xlab)) {
    xlab <- "Cycle Number"
  }
  
  # If ylim is null, use a pretty range to determine it
  if (is.null(ylim)) {
    ylim <- sapply(data.parameters, function(sublist) sublist$mu) %>%
      range(na.rm = T) %>%
      pretty() %>%
      range()
  }
  
  
  # Create the empty plot first
  plot(0, type = "n", xlim = xlim, ylim = ylim, main = main, xlab = xlab, ylab = ylab)
  
  # Add the legend
  legend(legend.pos, legend = list.names.factor, pch = list.names.numeric, col = list.names.numeric)
  
  # Iterate over the data.parameters list
  for (i in seq_along(data.parameters)) {
    
    # Add points to the plot
    # Create jitter
    x.jitter <- - (0.1 * scale(1:length(data.parameters), scale = F)[i])
    
    points(x = (seq(xlim[1], xlim[2]) + x.jitter), data.parameters[[i]]$mu, 
           pch = list.names.numeric[i], type = "b", col = list.names.numeric[i])
    
    # Add CI bars
    segments(x0 = (seq(xlim[1], xlim[2]) + x.jitter), 
             y0 = (data.parameters[[i]]$mu - data.parameters[[i]]$CI), 
             x1 = (seq(xlim[1], xlim[2]) + x.jitter), 
             y1 = (data.parameters[[i]]$mu + data.parameters[[i]]$CI))
    
  } # Close the for loop
  
} # Close the function