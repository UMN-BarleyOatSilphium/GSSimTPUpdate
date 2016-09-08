#' Calculate realized prediction accuracy
#' 
#' @description
#' Calculates the correlation between predicted genotypic values and known
#' genotypic values. Can also calculate a boostrapped 95% confidence interval.
#' 
#' @param predicted.values The predicted genotypic values (i.e. GEBVs from
#' genomic prediction).
#' @param observed.values The observed genotypic values (from using the 
#' \code{\link{GSsim.TPUpdate}[genotypic.value]}).
#' @param boot.reps The number of bootstrapping replicates when estimating
#' the confidence interval.
#' 
#' @import boot
#' 
#' @export
#' 
validate.predictions <- function(predicted.values, observed.values,
                                 boot.reps = NULL) {
  
  # Deal with input
  predicted.values <- as.matrix(predicted.values)
  observed.values <- as.matrix(observed.values)
  
  # If no boot.reps, do not perform bootstrapping
  if (is.null(boot.reps)) {
    
    pred.r <- cor(predicted.values, observed.values)
    pred.r.sd <- NULL
    
  } else { # Otherwise perform bootstrapping
    
    # Perform bootstrapping to measure the correlation between GEBVs and phenotypes
    # First combine the data
    data.to.bootstrap <- cbind(predicted.values, observed.values)
    # Define a function for performing a correlation of the sample given by the ith replication
    boot.cor <- function(input.data, i) {
      rep.data <- input.data[i,]
      rep.cor <- cor(rep.data[,1], rep.data[,2])
      return(rep.cor)
    }
    # Perform bootstrapping
    boot.results <- boot(data = data.to.bootstrap, statistic = boot.cor, R = boot.reps)
    
    # Parse the results
    pred.r <- boot.results$t0
    pred.r.sd <- sd(boot.results$t)
  }
  
  # Create an output list
  output.list <- list(pred.r = pred.r, pred.r.sd = pred.r.sd, boot.reps = boot.reps)
  # Return the data
  return(output.list)
} # Close the function