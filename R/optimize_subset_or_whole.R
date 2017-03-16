#' Optimize PEVmean or CDmean
#' 
#' @description 
#' Optimize the PEVmean or CDmean using a subset (for faster computation) or the
#' whole set of individuals.
#'  
#' @param method optimization procedure. Can be \code{"PEVmean"} or \code{"CDmean"}.
#' @param A The n x n relationship matrix.
#' @param n.training The target size of the training population.
#' @param phenotyped.index The index of the "phenotyped" entries in the matrix 
#' \code{A}. This can be interpreted as the entries from which you wish to
#' develop an optimized training population.
#' @param unphenotyped.index The index of "unphenotyped" entries
#' @param V_a The additive genetic variance
#' @param V_e The residual variance
#' @param max.iter The maximum number of iterations of the exhange algorithm
#' 
#' @importFrom MASS ginv
#' 
#' @export
#' 
optimize.subset <- function(method, A, n.training, phenotyped.index, 
                            unphenotyped.index, V_a, V_e, max.iter) {
    
  # Determine the total population size
  n.total = n.training + length(unphenotyped.index)
  
  # Create design matrix
  M <- design.M(n.training)
  
  # Randomly sample the "phenotyped" lines to start the TP
  ## These will serve as the first sample of phenotyped lines
  start.OTP <- sample(phenotyped.index, n.training) %>% 
    sort()
  # Save the sample
  save.OTP <- start.OTP
  
  # Find the "phenotyped" lines that were not included in the training set
  # but could be added later
  candidates.OTP <- setdiff(phenotyped.index, start.OTP)
  
  # This matrix will remain constant
  Z <- diag(n.total)[seq_along(start.OTP),] 
  
  # Create a vector of the indicies of the sampled training set and the 
  # unphenotyped lines
  subset.index <- c(start.OTP, unphenotyped.index)
  # Subset the A.matrix
  A1 <- A[subset.index, subset.index]
  
  # Find the inverse of the subset
  A1.inv <- invert.mat(A1, silent = T)
  
  # The index of the unphenotyped lines becomes the last n rows of the 
  # subsetted A matrix, where n is the number of unphenotyped lines
  unphenotyped.index.A1 <- setdiff( seq(nrow(A1)), seq_along(start.OTP))
  
  # Contrasts matrix
  ## The dimenstions should be row = n.total and column = n.unphenotyped 
  ## (i.e. the parents). This matrix will also remain constant
  c.mat <- make.contrast(unphenotyped.index = unphenotyped.index.A1, n.total = n.total)
  
  # Calculate the statistic
  if (method == "PEVmean") {
    statistic <- PEVmean(M = M, Z = Z, A = A1, A.inv = A1.inv, 
                         c = c.mat, V_e = V_e, V_a = V_a)
  }
  if (method == "CDmean") {
    statistic <- CDmean(M = M, Z = Z, A = A1, A.inv = A1.inv, 
                        c = c.mat, V_e = V_e, V_a = V_a)
  }
  
  statistic.save <- statistic
  
  # Create an empty vector to store stat values
  statistic.vector <- vector("numeric", length = max.iter)
  
  iter.counter = 1 # Iteration counter
  
  # While loop
  while(iter.counter <= max.iter) {
    
    # Add one to the count
    iter.counter = iter.counter + 1
    
    # Randomly remove one individual from the training set
    train.to.remove <- sample(start.OTP, 1)
    # Randomly choose one individual from the candidate set to add
    candidate.to.add <- sample(candidates.OTP, 1)
    
    # Create the new training set
    new.OTP <- sort(c(candidate.to.add, start.OTP[start.OTP != train.to.remove]))
    
    # Create a vector of the indicies of the sampled training set and the unphenotyped lines
    subset.index <- c(new.OTP, unphenotyped.index)
    # Subset the A.matrix
    A1 <- A[subset.index, subset.index]
    
    # Find the inverse of the subset
    A1.inv <- invert.mat(A1, silent = T)
    
    # Recalculate the statistic
    if (method == "PEVmean") {
      statistic <- PEVmean(M = M, Z = Z, A = A1, A.inv = A1.inv, 
                           c = c.mat, V_e = V_e, V_a = V_a)
      
      # If the statistic is lower, accept it, otherwise resample
      if (statistic < statistic.save) {
        start.OTP <- new.OTP
        candidates.OTP <- setdiff(phenotyped.index, new.OTP)
        statistic.save <- statistic
      }
      
    }
    if (method == "CDmean") {
      statistic <- CDmean(M = M, Z = Z, A = A1, A.inv = A1.inv, 
                          c = c.mat, V_e = V_e, V_a = V_a)
      
      # If the statistic is higher, accept it, otherwise resample
      if (statistic > statistic.save) {
        start.OTP <- new.OTP
        candidates.OTP <- setdiff(phenotyped.index, new.OTP)
        statistic.save <- statistic
      }
      
    }
    
    # Record the statistic 
    statistic.vector[iter.counter - 1] <- statistic.save
    
  } # Close the loop
  
  # Return the OTP, the current value of the statistic, etc.
  outlist <- list(OTP = start.OTP, statistic = statistic.save, 
                  statistic.vector = statistic.vector)
  return(outlist)

} # Close the function

#' 
#' Optimize PEVmean or CDmean
#' 
#' @describeIn optimize.subset
#' 
#' @importFrom MASS ginv
#' 
#' @export
#' 
optimize.whole <- function(method, A, n.training, phenotyped.index, 
                           unphenotyped.index, V_a, V_e, max.iter) {
  
  # The total number of entries for the contrasts matrix is the number of 
  # entries in the A matrix
  n.total <- nrow(A)
  
  # Create the M matrix
  M <- design.M(n.training)
  
  # Randomly sample the "phenotyped" lines to start the TP
  ## These will serve as the first sample of phenotyped lines
  start.OTP <- sort(sample(phenotyped.index, n.training))
  # Save the sample
  save.OTP <- start.OTP
  
  # Find the "phenotyped" lines that were not included, but could be included
  candidates.OTP <- setdiff(phenotyped.index, start.OTP)
  
  Z <- matrix(0, n.training, n.total)
  # Fill in the coordinates of the starting TP (rows = index in the start.OTP, column = index in the A.mat)
  Z[seq_along(start.OTP), start.OTP] <- 1
  
  # Use the whole A.mat
  A1 <- A
  # Get the inverse
  A1.inv <- invert.mat(A1, silent = T)
  
  # Contrasts matrix
  c.mat <- make.contrast(unphenotyped.index = unphenotyped.index, n.total = n.total)
  
  # Calculate the statistic
  if (method == "PEVmean") {
    statistic <- PEVmean(M = M, Z = Z, A = A1, A.inv = A1.inv, 
                         c = c.mat, V_e = V_e, V_a = V_a)
  }
  if (method == "CDmean") {
    statistic <- CDmean(M = M, Z = Z, A = A1, A.inv = A1.inv,
                        c = c.mat, V_e = V_e, V_a = V_a)
  }
  
  statistic.save <- statistic
  
  # Create an empty vector to store stat values
  statistic.vector <- vector("numeric", length = max.iter)
  
  iter.counter = 1 # Iteration counter

  # While loop
  while(iter.counter <= max.iter) {
    
    # Add one to the count
    iter.counter = iter.counter + 1
    
    # Randomly remove one individual from the training set
    train.to.remove <- sample(start.OTP, 1)
    # Randomly choose one individual from the candidate set to add
    candidate.to.add <- sample(candidates.OTP, 1)
    
    # Create the new training set
    new.OTP <- sort(c(candidate.to.add, start.OTP[start.OTP != train.to.remove]))
    
    Z <- matrix(0, n.training, n.total)
    # Fill in the coordinates of the starting TP 
    # (rows = index in the start.OTP, column = index in the A.mat)
    Z[seq_along(new.OTP), new.OTP] <- 1
    
    # No need to remake the contrast matrix since the index of the 
    # "unphenotyped" lines in the A matrix remains the same.
    
    # Recalculate the statistic
    if (method == "PEVmean") {
      statistic <- PEVmean(M = M, Z = Z, A = A1, A.inv = A1.inv, 
                           c = c.mat, V_e = V_e, V_a = V_a)
      
      # If the statistic is lower, accept it, otherwise resample
      if (statistic < statistic.save) {
        start.OTP <- new.OTP
        candidates.OTP <- setdiff(phenotyped.index, new.OTP)
        statistic.save <- statistic
      }
      
    }
    if (method == "CDmean") {
      statistic <- CDmean(M = M, Z = Z, A = A1, A.inv = A1.inv, 
                          c = c.mat, V_e = V_e, V_a = V_a)
      
      # If the statistic is higher, accept it, otherwise resample
      if (statistic > statistic.save) {
        start.OTP <- new.OTP
        candidates.OTP <- setdiff(phenotyped.index, new.OTP)
        statistic.save <- statistic
      }
      
    }
    
    # Record the statistic 
    statistic.vector[iter.counter - 1] <- statistic.save
    
  } # Close the loop
  
  # Return the OTP, the current value of the statistic, etc.
  outlist <- list(OTP = start.OTP, statistic = statistic.save, 
                  statistic.vector = statistic.vector)
  return(outlist)
  
} # Close the function
    



