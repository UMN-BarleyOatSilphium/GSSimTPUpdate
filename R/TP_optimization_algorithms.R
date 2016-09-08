#' Optimization algorithms
#' 
#' @export
#' 
TP.optimization.algorithms <- function(A, # relationship matrix of all candidate individuals
                                       phenotyped.lines, # Character of line name to be considered "phenotyped"
                                       unphenotyped.lines, # Character of line name to be considered "unphenotyped"
                                       n.TP, # The size of the desired TP
                                       V_e, # V_e for calculating lambda
                                       V_a, # V_a for calculating lambda
                                       optimization.method = c("PEVmean", "CDmean"),
                                       max.iter = 800,
                                       use.subset = FALSE # Logical indicating whether a subset of the whole should be used for measuring the contrats. This saves time but may be innacurate
) {
  
  # Define a function to create a matrix of contrasts between the individuals not
  ## in the training population (the candidates) and the mean of the
  ## whole population. This is taken from Rincent et al 2012
  make.contrast <- function(unphenotyped.index, # Index of the unphenotyped individuals in the A.mat
                            n.total) {
    # Number of unphenotyped individuals
    n.unphenotyped <- length(unphenotyped.index)
    
    # The positive value of the contrast
    pos.contrast = 1 - (1 / n.total)
    # The negative value of the contrast
    neg.contrast = -1 / n.total
    
    # Make the total constrast matrix
    contrast.mat <- matrix(neg.contrast, n.total, n.unphenotyped)
    # Fill in the positive values
    for (i in 1:n.unphenotyped) {
      contrast.mat[unphenotyped.index[i], i] <- pos.contrast
    }
    
    return(contrast.mat)
  }
  
  # Find the index of the phenotyped and unphenotyped lines in the A matrix
  phenotyped.index <- which(row.names(A) %in% phenotyped.lines)
  unphenotyped.index <- which(row.names(A) %in% unphenotyped.lines)
  
  # Set the total
  # If using a subset, the total is size of the training set + the number of 
  ## unphenotyped lines
  if (use.subset) {
    n.total = n.TP + length(unphenotyped.index)
  } else {
    # Otherwise it is the number of lines in the A matrix
    n.total = nrow(A)
  }
  
  # Calculate lambda using the provided variances
  lambda = V_e / V_a
  
  ## Begin the algorithm
  # The algorithm is the same for both PEVmean and CDmean (exchange), so we have two paths
  ## depending on which is desired
  if (optimization.method == "PEVmean") {
    
    # Set the CDmean to NA
    CDmean.save = NA
    CD.vector = NA
    
    ### Initialize the optimization algorithm
    ## Separate paths based on whether a subset will be used
    if (use.subset) {
      
      # Create design matricies
      I <- diag(n.TP) # Identity matrix of TP size
      X <- matrix(1, n.TP, 1) # X matrix
      M <- I - ( X %*% solve(t(X) %*% X) %*% t(X) )
      
      # Randomly sample the "phenotyped" lines to start the TP
      ## These will serve as the first sample of phenotyped lines
      start.OTP <- sort(sample(phenotyped.index, n.TP))
      # Save the sample
      save.OTP <- start.OTP
      
      # Find the "phenotyped" lines that were not included, but could be included
      candidates.OTP <- setdiff(phenotyped.index, start.OTP)
      
      # When using a subset (to save time)
      Z <- diag(n.total)[1:length(start.OTP),] # This matrix will remain constant
      # Create a vector of the indicies of the sampled training set and the unphenotyped lines
      subset.index <- c(start.OTP, unphenotyped.index)
      # Subset the A.matrix
      A1 <- A[subset.index, subset.index]
      
      # The index of the unphenotyped lines becomes the last n rows of the subsetted A matrix, where n
      ## is the number of unphenotyped lines
      unphenotyped.index.A1 <- setdiff(1:nrow(A1), 1:length(start.OTP))
      
      # Try to find the inverse of candidate.A.mat, if not add ridge
      A1.inv <- try(solve(A1), silent = T)
      # Assess whether it was an error
      if (class(A1.inv) == "try-error") {
        r = 1e-6
        ridge <- diag(nrow(A1)) * r
        A1 = A1 + ridge
        # Recalculate the inverse
        A1.inv <- solve(A1)
      }
      
      # Contrasts matrix
      ## The dimenstions should be row = n.total and column = n.unphenotyped (aka the parents)
      c.mat <- make.contrast(unphenotyped.index = unphenotyped.index.A1, n.total = n.total) # This matrix will also remain constant
      
      # Calculate the inital PEVmean of the set
      C_22 = solve( (t(Z) %*% M %*% Z) + (lambda * A1.inv) )
      numerator = ( t(c.mat) %*% C_22 %*% c.mat )
      denominator = ( t(c.mat) %*% c.mat )
      PEV.mat = numerator / denominator
      
      PEV = diag(PEV.mat) * V_e
      # Calculate the PEVmean
      PEVmean.save = mean(PEV)
      
      # Create an empty vector to store PEVmean values
      PEV.vector <- vector("numeric", length = max.iter)
      
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
        
        # Try to find the inverse of candidate.A.mat, if not add ridge
        A1.inv <- try(solve(A1), silent = T)
        # Assess whether it was an error
        if (class(A1.inv) == "try-error") {
          r = 1e-6
          ridge <- diag(nrow(A1)) * r
          A1 = A1 + ridge
          # Recalculate the inverse
          A1.inv <- solve(A1)
        }
        
        C_22 = solve( (t(Z) %*% M %*% Z) + (lambda * A1.inv) )
        numerator = ( t(c.mat) %*% C_22 %*% c.mat )
        denominator = ( t(c.mat) %*% c.mat )
        PEV.mat = numerator / denominator
        
        PEV = diag(PEV.mat) * V_e
        PEVmean.new <- mean(PEV)
        
        # If the new PEVmean is lower, accept, if not reject
        if (PEVmean.new < PEVmean.save) {
          start.OTP <- new.OTP
          PEVmean.save <- PEVmean.new
          # Create the new candidate set
          candidates.OTP <- setdiff(phenotyped.index, new.OTP)
        }
        
        # Record the PEVmean of the new PEVmean (if accepted) or the previous (if not)
        PEV.vector[iter.counter - 1] <- PEVmean.save
        
      } # Close the loop
      
    } else { # Path for using the whole set
      
      # Create design matricies
      I <- diag(n.TP) # Identity matrix of TP size
      X <- matrix(1, n.TP, 1) # X matrix
      M <- I - ( X %*% solve(t(X) %*% X) %*% t(X) )
      
      # Randomly sample the "phenotyped" lines to start the TP
      ## These will serve as the first sample of phenotyped lines
      start.OTP <- sort(sample(phenotyped.index, n.TP))
      # Save the sample
      save.OTP <- start.OTP
      
      # Find the "phenotyped" lines that were not included, but could be included
      candidates.OTP <- setdiff(phenotyped.index, start.OTP)
      
      # When using the whole candidates
      Z <- matrix(0, n.TP, n.total)
      # Fill in the coordinates of the starting TP (rows = index in the start.OTP, column = index in the A.mat)
      for (i in 1:length(start.OTP)) { Z[i, start.OTP[i] ] = 1 } 
      # Use the whole A.mat
      A1 <- A
      
      # Try to find the inverse of candidate.A.mat, if not add ridge
      A1.inv <- try(solve(A1), silent = T)
      # Assess whether it was an error
      if (class(A1.inv) == "try-error") {
        r = 1e-6
        ridge <- diag(nrow(A1)) * r
        A1 = A1 + ridge
        # Recalculate the inverse
        A1.inv <- solve(A1)
      }
      
      # Contrasts matrix
      ## The dimenstions should be row = n.total and column = n.unphenotyped (aka the parents)
      c.mat <- make.contrast(unphenotyped.index = unphenotyped.index, n.total = n.total)
      
      # Calculate the inital PEVmean of the set
      C_22 = solve( (t(Z) %*% M %*% Z) + (lambda * A1.inv) )
      numerator = ( t(c.mat) %*% C_22 %*% c.mat )
      denominator = ( t(c.mat) %*% c.mat )
      PEV.mat = numerator / denominator
      
      PEV = diag(PEV.mat) * V_e
      # Calculate the PEVmean
      PEVmean.save = mean(PEV)
      
      # Create an empty vector to store PEVmean values
      PEV.vector <- vector("numeric", length = max.iter)
      
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
        
        # When using the whole candidates
        Z <- matrix(0, n.TP, n.total)
        # Fill in the coordinates of the starting TP (rows = index in the start.OTP, column = index in the A.mat)
        for (i in 1:length(new.OTP)) { Z[i, new.OTP[i] ] = 1 } 
        
        # Calculate the inital PEVmean of the set
        # No need to remake the contrast matrix because the unphenotyped lines remain unchanged
        C_22 = solve( (t(Z) %*% M %*% Z) + (lambda * A1.inv) )
        numerator = ( t(c.mat) %*% C_22 %*% c.mat )
        denominator = ( t(c.mat) %*% c.mat )
        PEV.mat = numerator / denominator
        
        PEV = diag(PEV.mat) * V_e
        PEVmean.new <- mean(PEV)
        
        # If the new PEVmean is lower, accept, if not reject
        if (PEVmean.new < PEVmean.save) {
          start.OTP <- new.OTP
          PEVmean.save <- PEVmean.new
          # Create the new candidate set
          candidates.OTP <- setdiff(phenotyped.index, new.OTP)
        }
        
        # Record the PEVmean of the new PEVmean (if accepted) or the previous (if not)
        PEV.vector[iter.counter - 1] <- PEVmean.save
        
      } # Close the loop
      
    } # Close the use.subset if statement
    
  } # Close the PEVmean optimization if statement
  
  # New if statement for CDmean
  if(optimization.method == "CDmean") {
    
    # Set the PEVmean to NA
    PEVmean.save = NA
    PEV.vector = NA
    
    ### Initialize the optimization algorithm
    ## Separate paths based on whether a subset will be used
    if (use.subset) {
      
      # Create design matricies
      I <- diag(n.TP) # Identity matrix of TP size
      X <- matrix(1, n.TP, 1) # X matrix
      M <- I - ( X %*% solve(t(X) %*% X) %*% t(X) )
      
      # Randomly sample the "phenotyped" lines to start the TP
      ## These will serve as the first sample of phenotyped lines
      start.OTP <- sort(sample(phenotyped.index, n.TP))
      # Save the sample
      save.OTP <- start.OTP
      
      # Find the "phenotyped" lines that were not included, but could be included
      candidates.OTP <- setdiff(phenotyped.index, start.OTP)
      
      # When using a subset (to save time)
      Z <- diag(n.total)[1:length(start.OTP),] # This matrix will remain constant
      # Create a vector of the indicies of the sampled training set and the unphenotyped lines
      subset.index <- c(start.OTP, unphenotyped.index)
      # Subset the A.matrix
      A1 <- A[subset.index, subset.index]
      
      # The index of the unphenotyped lines becomes the last n rows of the subsetted A matrix, where n
      ## is the number of unphenotyped lines
      unphenotyped.index.A1 <- setdiff(1:nrow(A1), 1:length(start.OTP))
      
      # Try to find the inverse of candidate.A.mat, if not add ridge
      A1.inv <- try(solve(A1), silent = T)
      # Assess whether it was an error
      if (class(A1.inv) == "try-error") {
        r = 1e-6
        ridge <- diag(nrow(A1)) * r
        A1 = A1 + ridge
        # Recalculate the inverse
        A1.inv <- solve(A1)
      }
      
      # Contrasts matrix
      ## The dimenstions should be row = n.total and column = n.unphenotyped (aka the parents)
      c.mat <- make.contrast(unphenotyped.index = unphenotyped.index.A1, n.total = n.total) # This matrix will also remain constant
      
      # Calculate the inital CDmean of the set
      C_22 = solve( (t(Z) %*% M %*% Z) + (lambda * A1.inv) )
      numerator = ( t(c.mat) %*% (A1 - (lambda * C_22)) %*% c.mat )
      denominator = t(c.mat) %*% A1 %*% c.mat
      CD.mat = numerator / denominator
      
      CD = diag(CD.mat)
      # Calculate the CDmean
      CDmean.save = mean(CD)
      
      # Create an empty vector to store PEVmean values
      CD.vector <- vector("numeric", length = max.iter)
      
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
        
        # Try to find the inverse of candidate.A.mat, if not add ridge
        A1.inv <- try(solve(A1), silent = T)
        # Assess whether it was an error
        if (class(A1.inv) == "try-error") {
          r = 1e-6
          ridge <- diag(nrow(A1)) * r
          A1 = A1 + ridge
          # Recalculate the inverse
          A1.inv <- solve(A1)
        }
        
        C_22 = solve( (t(Z) %*% M %*% Z) + (lambda * A1.inv) )
        numerator = ( t(c.mat) %*% (A1 - (lambda * C_22)) %*% c.mat )
        denominator = t(c.mat) %*% A1 %*% c.mat
        CD.mat = numerator / denominator
        
        CD = diag(CD.mat)
        # Find the new mean
        CDmean.new <- mean(CD)
        
        # If the new CDmean is lower, accept, if not reject
        if (CDmean.new > CDmean.save) {
          start.OTP <- new.OTP
          CDmean.save <- CDmean.new
          # Create the new candidate set
          candidates.OTP <- setdiff(phenotyped.index, new.OTP)
        }
        
        # Record the CDmean of the new CDmean (if accepted) or the previous (if not)
        CD.vector[iter.counter - 1] <- CDmean.save
        
      } # Close the loop
      
    } else { # Path for using the whole set
      
      # Create design matricies
      I <- diag(n.TP) # Identity matrix of TP size
      X <- matrix(1, n.TP, 1) # X matrix
      M <- I - ( X %*% solve(t(X) %*% X) %*% t(X) )
      
      # Randomly sample the "phenotyped" lines to start the TP
      ## These will serve as the first sample of phenotyped lines
      start.OTP <- sort(sample(phenotyped.index, n.TP))
      # Save the sample
      save.OTP <- start.OTP
      
      # Find the "phenotyped" lines that were not included, but could be included
      candidates.OTP <- setdiff(phenotyped.index, start.OTP)
      
      # When using the whole candidates
      Z <- matrix(0, n.TP, n.total)
      # Fill in the coordinates of the starting TP (rows = index in the start.OTP, column = index in the A.mat)
      for (i in 1:length(start.OTP)) { Z[i, start.OTP[i] ] = 1 } 
      # Use the whole A.mat
      A1 <- A
      
      # Try to find the inverse of candidate.A.mat, if not add ridge
      A1.inv <- try(solve(A1), silent = T)
      # Assess whether it was an error
      if (class(A1.inv) == "try-error") {
        r = 1e-6
        ridge <- diag(nrow(A1)) * r
        A1 = A1 + ridge
        # Recalculate the inverse
        A1.inv <- solve(A1)
      }
      
      # Contrasts matrix
      ## The dimenstions should be row = n.total and column = n.unphenotyped (aka the parents)
      c.mat <- make.contrast(unphenotyped.index = unphenotyped.index, n.total = n.total)
      
      # Calculate the inital PEVmean of the set
      C_22 = solve( (t(Z) %*% M %*% Z) + (lambda * A1.inv) )
      numerator = ( t(c.mat) %*% (A1 - (lambda * C_22)) %*% c.mat )
      denominator = t(c.mat) %*% A1 %*% c.mat
      CD.mat = numerator / denominator
      
      CD = diag(CD.mat)
      # Calculate the CDmean
      CDmean.save = mean(CD)
      
      # Create an empty vector to store CDmean values
      CD.vector <- vector("numeric", length = max.iter)
      
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
        
        # When using the whole candidates
        Z <- matrix(0, n.TP, n.total)
        # Fill in the coordinates of the starting TP (rows = index in the start.OTP, column = index in the A.mat)
        for (i in 1:length(new.OTP)) { Z[i, new.OTP[i] ] = 1 } 
        
        # Calculate the inital PEVmean of the set
        # No need to remake the contrast matrix because the unphenotyped lines remain unchanged
        C_22 = solve( (t(Z) %*% M %*% Z) + (lambda * A1.inv) )
        numerator = ( t(c.mat) %*% (A1 - (lambda * C_22)) %*% c.mat )
        denominator = t(c.mat) %*% A1 %*% c.mat
        CD.mat = numerator / denominator
        
        CD = diag(CD.mat)
        # Find the new mean
        CDmean.new <- mean(CD)
        
        # If the new CDmean is lower, accept, if not reject
        if (CDmean.new > CDmean.save) {
          start.OTP <- new.OTP
          CDmean.save <- CDmean.new
          # Create the new candidate set
          candidates.OTP <- setdiff(phenotyped.index, new.OTP)
        }
        
        # Record the CDmean of the new CDmean (if accepted) or the previous (if not)
        CD.vector[iter.counter - 1] <- CDmean.save
        
      } # Close the loop
      
    } # Close the use.subset if statement
    
  } # Close the CD optimization if statement
  
  ### Output the optimized TP
  # Extract line names from A.mat
  line.names <- row.names(A)
  
  OTP.lines <- line.names[start.OTP]
  
  output.list <- list(optimized.lines = OTP.lines,
                      PEVmean = list(PEVmean = PEVmean.save,
                                     PEV.vector = PEV.vector),
                      CDmean = list(CDmean = CDmean.save,
                                    CD.vector = CD.vector) )
  
  
  return(output.list)
  
} # Close the function