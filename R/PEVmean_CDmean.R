#' PEVmean
#' 
#' @description 
#' Calculates the PEV mean or the CD mean from a contrasts matrix.
#' 
#' @param M Orthagonal projector. See \code{Details}
#' @param Z Incidence matrix
#' @param A Covariance matrix (i.e. additive relationship matrix)
#' @param A.inv The inverse of A (this is supplied to reduce computation time)
#' @param c Matrix of contrasts
#' @param V_e Error variance
#' @param V_a Additive genetic variance
#' 
#' @export
#' 
PEVmean <- function(M, Z, A, A.inv, c, V_e, V_a) {
  
  # Calculate lambda
  lambda <- V_e / V_a
  
  # Calculate the bottom-right matrix
  br <- (crossprod(Z, M) %*% Z) + (lambda * A.inv)
  # Calculate the inverse (C_22)
  C_22 <- GSsim.TPUpdate:::invert.mat(br, silent = T)
  
  numerator <- crossprod(c, C_22 ) %*% c
  denominator <- crossprod(c, c)
  PEV.mat <- numerator / denominator
  
  # Calculate the PEVmean and return it
  PEV <- diag(PEV.mat) * V_e
  return(mean(PEV))
  
} # Close the function

#' 
#' CDmean
#' 
#' @describeIn PEVmean
#'
#' @export
#' 
CDmean <- function(M, Z, A, A.inv, c, V_e, V_a) {
  
  # Calculate lambda
  lambda <- V_e / V_a
  
  # Calculate the bottom-right matrix
  br <- (crossprod(Z, M) %*% Z) + (lambda * A.inv)
  # Calculate the inverse (C_22)
  C_22 <- GSsim.TPUpdate:::invert.mat(br, silent = T)
  
  # Calculate the CDmean of the set
  numerator = crossprod(c, (A - (lambda * C_22 )) ) %*% c
  denominator = crossprod(c, A) %*% c
  CD.mat = numerator / denominator
  
  CD = diag(CD.mat)
  # Find the mean and return
  return(mean(CD))
  
} # Close the function
  
  
  
  
  