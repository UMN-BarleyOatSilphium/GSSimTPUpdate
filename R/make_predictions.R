#' Make genomic predictions
#' 
#' @description 
#' Uses a RR-BLUP mixed model to predict marker effects given phenotypic and
#' genotypic data from a training population. Then uses the marker effects
#' and genotypic data from selection candidates to predict genotypic values for
#' those candidates. The REML method of estimating variances is used.
#' 
#' @param pheno.train A matrix of phenotypes for the training population.
#' @param geno.train An incidence matrix of genotypes for the training population.
#' @param geno.pred An incidence matrix of genotypes for the prediction population
#' (i.e. selection candidates).
#' 
#' @return 
#' A matrix of GEBVs and the solution list to the mixed model
#' 
#' @import dplyr
#' @import rrBLUP
#' 
#' @export
#' 
make.predictions <- function(pheno.train,
                             geno.train,
                             geno.pred,
                             model = "RRBLUP") {
  
  # Deal with input
  pheno.train <- as.matrix(pheno.train)
  geno.train <- as.matrix(geno.train)
  geno.pred <- as.matrix(geno.pred)
  
  # Solve the mixed model
  solve.out <- mixed.solve(y = pheno.train, Z = geno.train, method = "REML")
  marker.effects <- solve.out$u
  
  # Calculate GEBVs
  GEBV <- geno.pred %*% marker.effects
  
  row.names(GEBV) <- row.names(geno.pred)
  
  # Create an output list
  output.list <- list(GEBV = GEBV, 
                      solve.out = solve.out)
  # Return the data
  return(output.list)
  
} # Close the function