#' Subset phenotypic or genotypic values
#' 
#' @param pheno.values.list The list returned when using the function
#' \code{\link{GSsim.TPUpdate}[phenotype.population]}.
#' @param line.names Individuals for which to extract values.
#' 
#' @export
#' 
select.values <- function(pheno.values.list, line.names) {
  
  # Subset pheno and geno values
  subset.pheno.values <- as.matrix(pheno.values.list$mean.pheno.values[line.names,])
  subset.geno.values <- as.matrix(pheno.values.list$geno.values[line.names,])
  
  # Calculate mean pheno and geno values
  subset.mu.p <- mean(subset.pheno.values)
  subset.mu.g <- mean(subset.geno.values)
  
  # Calculate genetic variance of the subset
  subset.V_g <- var(subset.geno.values)
  
  # Package and return a list
  subset.list <- list(geno.values = subset.geno.values,
                      mean.pheno.values = subset.pheno.values,
                      mu.p = subset.mu.p,
                      mu.g = subset.mu.g,
                      V_g = subset.V_g)
  return(subset.list)
  
} # Close the function