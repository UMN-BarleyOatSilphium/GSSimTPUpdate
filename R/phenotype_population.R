#' Generate phenotypic values
#' 
#' @description 
#' Generates value measures on the population of interest. The values can be the 
#' genotypic values (sum of the a values at each qtl) or the phenoypic values 
#' (genotypic values + error)
#' 
#' @param genome The hypred genome.
#' @param haploid.genos A 2n x m matrix of haploid alleles at m loci for n
#' individuals.
#' @param V_E The variance of the distribution from which to draw environmental
#' effects.
#' @param V_e The variance of the distribution from which to draw residual
#' effects.
#' @param just.geno.values Logical whether to just return the genotypic values.
#' @param n.evn The number of environments to "phenotype" the population.
#' @param n.rep The number of replicates of each individual in each environment.
#' @param run.anova Logical whether to run an anova of the simulated phenotypes.
#' This is useful to verify that the target heritability is being achieved.
#' 
#' @import dplyr
#' @import stringr
#' 
#' @export
#' 
phenotype.population <- function(genome, haploid.genos, V_E, V_e, # Residual variance
                                 just.geno.values = FALSE, n.env, n.rep, 
                                 run.anova = FALSE ){
  
  # Deal with input
  if(!is.list(haploid.genos)) {
    haploid.mat <- as.matrix(haploid.genos)
  } else {
    haploid.mat <- do.call("rbind", haploid.genos)
  }
  
  # Pull out the line names
  line.names <- row.names(haploid.mat)
  # Condense the names
  line.names1 <- line.names %>%
    str_replace(pattern = "\\.[0-9]$", replacement = "") %>%
    unique()
  
  # First generate genotype values based on the genome and genotypes
  g <- genotypic.value(genome = genome, haploid.genos = haploid.mat)
  row.names(g) <- line.names1
  # Find the length (i.e. number of genotypes)
  n.geno = length(g)
  
  # Find the genotype variance
  V_g <- var(g)
  
  # If just genotypic values are requested, stop and return a list
  if (just.geno.values) {
    output.list <- list(geno.values = g,
                        var.components = list(V_g = V_g))
    return(output.list)
  }
  
  # Sample the environmental effects
  e = rnorm(n.env, 0, sqrt(V_E))
  
  # Sample the residuals
  # Remember this includes error and gxe effect
  residual <- matrix(rnorm(n.geno * n.env * n.rep, 0, sqrt(V_e)), n.geno, n.env * n.rep)
  
  # Calculate phenotypic values
  y <- matrix(g, n.geno, (n.env * n.rep)) + 
    matrix(e, n.geno, (n.env * n.rep), byrow = T) + 
    residual
  
  # Convert to vector
  y <- as.numeric(y)
  
  # Reform into data.frame
  phenos <- data.frame(line = rep(line.names1, n.env * n.rep),
                       env = rep( rep( paste("env", seq(n.env), sep = ""), each = n.geno ), n.rep ),
                       rep = rep( rep( paste("rep", seq(n.rep), sep = ""), each = n.geno ), each = n.rep),
                       value = y )
  
  # Run a analysis of variance, if desired
  if (run.anova) {
    
    ## For unreplicated, multi-environment
    if (n.env > 1 & n.rep == 1) {
      # Fit the model
      lm.fit <- lm(formula = value ~ line + env, data = phenos)
      lm.anova <- anova(lm.fit)
      
      # Mean squares
      MS_error <- lm.anova["Residuals",]$`Mean Sq`
      MS_line <- lm.anova["line",]$`Mean Sq`
      
      # Variance components
      V_error <- MS_error
      V_line <- (MS_line - MS_error) / n.env
      
      # Heritability
      H = V_line / (V_line + (V_error / n.env))
      
    }
    
    if (n.env > 1 & n.rep > 1) {
      # Fit the model
      lm.fit <- lm(formula = value ~ line:env, data = phenos)
      lm.anova <- anova(lm.fit)
      
      # Mean squares
      MS_error <- lm.anova["Residuals",]$`Mean Sq`
      MS_line <- lm.anova["line",]$`Mean Sq`
      MS_linebyenv <- lm.anova["line:env"]$`Mean Sq`
      
      # Variance components
      V_error <- MS_error
      V_linebyenv = (MS_linebyenv - MS_error) / n.rep
      V_line = (MS_line - MS_linebyenv) / (n.rep * n.env)
      
      # Heritability
      H = V_line / (V_line + (V_linebyenv / n.rep) + (V_error / (n.rep * n.env)))
    }
    
    # Else NAs for the observed variance components
  } else {
    V_line = NA
    V_error = NA
    H = NA
    
  } # Close the run anova if statement
  
  # Find the average across environments as the value to use
  mu_y <- tapply(X = phenos$value, INDEX = phenos$line, FUN = mean) %>%
    as.matrix()
  
  # Find the average across genotypes
  mu.p <- mean(mu_y) # Phenotypic value
  mu.g <- mean(g) # Genotypic values
  
  
  # Output list
  output.list <- list(full.matrix = phenos,
                      geno.values = g,
                      mean.pheno.values = mu_y,
                      mu.p = mu.p,
                      mu.g = mu.g,
                      true.var.components = list(V_g = V_g, V_E = V_E, V_e = V_e),
                      obs.var.components = list(V_g = V_line, V_e = V_error),
                      obs.H = H)
  
  return(output.list)
  
} # Close the function
