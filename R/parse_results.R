#' Parse the results of the simulation
#' 
#' @description 
#' This function prases the R Data files created in interations
#' of the simulations and extract useful data
#' 
#' @param files A \code{character} vector of file paths. These files should be
#' \code{.RData} files.
#' @param filename A filepath to save the parsed, collective results
#' 
#' @import dplyr
#' 
#' @export
#' 
#' 
parse.results <- function(files, filename) {
  
  # Create an empty list
  collective.abbreviated.results <- list()
  
  # Loop over every file in files
  for (f in files) {
    
    # Load the file. This will be filled in later
    load(f)
    
    # Check to make sure that all of the elements of the results list are lists, too
    check.lists <- sapply(experiment.sub.results, is.list)
    
    # If one of the elements is not a list
    if (any(!check.lists)) {
      cat("\nOne of the elements in the results list is not a list. An error 
          likely occured.\n")
      cat("\nFilename: ", f)
      cat("\nElement index: ", which(!check.lists))
      cat("\nContents of element:\n")
      print(experiment.sub.results[!check.lists])
      
      # Remove the element from the list
      experiment.sub.results[!check.lists] <- NULL
      
    } # Close the if statement
      
    
    # Find the number of cycles
    n.cycles <- experiment.sub.results[[1]][[1]]$sim.results %>%
      length()
    
    # Variance components
    candidate.variance.components.list <- lapply(X = experiment.sub.results, FUN = function(set) 
      lapply(X = set, FUN = function(rep) 
        lapply(X = rep$sim.result, FUN = function(cycle) 
          return(cycle$candidate.values$true.var.components$V_g) )))
    
    # Genotypic value of candidates
    candidate.genotypic.value.list <- lapply(X = experiment.sub.results, FUN = function(set) {
      lapply(X = set, FUN = function(rep) {
        lapply(X = rep$sim.result, FUN = function(cycle) {
          return(cycle$candidate.values$mu.g) })})})
    
    # Variance components of the selections
    selection.variance.components.list <- lapply(X = experiment.sub.results, FUN = function(set) {
      lapply(X = set, FUN = function(rep) {
        lapply(X = rep$sim.result, FUN = function(cycle) {
          return(cycle$selection.values$V_g) })})})
    
    # Genotypic value of the selections
    selection.genotypic.value.list <- lapply(X = experiment.sub.results, FUN = function(set) {
      lapply(X = set, FUN = function(rep) {
        lapply(X = rep$sim.result, FUN = function(cycle) {
          return(cycle$selection.values$mu.g) })})})
    
    # Allele frequencies
    candidate.allele.freq.list <- lapply(X = experiment.sub.results, FUN = function(set) {
      lapply(X = set, FUN = function(rep) {
        lapply(X = rep$sim.result, FUN = function(cycle) {
          return(cycle$geno.summary.stats$candidate.maf) })})})
    
    TP.allele.freq.list <- lapply(X = experiment.sub.results, FUN = function(set) {
      lapply(X = set, FUN = function(rep) {
        lapply(X = rep$sim.result, FUN = function(cycle) {
          return(cycle$geno.summary.stats$TP.maf) })})})
    
    # Pairwise LD
    qtl.marker.LD.list <- lapply(X = experiment.sub.results, FUN = function(set) 
      lapply(X = set, FUN = function(rep) 
        lapply(X = rep$sim.result, FUN = function(cycle)
          list(mean.window = cycle$geno.summary.stats$qtl.marker.LD$mean.window,
               mean.max.window = cycle$geno.summary.stats$qtl.marker.LD$mean.max.window,
               mean.genome = cycle$geno.summary.stats$qtl.marker.LD$mean.genome,
               mean.max.genome = cycle$geno.summary.stats$qtl.marker.LD$mean.max.genome, 
               persistence = cycle$geno.summary.stats$qtl.marker.LD$persistance.of.phase ))))
    
    
    # Relationship of TP to the candidates
    relationship.list <- lapply(X = experiment.sub.results, FUN = function(set) {
      lapply(X = set, FUN = function(rep) {
        lapply(X = rep$sim.result, FUN = function(cycle) {
          return(cycle$geno.summary.stats$mu.TP.candidate.rel) })})})
    
    # Marker effects
    marker.effects.list <- lapply(X = experiment.sub.results, FUN = function(set) 
      lapply(X = set, FUN = function(rep) 
        lapply(X = rep$sim.result, FUN = function(cycle) 
          return(cycle$MM.solve$u) )))
    
    # Prediction accuracy results
    validation.results.list <- lapply(X = experiment.sub.results, FUN = function(set) {
      lapply(X = set, FUN = function(rep) {
        lapply(X = rep$sim.result, FUN = function(cycle) {
          return(cycle$prediction.accuracy) })})})
    
    # Training population updating - expected heterozygosity
    tp.update.exp.het.list <- lapply(X = experiment.sub.results, FUN = function(set) {
      lapply(X = set, FUN = function(rep) {
        lapply(X = rep$sim.result, FUN = function(cycle) {
          return(cycle$tp.update$Exp.het) })})})
    
    # The genome
    genome.list <- lapply(X = experiment.sub.results, FUN = function(set) {
      lapply(X = set, FUN = function(rep) {
        return(rep$genome) })})
    
    # Create a shorter output list
    abbreviate.output.list <- list(
      candidate.variance.components.list = candidate.variance.components.list,
      candidate.genotypic.value.list = candidate.genotypic.value.list,
      selection.variance.components.list = selection.variance.components.list,
      selection.genotypic.value.list = selection.genotypic.value.list,
      candidate.allele.freq.list = candidate.allele.freq.list,
      TP.allele.freq.list = TP.allele.freq.list,
      qtl.marker.LD.list = qtl.marker.LD.list,
      relationship.list = relationship.list,
      marker.effects.list = marker.effects.list,
      validation.results.list = validation.results.list,
      tp.update.exp.het.list = tp.update.exp.het.list,
      genome.list = genome.list,
      metadata = metadata
    )
    
    # Build a list
    collective.abbreviated.results[[change]] <- abbreviate.output.list
    
    
  } # Close the for loop
  
  # filename <- sub(pattern = ".RData", replacement = "_collective.RData", x = f)
  # Save the file
  save("collective.abbreviated.results", file = filename)
  
} # Close the function
