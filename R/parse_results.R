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
    

    # Genetic variance of the candidates
    candidate.gen.var <- lapply(X = experiment.sub.results, FUN = function(set) 
      lapply(X = set, FUN = function(rep) 
        lapply(X = rep$sim.result, FUN = function(cycle) 
          return(cycle$candidate.values$true.var.components$V_g) ))) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)
    
    # Genotypic value of candidates
    candidate.gen.val <- lapply(X = experiment.sub.results, FUN = function(set) {
      lapply(X = set, FUN = function(rep) {
        lapply(X = rep$sim.result, FUN = function(cycle) {
          return(cycle$candidate.values$mu.g) })})}) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)
    
    
    # Variance components of the selections
    selection.gen.var <- lapply(X = experiment.sub.results, FUN = function(set) {
      lapply(X = set, FUN = function(rep) {
        lapply(X = rep$sim.result, FUN = function(cycle) {
          return(cycle$selection.values$V_g) })})}) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)
    
    # Genotypic value of the selections
    selection.gen.val <- lapply(X = experiment.sub.results, FUN = function(set) {
      lapply(X = set, FUN = function(rep) {
        lapply(X = rep$sim.result, FUN = function(cycle) {
          return(cycle$selection.values$mu.g) })})}) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)
    
    # Allele frequencies
    candidate.allele.freq <- lapply(X = experiment.sub.results, FUN = function(set) {
      lapply(X = set, FUN = function(rep) {
        lapply(X = rep$sim.result, FUN = function(cycle) {
          return(cycle$geno.summary.stats$candidate.maf) })})}) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)
    
    TP.allele.freq <- lapply(X = experiment.sub.results, FUN = function(set) {
      lapply(X = set, FUN = function(rep) {
        lapply(X = rep$sim.result, FUN = function(cycle) {
          return(cycle$geno.summary.stats$TP.maf) })})}) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)
    
    # Pairwise LD
    qtl.marker.LD <- lapply(X = experiment.sub.results, FUN = function(set) 
      lapply(X = set, FUN = function(rep) 
        lapply(X = rep$sim.result, FUN = function(cycle)
          list(mean_window = cycle$geno.summary.stats$qtl.marker.LD$mean.window,
               mean_max_window = cycle$geno.summary.stats$qtl.marker.LD$mean.max.window,
               mean_genome = cycle$geno.summary.stats$qtl.marker.LD$mean.genome,
               mean_max_genome = cycle$geno.summary.stats$qtl.marker.LD$mean.max.genome, 
               persistence = cycle$geno.summary.stats$qtl.marker.LD$persistance.of.phase )))) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change) 
    
    
    # Relationship of TP to the candidates
    relationship <- lapply(X = experiment.sub.results, FUN = function(set) {
      lapply(X = set, FUN = function(rep) {
        lapply(X = rep$sim.result, FUN = function(cycle) {
          return(cycle$geno.summary.stats$mu.TP.candidate.rel) })})}) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change) 
    
    # Marker effects
    marker.effects <- lapply(X = experiment.sub.results, FUN = function(set) 
      lapply(X = set, FUN = function(rep) 
        lapply(X = rep$sim.result, FUN = function(cycle) 
          return(cycle$MM.solve$u) ))) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change) 
    
    # Prediction accuracy results
    validation.results <- lapply(X = experiment.sub.results, FUN = function(set) {
      lapply(X = set, FUN = function(rep) {
        lapply(X = rep$sim.result, FUN = function(cycle) {
          return(cycle$prediction.accuracy$pred.r) })})}) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change) 
    
    # Training population updating - expected heterozygosity
    tp.update.exp.het <- lapply(X = experiment.sub.results, FUN = function(set) {
      lapply(X = set, FUN = function(rep) {
        lapply(X = rep$sim.result, FUN = function(cycle) {
          return(cycle$tp.update$Exp.het) })})}) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)
    
    # The genome
    genome <- lapply(X = experiment.sub.results, FUN = function(set) {
      lapply(X = set, FUN = function(rep) {
        return(rep$genome) })})
    
    # Gather the data.frames / tibbles
    save.list <- list(
      candidate.gen.var = candidate.gen.var,
      candidate.gen.val = candidate.gen.val,
      selection.gen.var = selection.gen.var,
      selection.gen.val = selection.gen.val,
      candidate.allele.freq = candidate.allele.freq,
      TP.allele.freq = TP.allele.freq,
      qtl.marker.LD = qtl.marker.LD,
      relationship = relationship,
      marker.effects = marker.effects,
      validation.results = validation.results,
      validation.results = validation.results,
      tp.update.exp.het = tp.update.exp.het,
      genome = genome
    )
    
    # Build a list
    collective.abbreviated.results[[change]] <- save.list
    
    cat("File: ", f, " parsed.\n\n")
    
    
  } # Close the for loop
  
  # filename <- sub(pattern = ".RData", replacement = "_collective.RData", x = f)
  # Save the file
  save("collective.abbreviated.results", file = filename)
  
} # Close the function
