#' Parse the results of the simulation
#' 
#' @description 
#' This function prases the R Data files created in interations
#' of the simulations and extract useful data
#' 
#' @param files A \code{character} vector of file paths. These files should be
#' \code{.RData} files.
#' @param filename A filepath to save the parsed, collective results
#' @param max.reps The maximum number of replications to keep. Takes the first
#' \code{max.reps} replications.
#' 
#' @import dplyr
#' @import tidyr
#' 
#' @export
#' 
#' 
parse.results <- function(files, filename, max.reps) {
  
  # Create an empty list
  collective.abbreviated.results <- list()
  
  # Loop over every file in files
  for (f in files) {
    
    # Load the file. This will be filled in later
    load(f)
    
    # Unlist the first layer of the results list
    experiment.sub.results <- experiment.sub.results %>% 
      unlist(recursive = F)
    
    # Remove any rep results that are not lists
    check.lists <- sapply(X = experiment.sub.results, FUN = is.list)
    experiment.sub.results[!check.lists] <- NULL
    
    if (sum(!check.lists) > 0) {
    
      # Report the number of reps that were discarded
      cat("\n\nReps were removed from the file.")
      cat("\n\nFilename: ", f)
      cat("\nNumber of reps removed: ", sum(!check.lists))
      
    }
    
    # Find the number of cycles
    n.cycles <- experiment.sub.results[[1]]$sim.results %>%
      length()
    
    # Determine if the results were the Cumulative or Window
    tp.formation <- f %>% 
      str_extract(pattern = 'window|cumulative') %>% 
      str_to_title()
      
    
    # Determine if the max.iters are greater than the total number of replications
    if (max.reps > length(experiment.sub.results))
      stop("The max.reps arguments is greater than the total number of replications.")
    
    # Subet the first 'max.reps' replication
    experiment.sub.results <- experiment.sub.results[seq_len(max.reps)]
    
    # Create an empty list to save
    save.list <- list()
    plot.list <- list()
    

    # Genetic variance of the candidates
    save.list[["sc.gen.var"]] <- lapply(X = experiment.sub.results, FUN = function(rep)
        lapply(X = rep$sim.result, FUN = function(cycle)
          return(cycle$candidate.values$true.var.components$V_g) )) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)

    # Genotypic value of candidates
    save.list[["sc.gen.val"]] <- lapply(X = experiment.sub.results, FUN = function(rep)
      lapply(X = rep$sim.result, FUN = function(cycle)
          return(cycle$candidate.values$mu.g) )) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)


    # Allele frequencies
    save.list[["sc.allele.freq"]] <- lapply(X = experiment.sub.results, FUN = function(rep)
      lapply(X = rep$sim.result, FUN = function(cycle)
          list(qtl = cycle$geno.summary.stats$candidate.af$qtl,
               markers = cycle$geno.summary.stats$candidate.af$snp) )) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)
    

    save.list[["TP.allele.freq"]] <- lapply(X = experiment.sub.results, FUN = function(rep)
      lapply(X = rep$sim.result, FUN = function(cycle)
          list(qtl = cycle$geno.summary.stats$TP.af$qtl,
               markers = cycle$geno.summary.stats$TP.af$snp) )) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)

    # Pairwise LD
    save.list[["qtl.marker.LD"]] <- lapply(X = experiment.sub.results, FUN = function(rep)
      lapply(X = rep$sim.result, FUN = function(cycle)
          list(sc_mean_max_genome = cycle$geno.summary.stats$qtl.marker.LD$sc.mean.max.genome,
               tp_mean_max_genome = cycle$geno.summary.stats$qtl.marker.LD$tp.mean.max.genome,
               persistence = cycle$geno.summary.stats$qtl.marker.LD$persistance.of.phase ))) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)
    
    # Number of qtl and markers used to measure LD
    save.list[["n.loci.LD"]] <- lapply(X = experiment.sub.results, FUN = function(rep)
      lapply(X = rep$sim.result, FUN = function(cycle)
          list(n_qtl = cycle$geno.summary.stats$qtl.marker.LD$n.qtl.LD,
               n_markers = cycle$geno.summary.stats$qtl.marker.LD$n.marker.LD) )) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)

    # Prediction accuracy results
    save.list[["validation.results"]] <- lapply(X = experiment.sub.results, FUN = function(rep)
      lapply(X = rep$sim.result, FUN = function(cycle)
          return(cycle$prediction.accuracy) )) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)

    # Training population updating - expected heterozygosity
    save.list[["tp.update.exp.het"]] <- lapply(X = experiment.sub.results, FUN = function(rep)
      lapply(X = rep$sim.result, FUN = function(cycle)
          return(cycle$tp.update$Exp.het) )) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)
    
    # Determine the effect of the 1 allele for each QTL over iterations
    save.list[["qtl.effects"]] <- lapply(X = experiment.sub.results, FUN = function(rep)
      qtl.eff(rep$genome)$add) %>%
      as.data.frame() %>%
      tbl_df() %>%
      gather(iter, value) %>%
      mutate(iter = iter %>% as.factor() %>% as.numeric() ) 
      
    ## Inbreeding
    save.list[["sc.inbreeding"]] <- lapply(X = experiment.sub.results, FUN = function(rep)
      lapply(X = rep$sim.result, FUN = function(cycle)
          return(cycle$inbreeding$candidates) )) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)
    
    save.list[["tp.additions.inbreeding"]] <- lapply(X = experiment.sub.results, FUN = function(rep)
      lapply(X = rep$sim.result, FUN = function(cycle)
          return(cycle$inbreeding$TP.additions) )) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)
    
    ## Relationships
    save.list[["tp.sc.relationship"]] <- lapply(X = experiment.sub.results, FUN = function(rep)
      lapply(X = rep$sim.result, FUN = function(cycle)
          return(cycle$relationship$TP.candidates) )) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)
    
    save.list[["tp.additions.relationship"]] <- lapply(X = experiment.sub.results, FUN = function(rep)
      lapply(X = rep$sim.result, FUN = function(cycle)
          return(cycle$relationship$TP.additions) )) %>%
      unlist() %>%
      GSsim.TPUpdate:::nv_df(change = change)
    
    
    # # Results of the persistence of LD permutation tests
    # save.list[["permutation.results"]] <- lapply(X = experiment.sub.results, FUN = function(rep)
    #   lapply(X = rep$sim.result, FUN = function(cycle)
    #       return(cycle$geno.summary.stats$qtl.marker.LD$persistence.perm.results) ))
  

    # Build a list
    collective.abbreviated.results[[change]] <- save.list
    
    cat("\n\nFile: ", f, " parsed.\n\n")
    
    
  } # Close the for loop

  # Save the files
  save("collective.abbreviated.results", file = filename)

} # Close the function
