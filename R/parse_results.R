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
      qtl.eff(rep$genome)$add ) %>%
      as.data.frame() %>%
      tbl_df() %>%
      gather(iter, value) %>%
      mutate(iter = iter %>% as.factor() %>% as.numeric())
      
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
  
  ## Create data.frames for plotting
  
  # Change in accuracy
  plot.list[["df.acc"]] <- lapply(X = total.names, FUN = function(coll.name) 
    lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$validation.results) %>%
      bind_rows() %>%
      mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
               str_to_title()) ) %>%
    bind_rows() %>%
    sim.summarize() %>%
    mutate(variable = "Prediction Accuracy")
  
  # Change in Genetic Variance
  plot.list[["df.genvar"]] <- lapply(X = total.names, FUN = function(coll.name) 
    lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$sc.gen.var) %>%
      bind_rows() %>%
      mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
               str_to_title()) ) %>%
    bind_rows() %>% 
    sim.summarize() %>%
    mutate(variable = "Genetic Variance")
  
  
  # Change in mean genotypic value of the selection candidates
  plot.list[["df.genval"]] <- lapply(X = total.names, FUN = function(coll.name) 
    lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$sc.gen.val) %>%
      bind_rows() %>%
      mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
               str_to_title()) ) %>%
    bind_rows() %>%
    sim.summarize() %>%
    mutate(variable = "Genotypic Value")
  
  
  # Response to selection
  plot.list[["df.resp"]] <- lapply(X = total.names, FUN = function(coll.name) 
    lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$sc.gen.val) %>%
      bind_rows() %>%
      mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
               str_to_title()) ) %>%
    bind_rows() %>% 
    group_by(exp_name, change, iter, cycle) %>% 
    filter(row_number() == 1) %>%
    group_by(exp_name, change, iter) %>%
    mutate(value = c(NA, diff(value))) %>% 
    na.omit() %>%
    ungroup() %>%
    sim.summarize() %>%
    # Add the variable designator
    mutate(variable = "Response to Selection")
  
  # LD
  df <- lapply(X = total.names, FUN = function(coll.name) 
    lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$qtl.marker.LD) %>%
      bind_rows() %>%
      mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
               str_to_title()) ) %>%
    bind_rows()
  
  # Mean max LD in TP
  plot.list[["df.tpmeanmax"]] <- sim.summarize(df %>% filter(extra1 == "tp_mean_max_genome")) %>%
    mutate(variable = "Mean Max LD in Training Population")
  
  # Mean max LD in SC
  plot.list[["df.scmeanmax"]] <- sim.summarize(df %>% filter(extra1 == "sc_mean_max_genome")) %>%
    mutate(variable = "Mean Max LD in Selection Candidates")
  
  # Persistence of phase
  plot.list[["df.pers"]] <- sim.summarize(df %>% filter(extra1 == "persistence")) %>%
    mutate(variable = "Persistence of\nLD Phase")
  
  
  # Genomic relationship
  plot.list[["df.rel"]] <- lapply(X = total.names, FUN = function(coll.name) 
    lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$tp.sc.relationship) %>%
      bind_rows() %>%
      mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
               str_to_title()) ) %>%
    bind_rows() %>%
    sim.summarize() %>%
    mutate(variable = "Average Relationship")
  
  # Inbreeding
  plot.list[["df.inbred"]] <- lapply(X = total.names, FUN = function(coll.name) 
    lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$sc.inbreeding) %>%
      bind_rows() %>%
      mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
               str_to_title()) ) %>%
    bind_rows() %>%
    sim.summarize() %>%
    mutate(variable = "Inbreeding")
  
  # Rate of inbreeding
  plot.list[["df.rateinbred"]] <- lapply(X = total.names, FUN = function(coll.name) 
    lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$sc.inbreeding) %>%
      bind_rows() %>%
      mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
               str_to_title()) ) %>%
    bind_rows() %>% 
    group_by(exp_name, change, iter, cycle) %>% 
    filter(row_number() == 1) %>%
    group_by(exp_name, change, iter) %>%
    mutate(value = c(NA, diff(value))) %>% 
    na.omit() %>%
    ungroup() %>%
    sim.summarize() %>%
    mutate(variable = "Rate of Inbreeding")
  
  
  # QTL fixation
  plot.list[["df.fixedqtl"]] <- lapply(X = total.names, FUN = function(coll.name) 
    lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) tpc$sc.allele.freq) %>%
      bind_rows() %>%
      filter(extra1 == "qtl") %>% 
      group_by(change, iter, cycle) %>% 
      summarize(value = sum(value == 0 | value == 1)) %>%
      mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
               str_to_title()) ) %>%
    bind_rows() %>%
    ungroup() %>%
    sim.summarize()
  
  # QTL fixed for favorable allele
  plot.list[["df.fixedqtlfav"]] <- lapply(X = total.names, FUN = function(coll.name) 
    lapply(X = total.collective.data[[coll.name]], FUN = function(tpc) {
      # Gather the allele frequencies
      freq <- tpc$sc.allele.freq %>%
        filter(extra1 == "qtl")
      # Gather the effects
      eff <- tpc$qtl.effects
      # Add the effects to the frequencies
      freq$effect <- eff$value
      return(freq) }) %>%
      bind_rows() %>%
      mutate(exp_name = str_extract(string = coll.name, pattern = 'window|cumulative') %>% 
               str_to_title()) ) %>%
    bind_rows() %>%
    # Find QTL fixed for favorable allele
    group_by(exp_name, change, iter, cycle) %>% 
    filter((value == 0 & effect < 0) | (value == 1 & effect > 0) ) %>%
    summarize(value = n())
  
  filename2 <- sub(pattern = ".RData", replacement = "_plotdata.RData", x = filename)
  # Save the files
  save("collective.abbreviated.results", file = filename)
  save("plot.list", file = filename2)
  
} # Close the function
