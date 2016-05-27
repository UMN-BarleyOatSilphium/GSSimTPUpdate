## Script to parse simulation results
# This script is meant to be run on MSI and deliver a saved ".RData" file

# Capture arguments. This will be a vector of filenames
args <- commandArgs(trailingOnly = TRUE)
# The first argument will be the output filename
filename <- args[1]

# Output list
collective.abbreviated.results <- list()


# Loop over every file in "args"
for (f in args[-1]) {
  
  # Load the file. This will be filled in later
  load(f)
  
  # Find the number of cycles
  n.cycles <- length(experiment.sub.results[[1]][[1]]$sim.results)
  
  # Variance components
  candidate.variance.components.list <- lapply(X = experiment.sub.results, FUN = function(set) {
    lapply(X = set, FUN = function(rep) {
      lapply(X = rep$sim.result, FUN = function(cycle) {
        return(cycle$candidate.values$var.components) })})})
  
  # Genotypic value of candidates
  candidate.genotypic.value.list <- lapply(X = experiment.sub.results, FUN = function(set) {
    lapply(X = set, FUN = function(rep) {
      lapply(X = rep$sim.result, FUN = function(cycle) {
        return(cycle$candidate.values$mu.g) })})})
  
  # Variance components of the selections
  selection.variance.components.list <- lapply(X = experiment.sub.results, FUN = function(set) {
    lapply(X = set, FUN = function(rep) {
      lapply(X = rep$sim.result, FUN = function(cycle) {
        return(cycle$selection.values$var.components) })})})
  
  # Genotypic value of the selections
  selection.genotypic.value.list <- lapply(X = experiment.sub.results, FUN = function(set) {
    lapply(X = set, FUN = function(rep) {
      lapply(X = rep$sim.result, FUN = function(cycle) {
        return(cycle$selection.values$mu.g) })})})
  
  # Allele frequencies
  allele.freq.list <- lapply(X = experiment.sub.results, FUN = function(set) {
    lapply(X = set, FUN = function(rep) {
      lapply(X = rep$sim.result, FUN = function(cycle) {
        return(cycle$geno.summary.stats$allele.freq) })})})
  
  # Pairwise LD
  qtl.marker.LD.list <- lapply(X = experiment.sub.results, FUN = function(set) 
    lapply(X = set, FUN = function(rep) 
      lapply(X = rep$sim.result, FUN = function(cycle)
        list(
          # Measure the average LD between each QTL and the marker in which it is in
          ## highest LD
          mean.max.LD = mean(apply(X = cycle$geno.summary.stats$qtl.marker.LD, MARGIN = 1, FUN = max)),
          # Measure the average LD between all QTL-marker pairs
          mean.LD = mean(cycle$geno.summary.stats$qtl.marker.LD) ))))
  
  
  # Relationship of TP to the candidates
  relationship.list <- lapply(X = experiment.sub.results, FUN = function(set) {
    lapply(X = set, FUN = function(rep) {
      lapply(X = rep$sim.result, FUN = function(cycle) {
        return(cycle$geno.summary.stats$mu.TP.candidate.rel) })})})
  
  # Heterzygosity 
  heterozygosity.list <- lapply(X = experiment.sub.results, FUN = function(set) {
    lapply(X = set, FUN = function(rep) {
      lapply(X = rep$sim.result, FUN = function(cycle) {
        return(cycle$geno.summary.stats$heterozygosity) })})})
  
  # Genomic prediction results
  prediction.results.list<- lapply(X = experiment.sub.results, FUN = function(set) {
    lapply(X = set, FUN = function(rep) {
      lapply(X = rep$sim.result, FUN = function(cycle) {
        return(cycle$prediction.results) })})})
  
  # Prediction accuracy results
  validation.results.list <- lapply(X = experiment.sub.results, FUN = function(set) {
    lapply(X = set, FUN = function(rep) {
      lapply(X = rep$sim.result, FUN = function(cycle) {
        return(cycle$prediction.accuracy) })})})
  
  # Training population updating
  tp.update.list <- lapply(X = experiment.sub.results, FUN = function(set) {
    lapply(X = set, FUN = function(rep) {
      lapply(X = rep$sim.result, FUN = function(cycle) {
        return(cycle$tp.update) })})})
  
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
    allele.freq.list = allele.freq.list,
    qtl.marker.LD.list = qtl.marker.LD.list,
    relationship.list = relationship.list,
    heterozygosity.list = heterozygosity.list,
    prediction.results.list = prediction.results.list,
    validation.results.list = validation.results.list,
    tp.update.list = tp.update.list,
    genome.list = genome.list
  )
  
  # Build a list
  collective.abbreviated.results[[change]] <- abbreviate.output.list
  # collective.abbreviated.results[["no.change_0.1"]] <- abbreviate.output.list
  
  
} # Close the for loop

# filename <- sub(pattern = ".RData", replacement = "_collective.RData", x = f)
# Save the file
save("collective.abbreviated.results", file = filename)
