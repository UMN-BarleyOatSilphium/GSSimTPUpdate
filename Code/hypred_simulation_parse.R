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
  
  # Allele frequencies
  allele.freq.list <- lapply(X = experiment.sub.results, FUN = function(set) {
    lapply(X = set, FUN = function(rep) {
      lapply(X = rep$sim.result, FUN = function(cycle) {
        return(cycle$geno.summary.stats$allele.freq) })})})
  
  pairwise.div.list <- lapply(X = experiment.sub.results, FUN = function(set) {
    lapply(X = set, FUN = function(rep) {
      lapply(X = rep$sim.result, FUN = function(cycle) {
        return(cycle$geno.summary.stats$pairwise.div) })})})
  
  heterozygosity.list <- lapply(X = experiment.sub.results, FUN = function(set) {
    lapply(X = set, FUN = function(rep) {
      lapply(X = rep$sim.result, FUN = function(cycle) {
        return(cycle$geno.summary.stats$heterozygosity) })})})
  
  genome.list <- lapply(X = experiment.sub.results, FUN = function(set) {
    lapply(X = set, FUN = function(rep) {
      return(rep$genome) })})
  
  # pheno.value.list <- lapply(X = experiment.sub.results, FUN = function(set) {
  #   lapply(X = set, FUN = function(rep) {
  #     lapply(X = rep$sim.result, FUN = function(cycle) {
  #       return(cycle$selection.pheno.values) })})})
  
  prediction.results.list<- lapply(X = experiment.sub.results, FUN = function(set) {
    lapply(X = set, FUN = function(rep) {
      lapply(X = rep$sim.result, FUN = function(cycle) {
        return(cycle$prediction.results) })})})
  
  validation.results.list <- lapply(X = experiment.sub.results, FUN = function(set) {
    lapply(X = set, FUN = function(rep) {
      lapply(X = rep$sim.result, FUN = function(cycle) {
        return(cycle$prediction.accuracy) })})})
  
  
  
  # Create a shorter output list
  abbreviate.output.list <- list(
    allele.freq.list = allele.freq.list,
    pairwise.div.list = pairwise.div.list,
    heterozygosity.list = heterozygosity.list,
    genome.list = genome.list,
    prediction.results.list = prediction.results.list,
    validation.results.list = validation.results.list,
    # pheno.value.list = pheno.value.list,
    candidate.variance.components.list = candidate.variance.components.list,
    candidate.genotypic.value.list = candidate.genotypic.value.list
  )
  
  # Build a list
  collective.abbreviated.results[[change]] <- abbreviate.output.list
  
} # Close the for loop

# Save the file
save("collective.abbreviated.results", file = filename)
