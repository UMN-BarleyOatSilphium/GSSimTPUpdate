## Genomic selection simulations
# This code runs iterations of long-term genomic selection simulations to 
## test the usefulness of updating a training population and how it impacts
## prediction accuracy, genetic variance, and selection potential

# Arguments
args <- commandArgs(trailingOnly = T)

# First and only argument is the pop makeup (MN, ND, or MNxND)
# Second argument is how the TP should be combined after each cycle (cumulative or window)
if (all(is.na(args))) {
  pop.makeup <- "MN"
  tp.formation <- "cumulative"
  parents.sel.intensity = 100
  n.crosses = 50
  MSI = F
} else {
  pop.makeup <- args[1]
  tp.formation <- args[2]
  parents.sel.intensity <- as.numeric(args[3])
  n.crosses = as.numeric(args[4])
  MSI = T
}

# Load the packages
library(hypred, quietly = T)
library(rrBLUP, quietly = T)
library(boot, quietly = T)
library(parallel, quietly = T)

# # Other tools
if (MSI) {
  setwd("/panfs/roc/groups/6/smithkp/neyhartj/Genomic_Selection/Simulations/BarleySimGS-TPUpdate")
  # source("/panfs/roc/groups/6/smithkp/neyhartj/GitHub_Repos/Quant-Gen-Scripts/genotype_matrix_utils.R")
  n.cores = 16
} else {
  setwd("C:/Users/Jeff/Google Drive/Barley Lab/Projects/Side Projects/Simulations/BarleySimGS-TPUpdate/")
  n.cores = 1
}

# Load already-curated gamete data
load(file = "Files/Barley_CAP_simuation_starting_material.RData")

# Source functions
source("Code/simulation_functions.R")


# Other simulation parameters
# Entry-mean heritability
h2 = 0.5
# How many cycles?
n.cycles = 15
# Number of QTL underlying trait
n.QTL = 100
# Number of phenotyping environments and reps
n.env = 3
n.rep = 1

# Minor allele frequency cut-off for markers
min.maf = 0.03

# Barley population genetics data
mutation.rate.snp = 7e-8
mutation.rate.qtl = 7e-8

tp.change = c("best", "worst", "random", "no.change", "PEVmean", "CDmean")
# The number of lines to add the TP after each cycle
tp.update.increment = 150
# Size of the TP to maintain - this is the same as the starting TP
tp.size <- nrow(CAP.haploids) / 2

# Parent selection and crossing parameters
ind.per.cross = 30
cycle.candidate.size = n.crosses * ind.per.cross

# Standardized selection intensity
std.sel.intensity = parents.sel.intensity / cycle.candidate.size

# Computation parameters
n.iterations = 100

date <- format(Sys.time(), "%d%m%y-%H%M%S")

# Save the metadata to a list
metadata <- list(h2 = h2,
                 n.cycles = n.cycles,
                 n.QTL = 100,
                 min.marker.maf = min.maf,
                 parents.sel.intensity = parents.sel.intensity,
                 n.env = n.env, 
                 n.rep = n.rep,
                 mutation.rate.qtl = mutation.rate.qtl,
                 mutation.rate.snp = mutation.rate.snp,
                 tp.update.increment = tp.update.increment,
                 tp.size = tp.size,
                 n.crosses = n.crosses,
                 ind.per.cross = ind.per.cross,
                 cycle.candidate.size = cycle.candidate.size,
                 std.sel.intensity = std.sel.intensity,
                 n.iterations = n.iterations,
                 pop.makeup = pop.makeup,
                 date = date)


#### Define genome characteristics ####
# Find the snps per chromsome
n.chr.snps <- tapply(X = CAP.markers$rs, INDEX = CAP.markers$chrom, length)

# Find the chromsome lengths
chr.len <- as.numeric(tapply(X = CAP.markers$pos, INDEX = CAP.markers$chrom, FUN = max))

# Create a list of loci positions
genetic.map.list <- tapply(X = CAP.markers$pos, INDEX = CAP.markers$chrom, FUN = function(chr) list(chr))

# Make the initial genome
hv.genome <- make.genome( n.chr = length(chr.len), 
                          chr.len = chr.len, 
                          n.chr.snps = n.chr.snps,
                          genetic.map = genetic.map.list)



# Iterate over the different tp.change types
# for (change in tp.change) {
change = "worst"
  
  # Split iterations into cores
  if (n.cores > 1) {
    iters.per.core <- split(x = 1:n.iterations, factor(cut(x = 1:n.iterations, breaks = n.cores)))
    names(iters.per.core) <- paste("set", 1:length(iters.per.core), sep = "")
  } else {
    iters.per.core <- 1:n.iterations
  }
  
  # # Apply the iterations over cores
  # experiment.sub.results <- mclapply(X = iters.per.core, FUN = function(iter.set) {
  #   
  #   # Loop over each iteration
  #   lapply(X = 1:length(iter.set), FUN = function(rep.i) {
      
  ### DEBUGGING ###
  
  for (rep.i in 1:n.iterations) {    
      
      # All code below this line is variable in each iteration of the simulation
      
      #### Define trait parameters ####
      hv.genome <- trait.architecture(genome = hv.genome,
                                      n.QTL = n.QTL, 
                                      qtl.index = NULL, 
                                      qtl.dom.index = NULL, 
                                      qtl.perf.index = NULL, 
                                      qtl.add.eff = "geometric", 
                                      qtl.dom.eff = NULL)
      
      TP.haploids.i <- CAP.haploids
      # Convert the gametes to genotypes
      TP.genos <- genotype.loci(haploid.genos = TP.haploids.i, genome = hv.genome)
      
      # Find the MAF of all marker snps
      marker.maf <- measure.maf(geno.mat = TP.genos)
      # Determine which are below 
      markers.below.maf <- which(marker.maf < min.maf)
      
      
      # Calculate the frequency of the 1 allele in the base training population
      p_i <- apply(X = TP.genos + 1, MARGIN = 2, FUN = mean) / 2
      # Calculate the P matrix
      P = matrix(2 * (p_i - 0.5))
      # Calculate the normalization constant
      c = 2 * sum(p_i * (1 - p_i))
      
      ## Select the TP lines for use as parents
      # First separate MN and ND lines
      line.names <- row.names(TP.genos)
      ND.lines <- grep(pattern = "^ND", x = line.names, value = T)
      MN.lines <- setdiff(x = line.names, y = c(ND.lines))
      
      ## Set the inital variances for the heritability
      # True genetic variance
      TP.V_g <- var(genotypic.value(genome = hv.genome, haploid.genos = TP.haploids.i))
      # Environmental variance (scale * 8 as in Bernardo 2015)
      V_E <- TP.V_g * 8
      # Residual variance
      V_e <- n.rep * n.env * ((TP.V_g / h2) - TP.V_g)
      
      
      # Phenotype the training population
      TP.values <- phenotype.population(genome = hv.genome,
                                        haploid.genos = TP.haploids.i,
                                        V_E = V_E,
                                        V_e = V_e,
                                        n.env = n.env,
                                        n.rep = n.rep)
      
      TP.phenos <- TP.values$mean.pheno.values
      
      # Next select the top MN and top ND
      top.MN.lines <- names(sort(TP.phenos[MN.lines,], decreasing = T)[1:parents.sel.intensity])
      top.ND.lines <- names(sort(TP.phenos[ND.lines,], decreasing = T)[1:parents.sel.intensity])
      
      # Create a list to store the line names
      pop.makeup.list <- list(
        MN = list(p1 = top.MN.lines, p2 = top.MN.lines),
        ND = list(p1 = top.ND.lines, p2 = top.ND.lines),
        MNxND = list(p1 = top.MN.lines[1:(parents.sel.intensity / 2)], p2 = top.ND.lines[1:(parents.sel.intensity / 2)])
      )
      
      # Set the parent gamete data input
      parent.lines.list <- pop.makeup.list[[pop.makeup]]
      parent.haploids <- subset.gametes(gametes = TP.haploids.i, line.names = unique(unlist(parent.lines.list)))
      
      # Set dummy variables for the phenos and genos
      TP.phenos.i <- TP.phenos
      TP.genos.i <- TP.genos
      
      # Create an initial data list
      simulation.results <- list()
      debug.results <- list()
      
      # Loop over the number of cycles
      for (breeding.cycle in 1:n.cycles) {

        ##### Start the Cycle Executions #####

        ##### Step 1 - Crossing and inbreeding
        # Make a crossing block
        crossing.block.i <- make.crossing.block(parent1.lines = parent.lines.list$p1, 
                                                parent2.lines = parent.lines.list$p2, 
                                                n.crosses = n.crosses, 
                                                use.parents.once = T)
        
        # According to the crossing block, parents are crossed to form F1s, then
        ## the progeny are inbred to the F3 generation. Since each F1 plant is inbred
        ## individually, the resulting families consist of F1:3 lines
        candidate.haploid.i <- make.population(genome = hv.genome, 
                                               parental.haploids = parent.haploids,
                                               crossing.block = crossing.block.i,
                                               N = ind.per.cross,
                                               cycle.number = breeding.cycle,
                                               generations = 2,
                                               pop.type = "inbred",
                                               mutation.rate.snp = mutation.rate.snp,
                                               mutation.rate.qtl = mutation.rate.qtl)
        
        ##### Step 2 - Genotype
        # Find the genotypes of the markers and QTL
        candidate.genos.i <- genotype.loci(haploid.genos = candidate.haploid.i, 
                                           genome = hv.genome, 
                                           include.QTL = T)
        # Just the marker genotypes
        candidate.marker.genos.i <- genotype.loci(haploid.genos = candidate.haploid.i, 
                                                  genome = hv.genome, 
                                                  include.QTL = F)
        
        
        ##### Step 3 - Genotypic Summary Statistics
        # Measure the minor allele frequency of all loci
        candidate.genos.maf.i <- measure.maf(geno.mat = candidate.genos.i)
        # Measure the maf of just the markers
        candidate.marker.maf.i <- measure.maf(geno.mat = candidate.marker.genos.i)
        TP.genos.maf.i <- measure.maf(geno.mat = TP.genos.i)
        
        ### LD measures
        # Candidates
        candidate.qtl.marker.LD.i <- measure.LD(genome = hv.genome, 
                                                genos = candidate.genos.i, 
                                                Morgan.window = 0.5)
        
        # Per QTL, find the LD of the marker (within the window) with which the 
        ## QTL has the highest LD then find the mean of those values across 
        ## polymorphic qtl
        candidate.mean.max.LD.i <- mean(unlist(lapply(X = candidate.qtl.marker.LD.i, FUN = function(chr)
          lapply(X = chr, FUN = function(qtl) max(qtl^2)) )), na.rm = T)
        
        # Now just find the mean of LD across the whole window
        candidate.mean.window.LD.i <- mean(unlist(lapply(X = candidate.qtl.marker.LD.i, FUN = function(chr)
          lapply(X = chr, FUN = function(qtl) mean(qtl^2)) )), na.rm = T)
        
        # Measure LD on the TP
        TP.qtl.marker.LD.i <- measure.LD(genome = hv.genome, 
                                         genos = TP.haploids.i, 
                                         Morgan.window = 0.5)
        
        # Apply over chromosomes and combine the data
        TP.candidate.LD.i <- do.call("rbind", lapply(X = 1:length(TP.qtl.marker.LD.i), FUN = function(i) {
          # Find the common polymorphic QTL
          common.poly.qtl <- intersect(names(TP.qtl.marker.LD.i[[i]]), names(candidate.qtl.marker.LD.i[[i]]))
          # Apply over the common QTL and combine the data
          do.call("rbind", lapply(X = common.poly.qtl, FUN = function(poly.qtl) {
            # Find the common polymorphic markers within the polymorphic qtl
            common.poly.markers <- intersect(names(TP.qtl.marker.LD.i[[i]][[poly.qtl]]), names(candidate.qtl.marker.LD.i[[i]][[poly.qtl]]))
            # Return r for those markers
            data.frame( TP.LD = TP.qtl.marker.LD.i[[i]][[poly.qtl]][common.poly.markers], cand.LD = candidate.qtl.marker.LD.i[[i]][[poly.qtl]][common.poly.markers])
          })) }) )
        
        # If the data.frame has no data, return NA
        if (nrow(TP.candidate.LD.i) == 0) {
          TP.candidate.persistance.of.phase.i <- NA
        } else {
          # Find the correlation of r across all correlations
          TP.candidate.persistance.of.phase.i <- cor(TP.candidate.LD.i$TP.LD, TP.candidate.LD.i$cand.LD)
        }
          
        # Create a list to save
        qtl.marker.LD.i <- list(candidate.i.qtl.marker.LD = candidate.qtl.marker.LD.i,
                                TP.i.qtl.marker.LD = TP.qtl.marker.LD.i,
                                mean.max = candidate.mean.max.LD.i,
                                mean.window = candidate.mean.window.LD.i,
                                persistance.of.phase = TP.candidate.persistance.of.phase.i)
        
        ### Measure the average relationship between the TP and the candidates
        # Assign M
        M <- rbind(TP.genos.i, candidate.marker.genos.i)
        # Subtract P to make Z (need to convert P into a repeated matrix)
        W = M - matrix(P, nrow(M), length(P), byrow = T)
        # Calculate the relationship matrix
        A = tcrossprod(W) / c
        
        # Calculate the mean relationship between the TP and the candidates
        mu.relationship <- mean( rowMeans(A[row.names(TP.genos.i), row.names(candidate.marker.genos.i)]) )
        
        
        ##### Step 4 - Prediction
        # Remove the markers with maf below the threshold (set at the start of the sim)
        ## also remove monomorphic markers
        markers.to.remove <- sort(unique(c(
          which(candidate.marker.maf.i == 0),
          which(TP.genos.maf.i == 0),
          markers.below.maf
        )))
        
        # Filter the TP and candidate marker matrices for those markers
        TP.genos.use <- TP.genos.i[,-markers.to.remove]
        candidate.genos.use <- candidate.marker.genos.i[,-markers.to.remove]
        
        
        # Estimate marker effects
        marker.effects.solve <- mixed.solve(y = TP.phenos.i, Z = TP.genos.use, method = "REML")
        # Calculate GEBVs
        candidate.GEBV.i <- candidate.genos.use %*% marker.effects.solve$u
        

        
        
        ##### Step 5 - Phenotype the population
        # We use the haploid genotypes from the F1:3 generation to measure the
        ## the phenotypes
        
        # Measure the phenotype and true genotypic values of all selection candidates
        candidate.values.i <- phenotype.population( genome = hv.genome,
                                                    haploid.genos = candidate.haploid.i,
                                                    V_E = V_E,
                                                    V_e = V_e,
                                                    n.env = n.env,
                                                    n.rep = n.rep )
        
        # Validate the predictions
        # Find the correlation between the GEBVs and the true genotypic value
        pred.validation.i <- cor( candidate.GEBV.i, candidate.values.i$geno.values )


        
        ##### Step 6 - Select the parents of the next generation
        
        # Make selections on the GEBVs
        # Select the top 100 based on GEBVs for parents of the next cycle
        parent.selections.i <- select.population(pheno.mat = candidate.GEBV.i, 
                                                 sel.intensity = parents.sel.intensity, 
                                                 selection = "best")
        
        parent.lines.list <- list(p1 = parent.selections.i$lines.sel, p2 = parent.selections.i$lines.sel)
        
        # The parents are selected and crossed at the F3 stage, so subset the haploid genotpyes from the F1:3
        parent.haploids <- subset.gametes(gametes = candidate.haploid.i,
                                          line.names = parent.selections.i$lines.sel)
        
        parent.values <- subset.values(values.list = candidate.values.i, 
                                       lines.to.subset = parent.selections.i$lines.sel)
        
        
        ##### Step 7 - Update the TP
        
        # Skip this step if not called
        if (change != "no.change") {
          
          # If the TP change is best, worst, or random, simply subset the population.
          if (change %in% c("best", "worst", "random")) {
            
            TP.addition.list <- list(TP.addition.lines = select.population(pheno.mat = candidate.GEBV.i,
                                                                           sel.intensity = tp.update.increment,
                                                                           selection = change)$lines.sel )
          }

          if (change %in% c("PEVmean", "CDmean")) {
            # Analyze using PEVmean or CDmean
            # We want to see what optimized TP is best for the parents, so we will optimize the training set based
            ## on the lines from the whole candidate set, including the parents?
            phenotyped.lines <- row.names(candidate.marker.genos.i)
            unphenotyped.lines <- parent.selections.i$lines.sel
            
            # V_e is estimated from maximum liklihood
            V_e.i <- marker.effects.solve$Ve
            # V_a is estimated as the variance among marker effects * the number of markers
            V_a.i <- marker.effects.solve$Vu * ncol(TP.genos.use)
            
            # Subset the relationship matrix among candidates
            optimized.TP.additions <- TP.optimization.algorithms(A = A,
                                                                 phenotyped.lines = phenotyped.lines,
                                                                 unphenotyped.lines = unphenotyped.lines,
                                                                 n.TP = tp.update.increment,
                                                                 V_e = V_e.i,
                                                                 V_a = V_a.i,
                                                                 optimization.method = change,
                                                                 max.iter = 500,
                                                                 use.subset = T)
            
            # The optimized TP lines become the TP additions
            TP.addition.list <- list(TP.addition.lines = optimized.TP.additions$optimized.lines, optimization = optimized.TP.additions[-1])
            
          } # Close the tp optimization algorithm if statement
          
          ### Collect information on the new TP lines
          
          # TP additions
          TP.addition.lines <- TP.addition.list$TP.addition.lines
          # Subset the haploid genotypes for these lines
          TP.addition.haploids <- subset.gametes(gametes = candidate.haploid.i,
                                                line.names = TP.addition.lines)
          
          # Subset the geno matrix for these lines
          TP.addtion.genos <- candidate.marker.genos.i[TP.addition.lines,]
          
          # Gather genotypic and phenotypic values of the TP additions
          TP.addition.values <- subset.values(values.list = candidate.values.i, TP.addition.lines)
          # Separate the phenotypes
          TP.addtion.phenos <- TP.addition.values$mean.pheno.values
          
          
          # Combine the new data to the TP
          TP.phenos.i <- rbind(TP.phenos.i, TP.addtion.phenos)
          TP.genos.i <- rbind(TP.genos.i, TP.addtion.genos)
          TP.haploids.i <- rbind(TP.haploids.i, TP.addition.haploids)
          
          
          # If the TP formation calls for a sliding window, use only the ~750 most recent training individuals
          if (tp.formation == "window") {
            
            # If the breeding cycle * tp addition size is less than the starting tp size, 
            ## include the most recent 150 additions, then randomly sample from the
            ## remaining TP
            if ((breeding.cycle * tp.update.increment) < tp.size) {
              
              # Find the index of the most recent additions
              tp.recent.index <- tail(1:nrow(TP.genos.i), (breeding.cycle * tp.update.increment))
              # Randomly select among the remaining index to maintain the tp.size
              tp.random.index <- sort( sample(setdiff(1:nrow(TP.genos.i), tp.recent.index), size = (tp.size - length(tp.recent.index)) ) )
              # Combine
              tp.keep.index <- sort(c(tp.recent.index, tp.random.index))
              
            } else { # Otherwise just take the last tp.size individuals added to the TP
              tp.keep.index <- tail(1:nrow(TP.genos.i), tp.size)
              
            }
              # Set the TP.pheno and TP.genos
              TP.phenos.i <- as.matrix(TP.phenos.i[tp.keep.index,])
              TP.genos.i <- as.matrix(TP.genos.i[tp.keep.index,])
              TP.haploids.i <- subset.gametes(gametes = TP.haploids.i, line.names = row.names(TP.genos.i))
          }
          
        } else {
          
          TP.addition.list <- NA
          
        } # Close the tp.change if statement
                              
        
        print( paste("Cycle", breeding.cycle, "complete.") )
        
        cycle.name <- paste("cycle", breeding.cycle, sep = "")
        
        # Gather data for analysis
        simulation.results[[cycle.name]] <- list(geno.summary.stats = list(candidate.maf = candidate.genos.maf.i,
                                                                           TP.maf = TP.genos.maf.i,
                                                                           qtl.marker.LD = qtl.marker.LD.i,
                                                                           mu.TP.candidate.rel = mu.relationship),
                                                 marker.effects.solve = marker.effects.solve,
                                                 candidate.GEBV = candidate.GEBV.i,
                                                 candidate.values = candidate.values.i,
                                                 selection.values = parent.values,
                                                 prediction.accuracy = pred.validation.i,
                                                 tp.update = TP.addition.list )

                                   
        
      } # Close the per-cycle loop
      
      # Return the simulation data
      # return( list(sim.results = simulation.results, genome = hv.genome) )
      debug.results[[rep.i]] <- simulation.results
      
      } # Replicate for debugging

# Save everything
save.image(file = paste("simulation_debugging_", date, ".RData", sep = ""))
savehistory(file = paste("simulation_debugging_", date, ".Rhistory", sep = ""))
      
  #   }) # Close the iteration lapply
  #   
  # }, mc.cores = n.cores)
  
#   # Save the tp.change data
#   filename = paste("Files/", "simulation_results_q", n.QTL, "_sel", parents.sel.intensity, "_popmakeup-", pop.makeup, "_tpchange-", change, "_tpformation-", tp.formation, "_", date, ".RData", sep = "")
#   save(list = c("experiment.sub.results", "change", "metadata"), file = filename)
#   
#   
# } # Close the tp.change for loop
