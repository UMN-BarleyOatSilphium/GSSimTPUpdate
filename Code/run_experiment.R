## Genomic selection simulations
# This code runs iterations of long-term genomic selection simulations to 
## test the usefulness of updating a training population and how it impacts
## prediction accuracy, genetic variance, and selection potential

## Update: April 26, 2016
# This update will test whether including the parents is sufficient for updating the TP and maintaining predictive ability

# Are we using MSI?
MSI = T

# Arguments
args <- commandArgs(trailingOnly = T)

# First and only argument is the pop makeup (MN, ND, or MNxND)
# Second argument is how the TP should be combined after each cycle (cumulative or window)
if (all(is.na(args))) {
  pop.makeup <- "MN"
  tp.formation <- "cumulative"
  parents.sel.intensity = 100
  n.crosses = 50
} else {
  pop.makeup <- args[1]
  tp.formation <- args[2]
  parents.sel.intensity <- as.numeric(args[3])
  n.crosses = as.numeric(args[4])
}

# Load the packages
library(hypred, quietly = T)
library(rrBLUP, quietly = T)
library(boot, quietly = T)
library(parallel, quietly = T)
library(EMMREML, quietly = T)

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
source("Code/hypred_simulation_FUNCTIONS.R")


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
                 parents.sel.intensity = parents.sel.intensity,
                 n.env = n.env, 
                 n.rep = n.rep,
                 mutation.rate.qtl = mutation.rate.qtl,
                 mutation.rate.snp = mutation.rate.snp,
                 tp.update.increment = tp.update.increment,
                 tp.size = tp.size,
                 n.crosses = n.crosses,
                 ind.per.cross = ind.per.cross,
                 std.sel.intensity = std.sel.intensity,
                 n.iterations = n.iterations,
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
for (change in tp.change) {
  
  # Split iterations into cores
  if (n.cores > 1) {
    iters.per.core <- split(x = 1:n.iterations, factor(cut(x = 1:n.iterations, breaks = n.cores)))
    names(iters.per.core) <- paste("set", 1:length(iters.per.core), sep = "")
  } else {
    iters.per.core <- 1:n.iterations
  }
  
  # Apply the iterations over cores
  experiment.sub.results <- mclapply(X = iters.per.core, FUN = function(iter.set) {
    
    # Loop over each iteration
    set.results <- lapply(X = 1:length(iter.set), FUN = function(i) {
      
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
      
      # Find the MAF of all sites
      marker.maf <- sapply(X = colMeans(TP.genos+1)/2, FUN = function(freq) min(freq, 1-freq))
      # Determine those below the threshold
      snps.below.maf <- which(marker.maf < min.maf)
      
      
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
      TP.values <- evaluate.population2(genome = hv.genome,
                                        haploid.genos = TP.haploids.i,
                                        V_E = V_E,
                                        V_e = V_e,
                                        n.env = n.env,
                                        n.rep = n.rep)
      
      TP.phenos <- TP.values$mean.pheno.values
      
      # Next select the top 80 MN and top 80 ND
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
      parent.haploids <- TP.haploids.i
      
      # Set dummy variables for the phenos and genos
      TP.phenos.i <- TP.phenos
      TP.genos.i <- TP.genos
      
      # Create an initial data list
      simulation.results <- list()
      
      # Loop over the number of cycles
      for (breeding.cycle in 1:n.cycles) {
      # for (breeding.cycle in 1:5) {
        
        ##### Start the Cycle Executions #####

        ##### Step 1 - Crossing
        # Make a crossing block
        crossing.block.i <- make.crossing.block(parent1.lines = parent.lines.list$p1, 
                                                parent2.lines = parent.lines.list$p2, 
                                                n.crosses = n.crosses, 
                                                use.parents.once = T)
        
        ##### Step 2 - Make the crosses and inbreed to genotyping
        candidate.haploid.i <- make.population(genome = hv.genome, 
                                                           parental.haploids = parent.haploids,
                                                           crossing.block = crossing.block.i,
                                                           N = ind.per.cross,
                                                           cycle.number = breeding.cycle,
                                                           generations = 2,
                                                           pop.type = "inbred",
                                                           mutation.rate.snp = mutation.rate.snp,
                                                           mutation.rate.qtl = mutation.rate.qtl)
        
        # Find the genotypes of the markers and QTL
        candidate.i.genos <- genotype.loci(haploid.genos = candidate.haploid.i, 
                                           genome = hv.genome, 
                                           include.QTL = T)
        # Just the marker genotypes
        candidate.i.marker.genos <- genotype.loci(haploid.genos = candidate.haploid.i, 
                                                  genome = hv.genome, 
                                                  include.QTL = F)
        
      
        # Measure summary statistics
        candidate.i.genos.allele.freq <- calculate.allele.freq(geno.mat = candidate.i.genos)

        # Measure LD
        candidate.i.qtl.marker.LD <- measure.LD(genome = hv.genome, genos = candidate.i.genos, Morgan.window = 0.5)
        
        # Per QTL, find the LD of the marker with which the QTL has the highest LD
        candidate.i.mean.max.LD <- mean(unlist(
          sapply(X = candidate.i.qtl.marker.LD, FUN = function(chr) 
            sapply(X = chr, FUN = function(r) 
              ifelse(test = all(is.na(r)), yes = max((r^2)), no = max((r^2), na.rm = T))) )
          ), na.rm = T)
        
        candidate.i.mean.window.LD <- mean(unlist(
          sapply(X = candidate.i.qtl.marker.LD, FUN = function(chr) 
            sapply(X = chr, FUN = function(r) 
              mean(r^2, na.rm = T)))
        ), na.rm = T)
        
        # Measure LD on the TP
        TP.i.qtl.marker.LD <- measure.LD(genome = hv.genome, genos = TP.haploids.i, Morgan.window = 0.5)
        
        # Apply over chromosomes and combine the data
        TP.candidate.LD <- do.call("rbind", sapply(X = 1:length(TP.i.qtl.marker.LD), FUN = function(i) {
          # Find the common polymorphic QTL
          common.poly.qtl <- intersect(names(TP.i.qtl.marker.LD[[i]]), names(candidate.i.qtl.marker.LD[[i]]))
          # Apply over the common QTL and combine the data
          do.call("rbind", sapply(X = common.poly.qtl, FUN = function(poly.qtl) {
            # Find the common polymorphic markers
            common.poly.markers <- intersect(names(TP.i.qtl.marker.LD[[i]][[poly.qtl]]), names(candidate.i.qtl.marker.LD[[i]][[poly.qtl]]))
            # Return r for those markers
            cbind( TP.i.qtl.marker.LD[[i]][[poly.qtl]][common.poly.markers], candidate.i.qtl.marker.LD[[i]][[poly.qtl]][common.poly.markers])
          })) }))
        
        # Find the correlation of r across all correlations
        TP.candidate.persistance.of.phase <- cor(TP.candidate.LD, use = "complete.obs")[1,2]
        
        # Create a list to save
        qtl.marker.LD <- list(candidate.i.qtl.marker.LD = candidate.i.qtl.marker.LD,
                              TP.i.qtl.marker.LD = TP.i.qtl.marker.LD,
                              mean.max = candidate.i.mean.max.LD,
                              mean.window = candidate.i.mean.window.LD,
                              persistance.of.phase = TP.candidate.persistance.of.phase)
        
        
        
        ### #3 Prediction ###
        # Only use polymorphic markers and those with a MAF above the threshold
        # Determine monomorphic markers in the candidates
        mono.snps.candidates <- which(abs(colMeans(candidate.i.marker.genos)) == 1)
        # Determine monomorphic markers in the training population
        mono.snps.training <- which(abs(colMeans(TP.genos.i)) == 1)
        # Determine the common monomorphic markers and those below the maf freq
        markers.to.remove <- sort(unique(c(mono.snps.candidates, mono.snps.training, snps.below.maf)))
      
        
        # Filter the TP and candidate marker matrices for those markers
        TP.genos.use <- TP.genos.i[,-markers.to.remove]
        candidate.genos.use <- candidate.i.marker.genos[,-markers.to.remove]
        
        ### Calculate Relationship
        # Assign M
        M <- rbind(TP.genos.i, candidate.i.marker.genos)
        # Subtract P to make Z (need to convert P into a repeated matrix)
        Z = M - matrix(P, nrow(M), length(P), byrow = T)
        # Calculate the relationship matrix
        A = tcrossprod(Z) / c
        
        # Estimate marker effects
        marker.effects.solve <- mixed.solve(y = TP.phenos.i, Z = TP.genos.use, method = "REML")
        # Calculate GEBVs
        candidate.i.GEBV <- candidate.genos.use %*% marker.effects.solve$u
        
        # Measure the phenotype and true genotypic values of all selection candidates
        candidate.i.values <- evaluate.population2( genome = hv.genome,
                                                    haploid.genos = candidate.haploid.i,
                                                    V_E = V_E,
                                                    V_e = V_e,
                                                    n.env = n.env,
                                                    n.rep = n.rep )
        
        # Validate the predictions
        # Find the correlation between the GEBVs and the true genotypic value
        pred.validation.i <- cor( candidate.i.GEBV, candidate.i.values$geno.values )

        # Calculate the mean relationship between the TP and the candidates
        mu.relationship <- mean(rowMeans(A[row.names(TP.genos.use), row.names(candidate.genos.use)]))
        
        ##### Step 4 - Select the Next Parents #####
        
        # Make selections on the GEBVs
        # Select the top 100 based on GEBVs for parents of the next cycle
        parent.selections.i <- select.population(pheno.mat = candidate.i.GEBV, 
                                                 sel.intensity = parents.sel.intensity, 
                                                 selection = "best")
        
        parent.lines <- parent.selections.i$lines.sel
        parent.lines.list <- list(p1 = parent.lines, p2 = parent.lines)
        # The parents are selected and crossed at the F3 stage, so subset the gametes from the F3
        parent.haploids <- subset.gametes(gametes = candidate.haploid.i,
                                         line.names = parent.lines)
        
        parent.values <- subset.values(values.list = candidate.i.values, lines.to.subset = parent.lines)
        
        
        ##### Step 5 - Update the TP #####
        # Skip this step if not called
        if (change != "no.change") {
          
          # If the TP change is best, worst, or random
          if (change %in% c("best", "worst", "random")) {
            
            TP.addition.list <- list(TP.addition.lines = select.population(pheno.mat = candidate.i.GEBV,
                                                         sel.intensity = tp.update.increment,
                                                         selection = change)$lines.sel )
          }

          if (change %in% c("PEVmean", "CDmean")) {
            # Analyze using PEVmean or CDmean
            # We want to see what optimized TP is best for the parents, so we will optimize the training set based
            ## on the lines from the whole candidate set, including the parents?
            phenotyped.lines <- row.names(candidate.i.marker.genos)
            unphenotyped.lines <- parent.lines
            
            # V_e is estimated from maximum liklihood
            V_e.i <- marker.effects.solve$Ve
            # V_a is estimated as the variance among marker effects * the number of markers
            V_a.i <- marker.effects.solve$Vu * ncol(TP.genos.use)
            
            # Subset the relationship matrix among candidates
            optimized.TP.additions <- try(TP.optimization.algorithms(A = A,
                                                                     phenotyped.lines = phenotyped.lines,
                                                                     unphenotyped.lines = unphenotyped.lines,
                                                                     n.TP = tp.update.increment,
                                                                     V_e = V_e.i,
                                                                     V_a = V_a.i,
                                                                     optimization.method = change,
                                                                     max.iter = 500,
                                                                     use.subset = T))
            
            # If an error is found, just try again
            # Use a counter to make sure this doesn't implode
            counter = 0
            while (all(class(optimized.TP.additions) == "try-error", counter < 100) ) {
              optimized.TP.additions <- try(TP.optimization.algorithms(A = A,
                                                                       phenotyped.lines = phenotyped.lines,
                                                                       unphenotyped.lines = unphenotyped.lines,
                                                                       n.TP = tp.update.increment,
                                                                       V_e = V_e.i,
                                                                       V_a = V_a.i,
                                                                       optimization.method = change,
                                                                       max.iter = 500,
                                                                       use.subset = T))
              counter = counter + 1
            }
            
            # The optimized TP lines become the TP additions
            TP.addition.list <- list(TP.addition.lines = optimized.TP.additions$optimized.lines, optimization = optimized.TP.additions[-1])
            
          } # Close the tp optimization algorithm if statement
          
          # TP additions
          TP.addition.lines <- TP.addition.list$TP.addition.lines
          # Subset the gametes for these lines
          TP.addition.haploids <- subset.gametes(gametes = candidate.haploid.i,
                                                line.names = TP.addition.lines)
          
          # Subset the geno matrix for these lines
          TP.addtion.genos <- candidate.i.marker.genos[TP.addition.lines,]
          
          # Gather genotypic and phenotypic values of the TP additions
          TP.addition.values <- subset.values(values.list = candidate.i.values, TP.addition.lines)
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
        simulation.results[[cycle.name]] <- list(geno.summary.stats = list(allele.freq = candidate.i.genos.allele.freq,
                                                                           qtl.marker.LD = qtl.marker.LD,
                                                                           mu.TP.candidate.rel = mu.relationship),
                                                 marker.effects.solve = marker.effects.solve,
                                                 candidate.values = candidate.i.values,
                                                 selection.values = parent.values,
                                                 prediction.accuracy = pred.validation.i,
                                                 tp.update = TP.addition.list )

                                   
        
      } # Close the per-cycle loop
      
      # Return the simulation data
      return( list(sim.results = simulation.results, genome = hv.genome) )
      
    }) # Close the iteration lapply
  
    # Return the set data
    return(set.results)
    
  }, mc.cores = n.cores)
  
  # Save the tp.change data
  filename = paste("Files/", "simulation_results_q", n.QTL, "_sel", parents.sel.intensity, "_popmakeup-", pop.makeup, "_tpchange-", change, "_tpformation-", tp.formation, "_", date, ".RData", sep = "")
  save(list = c("experiment.sub.results", "change", "metadata"), file = filename)
  
  
} # Close the tp.change for loop
