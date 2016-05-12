## Genomic selection simulations
# This code runs iterations of long-term genomic selection simulations to 
## test the usefulness of updating a training population and how it impacts
## prediction accuracy, genetic variance, and selection potential

## Update: April 26, 2016
# This update will test whether including the parents is sufficient for updating the TP and maintaining predictive ability

# Arguments
args <- commandArgs(trailingOnly = T)

# First and only argument is the pop makeup (MN, ND, or MNxND)
# Second argument is how the TP should be combined after each cycle (cumulative or window)
if (all(is.na(args))) {
  pop.makeup <- "MN"
  tp.formation <- "cumulative"
} else {
  pop.makeup <- args[1]
  tp.formation <- args[2]
}

# Load the packages
library(hypred)
library(rrBLUP)
library(boot)
library(parallel)

# Load already-curated gamete data
load(file = "Files/Barley_CAP_simuation_starting_material.RData")

# Are we using MSI?
MSI = FALSE

# # Other tools
if (MSI) {
  setwd("/panfs/roc/groups/6/smithkp/neyhartj/Genomic_Selection/Simulations/Barley_GS_PEV/Barley_GS_Simulations")
  # source("/panfs/roc/groups/6/smithkp/neyhartj/GitHub_Repos/Quant-Gen-Scripts/genotype_matrix_utils.R")
  n.cores = 16
} else {
  setwd("C:/Users/Jeff/Google Drive/Barley Lab/Projects/Side Projects/Simulations/Barley_GS_Simulations/")
  n.cores = 1
}

# Source functions
source("Code/hypred_simulation_FUNCTIONS.R")


# Other simulation parameters
# Entry-mean heritability
h2 = 0.5
# How many cycles?
n.cycles = 15
# Number of QTL underlying trait
n.QTL = 100
# Selection intensity
GEBV.sel.intensity = 0.1
# Number of phenotyping environments and reps
n.env = 5
n.rep = 1

# Barley population genetics data
mutation.rate.snp = 7e-8
mutation.rate.qtl = 2.5e-5

# tp.change = c("best", "worst", "random", "parents", "no.change", "PEVmean", "CDmean")
tp.change = "no.change"
# The number of lines to add the TP after each cycle
tp.update.increment = 150
# Size of the TP to maintain - this is the same as the starting TP
tp.size <- nrow(CAP.gametes.0.03) / 2

# Parent selection and crossing parameters
n.crosses = 40
ind.per.cross = 30

# Computation parameters
n.iterations = 100

date <- format(Sys.time(), "%d%m%y-%H%M%S")

# Create the final experiment output list
experiment.results <- list()

#### Define genome characteristics ####

### 05.12.16 Filtering markers based on minor allele prior to running the simulation
# The idea here is that minor alleles may have a inaccurate marker effect estimates,
# so if we prime the simulation by filtering out the markers with lower MAF,
# we might get more accurate marker estimates earlier on

maf.filt.results.list <- list()

for (min.maf in gsub(pattern = "CAP.gametes.", apropos("CAP.gametes"), replacement = "")) {
  
  assign(value = get(apropos(paste("CAP.gametes", min.maf, sep = "."))), x = "CAP.gametes")
  assign(value = get(apropos(paste("sampled.markers", min.maf, sep = "."))), x = "sampled.markers")
  

  # Find the snps per chromsome
  n.chr.snps = nrow(sampled.markers) / length(unique(sampled.markers$chrom))
  
  # Find the chromsome lengths
  chr.len <- as.numeric(tapply(X = sampled.markers$pos, INDEX = sampled.markers$chrom, FUN = max))
  
  # Make the initial genome
  hv.genome <- make.genome( n.chr = 7, 
                            chr.len = chr.len, 
                            n.chr.snps = n.chr.snps,
                            genetic.map = sampled.markers$pos)
  
  
  
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
                                        qtl.ids = NULL, 
                                        qtl.dom.ids = NULL, 
                                        qtl.per.ids = NULL, 
                                        qtl.add.eff = "geometric", 
                                        qtl.dom.eff = NULL,
                                        keep.all.snps = FALSE)
        
        # Since more loci have been added to the genome, we need to add columns to the gamete matricies to make up for this
        # TP.gametes <- augment.gamete.mat(genome = hv.genome, 
        #                                  gamete.mat = CAP.gametes)
        
        TP.gametes <- CAP.gametes
        # Convert the gametes to genotypes
        TP.genos <- genotype.markers(gametes = TP.gametes, genome = hv.genome, DH = F)
        
        ## Select the TP lines for use as parents
        # First separate MN and ND lines
        line.names <- row.names(TP.genos)
        ND.lines <- grep(pattern = "^ND", x = line.names, value = T)
        MN.lines <- setdiff(x = line.names, y = c(ND.lines))
        
        # Phenotype the training population
        TP.values <- measure.values(genome = hv.genome, gametes = TP.gametes, h2 = h2, n.env = n.env, n.rep = n.rep, model = NULL)
        TP.phenos <- TP.values$mean.pheno.values
        
        # Next select the top 80 MN and top 80 ND
        top.MN.lines <- names(sort(TP.phenos[MN.lines,], decreasing = T)[1:80])
        top.ND.lines <- names(sort(TP.phenos[ND.lines,], decreasing = T)[1:80])
        
        # Create a list to store the line names
        pop.makeup.list <- list(
          MN = list(p1 = top.MN.lines, p2 = top.MN.lines),
          ND = list(p1 = top.ND.lines, p2 = top.ND.lines),
          MNxND = list(p1 = top.MN.lines[1:40], p2 = top.ND.lines[1:40])
        )
        
        # Set the parent gamete data input
        parent.lines.list <- pop.makeup.list[[pop.makeup]]
        parent.gametes <- TP.gametes
        
        # Set dummy variables for the phenos and genos
        TP.phenos.i <- TP.phenos
        TP.genos.i <- TP.genos
  
        
        # Create an initial data list
        simulation.results <- list()
        # Loop over the number of cycles
        for (breeding.cycle in 1:n.cycles) {
          
          ##### Start the Cycle Executions #####
  
          ##### Step 1 - Crossing
          # Make a crossing block
          crossing.block.i <- make.crossing.block(parent1.lines = parent.lines.list$p1, 
                                                  parent2.lines = parent.lines.list$p2, 
                                                  n.crosses = n.crosses, 
                                                  use.parents.once = T)
          
          ##### Step 2 - Make the crosses and inbreed to genotyping
          system.time(candidate.gametes.i <- make.population(genome = hv.genome, 
                                                             named.parental.gametes = parent.gametes,
                                                             crossing.block = crossing.block.i,
                                                             N = ind.per.cross,
                                                             cycle.number = breeding.cycle,
                                                             generations = 2,
                                                             pop.type = "inbred",
                                                             mutation.rate.snp = mutation.rate.snp,
                                                             mutation.rate.qtl = mutation.rate.qtl))
          
          # Genotype the population
          candidate.i.genos <- genotype.markers(gametes = candidate.gametes.i, 
                                                genome = hv.genome, 
                                                DH = F)
        
          # Measure summary statistics
          candidate.i.genos.allele.freq <- calculate.allele.freq(geno.mat = candidate.i.genos)
          candidate.i.genos.pairwise.div <- SNP.pairwise.div(geno.mat = candidate.i.genos)
          # Measure heterozygosity of the entries
          candidate.i.genos.het <- apply(X = candidate.i.genos, MARGIN = 1, FUN = function(geno) sum(geno == 0) / length(geno))
        
          # Remove monomorphic markers before making predictions
          # Determine monomorphic markers
          poly.snps <- apply(X = candidate.i.genos, MARGIN = 2, FUN = function(snp) length(unique(snp)) > 1)
          # Filter the TP and candidate marker matrices for those markers
          TP.genos.use <- TP.genos.i[,poly.snps]
          candidate.genos.use <- candidate.i.genos[,poly.snps]
          
          ##### Step 3 - Genomic prediction
          system.time(candidate.i.prediction <- make.predictions(pheno.train = TP.phenos.i, 
                                                                 geno.train = TP.genos.use, 
                                                                 geno.pred = candidate.genos.use, 
                                                                 model = "RRBLUP"))
          
          # Retrieve GEBVs
          candidate.i.GEBV <- candidate.i.prediction$GEBV
          
          # Measure the phenotype and true genotypic values of all selection candidates
          candidate.i.values <- measure.values( genome = hv.genome,
                                                gametes = candidate.gametes.i,
                                                h2 = h2,
                                                n.env = n.env,
                                                n.rep = n.rep,
                                                V_e.scale = 8 )
          
          # Validate the predictions
          # Find the correlation between the GEBVs and the true genotypic value
          pred.validation.i <- validate.predictions(predicted.GEBVs = candidate.i.GEBV,
                                                    observed.values = candidate.i.values$geno.values,
                                                    boot.reps = 1000)
          
          ##### Step 4 - Select the Next Parents #####
          
          # Make selections on the GEBVs
          # Select the top 100 based on GEBVs for parents of the next cycle
          parent.selections.i <- select.population(pheno.mat = candidate.i.GEBV, 
                                                   sel.intensity = 100, 
                                                   selection = "best")
          
          parent.lines <- parent.selections.i$lines.sel
          parent.lines.list <- list(p1 = parent.lines, p2 = parent.lines)
          # The parents are selected and crossed at the F3 stage, so subset the gametes from the F3
          parent.gametes <- subset.gametes(gametes = candidate.gametes.i,
                                           line.names = parent.lines)
          
          parent.values <- subset.values(values.list = candidate.i.values, lines.to.subset = parent.lines)
          
          
          ##### Step 5 - Update the TP #####
          # Skip this step if not called
          if (change != "no.change") {
            
            # If the change is just the parents, just add the parent candidates
            if (change == "parents") {
              
              TP.addition.list <- list(TP.addition.lines = parent.selections.i$lines.sel)
            }
            
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
              phenotyped.lines <- row.names(candidate.i.genos)
              unphenotyped.lines <- parent.lines
              
              # The V_e and V_a will be taken from the REML estimates of the previous mixed model
              V_e.i <- candidate.i.prediction$solve.out$Ve
              V_a.i <- candidate.i.prediction$solve.out$Vu * ncol(candidate.i.genos) # genetic variance is V_u * n.markers
              
              optimized.TP.additions <- try(TP.optimization(genos = candidate.genos.use,
                                                            phenotyped.lines = phenotyped.lines,
                                                            unphenotyped.lines = unphenotyped.lines,
                                                            n.TP = 150,
                                                            V_e = V_e.i,
                                                            V_a = V_a.i,
                                                            optimization = change,
                                                            max.iter = 500,
                                                            use.subset = T))
              
              # If an error is found, just try again
              # Use a counter to make sure this doesn't implode
              counter = 0
              while (all(class(optimized.TP.additions) == "try-error", counter < 50) ) {
                optimized.TP.additions <- try(TP.optimization(genos = candidate.genos.use,
                                                              phenotyped.lines = phenotyped.lines,
                                                              unphenotyped.lines = unphenotyped.lines,
                                                              n.TP = 150,
                                                              V_e = V_e.i,
                                                              V_a = V_a.i,
                                                              optimization = change,
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
            TP.addition.gametes <- subset.gametes(gametes = candidate.gametes.i,
                                                  line.names = TP.addition.lines)
            
            # Subset the geno matrix for these lines
            TP.addtion.genos <- candidate.i.genos[TP.addition.lines,]
            
            # Gather genotypic and phenotypic values of the TP additions
            TP.addition.values <- subset.values(values.list = candidate.i.values, TP.addition.lines)
            # Separate the phenotypes
            TP.addtion.phenos <- TP.addition.values$mean.pheno.values
            
            
            # Combine the new data to the TP
            TP.phenos.i <- rbind(TP.phenos.i, TP.addtion.phenos)
            TP.genos.i <- rbind(TP.genos.i, TP.addtion.genos)
            
            
            # If the TP formation calls for a sliding window, use only the ~750 most recent training individuals
            if (tp.formation == "window") {
              
              # If the breeding cycle * tp addition size is less than 750, select the last 750, then randomly sample the TP
              if ((breeding.cycle * tp.update.increment) < tp.size) {
                
                tp.keep.index <- tail(1:nrow(TP.genos.i), (breeding.cycle * tp.update.increment))
                tp.random.index <- sort(sample(setdiff(1:nrow(TP.genos.i), tp.keep.index), (tp.size - length(tp.keep.index))))
                tp.all.index <- sort(c(tp.keep.index, tp.random.index))
                
                # Set the TP.pheno and TP.genos
                TP.phenos.i <- as.matrix(TP.phenos.i[tp.all.index,])
                TP.genos.i <- as.matrix(TP.genos.i[tp.all.index,])
                
              } else { # Otherwise just take the tail
                
                tp.keep.index <- tail(1:nrow(TP.genos.i), (breeding.cycle * tp.update.increment))
                
                TP.phenos.i <- as.matrix(TP.phenos.i[tp.keep.index,])
                TP.genos.i <- as.matrix(TP.genos.i[tp.keep.index,])
              }
            }
            
          } else {
            TP.addition.list <- NA
          } # Close the tp.change if statement
                                
          
          print( paste("Cycle", breeding.cycle, "complete.") )
          
          cycle.name <- paste("cycle", breeding.cycle, sep = "")
          
          # Gather data for analysis
          simulation.results[[cycle.name]] <- list(geno.summary.stats = list(pairwise.div = candidate.i.genos.pairwise.div,
                                                                             allele.freq = candidate.i.genos.allele.freq,
                                                                             heterozygosity = candidate.i.genos.het),
                                                   prediction.results = list(GEBV = candidate.i.prediction$GEBV,
                                                                             parameters = candidate.i.prediction$parameters),
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
    filename = paste("Files/", "simulation_results_q", n.QTL, "_sel", GEBV.sel.intensity, "_popmakeup-", pop.makeup, "_tpchange-", change, "_tpformation-", tp.formation, "_", date, ".RData", sep = "")
    save(list = c("experiment.sub.results", "change"), file = filename)
    
    
  } # Close the tp.change for loop
  
  maf.filt.results.list[[as.character(min.maf)]] <- experiment.results
  
} # Close the maf.filter loop