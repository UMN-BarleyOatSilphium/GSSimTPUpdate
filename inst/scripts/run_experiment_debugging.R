## Genomic selection simulations
# This code runs iterations of long-term genomic selection simulations to 
## test the usefulness of updating a training population and how it impacts
## prediction accuracy, genetic variance, and selection potential

# Arguments
args <- commandArgs(trailingOnly = T)

# First and only argument is the pop makeup (MN, ND, or MNxND)
# Second argument is how the TP should be combined after each cycle (cumulative or window)
# If no arguments are given (i.e. it's being run locally), set defaults for testing
if (all(is.na(args))) {
  pop.makeup <- "MNxND"
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

# Other information and pre-processing
if (MSI) {
  # Directory to save the files
  save.dir <- "/panfs/roc/groups/6/smithkp/neyhartj/Genomic_Selection/Simulations/GSsim.TPUpdate/inst/output/"
  
  # Set the directory of the R packages
  package.dir <- "/panfs/roc/groups/6/smithkp/neyhartj/R/x86_64-pc-linux-gnu-library/3.3/"
  
  # Load the packages
  library(GSsim.TPUpdate, quietly = T, lib.loc = package.dir)
  library(parallel, quietly = T, package.dir)
  library(stringr, quietly = T, package.dir)
  
  # If we are not running MSI
} else {
  
  save.dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Side Projects/Simulations/GSsim.TPUpdate/inst/output/"
  
  # Load the packages
  library(GSsim.TPUpdate)
  library(parallel)
  library(stringr)
  
}

# Set the number of cores by detection
n.cores <- detectCores()

# Other simulation parameters
# Entry-mean heritability in the base population
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

# Barley population genetics data (citation!)
mutation.rate.snp = 7e-8
mutation.rate.qtl = 7e-8

tp.change = c("best", "worst", "random", "nochange", "PEVmean", "CDmean")

# The number of lines to add the TP after each cycle
tp.update.increment = 150
# Size of the TP to maintain - this is the same as the starting TP
tp.size <- nrow(CAP.haploids) / 2

# Parent selection and crossing parameters
ind.per.cross = 20
cycle.candidate.size = n.crosses * ind.per.cross

# Standardized selection intensity
std.sel.intensity = parents.sel.intensity / cycle.candidate.size

# Computation parameters
n.iterations = 500

date <- format(Sys.time(), "%d%m%y-%H%M%S")

# Save the metadata to a list
metadata <- list(h2 = h2,
                 n.cycles = n.cycles,
                 n.QTL = n.QTL,
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
for (change in tp.change) {
  
  #   # Split iterations into cores
  # if (n.cores > 1) {
  #   iters.per.core <- split(x = 1:n.iterations, factor(cut(x = 1:n.iterations, breaks = n.cores)))
  #   names(iters.per.core) <- paste("set", seq_along(iters.per.core), sep = "")
  # } else {
  #   iters.per.core <- 1:n.iterations
  # }    

  # Apply the iterations over cores
  # experiment.sub.results <- mclapply(X = iters.per.core,  FUN = function(iter.set) {
    
  for (iter in n.iterations) {

    # # Create an empty list to store repetition results
    # rep.results <- list()
    # 
    # # Loop over each iteration
    # for (r in seq_along(iter.set)) {
      
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
      marker.maf <- measure.maf(TP.genos)
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
      ND.lines <- str_subset(string = line.names, pattern = "^ND")
      MN.lines <- setdiff(x = line.names, y = ND.lines)
      
      ## Set the inital variances for the heritability
      # True genetic variance
      TP.V_g <- genotypic.value(genome = hv.genome, haploid.genos = TP.haploids.i) %>%
        var()
      # Environmental variance (scale * 8 as in Bernardo 2015)
      V_E = TP.V_g * 8
      # Residual variance scaled to achieve desired h2
      V_e = n.rep * n.env * ((TP.V_g / h2) - TP.V_g)
      
      
      # Phenotype the training population
      TP.values <- phenotype.population(genome = hv.genome,
                                        haploid.genos = TP.haploids.i,
                                        V_E = V_E,
                                        V_e = V_e,
                                        n.env = n.env,
                                        n.rep = n.rep,
                                        run.anova = T)
      
      TP.phenos <- TP.values$mean.pheno.values
      
      # Next select the top MN and top ND
      top.MN.lines <- names(sort(TP.phenos[MN.lines,], decreasing = T)[1:parents.sel.intensity])
      top.ND.lines <- names(sort(TP.phenos[ND.lines,], decreasing = T)[1:parents.sel.intensity])
      
      # Create a list to store the line names
      pop.makeup.list <- list(
        MN = list(p1 = top.MN.lines, p2 = top.MN.lines),
        ND = list(p1 = top.ND.lines, p2 = top.ND.lines),
        MNxND = list(p1 = top.MN.lines[seq(parents.sel.intensity / 2)], p2 = top.ND.lines[seq(parents.sel.intensity / 2)])
      )
      
      # Set the parent gamete data input
      parent.lines.list <- pop.makeup.list[[pop.makeup]]
      parent.haploids <- select.haploids(haploid.genos = TP.haploids.i, line.names = unique(unlist(parent.lines.list)))
      
      # Set dummy variables for the phenos and genos
      TP.phenos.i <- TP.phenos
      TP.genos.i <- TP.genos
      
      # Create an initial data list
      simulation.results <- list()
      
      # Loop over the number of cycles
      for (breeding.cycle in seq(n.cycles)) {

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
        candidate.genos.maf.i <- measure.maf(genos = candidate.genos.i)
        # Measure the maf of just the markers
        candidate.marker.maf.i <- measure.maf(genos = candidate.marker.genos.i)
        TP.genos.maf.i <- measure.maf(genos = TP.genos.i)
        
        ### LD measures
        # Candidates
        candidate.LD.window <- measure.LD(genome = hv.genome, 
                                          genos = candidate.genos.i, 
                                          Morgan.window = 0.25)
        
        # No window
        candidate.LD.genome <- measure.LD(genome = hv.genome,
                                          genos = candidate.genos.i)
        
        # Per QTL, find the LD of the marker (within the window) with which the 
        ## QTL has the highest LD then find the mean of those values across 
        ## polymorphic qtl
        candidate.mean.max.LD.window <- lapply(X = candidate.LD.window, FUN = function(chr)
          lapply(X = chr, FUN = function(qtl) max(qtl^2)) ) %>%
          unlist() %>%
          mean(na.rm = T)
        
        # Now just find the mean of LD across the whole window
        candidate.mean.LD.window <- candidate.LD.window %>% 
          unlist() %>% 
          .^2 %>% 
          mean(na.rm = T)
        
        # For the whole genome, find the mean LD value across those
        ## max LD values per QTL
        candidate.mean.max.LD.genome <- apply(X = candidate.LD.genome, MARGIN = 1, FUN = function(qtl)
          max(qtl^2) ) %>%
          mean(na.rm = T)
        
        # For the whole, genome, find the mean LD across all QTL-marker pairs
        candidate.mean.LD.genome <- mean(candidate.LD.genome ^ 2, na.rm = T)
        
        # Measure genomic LD on the TP
        TP.LD.genome <- measure.LD(genome = hv.genome, 
                                   genos = genotype.loci(haploid.genos = TP.haploids.i,
                                                         genome = hv.genome,
                                                         include.QTL = T) )
        
        ## Persistance of LD phase
        # First find the common polymorphic QTL
        common.poly.QTL <- intersect( row.names(TP.LD.genome), row.names(candidate.LD.genome) )
        common.poly.markers <- intersect( colnames(TP.LD.genome), colnames(candidate.LD.genome) )
        
        # Subset the TP and candidates for those markers and QTL, then vectorize
        TP.LD.vector <- TP.LD.genome[common.poly.QTL, common.poly.markers] %>%
          as.vector()
        candidate.LD.vector <- candidate.LD.genome[common.poly.QTL, common.poly.markers] %>%
          as.vector()
        
        # Correlate
        TP.candidate.persistance.of.phase <- cor(TP.LD.vector, candidate.LD.vector)
          
        # Create a list to save
        qtl.marker.LD.i <- list(mean.max.window = candidate.mean.max.LD.window,
                                mean.window = candidate.mean.LD.window,
                                mean.max.genome = candidate.mean.max.LD.genome,
                                mean.genome = candidate.mean.LD.genome,
                                persistance.of.phase = TP.candidate.persistance.of.phase)
        
        ### Measure the average relationship between the TP and the candidates
        # Assign M
        M <- rbind(TP.genos.i, candidate.marker.genos.i)
        # Subtract P to make Z (need to convert P into a repeated matrix)
        W = M - matrix(P, nrow(M), length(P), byrow = T)
        # Calculate the relationship matrix
        A = tcrossprod(W) / c
        
        # Subset the relationship matrix for the selection candidates
        A.sc <- A[row.names(candidate.marker.genos.i), row.names(candidate.marker.genos.i)]
        
        # Calculate the mean relationship between the TP and the candidates
        mu.relationship <- A[row.names(TP.genos.i), row.names(candidate.marker.genos.i)] %>%
          mean()
        
        ## Find the average inbreeding coefficient among the selection candidates
        candidate.inbreeding <- A.sc %>%
          diag() %>%
          - 1 %>%
          mean()
        
        
        ##### Step 4 - Prediction
        # Remove the markers with maf below the threshold (set at the start of the sim)
        ## also remove monomorphic markers
        markers.to.remove <- c( which(candidate.marker.maf.i == 0), 
                                which(TP.genos.maf.i == 0),
                                markers.below.maf ) %>%
          unique() %>%
          sort()

        
        # Filter the TP and candidate marker matrices for those markers
        TP.genos.use <- TP.genos.i[,-markers.to.remove]
        candidate.genos.use <- candidate.marker.genos.i[,-markers.to.remove]
        
        
        # Estimate marker effects
        predictions.out <- make.predictions(pheno.train = TP.phenos.i,
                                            geno.train = TP.genos.use,
                                            geno.pred = candidate.genos.use)
        
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
        pred.validation.i <- validate.predictions(predicted.values = predictions.out$GEBV,
                                                  observed.values = candidate.values.i$geno.values)


        
        ##### Step 6 - Select the parents of the next generation
        
        # Make selections on the GEBVs
        # Select the top 100 based on GEBVs for parents of the next cycle
        parent.selections.i <- select.population(value.mat = predictions.out$GEBV, 
                                                 sel.intensity = parents.sel.intensity, 
                                                 selection = "best")
        
        parent.lines.list <- list(p1 = parent.selections.i$lines.sel, p2 = parent.selections.i$lines.sel)
        
        # The parents are selected and crossed at the F3 stage, so subset the haploid genotpyes from the F1:3
        parent.haploids <- select.haploids(haploid.genos = candidate.haploid.i,
                                          line.names = parent.selections.i$lines.sel)
        
        parent.values <- select.values(pheno.values.list = candidate.values.i, 
                                       line.names = parent.selections.i$lines.sel)
        
        # Calculate the mean relationship among the parents
        parents.mu.relationship <- A[parent.lines.list$p1, parent.lines.list$p1] %>%
          .[upper.tri(.)] %>% 
          mean()
        
        # Calculat the mean inbreeding coefficient among the parents
        parents.inbreeding <- A[parent.lines.list$p1, parent.lines.list$p1] %>%
          diag() %>%
          - 1 %>%
          mean()
        
        
        ##### Step 7 - Update the TP
        
        # Skip this step if not called
        if (change != "nochange") {
          
          # If the TP change is best, worst, or random, simply subset the population.
          if (change %in% c("best", "worst", "random")) {
            
            TP.addition.list <- 
              list(TP.addition.lines = select.population(value.mat = predictions.out$GEBV,
                                                         sel.intensity = tp.update.increment,
                                                         selection = change)$lines.sel )
          }

          if (change %in% c("PEVmean", "CDmean")) {
            
            # Analyze using PEVmean or CDmean
            # We want to see what optimized TP is best for the parents, so we will optimize the training set based
            ## on the lines from the whole candidate set, including the parents
            phenotyped.index <- which( row.names(candidate.marker.genos.i) %in%
                                         row.names(A.sc) )
            unphenotyped.index <- which( parent.selections.i$lines.sel %in%
                                           row.names(A.sc) )

            # V_e is estimated from maximum likelihood
            V_e.i <- predictions.out$solve.out$Ve
            # V_a is estimated as the variance among marker effects * the number of markers
            V_a.i <- predictions.out$solve.out$Vu * ncol(TP.genos.use)
            
            # Run the optimization algorithm
            optimized.TP <- optimize.subset(method = change, A = A.sc,
                                            n.training = tp.update.increment,
                                            phenotyped.index = phenotyped.index,
                                            unphenotyped.index = unphenotyped.index,
                                            V_a = V_a.i, V_e = V_e.i, max.iter = 200)
            
            # Using the optimized index of "phenotyped" entries, determine
            # which candidates should be added to the TP
            TP.addition.list <- list(TP.addition.lines = row.names(A.sc)[optimized.TP$OTP],
                                     optimization.info = optimized.TP)
            
          } # Close the tp optimization algorithm if statement
          
          ### Collect information on the new TP lines
          
          # TP additions
          TP.addition.lines <- TP.addition.list$TP.addition.lines
          
          TP.addition.haploids <- select.haploids(haploid.genos = candidate.haploid.i,
                                                  line.names = TP.addition.lines)
          
          # Subset the geno matrix for these lines
          TP.addition.genos <- candidate.marker.genos.i[TP.addition.lines,]
          TP.addition.genos.qtl <- candidate.genos.i[TP.addition.lines,]
          
          # Gather genotypic and phenotypic values of the TP additions
          TP.addition.values <- select.values(pheno.values.list = candidate.values.i, 
                                              line.names = TP.addition.lines)
          # Separate the phenotypes
          TP.addition.phenos <- TP.addition.values$mean.pheno.values
          
          # Measure the average relationship among these lines
          TP.addition.mu.relationship <- A[TP.addition.lines, TP.addition.lines] %>%
            .[upper.tri(.)] %>% 
            mean()
          
          # Measure inbreeding among the additions
          TP.addition.inbreeding <- A[TP.addition.lines, TP.addition.lines] %>%
            diag() %>% 
            - 1 %>%
            mean()
          
          # Combine the new data to the TP
          TP.phenos.i <- rbind(TP.phenos.i, TP.addition.phenos)
          TP.genos.i <- rbind(TP.genos.i, TP.addition.genos)
          TP.haploids.i <- rbind(TP.haploids.i, TP.addition.haploids)
          
          ## Measure the expected heterozygosity of the additions, using all
          ## loci including QTL
          TP.addition.list[["Exp.het"]] <- 
            measure.expected.het(genos = TP.addition.genos.qtl)
          
          
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
              TP.haploids.i <- select.haploids(haploid.genos = TP.haploids.i,
                                               line.names = row.names(TP.genos.i))
          }
          
        } else {
          
          TP.addition.list <- list(TP.addition.lines = NA,
                                   Exp.het = NA)
          
          TP.addition.inbreeding <- NA
          TP.addition.mu.relationship <- NA
          
          
        } # Close the tp.change if statement
                              
        
        print( paste("Cycle", breeding.cycle, "complete.") )
        
        cycle.name <- paste("cycle", breeding.cycle, sep = "")
        
        # Gather data for analysis
        simulation.results[[cycle.name]] <- 
          list(geno.summary.stats = list(candidate.maf = candidate.genos.maf.i,
                                         TP.maf = TP.genos.maf.i,
                                         qtl.marker.LD = qtl.marker.LD.i,
                                         mu.TP.candidate.rel = mu.relationship),
               MM.solve = predictions.out$solve.out,
               candidate.GEBV = predictions.out$GEBV,
               candidate.values = candidate.values.i,
               selection.values = parent.values,
               prediction.accuracy = pred.validation.i,
               parents <- parent.lines.list,
               tp.update = TP.addition.list,
               inbreeding = list(candidates = candidate.inbreeding,
                                 TP.additions = TP.addition.inbreeding,
                                 parents = parents.inbreeding),
               relationship = list(TP.candidates = mu.relationship,
                                   TP.additions = TP.addition.mu.relationship,
                                   parents = parents.mu.relationship))

                                   
        
      } # Close the per-cycle loop
      
      # Rep number
      # rep.name <- paste("rep", r, sep = "")
      
      # Add the simulation results to the set results
      # rep.results[[rep.name]] <- list(sim.results = simulation.results, genome = hv.genome)
      
    } # Close the iteration loop
    
    # Return the rep list
    # return(rep.results)
    
  # End parlapply
  # }, mc.cores = n.cores)
  
  # Save the tp.change data
  filename <- file.path(save.dir, paste("simulation_results_", pop.makeup, "_", change, "_", 
                                        tp.formation, "_", date, ".RData", sep = "") )
                        
  save(list = c("experiment.sub.results", "change", "metadata"), file = filename)
  
  
} # Close the tp.change for loop


## Parse the results
# Create a filename to save
filename <- file.path(save.dir, paste("simulation_results_", pop.makeup, "_", 
                                      tp.formation, "_collective.RData", sep = ""))

# Create a vector of file paths
files <- list.files(save.dir, full.names = T) %>% 
  str_subset(pattern = paste("simulation_results_", pop.makeup, "_[A-Za-z]*_", 
                             tp.formation, "_[0-9]*-[0-9]*.RData", sep = ""))
  
# Parse
parse.results(files = files, filename = filename)
