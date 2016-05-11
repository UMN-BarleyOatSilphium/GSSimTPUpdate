## Genomic selection simulations
# Using hypred to simulate 10 cycles of selection using genomic prediction
# February 26, 2016

# Important information to save in the subset
# 1. Variance components
# 2. Allele frequencies
# 3. Pairwise diversity
# 4. The genome
# 5. Phenotypic mean
# 6. Genotypic mean
# 7. Prediction results and accuracy

# Arguments
args <- commandArgs(trailingOnly = T)

# First and only argument is the pop makeup
# If all the arguments are NA, stop
if (all(is.na(args))) {
  pop.makeup <- "MN"
  tp.formation <- "cumulative"
} else {
  pop.makeup <- args[1]
  tp.formation <- args[2]
}

# Are we using MSI?
MSI = FALSE

# Other simulation parameters
# Entry-mean heritability
h2 = 0.5
# How many cycles?
n.cycles = 25
# Number of QTL underlying trait
n.QTL = 100
# Selection intensity
GEBV.sel.intensity = 0.1
# Number of environments from which to draw phenotypes
n.env = 3

# How should the training population be formed in the next cycle?
# Questions! 
## 1. Should the top 150 lines be chosen based on GEBVs before phenotyping?
## 2. The parents for the next generation are selected based on the top x GEBVs. Is this suitable
# "best" means that the top 50 cycle x lines are added to the TP (based on phenotype)
# "worst" means that the bottom 50 cycle x lines are added to the TP (based on phenotype)
# "random" means that a random 50 of cycle x lines are added to the TP (based on phenotype)
# "best+random" means that the top 25 lines based on phenotype are combined with a random 25 (out of 120) to add to the TP
# "best+worst" means that the top 25 and the bottom 25 lines based on phenotype are added to the TP
# "no.change" # The same TP is used each cycle
tp.change = c("best", "worst", "random", "best+random", "no.change", "best+worst")
tp.size = 200
tp.update.increment = 50

# How to form the TP?
# "cumulative" means that a TP will grow with each addition of lines
# "random" means the TP candidates will be pooled and randomly sampled
# tp.formation = "cumulative"

# What combination of parents should we work with
# Choices are "MN", "ND", or "MNxND"
# pop.makeup = "MN"


# Computation parameters
n.iterations = 100

# Load the packages
library(hypred)
library(EMMREML)
library(rrBLUP)
library(boot)

# # Other tools
if (MSI) {
  setwd("/panfs/roc/groups/6/smithkp/neyhartj/Genomic_Selection/Simulations/Barley_GS_PEV/Hypred")
  # source("/panfs/roc/groups/6/smithkp/neyhartj/GitHub_Repos/Quant-Gen-Scripts/genotype_matrix_utils.R")
  n.cores = 16
} else {
  setwd("C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomic Selection/Side Projects/Simulations/Hypred simulations/")

  # devtools::source_url("https://raw.githubusercontent.com/neyhartj/Quant-Gen-Scripts/master/genotype_matrix_utils.R")
  n.cores = 1
}


# Source functions
source("Code/hypred_simulation_FUNCTIONS.R")
##

# ##### Manage input genotype data #####
# 
# # Load S6 genotype data
# parents.hmp <- read.table(file = "Files/21_parents_genotype_hmp.txt", header = T, check.names = F, as.is = T)
# CAP.hmp <- read.table(file = "Files/MN_ND_CAP_genotypes_hmp.txt", header = T, check.names = F, as.is = T)
# # Extract marker information
# marker.info <- parents.hmp[,c(1:4)]
# # Filter CAP markers that are not in the parent list
# CAP.hmp <- CAP.hmp[CAP.hmp$rs %in% parents.hmp$rs,]
# # Add marker information to the CAP.hmp data.frame
# CAP.hmp[,c(1:4)] <- parents.hmp[parents.hmp$rs %in% CAP.hmp$rs, c(1:4)]
# # Record demographic information for the parents and the CAP
# 
# 
# # Remove markers without map positions
# parents.hmp <- parents.hmp[!is.na(parents.hmp$pos),]
# CAP.hmp <- CAP.hmp[!is.na(CAP.hmp$pos),]
# # Order on chromosome, then position
# parents.hmp <- parents.hmp[order(parents.hmp$chrom, parents.hmp$pos),]
# CAP.hmp <- CAP.hmp[order(CAP.hmp$chrom, CAP.hmp$pos),]
# 
# # Convert to matrix
# parents.mat <- hmp2mat(parents.hmp)
# CAP.mat <- hmp2mat(CAP.hmp)
# 
# # Filter
# parents.filt <- filter.genos(geno.mat = parents.mat, min.MAF = 0.05, max.marker.missing = 0.1, max.entry.missing = 0.1, max.entry.het = 0.03, max.marker.het = 1)
# CAP.filt <- filter.genos(geno.mat = CAP.mat, min.MAF = 0.03, max.marker.missing = 0.1, max.marker.het = 1, max.entry.het = 0.25, max.entry.missing = 0.1)
# # Impute and Relationship matrix
# parents.impute <- A.mat(X = parents.filt, min.MAF = 0, max.missing = 1, impute.method = "mean", return.imputed = T)$imputed
# CAP.impute <- A.mat(X = CAP.filt, min.MAF = 0, max.missing = 1, impute.method = "mean", return.imputed = T)$imputed
# 
# # Demographic information
# CAP.demographics <- data.frame(Lines = c(row.names(CAP.impute)[!grepl(pattern = "^ND", x = row.names(CAP.impute))], row.names(CAP.impute)[grepl(pattern = "^ND", x = row.names(CAP.impute))]), Program = c( rep("MN", length(row.names(CAP.impute)[!grepl(pattern = "^ND", x = row.names(CAP.impute))])), rep("ND", length(row.names(CAP.impute)[grepl(pattern = "^ND", x = row.names(CAP.impute))])) ))
# # Sort on line name
# CAP.demographics <- CAP.demographics[order(CAP.demographics$Lines),]
# CAP.impute <- CAP.impute[order(row.names(CAP.impute)),]
# 
# # Round the marker codes
# parents.impute <- round(parents.impute)
# CAP.impute <- round(CAP.impute)
# 
# # Collect information on retained markers
# parents.markers <- colnames(parents.impute)
# parents.marker.info <- parents.hmp[parents.hmp$rs %in% parents.markers, c(1,3,4)]
# CAP.markers <- colnames(CAP.impute)
# CAP.marker.info <- CAP.hmp[CAP.hmp$rs %in% CAP.markers, c(1,3,4)]
# 
# # Change heterozygous markers to homozygous based on flanking markers
# set.seed(858)
# parents.impute <- het.impute(geno.mat = parents.impute, marker.map = parents.marker.info)
# set.seed(858)
# CAP.impute <- het.impute(geno.mat = CAP.impute, marker.map = CAP.marker.info)
# 
# # Relationship matrix
# parents.Amat <- A.mat(parents.impute, min.MAF = 0, max.missing = 1)
# CAP.Amat <- A.mat(CAP.impute, min.MAF = 0, max.missing = 1)
# # PCA plots
# plot.PCA(Gmat = parents.Amat, x.PC = 1, y.PC = 2)
# plot.PCA(Gmat = CAP.Amat, x.PC = 1, y.PC = 2, color.factor = factor(CAP.demographics$Program))
# 
# # Pairwise diversity
# parents.div <- SNP.pairwise.div(geno.mat = parents.impute)
# CAP.div <- SNP.pairwise.div(geno.mat = CAP.impute)
# 
# 
# 
# # Remove markers with duplicated positions
# parents.marker.info <- parents.marker.info[!duplicated(parents.marker.info$pos),]
# max.parent.markers <- min(table(parents.marker.info$chrom))
# CAP.marker.info <- CAP.marker.info[!duplicated(CAP.marker.info$pos),]
# max.CAP.markers <- min(table(CAP.marker.info$chrom))
# 
# # Sample the same number of markers per chromosome
# set.seed(12192)
# parent.rs.sampled <- unlist(tapply(X = parents.marker.info$rs, INDEX = parents.marker.info$chrom, FUN = function(rs) {
#   rs.sampled <- sample(x = rs, size = max.parent.markers)
# }))
# set.seed(12192)
# CAP.rs.sampled <- unlist(tapply(X = CAP.marker.info$rs, INDEX = CAP.marker.info$chrom, FUN = function(rs) {
#   rs.sampled <- sample(x = rs, size = max.CAP.markers)
# }))
# 
# # Find marker information for these sampled markers
# parents.marker.info.sampled <- parents.marker.info[parents.marker.info$rs %in% parent.rs.sampled,]
# CAP.marker.info.sampled <- CAP.marker.info[CAP.marker.info$rs %in% CAP.rs.sampled,]
# # Create maps for the genome
# parents.map.sampled <- as.numeric(parents.marker.info.sampled$pos) / 100
# parents.markers.sampled <- parents.marker.info.sampled$rs
# CAP.map.sampled <- as.numeric(CAP.marker.info.sampled$pos) / 100
# CAP.markers.sampled <- CAP.marker.info.sampled$rs
# 
# # Subsample the genotype matrices for the sampled markers
# parents.sampled <- parents.impute[, parents.markers.sampled]
# CAP.sampled <- CAP.impute[,CAP.markers.sampled]
# 
# # Change the homozogous parents into gametes
# parents.gametes <- geno.reverse.code(geno.mat = parents.sampled)
# CAP.gametes <- geno.reverse.code(geno.mat = CAP.sampled)
# # Add line names back in
# row.names(parents.gametes) <- paste(rep(row.names(parents.sampled), each = 2), 1:2, sep = ".")
# row.names(CAP.gametes) <- paste(rep(row.names(CAP.sampled), each = 2), 1:2, sep = ".")
# 

########## Start the simuation rounds ##########

# Load already-curated gamete data
# save(list = c("parents.gametes", "CAP.gametes", "max.CAP.markers", "CAP.map.sampled", "max.parent.markers", "parents.map.sampled", "CAP.sampled", "parents.sampled"), file = "Files/simulation_foundation_data_031116.RData")
load("Files/simulation_foundation_data_031116.RData")

hv.genome <- make.genome( n.chr = 7, 
                          chr.len = c(1.432, 1.7292, 1.801, 1.465, 1.899, 1.422, 1.625), 
                          n.chr.snps = max.CAP.markers,
                          genetic.map = CAP.map.sampled)

# Other input data
input.geno <- CAP.sampled
input.gametes <- CAP.gametes

# Global setup (this is used across all iterations)
# First we need to create the founding families. We will randomly mate 
## the MN x MN varieties, the MN x ND varieties, and the ND x ND varieties
line.names <- row.names(input.geno)
ND.lines <- grep(pattern = "^ND", x = line.names, value = T)
MN.lines <- setdiff(x = line.names, y = c(ND.lines))
# Create a list to store the line names
pop.makeup.list <- list(
  MN = list(p1 = MN.lines, p2 = MN.lines),
  ND = list(p1 = ND.lines, p2 = ND.lines),
  MNxND = list(p1 = MN.lines, p2 = ND.lines)
)

# Set the parent gamete data input
parent.lines.list <- pop.makeup.list[[pop.makeup]]
parent.gametes <- input.gametes

date <- format(Sys.time(), "%d%m%y-%H%M%S")

# Create the final experiment output list
experiment.results <- list()


library(parallel)

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
  
      #### Define trait parameters ####
      
      hv.genome <- trait.architecture(genome = hv.genome,
                                      n.QTL = n.QTL, 
                                      qtl.ids = "random", 
                                      qtl.dom.ids = NULL, 
                                      qtl.per.ids = NULL, 
                                      qtl.add.eff = "geometric", 
                                      qtl.dom.eff = NULL)
      
      # Summary of the genome
      # summary(hv.genome)
      
      # Create the genotype and gamete matrices for the TP
      TP.gametes <- subset.gametes(input.gametes, line.names = line.names)
      TP.genos <- genotype.markers(gametes = TP.gametes, genome = hv.genome, DH = F)
      # Phenotype the training population
      TP.phenos <- measure.values(genome = hv.genome, gametes = input.gametes, h2 = h2, n.env = 3)$pheno.values
      
      TP.phenos.i <- TP.phenos
      TP.genos.i <- TP.genos
      
      # If the TP format is random, create an initial random TP
      if (tp.formation == "random") {
        # Update the TP
        updated.TP <- update.training.population(current.TP.genos = TP.genos.i,
                                                 current.TP.phenos = TP.phenos.i,
                                                 update.method = "random",
                                                 TP.size = tp.size)
        TP.genos.i <- updated.TP$updated.TP.genos
        TP.phenos.i <- updated.TP$updated.TP.phenos
      }
      
      # The timing is about 180 seconds per simulation repetition for 20 cycles
      # Create an initial data list
      simulation.results <- list()
      # Loop over the number of cycles
      for (bc in 1:n.cycles) {

        # Complete a cycle
        cycle.results <- complete.cycle(genome = hv.genome, 
                                        h2 = h2, 
                                        female.parents = parent.lines.list$p1, 
                                        male.parents = parent.lines.list$p2, 
                                        parent.gamete.matrix = parent.gametes, 
                                        n.crosses = 40, 
                                        family.size = 30, 
                                        cycle.number = bc, 
                                        TP.genos = TP.genos.i, 
                                        TP.phenos = TP.phenos.i, 
                                        genotyping.S.gen = 2, 
                                        GEBV.sel.intensity = GEBV.sel.intensity, 
                                        GEBV.sel.type = "best", 
                                        phenotyping.S.gen = 4,
                                        n.env = n.env,
                                        mutation.rate.snp = 7e-8, 
                                        mutation.rate.qtl = 2.5e-5, 
                                        filter.when.genotyping = F, 
                                        validation = T,
                                        geno.summary.stats = T)
  
        
        # Designate potential new parents
        # This should be the selections made on GEBVs
        parent.lines <- cycle.results$genomic.selection.info$GEBV.selection.info$lines.sel
        parent.lines.list <- list(p1 = parent.lines, p2 = parent.lines)
        # The parents are selected and crossed at the F3 stage, so subset the gametes from the F3
        parent.gametes <- cycle.results$genomic.selection.info$selection.gametes
        
        # Update the training population
        updated.TP <- update.training.population(cycle.results = cycle.results,
                                                 current.TP.genos = TP.genos.i,
                                                 current.TP.phenos = TP.phenos.i,
                                                 n.additions = tp.update.increment,
                                                 update.method = tp.formation,
                                                 TP.size = tp.size,
                                                 selection.method = change)
        
        # Set the new TP genos and phenos
        TP.genos.i <- updated.TP$updated.TP.genos
        TP.phenos.i <- updated.TP$updated.TP.phenos
        
        # Parse and add data to the list
        # Relevant data includes:
        ## Crossing block
        ## Genotypic summary stats
        ## The results of the mixed model solution
        ## The phenotypic values of the selection set from genomic selection (after inbreeding)
        ## The results of validation (i.e. the accuracy of the GEBVs to the true genetic value)
        ## Variance components for the selection candidates and for the phenotyped individuals
        cycle.data <- list(
          crossing.block = cycle.results$crossing.info$crossing.block,
          geno.summary.stats = cycle.results$inbreeding.and.genotyping.info$summary.stats,
          model.fit = cycle.results$genomic.selection.info$prediction.info,
          selection.pheno.values = cycle.results$phenotyping.info$selection.values$pheno.values,
          validation.results = cycle.results$phenotyping.info$validation.results,
          var.components = list(selection.var.components = cycle.results$phenotyping.info$selection.values$var.components,
                                candidate.var.components = cycle.results$genomic.selection.info$candidate.var.components)
        )
        
        print( paste("Cycle", bc, "complete.") )
        
        cycle.name <- paste("cycle", bc, sep = "")
        simulation.results[[cycle.name]] <- cycle.data
        
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

  
  

################ PLOTTING ##################


# Checking function correctness
# Create founders
# founders <- hypredFounder(object = hv.genome, prob.snp = 0.5)
# # Make doubled haploid
# parent1 <- rbind(founders[1,], founders[1,])
# parent2 <- rbind(founders[2,], founders[2,])
# # Create F1
# F_1 <- make.family(genome = hv.genome, parent1.genome = parent1, parent2.genome = parent2, 
#                    N = 1000, generations = 0, pop.type = "inbred", cycle.number = 1, family.number = 1,
#                    mutation.rate.snp = 7 * 10^-8, mutation.rate.qtl = 2.5 * 10^-5)
# 
# # Get the genotypes
# F_1.genos <- hypredCode(hv.genome, F_1, DH = F, type = "-101")
# # Proportion of hets
# mean(apply(X = F_1.genos, MARGIN = 1, FUN = function(geno) sum(geno == 0) / length(geno)))
# # Genetic variance
# F_1.GVs <- hypredTruePerformance(hv.genome, F_1, DH = F)
# var(F_1.GVs)
# 
# # Create an F_3
# F_3 <- make.family(genome = hv.genome, parent1.genome = parent1, parent2.genome = parent2, 
#                    N = 1000, generations = 2, pop.type = "inbred", cycle.number = 1, family.number = 1,
#                    mutation.rate.snp = 7 * 10^-8, mutation.rate.qtl = 2.5 * 10^-5)
# 
# F_3.genos <- hypredCode(hv.genome, F_3, DH = F, type = "-101")
# mean(apply(X = F_3.genos, MARGIN = 1, FUN = function(geno) sum(geno == 0) / length(geno)))
# # Genetic variance
# F_3.GVs <- hypredTruePerformance(hv.genome, F_3, DH = F)
# var(F_3.GVs)
# 
# # Make an F_8
# F_8 <- make.family(genome = hv.genome, parent1.genome = parent1, parent2.genome = parent2, 
#                    N = 1000, generations = 7, pop.type = "inbred", cycle.number = 1, family.number = 1,
#                    mutation.rate.snp = 7 * 10^-8, mutation.rate.qtl = 2.5 * 10^-5)
# 
# F_8.genos <- hypredCode(hv.genome, F_8, DH = F, type = "-101")
# mean(apply(X = F_8.genos, MARGIN = 1, FUN = function(geno) sum(geno == 0) / length(geno)))
# # Genetic variance
# F_8.GVs <- hypredTruePerformance(hv.genome, F_8, DH = F)
# var(F_8.GVs)



# ##### SANDBOX #####
# 
# # Select the top 500 lines
# top500.phenos <- select.population(pheno.mat = MN.MN.phenos, sel.intensity = (500/2500), selection = "best")
# # Find the line names
# top500.lines <- top500.phenos$lines.sel
# # Randomly select 200 for the TP
# set.seed(424)
# TP.phenos <- select.population(pheno.mat = top500.phenos$value.sel, sel.intensity = (200/500), selection = "random")
# TP.lines <- TP.phenos$lines.sel
# TP.phenos <- TP.phenos$value.sel
# # Subset the TP gametes
# TP.gametes <- subset.gametes(gametes = MN.MN.pop, line.names = TP.lines)
# # Create genotype data
# TP.genos <- genotype.markers(gametes = TP.gametes, genome = hv.genome)
# # Relationship matrix
# TP.Amat <- A.mat(X = TP.genos, min.MAF = 0, max.missing = 1)
# 
# 
# # # Create a crossing block of the TP
# # set.seed(436)
# # TP.crossing.block <- make.crossing.block(parent1.lines = TP.lines, parent2.lines = TP.lines, n.crosses = 40, reciprocal = F)
# # # Make a population
# # C1.pop.F3 <- make.population(genome = hv.genome, 
# #                              named.parental.gametes = TP.gametes, 
# #                              crossing.block = TP.crossing.block, 
# #                              N = 30, 
# #                              cycle.number = 1, 
# #                              generations = 2, 
# #                              pop.type = "inbred", 
# #                              mutation.rate.snp = SNP.mu, 
# #                              mutation.rate.qtl = QTL.mu)
# # 
# # # Genetic variance of the C1.F3
# # C1.var.components <- phenotype.population(genome = hv.genome, gametes = C1.pop.F3, h2 = h2, return.genotype.values = T)$var.components
# # 
# # # Create genos of the C1
# # C1.genos <- genotype.markers(gametes = C1.pop.F3, genome = hv.genome)
# # 
# # # Make predictions
# # C1.solve <- make.predictions(pheno.train = TP.phenos, geno.train = TP.genos, geno.pred = C1.genos, model = "RRBLUP")
# # # Extract GEBVs
# # C1.GEBV <- C1.solve$GEBV
# # # Select the top 120 lines on GEBV
# # C1.sel <- select.population(pheno.mat = C1.GEBV, sel.intensity = 0.1, selection = "best")
# # C1.sel.values <- C1.sel$value.sel
# # C1.sel.lines <- C1.sel$lines.sel
# # 
# # # Subset C1.sel gametes
# # C1.sel.gametes <- subset.gametes(gametes = C1.pop.F3, line.names = C1.sel.lines)
# # # Advance the selected C1
# # C1.pop.F5 <- advance.family(genome = hv.genome,
# #                             pop.mat = C1.sel.gametes, 
# #                             pop.type = "inbred",
# #                             cycle.number = 1,
# #                             starting.generation = 3,
# #                             generations = 2, 
# #                             mutation.rate.snp = SNP.mu,
# #                             mutation.rate.qtl = QTL.mu)
# # 
# # # Variance components of the F5
# # C1.F5.var.components <- phenotype.population(genome = hv.genome, gametes = C1.pop.F5, h2 = h2, return.genotype.values = T)$var.components
# # 
# # # Phenotype the population
# # C1.sel.phenos <- phenotype.population(genome = hv.genome,
# #                                       gametes = C1.pop.F5, 
# #                                       h2 = h2, 
# #                                       return.genotype.values = T,
# #                                       pheno.dist = "normal")$pheno.values
# # 
# # C1.sel.genos <- genotype.markers(gametes = C1.sel.gametes, genome = hv.genome)
# # TP.allele.freq <- calculate.allele.freq(geno.mat = TP.genos, plot.sfs = T)
# # C1.allele.freq <- calculate.allele.freq(geno.mat = C1.sel.genos, plot.sfs = T)
# 
# # Make the crossing block for each scenario
# # set.seed(1910)
# MN.MN.crosses <- make.crossing.block(parent1.lines = MN.lines, parent2.lines = MN.lines, n.crosses = 10, reciprocal = F)
# # set.seed(0304)
# ND.ND.crosses <- make.crossing.block(parent1.lines = ND.lines, parent2.lines = ND.lines, n.crosses = 10, reciprocal = F)
# # set.seed(13300)
# MN.ND.crosses <- make.crossing.block(parent1.lines = MN.lines, parent2.lines = ND.lines, n.crosses = 10, reciprocal = F)
# 
# # # Make the initial populations based on each cross
# # MN.MN.pop <- make.population(genome = hv.genome, 
# #                              named.parental.gametes = parents.gametes, 
# #                              crossing.block = MN.MN.crosses, 
# #                              N = 250, 
# #                              generations = 5, 
# #                              pop.type = "inbred",
# #                              cycle.number = 0,
# #                              mutation.rate.snp = SNP.mu, 
# #                              mutation.rate.qtl = QTL.mu)
# 
# # ND.ND.pop <- make.population(genome = hv.genome, 
# #                            named.parental.gametes = parents.gametes, 
# #                            crossing.block = ND.ND.crosses, 
# #                            N = 250, 
# #                            generations = 5, 
# #                            pop.type = "inbred",
# #                            cycle.number = 0,
# #                            mutation.rate.snp = SNP.mu, 
# #                            mutation.rate.qtl = QTL.mu)
# # 
# # MN.ND.pop <- make.population(genome = hv.genome, 
# #                            named.parental.gametes = parents.gametes, 
# #                            crossing.block = MN.ND.crosses, 
# #                            N = 250, 
# #                            generations = 5, 
# #                            cycle.number = 0,
# #                            pop.type = "inbred",
# #                            mutation.rate.snp = SNP.mu, 
# #                            mutation.rate.qtl = QTL.mu)
# 
# 
# # Exploring ms for the simulation
# library(phyclust)
# 
# 
# 
# 
# 
# 
# # First create a randomly mating population of 10000 individuals over 50 generations
# # random.pool <- make.population(genome = hv.genome, founder1 = founder1, founder2 = founder2, N = 10000, pop.type = "random", gens = 50, mutation.rate.snp = SNP.mu, mutation.rate.qtl = QTL.mu)
# 
# # Select the top 500
# # Extract genotypic values
# random.pool.geno.values <- hypredTruePerformance(hv.genome, genotypes = random.pool, DH = F)
# # Generate phenotypic values
# V_g = var(random.pool.geno.values)
# V_e <- (V_g / h2) - V_g
# random.pool.pheno.values <- random.pool.geno.values + rnorm(N, 0, sqrt(V_e))
# # Select the top 500
# breeding.pool.index <- order(random.pool.pheno.values, decreasing = T)[1:500]
# 
# # Subset the random pool
# breeding.pool.index <- (breeding.pool.index * 2) - 1
# breeding.pool.index <- sort(c( sort(breeding.pool.index), sort(breeding.pool.index + 1) ))
# breeding.pool <- random.pool[breeding.pool.index,]
# 
# # Randomly mate the breeding pool for 20 generations
# breeding.pool.rand <- change.population(genome = hv.genome, population = breeding.pool, pop.type = "random", gens = 20, mutation.rate.snp = SNP.mu, mutation.rate.qtl = QTL.mu)
# 
# 
# 
# # Inbreed the whole population
# inbred.pool <- change.population(genome = hv.genome, population = random.pool, pop.type = "inbred", gens = 10, mutation.rate.snp = SNP.mu, mutation.rate.qtl = QTL.mu)
# 
# 
# ## Create the training population
# # Randomly select 200 lines from the top 
# 
# 
# 
# 
# # Make an inbred population
# inbred.pop <- make.population(genome = hv.genome, founder1 = founder1, founder2 = founder2, N = N, pop.type = "inbred", gens = 8, mutation.rate.snp = SNP.mu, mutation.rate.qtl = QTL.mu)
# 
# 
# # Obtain genotypic values
# geno.values <- hypredTruePerformance(hv.genome, inbred.pop, DH = F)
# 
# # Genetic variance
# V_g <- var(geno.values)
# # Error variance
# V_e <- (V_g / h2) - V_g
# 
# # Generate phenotypes
# pheno.values <- geno.values + rnorm(N, 0, sqrt(V_e))
# 
# system.time(inbred.pop <- make.population(genome = hv.genome, founder1 = founder1, founder2 = founder2, N = N, pop.type = "inbred", gens = 8, mutation.rate.snp = SNP.mu, mutation.rate.qtl = QTL.mu))
# 
# MAF <- apply(X = hypredCode(hv.genome, inbred.pop, DH = F, type = "-101") + 1, MARGIN = 2, FUN = function(m) {
#   p <- mean(m) / 2
#   MAF <- min(p, 1-p)
#   return(MAF)
# })
# 
# #################
