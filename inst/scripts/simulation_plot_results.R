## Script to graph results of GS simulations

# Load packages
library(dplyr)
library(stringr)
library(GSsim.TPUpdate)

# Designate the directory with the files
results.dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Side Projects/Simulations/GSsim.TPUpdate/Results/Base Experiment/"

# Load data from the allele frequency experiment
all.files <- list.files(results.dir, full.names = T) %>%
  str_subset(pattern = "collective")

filename <- all.files[1]

load(filename)

# Experimental parameters
pop.makeup <- str_extract(string = filename, pattern = 'results_[A-Za-z]{2,5}') %>%
  str_extract(pattern = '[A-Za-z]{2,5}$')

tp.formation <- 
  str_extract(string = filename, pattern = paste0(c("cumulative", "window"), collapse = "|"))

n.cycles = collective.abbreviated.results[[1]][[1]][[1]][[1]] %>%
  length()
n.reps = lapply(X = collective.abbreviated.results[[1]][[1]], FUN = length) %>%
  unlist() %>%
  sum()

tp.change.factors <- collective.abbreviated.results %>%
  names() %>%
  as.factor()



    


# Change in genetic variance over cycles  
V_g.list <- lapply(X = collective.abbreviated.results, function(tpc)
  do.call("cbind", lapply(X = tpc$candidate.variance.components.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) return(cycle) )))))

# Plot
plot.result(data.list = V_g.list, 
            ylab = "Genetic Variance", 
            main = paste("Genetic Variance", paste("Population:", pop.makeup, 
                                                   ", TP formation:", tp.formation), sep = "\n"), 
            legend.pos = "bottomleft")
  



# Change in genotypic value over cycles
gen.mu.list <- lapply(X = collective.abbreviated.results, function(tpc)
  do.call("cbind", lapply(X = tpc$candidate.genotypic.value.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) cycle )))))


# Plot
plot.result(data.list = gen.mu.list, 
         ylab = "Genotypic Value", 
         main = paste("Genotypic Value", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
         legend.pos = "topleft")




# Change in the prediction accuracy over cycles
val.pred.list <- lapply(X = collective.abbreviated.results, function(tpc)
  do.call("cbind", lapply(X = tpc$validation.results.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) cycle$pred.r ) ))))

# Plot
plot.result(data.list = val.pred.list,
         ylab = "Realized Prediction Accuracy", 
         main = paste("Realized Prediction Accuracy", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
         legend.pos = "bottomleft")





# Change in QTL-marker LD over cycles
# Each QTL's max LD
qtl.marker.max.LD.list <- lapply(X = collective.abbreviated.results, function(tpc)
  do.call("cbind", lapply(X = tpc$qtl.marker.LD.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) cycle$mean.max.genome ) ))))

# Plot
plot.result(data.list = qtl.marker.max.LD.list,
         ylab = "Linkage Disequilibrium (r)", 
         main = paste("Mean LD of Polymorphic QTL with Marker in Highest LD", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
         legend.pos = "bottomleft")

# Mean LD across whole genome
qtl.marker.mean.LD.list <- lapply(X = collective.abbreviated.results, function(tpc)
  do.call("cbind", lapply(X = tpc$qtl.marker.LD.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) cycle$mean.genome ) ))))

plot.result(data.list = qtl.marker.mean.LD.list,
         ylab = "Linkage Disequilibrium (r)",
         main = paste("LD of Polymorphic QTL with Markers Within 50 cM", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
         legend.pos = "bottomleft")

# Persistence of LD phase
persistence.of.phase.list <- lapply(X = collective.abbreviated.results, function(tpc)
  do.call("cbind", lapply(X = tpc$qtl.marker.LD.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) cycle$persistence ) ))))

plot.result(data.list = persistence.of.phase.list,
         ylab = "Correlatio of r",
         main = paste("Persistance of LD Phase Between TP and Selection Candidates", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
         legend.pos = "bottomleft")


# Change in average relationship between TP and candidates
relationship.list <- lapply(X = collective.abbreviated.results, function(tpc)
  do.call("cbind", lapply(X = tpc$relationship.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) cycle )))) )

# Plot
plot.result(data.list = relationship.list,
         ylab = "Additive Genetic Relationship\n(With Respect to the Base Population)", 
         main = paste("Scaled Additive Relationship Between Training Set and Candidates", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
         legend.pos = "topleft")


# Expected heterozygosity
exp.het.list <- lapply(X = collective.abbreviated.results, function(tpc)
  do.call("cbind", lapply(X = tpc$tp.update.exp.het.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) cycle )))) )

# Plot
plot.result(data.list = exp.het.list,
         ylab = "Additive Genetic Relationship\n(With Respect to the Base Population)", 
         main = paste("Expected Heterozygosity of TP Additions", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
         legend.pos = "bottomleft")







# 
# 
# ##### Calculate response to selection
# # The predicted response to selection will be calculated as R = k_p * sqrt(V_a) * r_MG where:
# ## R is response to selection 
# ## k_p is the standardized selection coefficient (for 100 / 1200, k_p = 1.839)
# ## sqrt(V_a) is the additive genetic standard deviation (V_g = V_a in our simulation)
# ## r_MG is the correlation between predicted and observed genotypic values
# # The observed reseponse to selection will be calculated as R = mu_Ck+1 - mu_Ck where
# ## S is the selection differential, calculated as mu_Ck.selected - mu_Ck where
# ### mu_Ck+1 is the mean of all candidates at cycle k + 1
# ### mu_Ck is the mean of all candidates at cycle k
# 
# k_p = 1.839 # Set k_p
# 
# # Predicted response
# pred.R.list <- lapply(X = collective.abbreviated.results, function(tpc) {
#   # Make a k x r data.frame for each parameter where k is the number of cycle and r is the number of iterations
#   # Predicted response
#   V_g.df <- do.call("cbind", lapply(X = tpc$candidate.variance.components.list, FUN = function(set) 
#     sapply(set, function(rep) 
#       sapply(rep, function(cycle) cycle$V_g ))))
#   
#   r_MG.df <- do.call("cbind", lapply(X = tpc$validation.results.list, FUN = function(set) 
#     sapply(set, function(rep) 
#       sapply(rep, function(cycle) cycle$pred.r ))))
#   
#   # Calculate predicted response
#   R_exp.df = k_p * sqrt(V_g.df) * r_MG.df
#   # Trim the last row and rename rows
#   R_exp.df.i <- R_exp.df[-nrow(R_exp.df),]; row.names(R_exp.df.i) <- row.names(R_exp.df)[-1]
#   
#   return(R_exp.df.i) })
# 
# obs.R.list <- lapply(X = collective.abbreviated.results, function(tpc) {
#   # Observed response
#   mu.df <- do.call("cbind", lapply(X = tpc$candidate.genotypic.value.list, FUN = function(set) 
#     sapply(set, function(rep) 
#       sapply(rep, function(cycle) cycle ))))
#   
#   # Calculate the difference between adjacent rows (cycles)
#   R_obs.df = diff(mu.df)
#   
#   return(R_obs.df) })
# 
# # Plot
# sim.plot(data.list = obs.R.list,
#          xlim = c(2, n.cycles),
#          ylab = "Response to Selection", 
#          main = paste("Observed Response to Selection", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
#          legend.pos = "bottomleft")
# 



# 
# 
# 
# ## Find the number of polymorphic markers used in each cycle
# poly.marker.list <- lapply(X = collective.abbreviated.results, FUN = function(tpc) 
#   do.call("cbind", sapply(tpc$allele.freq.list, FUN = function(set) 
#     sapply(set, function(rep) 
#       sapply(rep, FUN = function(cycle) {
#         1 - (sum(cycle %in% c(0,1)) / length(cycle))
#       }) ))))
# 
# # Plot
# sim.plot(data.list = poly.marker.list,
#          ylab = "Proportion of Markers that Are Polymorphic", 
#          main = paste("Polymorphic Marker Proportion", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
#          legend.pos = "bottomleft")
# 

# # Size of training population
# tp.size.list <- lapply(X = collective.abbreviated.results, FUN = function(tpc) 
#   do.call("cbind", sapply(tpc$prediction.results.list, FUN = function(set) 
#     sapply(set, function(rep) 
#       sapply(rep, FUN = function(cycle) 
#         return(cycle$parameters$n.TP))))) )
# 
# # Plot
# sim.plot(data.list = tp.size.list,
#          ylab = "Training Population Size", 
#          main = paste("Training Population Size Across Cycles", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
#          legend.pos = "topleft")

# 
# 
# # Plot site frequency spectra over the reps for each cycle
# library(plyr)
# # First parse the results to get the marker count at the designated maf breaks
# sfs.count.list <- lapply(collective.abbreviated.results, function(tpc)
#   unlist(lapply(tpc$allele.freq.list, function(set)
#     lapply(set, function(rep)
#       sapply(rep, function(cycle) {
#         maf <- sapply(cycle, function(freq) min(freq, 1 - freq))
#         hist(maf, breaks = seq(0, 0.5, 0.1), plot = F)$counts }))),
#     recursive = F))
# 
# # Iterate over the different treatments in the list
# for (i in 1:length(sfs.count.list)) {
#   
#   # Create a single data.frame of the mean count in each break at each cycle across reps
#   count.cycle.mu <- aaply(laply(sfs.count.list[[i]], as.matrix), c(2,3), mean)
#   count.cycle.mu.CI <- aaply(laply(sfs.count.list[[i]], as.matrix), c(2,3), function(data) {
#     t.per <- qt(p = (1 - (0.05 / 2)), df = length(data) - 1)
#     t.per * ( sd(data) / sqrt(length(data)))
#   })
#   
#   # Set the mapping margins
#   frame()
#   par(mfrow = c(2,8))
#   
#   # Iterate over each cycle
#   for (y in 1:n.cycles) {
#   
#     # Empty plot
#     plot(0, 
#          type = "n",
#          xlim = c(0, 0.5),
#          xlab = "",
#          ylim = c(0, 600),
#          ylab = ""
#     )
#     
#     # Sequence of breaks
#     breaks = seq(0, 0.5, 0.1)
#     
#     # Iterate over the counts in each bin
#     for (t in 1:length(count.cycle.mu[,y])) {
#       xleft = breaks[t]
#       xright = breaks[t+1]
#       ytop = count.cycle.mu[t,y]
#       rect(xleft = xleft, ybottom = 0, xright = xright, ytop = ytop)
#     }
#     
#   } # Close per-cycle loop
#   
# } # Close the treatment loop
# 

