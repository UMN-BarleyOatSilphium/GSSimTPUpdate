## Script to graph results of GS simulations

# Load packages
library(plyr)

# Set working directory
# Base experiment
setwd("C:/Users/Jeff/Google Drive/Barley Lab/Projects/Side Projects/Simulations/BarleySimGS-TPUpdate/Results/Base Experiment/")

# Load data
all.files <- list.files()
filename <- all.files[3]

# # Allele freq experiment
# setwd("C:/Users/Jeff/Google Drive/Barley Lab/Projects/Side Projects/Simulations/BarleySimGS-TPUpdate/Results/Allele Freq Experiment/")
# 
# # Load data
# all.files <- list.files()
# filename <- all.files[1]

load(filename)

pop.makeup = strsplit(x = substring(text = filename, first = regexpr(pattern = "popmakeup", text = filename)[1] + 10), split = "_")[[1]][1]
tp.formation = strsplit(x = substring(text = filename, first = regexpr(pattern = "tpformation", text = filename)[1] + 12), split = "_")[[1]][1]

n.cycles = length(collective.abbreviated.results[[1]][[1]][[1]][[1]])
n.reps = sum(unlist(lapply(X = collective.abbreviated.results[[1]][[1]], FUN = length)))

tp.change.factors <- as.factor(names(collective.abbreviated.results))






# Define a function to plot
sim.plot <- function(data.list, 
                     xlim = c(1, n.cycles), 
                     xlab = "Cycle Number", 
                     ylim = NULL, 
                     ylab, 
                     main,
                     legend.pos,
                     just.data = FALSE # If true, return the estimates of the mean and CI and stop the function.
){
  
  # Extract the names from the data.list for use as a factor
  list.names <- names(data.list)
  list.names.factor <- factor(list.names)
  list.names.numeric <- as.numeric(list.names.factor)
  
  # Apply a function over the data.list to calculate the mean and a confidence
  ## interval at every cycle
  data.parameters <- lapply(X = data.list, FUN = function(data) {
    
    # Find the mean across iterations
    mu <- apply(X = data, MARGIN = 1, FUN = mean)
    
    # Calculate a confidence interval based on a t-distribution
    CI <- apply(X = data, MARGIN = 1, FUN = function(cycle) {
      t.per <- qt(p = (1 - (0.05 / 2)), df = length(cycle) - 1)
      t.per * ( sd(cycle) / sqrt(length(cycle)) ) }) 
    
    list(mu = mu, CI = CI) })
  
  # If just the data is requested, return it
  if (just.data) {
    return(data.parameters)
  }
  
  
  # If ylim is null, use a pretty range to determine it
  if (is.null(ylim)) {
    ylim <- range(pretty(range(sapply(data.parameters, function(sublist) sublist$mu), na.rm = T)))
  }
  
  # Create the empty plot first
  plot(0, type = "n", xlim = xlim, xlab = xlab, ylim = ylim, ylab = ylab, main = main)
  
  # Add the legend
  legend(legend.pos, legend = list.names.factor, pch = list.names.numeric, col = list.names.numeric)
  
  # Iterate over the data.parameters list
  for (i in 1:length(data.parameters)) {
    
    # Add points to the plot
    # Create jitter
    x.jitter <- - (0.1 * scale(1:length(data.parameters), scale = F)[i])
    points(x = (seq(xlim[1], xlim[2]) + x.jitter), data.parameters[[i]]$mu, pch = list.names.numeric[i], type = "b", col = list.names.numeric[i])

    # Add CI bars
    segments(x0 = (seq(xlim[1], xlim[2]) + x.jitter), 
             y0 = (data.parameters[[i]]$mu - data.parameters[[i]]$CI), 
             x1 = (seq(xlim[1], xlim[2]) + x.jitter), 
             y1 = (data.parameters[[i]]$mu + data.parameters[[i]]$CI))
  } # Close the for loop
  
} # Close the function
    


# Change in genetic variance over cycles  
V_g.list <- lapply(X = collective.abbreviated.results, function(tpc)
  do.call("cbind", lapply(X = tpc$candidate.variance.components.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) return(cycle$V_g) )))))

# Plot
sim.plot(data.list = V_g.list, 
         ylab = "Genetic Variance", 
         main = paste("Genetic Variance", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
         legend.pos = "bottomleft")
  



# Change in genotypic value over cycles
gen.mu.list <- lapply(X = collective.abbreviated.results, function(tpc)
  do.call("cbind", lapply(X = tpc$candidate.genotypic.value.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) cycle )))))


# Plot
sim.plot(data.list = gen.mu.list, 
         ylab = "Genotypic Value", 
         main = paste("Genotypic Value", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
         legend.pos = "topleft")





##### Calculate response to selection
# The predicted response to selection will be calculated as R = k_p * sqrt(V_a) * r_MG where:
## R is response to selection 
## k_p is the standardized selection coefficient (for 100 / 1200, k_p = 1.839)
## sqrt(V_a) is the additive genetic standard deviation (V_g = V_a in our simulation)
## r_MG is the correlation between predicted and observed genotypic values
# The observed reseponse to selection will be calculated as R = mu_Ck+1 - mu_Ck where
## S is the selection differential, calculated as mu_Ck.selected - mu_Ck where
### mu_Ck+1 is the mean of all candidates at cycle k + 1
### mu_Ck is the mean of all candidates at cycle k

k_p = 1.839 # Set k_p

# Predicted response
pred.R.list <- lapply(X = collective.abbreviated.results, function(tpc) {
  # Make a k x r data.frame for each parameter where k is the number of cycle and r is the number of iterations
  # Predicted response
  V_g.df <- do.call("cbind", lapply(X = tpc$candidate.variance.components.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) cycle$V_g ))))
  
  r_MG.df <- do.call("cbind", lapply(X = tpc$validation.results.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) cycle$pred.r ))))
  
  # Calculate predicted response
  R_exp.df = k_p * sqrt(V_g.df) * r_MG.df
  # Trim the last row and rename rows
  R_exp.df.i <- R_exp.df[-nrow(R_exp.df),]; row.names(R_exp.df.i) <- row.names(R_exp.df)[-1]
  
  return(R_exp.df.i) })

obs.R.list <- lapply(X = collective.abbreviated.results, function(tpc) {
  # Observed response
  mu.df <- do.call("cbind", lapply(X = tpc$candidate.genotypic.value.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) cycle ))))
  
  # Calculate the difference between adjacent rows (cycles)
  R_obs.df = diff(mu.df)
  
  return(R_obs.df) })

# Plot
sim.plot(data.list = obs.R.list,
         xlim = c(2, n.cycles),
         ylab = "Response to Selection", 
         main = paste("Observed Response to Selection", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
         legend.pos = "bottomleft")




# Change in the prediction accuracy over cycles
val.pred.list <- lapply(X = collective.abbreviated.results, function(tpc)
  do.call("cbind", lapply(X = tpc$validation.results.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) cycle$pred.r ) ))))

# Plot
sim.plot(data.list = val.pred.list,
         ylab = "Realized Prediction Accuracy", 
         main = paste("Realized Prediction Accuracy", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
         legend.pos = "bottomleft")





## Find the number of polymorphic markers used in each cycle
poly.marker.list <- lapply(X = collective.abbreviated.results, FUN = function(tpc) 
  do.call("cbind", sapply(tpc$allele.freq.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, FUN = function(cycle) {
        1 - (sum(cycle %in% c(0,1)) / length(cycle))
      }) ))))

# Plot
sim.plot(data.list = poly.marker.list,
         ylab = "Proportion of Markers that Are Polymorphic", 
         main = paste("Polymorphic Marker Proportion", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
         legend.pos = "bottomleft")


# Size of training population
tp.size.list <- lapply(X = collective.abbreviated.results, FUN = function(tpc) 
  do.call("cbind", sapply(tpc$prediction.results.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, FUN = function(cycle) 
        return(cycle$parameters$n.TP))))) )

# Plot
sim.plot(data.list = tp.size.list,
         ylab = "Training Population Size", 
         main = paste("Training Population Size Across Cycles", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
         legend.pos = "topleft")



# Plot site frequency spectra over the reps for each cycle
library(plyr)
# First parse the results to get the marker count at the designated maf breaks
sfs.count.list <- lapply(collective.abbreviated.results, function(tpc)
  unlist(lapply(tpc$allele.freq.list, function(set)
    lapply(set, function(rep)
      sapply(rep, function(cycle) {
        maf <- sapply(cycle, function(freq) min(freq, 1 - freq))
        hist(maf, breaks = seq(0, 0.5, 0.1), plot = F)$counts }))),
    recursive = F))

# Iterate over the different treatments in the list
for (i in 1:length(sfs.count.list)) {
  
  # Create a single data.frame of the mean count in each break at each cycle across reps
  count.cycle.mu <- aaply(laply(sfs.count.list[[i]], as.matrix), c(2,3), mean)
  count.cycle.mu.CI <- aaply(laply(sfs.count.list[[i]], as.matrix), c(2,3), function(data) {
    t.per <- qt(p = (1 - (0.05 / 2)), df = length(data) - 1)
    t.per * ( sd(data) / sqrt(length(data)))
  })
  
  # Set the mapping margins
  frame()
  par(mfrow = c(2,8))
  
  # Iterate over each cycle
  for (y in 1:n.cycles) {
  
    # Empty plot
    plot(0, 
         type = "n",
         xlim = c(0, 0.5),
         xlab = "",
         ylim = c(0, 600),
         ylab = ""
    )
    
    # Sequence of breaks
    breaks = seq(0, 0.5, 0.1)
    
    # Iterate over the counts in each bin
    for (t in 1:length(count.cycle.mu[,y])) {
      xleft = breaks[t]
      xright = breaks[t+1]
      ytop = count.cycle.mu[t,y]
      rect(xleft = xleft, ybottom = 0, xright = xright, ytop = ytop)
    }
    
  } # Close per-cycle loop
  
} # Close the treatment loop




# Change in QTL-marker LD over cycles
# Each QTL's max LD
qtl.marker.max.LD.list <- lapply(X = collective.abbreviated.results, function(tpc)
  do.call("cbind", lapply(X = tpc$qtl.marker.LD.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) cycle$mean.max ) ))))

# Plot
sim.plot(data.list = qtl.marker.max.LD.list,
         ylab = "Linkage Disequilibrium (r)", 
         main = paste("Mean LD of Polymorphic QTL with Marker in Highest LD", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
         legend.pos = "bottomleft")

# Mean LD within 50 cM window
qtl.marker.mean.LD.list <- lapply(X = collective.abbreviated.results, function(tpc)
  do.call("cbind", lapply(X = tpc$qtl.marker.LD.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) cycle$mean.window ) ))))

sim.plot(data.list = qtl.marker.mean.LD.list,
         ylab = "Linkage Disequilibrium (r)",
         main = paste("LD of Polymorphic QTL with Markers Within 50 cM", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
         legend.pos = "bottomleft")



# Change in average relationship between TP and candidates
relationship.list <- lapply(X = collective.abbreviated.results, function(tpc)
  do.call("cbind", lapply(X = tpc$relationship.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) cycle )))) )

# Plot
sim.plot(data.list = relationship.list,
         ylab = "Additive Genetic Relationship", 
         main = paste("Scaled Additive Relationship Between Training Set and Candidates", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n"), 
         legend.pos = "bottomleft")
  
