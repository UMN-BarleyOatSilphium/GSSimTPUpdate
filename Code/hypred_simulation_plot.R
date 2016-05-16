## Script to graph results of GS simulations

# Load packages
library(plyr)

# Set working directory
setwd("C:/Users/Jeff/Google Drive/Barley Lab/Projects/Side Projects/Simulations/Barley_GS_Simulations/Results/")


# Load data
# filename <- "simulation_results_q100_sel0.1_popmakeup-MN_tpformation-window_collective_300416.RData"
# filename <- "simulation_results_q100_sel0.1_popmakeup-ND_tpformation-window_collective_300416.RData"
# filename <- "simulation_results_q100_sel0.1_popmakeup-MNxND_tpformation-window_collective_part3.RData"
# filename <- "simulation_results_q100_sel0.1_popmakeup-MN_tpformation-cumulative_collective_part3.RData"
# filename <- "simulation_results_q100_sel0.1_popmakeup-ND_tpformation-cumulative_collective_part3.RData"
# filename <- "simulation_results_q100_sel0.1_popmakeup-MNxND_tpformation-cumulative_collective_part3.RData"

# Allele freq experiment
setwd("C:/Users/Jeff/Google Drive/Barley Lab/Projects/Side Projects/Simulations/Barley_GS_Simulations/Results/Allele Freq Experiment/")

all.files <- list.files()
filename <- all.files[1]

load(filename)

pop.makeup = strsplit(x = substring(text = filename, first = regexpr(pattern = "popmakeup", text = filename)[1] + 10), split = "_")[[1]][1]
tp.formation = strsplit(x = substring(text = filename, first = regexpr(pattern = "tpformation", text = filename)[1] + 12), split = "_")[[1]][1]

n.cycles = length(collective.abbreviated.results[[1]][[1]][[1]][[1]])
n.reps = sum(unlist(lapply(X = collective.abbreviated.results[[1]][[1]], FUN = length)))

# Change in genetic variance over cycles  
V_g.list <- lapply(X = collective.abbreviated.results, function(tpc)
  do.call("cbind", lapply(X = tpc$candidate.variance.components.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) return(cycle$true$V_g) )))))

# Empty plot
plot(0, 
     type = "n",
     xlim = c(0, n.cycles),
     xlab = "Cycle Number",
     # ylim = range(pretty(range(V_g.list)))
     ylim = c(0,10),
     ylab = "Mean V_g",
     main = paste("Mean Genetic Variance Across Cycles", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n")
)

# Plotting shape factors
plot.shapes.factor <- factor(names(V_g.list))

# Add legend
legend("topright", legend = names(V_g.list), pch = as.numeric(factor(names(V_g.list))))

for (i in 1:length(V_g.list)) {
  
  # Find the mean and sd
  V_g.mu <- apply(X = V_g.list[[i]], MARGIN = 1, FUN = mean, na.rm = T)
  # Determine 95% confidence interval based on t distribution
  V_g.CI <- apply(X = V_g.list[[i]], MARGIN = 1, FUN = function(cycle) {
    t.per <- qt(p = (1 - (0.05 / 2)), df = length(cycle) - 1)
    t.per * ( sd(cycle) / sqrt(length(cycle)) )
  })
  
  # Add points to the plot
  x.jitter <- - (0.1 * scale(1:length(V_g.list), scale = F)[i])
  points(x = (1:n.cycles + x.jitter), V_g.mu, pch = as.numeric(plot.shapes.factor[i]), type = "p")
  # points(x = 1:n.cycles, scale = F)[i])), V_g.mu, pch = as.numeric(plot.shapes.factor[i]))


  # Add standard deviation bars
  segments(x0 = (1:n.cycles + x.jitter), y0 = (V_g.mu - V_g.CI), x1 = (1:n.cycles + x.jitter), y1 = (V_g.mu + V_g.CI))
  
}
  


# Change in pairwise diversity over cycles
div.list <- lapply(X = collective.abbreviated.results, function(tpc)
  do.call("cbind", lapply(X = tpc$pairwise.div.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) return(mean(cycle) ))))))

# Empty plot
plot(0, 
     type = "n",
     xlim = c(0, n.cycles),
     xlab = "Cycle Number",
     # ylim = range(pretty(range(V_g.list)))
     ylim = c(0, 0.3),
     ylab = "Mean Average Pairwise Diversity",
     main = paste("Mean Average Pairwise Diversity Across Cycles", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n")
)

# Plotting shape factors
plot.shapes.factor <- factor(names(div.list))

# Add legend
legend("bottomleft", legend = names(div.list), pch = as.numeric(factor(names(div.list))))

for (i in 1:length(div.list)) {
  
  # Find the mean and sd
  div.mu <- apply(X = div.list[[i]], MARGIN = 1, FUN = mean, na.rm = T)
  # Determine 95% confidence interval based on t distribution
  div.mu.CI <- apply(X = div.list[[i]], MARGIN = 1, FUN = function(cycle) {
    t.per <- qt(p = (1 - (0.05 / 2)), df = length(cycle) - 1)
    t.per * ( sd(cycle) / sqrt(length(cycle)) )
  })
  
  # Add points to the plot
  points(x = 1:n.cycles, div.mu, pch = as.numeric(plot.shapes.factor[i]))
  
  # Add standard deviation bars
  segments(x0 = 1:n.cycles, y0 = (div.mu - div.mu.CI), x1 = 1:n.cycles, y1 = (div.mu + div.mu.CI))
  
}



# Change in genotypic value over cycles
gen.mu.list <- lapply(X = collective.abbreviated.results, function(tpc)
  do.call("cbind", lapply(X = tpc$candidate.genotypic.value.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) return(mean(cycle) ))))))

# Empty plot
plot(0, 
     type = "n",
     xlim = c(0, n.cycles),
     xlab = "Cycle Number",
     # ylim = range(pretty(range(V_g.list)))
     ylim = c(-5, 30),
     ylab = "Mean Average Genotypic Value",
     main = paste("Mean Average Genotypic Values Across Cycles", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n")
)

# Plotting shape factors
plot.shapes.factor <- factor(names(gen.mu.list))

# Add legend
legend("bottomright", legend = names(gen.mu.list), pch = as.numeric(factor(names(gen.mu.list))))

for (i in 1:length(gen.mu.list)) {
  
  # Find the mean and sd
  gen.mu <- apply(X = gen.mu.list[[i]], MARGIN = 1, FUN = mean, na.rm = T)
  gen.mu.CI <- apply(X = gen.mu.list[[i]], MARGIN = 1, FUN = function(cycle) {
    t.per <- qt(p = (1 - (0.05 / 2)), df = length(cycle) - 1)
    t.per * ( sd(cycle) / sqrt(length(cycle)) )
  })
  
  # Add points to the plot
  x.jitter <- - (0.1 * scale(1:length(V_g.list), scale = F)[i])
  points(x = 1:n.cycles + x.jitter, gen.mu, pch = as.numeric(plot.shapes.factor[i]))
  
  # Add standard deviation bars
  segments(x0 = 1:n.cycles + x.jitter, y0 = (gen.mu - gen.mu.CI), x1 = 1:n.cycles + x.jitter, y1 = (gen.mu + gen.mu.CI))
  
}




# Change in the prediction accuracy over cycles
val.pred.list <- lapply(X = collective.abbreviated.results, function(tpc)
  list(pred_r = 
        do.call("cbind", lapply(X = tpc$validation.results.list, FUN = function(set) 
          sapply(set, function(rep) 
            sapply(rep, function(cycle) return(mean(cycle$pred.r) ))))),
       pred_r_sd = 
         do.call("cbind", lapply(X = tpc$validation.results.list, FUN = function(set) 
           sapply(set, function(rep) 
             sapply(rep, function(cycle) return(mean(cycle$pred.r.sd) )))))
  ))

# Empty plot
plot(0, 
     type = "n",
     xlim = c(0, n.cycles),
     xlab = "Cycle Number",
     # ylim = range(pretty(range(V_g.list)))
     ylim = c(-0.1, 1),
     ylab = "Mean Prediction Accuracy of True Genotypic Value (r)",
     main = paste("Mean Prediction Accuracy of True Genotypic Value Across Cycles", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n")
)

# Plotting shape factors
plot.shapes.factor <- factor(names(val.pred.list))

# Add legend
legend("topright", legend = names(val.pred.list), pch = as.numeric(factor(names(val.pred.list))))

for (i in 1:length(val.pred.list)) {
  
  # Find the mean and CI of the mean
  pred_r.mu <- apply(X = val.pred.list[[i]]$pred_r, MARGIN = 1, FUN = mean, na.rm = T)
  # Determine 95% confidence interval based on t distribution
  pred_r.mu.CI <- apply(X = val.pred.list[[i]]$pred_r, MARGIN = 1, FUN = function(cycle) {
    t.per <- qt(p = (1 - (0.05 / 2)), df = length(cycle) - 1)
    t.per * ( sd(cycle) / sqrt(length(cycle)) )
  })
  
  # Add points to the plot
  x.jitter <- - (0.1 * scale(1:length(V_g.list), scale = F)[i])
  points(x = 1:n.cycles + x.jitter, pred_r.mu, pch = as.numeric(plot.shapes.factor[i]))
  
  # Add standard deviation bars
  segments(x0 = 1:n.cycles + x.jitter, y0 = (pred_r.mu - pred_r.mu.CI), x1 = 1:n.cycles + x.jitter, y1 = (pred_r.mu + pred_r.mu.CI))
  
}


## Find the number of polymorphic markers used in each cycle
poly.marker.list <- lapply(X = collective.abbreviated.results, FUN = function(tpc) 
  do.call("cbind", sapply(tpc$allele.freq.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, FUN = function(cycle) {
        1 - (sum(cycle %in% c(0,1)) / length(cycle))
      }) ))))


# Empty plot
plot(0, 
     type = "n",
     xlim = c(0, n.cycles),
     xlab = "Cycle Number",
     # ylim = range(pretty(range(V_g.list)))
     ylim = c(0, 1),
     ylab = "Proportion of Markers That are Polymorphic",
     main = paste("Proportion of Markers That are Polymorphic", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n")
)

# Plotting shape factors
plot.shapes.factor <- factor(names(poly.marker.list))

# Add legend
legend("topright", legend = names(poly.marker.list), pch = as.numeric(factor(names(poly.marker.list))))

for (i in 1:length(poly.marker.list)) {
  
  # Find the mean and sd
  marker.mu <- apply(X = poly.marker.list[[i]], MARGIN = 1, FUN = mean, na.rm = T)
  marker.mu.CI <- apply(X = poly.marker.list[[i]], MARGIN = 1, FUN = function(cycle) {
    t.per <- qt(p = (1 - (0.05 / 2)), df = length(cycle) - 1)
    t.per * ( sd(cycle) / sqrt(length(cycle)) )
  })
  
  # Add points to the plot
  points(x = 1:n.cycles, marker.mu, pch = as.numeric(plot.shapes.factor[i]))
  
  # Add standard deviation bars
  segments(x0 = 1:n.cycles, y0 = (marker.mu - marker.mu.CI), x1 = 1:n.cycles, y1 = (marker.mu + marker.mu.CI))
  
}


# Size of training population
tp.size.list <- lapply(X = collective.abbreviated.results, FUN = function(tpc) 
  do.call("cbind", sapply(tpc$prediction.results.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, FUN = function(cycle) 
        return(cycle$parameters$n.TP))))) )

# Empty plot
plot(0, 
     type = "n",
     xlim = c(0, n.cycles),
     xlab = "Cycle Number",
     # ylim = range(pretty(range(V_g.list)))
     ylim = range(pretty(range(0, max(unlist(tp.size.list))))),
     ylab = "Training Population Size Used in Prediction",
     main = paste("Training Population Size", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n")
)

# Plotting shape factors
plot.shapes.factor <- factor(names(tp.size.list))

# Add legend
legend("topright", legend = names(tp.size.list), pch = as.numeric(factor(names(tp.size.list))))

for (i in 1:length(tp.size.list)) {
  
  # Find the mean and sd
  n.tp.mu <- apply(X = tp.size.list[[i]], MARGIN = 1, FUN = mean, na.rm = T)
  n.tp.mu.CI <- apply(X = tp.size.list[[i]], MARGIN = 1, FUN = function(cycle) {
    t.per <- qt(p = (1 - (0.05 / 2)), df = length(cycle) - 1)
    t.per * ( sd(cycle) / sqrt(length(cycle)) )
  })
  
  # Add points to the plot
  points(x = 1:n.cycles, n.tp.mu, pch = as.numeric(plot.shapes.factor[i]))
  
  # Add standard deviation bars
  segments(x0 = 1:n.cycles, y0 = (n.tp.mu - n.tp.mu.CI), x1 = 1:n.cycles, y1 = (n.tp.mu + n.tp.mu.CI))
  
}


# Plot site frequency spectra over the reps for each cycle
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
  

  
  

# Plotting shape factors
plot.shapes.factor <- factor(names(tp.size.list))

# Add legend
legend("topright", legend = names(tp.size.list), pch = as.numeric(factor(names(tp.size.list))))

for (i in 1:length(tp.size.list)) {
  
  # Find the mean and sd
  n.tp.mu <- apply(X = tp.size.list[[i]], MARGIN = 1, FUN = mean, na.rm = T)
  n.tp.mu.CI <- apply(X = tp.size.list[[i]], MARGIN = 1, FUN = function(cycle) {
    t.per <- qt(p = (1 - (0.05 / 2)), df = length(cycle) - 1)
    t.per * ( sd(cycle) / sqrt(length(cycle)) )
  })
  
  # Add points to the plot
  points(x = 1:n.cycles, n.tp.mu, pch = as.numeric(plot.shapes.factor[i]))
  
  # Add standard deviation bars
  segments(x0 = 1:n.cycles, y0 = (n.tp.mu - n.tp.mu.CI), x1 = 1:n.cycles, y1 = (n.tp.mu + n.tp.mu.CI))
  
}


# 
# # Look at the proportion of QTL that are fixed for one allele or the other at eahc
# ## cycle
# ## Find the number of polymorphic markers used in each cycle
# poly.marker.list <- lapply(X = collective.abbreviated.results, FUN = function(tpc) 
#   do.call("cbind", sapply(tpc$allele.freq.list, FUN = function(set) 
#     sapply(set, function(rep) 
#       sapply(rep, FUN = function(cycle) {
#         1 - (sum(cycle %in% c(0,1)) / length(cycle))
#       }) ))))
# 
# 
# # Empty plot
# plot(0, 
#      type = "n",
#      xlim = c(0, n.cycles),
#      xlab = "Cycle Number",
#      # ylim = range(pretty(range(V_g.list)))
#      ylim = c(0, 1),
#      ylab = "Proportion of Markers That are Polymorphic",
#      main = paste("Proportion of Markers That are Polymorphic", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n")
# )
# 
# # Plotting shape factors
# plot.shapes.factor <- factor(names(poly.marker.list))
# 
# # Add legend
# legend("topright", legend = names(poly.marker.list), pch = as.numeric(factor(names(poly.marker.list))))
# 
# for (i in 1:length(poly.marker.list)) {
#   
#   # Find the mean and sd
#   marker.mu <- apply(X = poly.marker.list[[i]], MARGIN = 1, FUN = mean, na.rm = T)
#   marker.mu.CI <- apply(X = poly.marker.list[[i]], MARGIN = 1, FUN = function(cycle) {
#     t.per <- qt(p = (1 - (0.05 / 2)), df = length(cycle) - 1)
#     t.per * ( sd(cycle) / sqrt(length(cycle)) )
#   })
#   
#   # Add points to the plot
#   points(x = 1:n.cycles, marker.mu, pch = as.numeric(plot.shapes.factor[i]))
#   
#   # Add standard deviation bars
#   segments(x0 = 1:n.cycles, y0 = (marker.mu - marker.mu.CI), x1 = 1:n.cycles, y1 = (marker.mu + marker.mu.CI))
#   
# }
