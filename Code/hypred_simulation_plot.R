## Script to graph results of GS simulations

# Set working directory
# setwd("C:/Users/Jeff/Google Drive/Barley Lab/Projects/Side Projects/Simulations/Barley_GS_Simulations/Results/")
setwd("C:/Users/Jeff/Google Drive/Barley Lab/Projects/Side Projects/Simulations/Barley_GS_Simulations/Results/Allele Freq Experiment/")



# Load data
# filename <- "simulation_results_q100_sel0.1_popmakeup-MN_tpformation-window_collective_part3.RData"
# filename <- "simulation_results_q100_sel0.1_popmakeup-ND_tpformation-window_collective_part3.RData"
# filename <- "simulation_results_q100_sel0.1_popmakeup-MNxND_tpformation-window_collective_part3.RData"
# filename <- "simulation_results_q100_sel0.1_popmakeup-MN_tpformation-cumulative_collective_part3.RData"
# filename <- "simulation_results_q100_sel0.1_popmakeup-ND_tpformation-cumulative_collective_part3.RData"
filename <- "simulation_results_q100_sel0.1_popmakeup-MNxND_tpformation-cumulative_collective_part3.RData"

# Allele freq experiment
all.files <- list.files()
filename <- all.files[5]

load(filename)

pop.makeup = strsplit(x = substring(text = filename, first = regexpr(pattern = "popmakeup", text = filename)[1] + 10), split = "_")[[1]][1]
tp.formation = strsplit(x = substring(text = filename, first = regexpr(pattern = "tpformation", text = filename)[1] + 12), split = "_")[[1]][1]

n.cycles = length(collective.abbreviated.results[[1]][[1]][[1]][[1]])


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
     ylim = c(0,15),
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
  V_g.sd <- apply(X = V_g.list[[i]], MARGIN = 1, FUN = sd, na.rm = T)
  
  # Add points to the plot
  x.jitter <- - (0.1 * scale(1:length(V_g.list), scale = F)[i])
  points(x = (1:n.cycles + x.jitter), V_g.mu, pch = as.numeric(plot.shapes.factor[i]), type = "p")
  # points(x = 1:n.cycles, scale = F)[i])), V_g.mu, pch = as.numeric(plot.shapes.factor[i]))


  # Add standard deviation bars
  segments(x0 = (1:n.cycles + x.jitter), y0 = (V_g.mu - V_g.sd), x1 = (1:n.cycles + x.jitter), y1 = (V_g.mu + V_g.sd))
  
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
  div.sd <- apply(X = div.list[[i]], MARGIN = 1, FUN = sd, na.rm = T)
  
  # Add points to the plot
  points(x = 1:n.cycles, div.mu, pch = as.numeric(plot.shapes.factor[i]))
  
  # Add standard deviation bars
  segments(x0 = 1:n.cycles, y0 = (div.mu - div.sd), x1 = 1:n.cycles, y1 = (div.mu + div.sd))
  
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
  gen.sd <- apply(X = gen.mu.list[[i]], MARGIN = 1, FUN = sd, na.rm = T)
  
  # Add points to the plot
  x.jitter <- - (0.1 * scale(1:length(V_g.list), scale = F)[i])
  points(x = 1:n.cycles + x.jitter, gen.mu, pch = as.numeric(plot.shapes.factor[i]))
  
  # Add standard deviation bars
  segments(x0 = 1:n.cycles + x.jitter, y0 = (gen.mu - gen.sd), x1 = 1:n.cycles + x.jitter, y1 = (gen.mu + gen.sd))
  
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
  
  # Find the mean and sd
  pred_r.mu <- apply(X = val.pred.list[[i]]$pred_r, MARGIN = 1, FUN = mean, na.rm = T)
  phen_r_sd.mu <- apply(X = val.pred.list[[i]]$pred_r_sd, MARGIN = 1, FUN = mean, na.rm = T)
  
  # Add points to the plot
  x.jitter <- - (0.1 * scale(1:length(V_g.list), scale = F)[i])
  points(x = 1:n.cycles + x.jitter, pred_r.mu, pch = as.numeric(plot.shapes.factor[i]))
  
  # Add standard deviation bars
  segments(x0 = 1:n.cycles + x.jitter, y0 = (pred_r.mu - phen_r_sd.mu), x1 = 1:n.cycles + x.jitter, y1 = (pred_r.mu + phen_r_sd.mu))
  
}


## Find the number of polymorphic markers used in each cycle
poly.marker.list <- lapply(X = collective.abbreviated.results, FUN = function(tpc) 
  do.call("cbind", sapply(tpc$prediction.results.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, FUN = function(cycle) 
        return(cycle$parameters$n.markers))))) )

# Empty plot
plot(0, 
     type = "n",
     xlim = c(0, n.cycles),
     xlab = "Cycle Number",
     # ylim = range(pretty(range(V_g.list)))
     ylim = c(0, 700),
     ylab = "Mean Number of Polymorphic Markers Used in Prediction",
     main = paste("Mean Number of Polymorphic Markers", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n")
)

# Plotting shape factors
plot.shapes.factor <- factor(names(poly.marker.list))

# Add legend
legend("topright", legend = names(poly.marker.list), pch = as.numeric(factor(names(poly.marker.list))))

for (i in 1:length(poly.marker.list)) {
  
  # Find the mean and sd
  marker.mu <- apply(X = poly.marker.list[[i]], MARGIN = 1, FUN = mean, na.rm = T)
  marker.sd <- apply(X = poly.marker.list[[i]], MARGIN = 1, FUN = sd, na.rm = T)
  
  # Add points to the plot
  points(x = 1:n.cycles, marker.mu, pch = as.numeric(plot.shapes.factor[i]))
  
  # Add standard deviation bars
  segments(x0 = 1:n.cycles, y0 = (marker.mu - marker.sd), x1 = 1:n.cycles, y1 = (marker.mu + marker.sd))
  
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
  n.tp.sd <- apply(X = tp.size.list[[i]], MARGIN = 1, FUN = sd, na.rm = T)
  
  # Add points to the plot
  points(x = 1:n.cycles, n.tp.mu, pch = as.numeric(plot.shapes.factor[i]))
  
  # Add standard deviation bars
  segments(x0 = 1:n.cycles, y0 = (n.tp.mu - n.tp.sd), x1 = 1:n.cycles, y1 = (n.tp.mu + n.tp.sd))
  
}


## Calculate the proportion of loci that are fixed for one allele or the other
# Patience. This takes some time
rare.var.list <- lapply(X = collective.abbreviated.results, function(tpc)
  do.call("cbind", lapply(X = tpc$allele.freq.list, FUN = function(set) 
    sapply(set, function(rep) 
      sapply(rep, function(cycle) {
        MAF <- sapply(cycle, function(freq) min(freq, 1 - freq))
        prop.rare <- sum(MAF == 0) / length(MAF)
        return(prop.rare) })))))

# Empty plot
plot(0, 
     type = "n",
     xlim = c(0, n.cycles),
     xlab = "Cycle Number",
     # ylim = range(pretty(range(V_g.list)))
     ylim = c(0, 1),
     ylab = "Mean Proportion of Loci With Fixed Allele",
     main = paste("Mean Proportion of Loci With Fixed Allele", paste("Population:", pop.makeup, ", TP formation:", tp.formation), sep = "\n")
)

# Plotting shape factors
plot.shapes.factor <- factor(names(rare.var.list))

# Add legend
legend("bottomright", legend = names(rare.var.list), pch = as.numeric(factor(names(rare.var.list))))

for (i in 1:length(val.pred.list)) {
  
  # Find the mean and sd
  ra.mu <- apply(X = rare.var.list[[i]], MARGIN = 1, FUN = mean, na.rm = T)
  ra.sd <- apply(X = rare.var.list[[i]], MARGIN = 1, FUN = sd, na.rm = T)
  
  # Add points to the plot
  points(x = 1:n.cycles, ra.mu, pch = as.numeric(plot.shapes.factor[i]))
  
  # Add standard deviation bars
  segments(x0 = 1:n.cycles, y0 = (ra.mu - ra.sd), x1 = 1:n.cycles, y1 = (ra.mu + ra.sd))
  
}

