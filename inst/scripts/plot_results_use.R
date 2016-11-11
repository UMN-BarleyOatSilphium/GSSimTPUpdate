## Plot the parsed results of the simulation

# Load packages
library(tidyverse, quietly = T)
library(stringr, quietly = T)
library(GSSimTPUpdate)

data("plotting_data")

tp.change.factors <- c(best = "Top", CDmean = "CDmean", nochange = "No Change", 
                       PEVmean = "PEVmean", random = "Random", worst = "Bottom") %>%
  as.factor()

# Number of cycles
n.cycles <- 15
  

### Figure 3 is the change in prediction accuracy

df1.acc <- plotting_data$accuracy %>%
  sim.summarize() %>%
  mutate(variable = "Prediction Accuracy") %>%
  filter(heritability == "0.5") %>%
  select(-heritability)

sim.ggplot(df.summary = df1.acc, main = "", 
           ylab = "Prediction Accuracy", 
           col.factors = tp.change.factors, facet.vars = "exp_name")

### Figure 4 is the change in genetic variance and genetic value

# Calculate mean and CI for each change-cycle combination over iterations
df1.genvar <- plotting_data$genvar %>%
  sim.summarize() %>%
  mutate(variable = "Genetic Variance") %>%
  filter(heritability == "0.5") %>%
  select(-heritability)

# Calculate mean and CI for each change-cycle combination over iterations
df1.genval <- plotting_data$genval %>%
  sim.summarize() %>%
  mutate(variable = "Genotypic Value") %>%
  filter(heritability == "0.5") %>%
  select(-heritability)


### Create a figure with both plots
df1 <- bind_rows(df1.genval, df1.genvar)

# Plot
sim.ggplot(df.summary = df1, main = "", ylab = "", col.factors = tp.change.factors)






### Figure 5 will have the genomic relationship, persistence of phase, inbreeding,
### and number of QTL fixed for an allele

# Persistence of phase

df1.persistence <- plotting_data$qtl_marker_LD %>%
  filter(extra1 == "persistence") %>%
  sim.summarize() %>%
  mutate(variable = "Persistence of LD Phase") %>%
  filter(heritability == "0.5") %>%
  select(-heritability)


### Relationship

df1.relatioship <- plotting_data$rel %>%
  sim.summarize() %>%
  mutate(variable = "Average Relationship\n(Scaled to Base Population)") %>%
  filter(heritability == "0.5") %>%
  select(-heritability)


# Inbreeding

df1.inbreeding <- plotting_data$inbreeding %>%
  sim.summarize() %>%
  mutate(variable = "Inbreeding\n(Scaled to Base Population)") %>%
  filter(heritability == "0.5") %>%
  select(-heritability)

# Fixed QTL

df1.qtlfreq <- plotting_data$fixedqtl %>%
  ungroup() %>%
  mutate(value = value / 100) %>%
  sim.summarize() %>%
  mutate(variable = "Proportion of Fixed QTL") %>%
  filter(heritability == "0.5") %>%
  select(-heritability)


df1 <- bind_rows(df1.persistence, df1.relatioship, df1.inbreeding, df1.qtlfreq) %>%
  mutate(variable = factor(variable, 
                           levels = c("Average Relationship\n(Scaled to Base Population)", 
                                      "Persistence of LD Phase",
                                      "Inbreeding\n(Scaled to Base Population)", 
                                      "Proportion of Fixed QTL")))

gp.combined <- sim.ggplot(df.summary = df1, main = "", ylab = "", 
                          col.factors = tp.change.factors, text.y.scaling = 1.01)
