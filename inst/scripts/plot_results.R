## Plot the parsed results of the simulation

# Load packages
library(tidyverse, quietly = T)
library(stringr, quietly = T)
library(GSsim.TPUpdate)


# Set directory and grab the files
results.dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Side Projects/Simulations/GSsim.TPUpdate/Results"

# File path for figures
figures.dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Side Projects/Simulations/GSsim.TPUpdate/figures/"

plot.data.file <- list.files(results.dir, full.names = T)
load(plot.data.file)

tp.change.factors <- c(best = "Top", CDmean = "CDmean", nochange = "No Change", 
                       PEVmean = "PEVmean", random = "Random", worst = "Bottom") %>%
  as.factor()

# Number of cycles
n.cycles <- df.list$accuracy$cycle %>% unique() %>% na.omit() %>% length()
  

### Figure 3 is the change in prediction accuracy

df1.acc <- df.list$accuracy %>%
  sim.summarize() %>%
  mutate(variable = "Prediction Accuracy") %>%
  filter(heritability == "0.5") %>%
  select(-heritability)

sim.ggplot(df.summary = df1.acc, main = "", 
           ylab = "Prediction Accuracy", 
           col.factors = tp.change.factors, facet.vars = "exp_name")

ggsave(filename = file.path(figures.dir, "pred_acc_combined.jpg"),
       height = 5, width = 9)


### Figure 4 is the change in genetic variance and genetic value

# Calculate mean and CI for each change-cycle combination over iterations
df1.genvar <- df.list$genvar %>%
  sim.summarize() %>%
  mutate(variable = "Genetic Variance") %>%
  filter(heritability == "0.5") %>%
  select(-heritability)

# Calculate mean and CI for each change-cycle combination over iterations
df1.genval <- df.list$genval %>%
  sim.summarize() %>%
  mutate(variable = "Genotypic Value") %>%
  filter(heritability == "0.5") %>%
  select(-heritability)


### Create a figure with both plots
df1 <- bind_rows(df1.genval, df1.genvar)

# Plot
sim.ggplot(df.summary = df1, main = "", ylab = "", col.factors = tp.change.factors)

ggsave(filename = file.path(figures.dir, "gen_val_gen_var_combined.jpg"),
       height = 8, width = 9)



## Alternatively, the true response to selection

df1.resp <- df.list$genval %>% 
  group_by(heritability, exp_name, change, iter, cycle) %>% 
  filter(row_number() == 1) %>%
  group_by(heritability, exp_name, change, iter) %>%
  mutate(value = c(NA, diff(value))) %>% 
  na.omit() %>%
  ungroup() %>%
  sim.summarize() %>%
  # Add the variable designator
  mutate(variable = "Response to Selection") %>%
  filter(heritability == "0.5") %>%
  select(-heritability)

sim.ggplot(df.summary = df1.resp, main = "Response to Selection Across Cycles", 
           ylab = "Response to Selection", 
           col.factors = tp.change.factors, facet.vars = "exp_name")

ggsave(filename = file.path(figures.dir, str_c(pop.type, "_resp_sel.jpg")),
       height = 5, width = 9)



### Figure 5 will have the genomic relationship, persistence of phase, inbreeding,
### and number of QTL fixed for an allele

# Persistence of phase

df1.persistence <- df.list$qtl_marker_LD %>%
  filter(extra1 == "persistence") %>%
  sim.summarize() %>%
  mutate(variable = "Persistence of LD Phase") %>%
  filter(heritability == "0.5") %>%
  select(-heritability)


### Relationship

df1.relatioship <- df.list$rel %>%
  sim.summarize() %>%
  mutate(variable = "Average Relationship\n(Scaled to Base Population)") %>%
  filter(heritability == "0.5") %>%
  select(-heritability)


# Inbreeding

df1.inbreeding <- df.list$inbreeding %>%
  sim.summarize() %>%
  mutate(variable = "Inbreeding\n(Scaled to Base Population)") %>%
  filter(heritability == "0.5") %>%
  select(-heritability)

# Fixed QTL

df1.qtlfreq <- df.list$fixedqtl %>%
  ungroup() %>%
  mutate(value = value / 100) %>%
  sim.summarize() %>%
  mutate(variable = "Proportion of Fixed QTL") %>%
  filter(heritability == "0.5") %>%
  select(-heritability)

# Fixed markers
# df1.marfreq <- df.list$fixedmar %>%
#   ungroup() %>%
#   mutate(value = value / 1490) %>%
#   sim.summarize() %>%
#   mutate(variable = "Proportion of Fixed Markers") %>%
#   filter(heritability == "0.5") %>%
#   select(-heritability)

# df1.qtlfreq <- df.list$fixedqtl %>%
#   ungroup() %>%
#   sim.summarize() %>%
#   mutate(variable = "Number of Fixed QTL") %>%
#   filter(heritability == "0.5") %>%
#   select(-heritability)


df1 <- bind_rows(df1.persistence, df1.relatioship, df1.inbreeding, df1.qtlfreq) %>%
  mutate(variable = factor(variable, 
                           levels = c("Average Relationship\n(Scaled to Base Population)", 
                                      "Persistence of LD Phase",
                                      "Inbreeding\n(Scaled to Base Population)", 
                                      "Proportion of Fixed QTL")))

gp.combined <- sim.ggplot(df.summary = df1, main = "", ylab = "", 
                          col.factors = tp.change.factors, text.y.scaling = 1.01)

ggsave(filename = file.path(figures.dir, "explanatory_vars.jpg"),
       plot = gp.combined, height = 15, width = 9)


