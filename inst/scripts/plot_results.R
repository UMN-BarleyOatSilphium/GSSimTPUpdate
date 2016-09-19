## Plot the parsed results of the simulation

# Load packages
library(tidyverse, quietly = T)
library(stringr, quietly = T)
library(GSsim.TPUpdate, quietly = T, verbose = F)


# Set directory and grab the files
results.dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Side Projects/Simulations/GSsim.TPUpdate/Results/Base Experiment/"

# File path for figures
figures.dir <- file.path(results.dir, "Figures")

# Load data from the allele frequency experiment
all.files <- list.files(results.dir, full.names = T) %>%
  str_subset(pattern = "simulation_results") %>%
  str_subset(pattern = "collective")

# Just MNxND files
filename <- all.files %>%
  str_subset(pattern = "MNxND") %>%
  .[1]

load(filename)

### Initial experimen settings
# Experimental parameters
pop.makeup <- str_extract(string = filename, pattern = 'results_[A-Za-z]{2,5}') %>%
  str_extract(pattern = '[A-Za-z]{2,5}$')

tp.formation <- 
  str_extract(string = filename, pattern = paste0(c("cumulative", "window"), collapse = "|"))

n.cycles <- collective.abbreviated.results[[1]][[1]] %>% 
  select(cycle) %>% 
  unique() %>% 
  nrow()

n.iter = collective.abbreviated.results[[1]][[1]] %>% 
  select(iter) %>% 
  unique() %>% 
  nrow()

# The different updating methods
collective.abbreviated.results %>%
  names()

tp.change.factors <- c(best = "Best", CDmean = "CDmean", nochange = "No Change", 
                       PEVmean = "PEVmean", random = "Random", worst = "Worst") %>%
  as.factor()


## Define a function for summarization and calculation of mean and confidence
## interval
# df - a data.frame of the update method, cycle number, and parameter value.
sim.summarize <- function(df) {
  
  df %>%
    group_by(change, cycle) %>%
    summarize(mean = mean(value), 
              sd = sd(value),
              n = n()) %>%
    mutate(se = sd / sqrt(n),
           ci = qt(p = 1 - (0.05 / 2), df = n - 1) * se)
  
}

## Define a function to plot the summary data.frame using ggplot2

sim.ggplot <- function(df.summary, main, ylab, xlab = "Breeding Cycle", ylim, 
                       xlim = c(1, n.cycles), tp.change.factors) {
  
  gp <- ggplot(data = df.summary, aes(x = cycle, y = mean, col = change, shape = change)) +
    geom_point(size = 2) +
    geom_line(size = 1) +
    geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), col = "black", width = 0.10) +
    ggtitle(main) + 
    ylab(ylab) +
    xlab(xlab) +
    ylim(ylim) +
    xlim(xlim) +
    scale_color_discrete(name = "Update Method",
                         labels = as.character(tp.change.factors)) +
    scale_shape_discrete(name = "Update Method",
                         labels = as.character(tp.change.factors))
  
  print(gp)
  
}
  


### Change in Genetic Variance
df <- lapply(X = collective.abbreviated.results, FUN = function(tpc) tpc$candidate.gen.var) %>%
  bind_rows()

# Calculate mean and CI for each change-cycle combination over iterations
df1 <- sim.summarize(df)

sim.ggplot(df.summary = df1, main = "Genetic Variance", ylab = "Genetic Variance",
           ylim = c(0,8), tp.change.factors = tp.change.factors)

ggsave(filename = file.path(figures.dir, paste("figure_", pop.makeup, "_", tp.formation, 
                                               "_gen_var.jpg", sep = "")),
       height = 5, width = 6.5)


```

This is a little different from before, but follows the same trend. Not the most
convincing piece of evidence.



### Change in mean genotypic value of the selection candidates

``` {r}

df <- lapply(X = collective.abbreviated.results, FUN = function(tpc) tpc$candidate.gen.val) %>%
  bind_rows()

# Calculate mean and CI for each change-cycle combination over iterations
df1 <- df %>% 
  group_by(change, cycle) %>% 
  summarize(mean = mean(value), 
            sd = sd(value),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         ci = qt(p = 1 - (0.05 / 2), df = n - 1) * se)

ggplot(data = df1, aes(x = cycle, y = mean, col = change, shape = change)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), col = "black", width = 0.10) +
  ggtitle("Genotypic Value of the Selection Candidates") + 
  ylab("Genotypic Value") +
  xlab("Breeding Cycle") +
  ylim(c(0,25)) +
  scale_color_discrete(name = "Update Method",
                       labels = as.character(tp.change.factors)) +
  scale_shape_discrete(name = "Update Method",
                       labels = as.character(tp.change.factors))

ggsave(filename = file.path(figures.dir, paste("figure_", pop.makeup, "_", tp.formation, 
                                               "_geno_value.jpg", sep = "")),
       height = 5, width = 6.5)


```

This is an important figure to demonstrate practicality in a breeding program. 
All methods are pretty similar for the first few cycles, but using the "best"
method results in the highest genotypic values in the population after 4 cycles
and holds onto that lead for the remaining cycles.


### Change in the prediction accuracy over cycles

``` {r}
# Plot
df <- lapply(X = collective.abbreviated.results, FUN = function(tpc) tpc$validation.results) %>%
  bind_rows()

# Calculate mean and CI for each change-cycle combination over iterations
df1 <- df %>% 
  group_by(change, cycle) %>% 
  summarize(mean = mean(value), 
            sd = sd(value),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         ci = qt(p = 1 - (0.05 / 2), df = n - 1) * se)

ggplot(data = df1, aes(x = cycle, y = mean, col = change, shape = change)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), col = "black", width = 0.10) +
  ggtitle("Realized Prediction Accuracy") + 
  ylab("Prediction Accuracy (r)") +
  xlab("Breeding Cycle") +
  ylim(c(0,0.75)) +
  scale_color_discrete(name = "Update Method",
                       labels = as.character(tp.change.factors)) +
  scale_shape_discrete(name = "Update Method",
                       labels = as.character(tp.change.factors))

ggsave(filename = file.path(figures.dir, paste("figure_", pop.makeup, "_", tp.formation, 
                                               "_prediction_accuracy.jpg", sep = "")),
       height = 5, width = 6.5)


```

This is a heavy-hitter for the "why should I care" factor of this experiment. Clearly
not updating is unfavorable, but "best" maintains a nice advantage during the more
critical cycles. As usual all methods (except "no change") converge after 10
cycles or so. CDmean and PEVmean swap places, but generally do similarly.


## Potential Explanatory Variables


### Change in QTL-marker LD over cycles

``` {r}
# Each QTL's max LD
# Plot
df <- lapply(X = collective.abbreviated.results, FUN = function(tpc) tpc$qtl.marker.LD) %>%
  bind_rows()

# Calculate mean and CI for each change-cycle combination over iterations
df1 <- df %>% 
  filter(extra == "mean_max_genome") %>%
  group_by(change, cycle) %>% 
  summarize(mean = mean(value), 
            sd = sd(value),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         ci = qt(p = 1 - (0.05 / 2), df = n - 1) * se)

ggplot(data = df1, aes(x = cycle, y = mean, col = change, shape = change)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  ylim(c(0.6,0.9)) +
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), col = "black", width = 0.10) +
  ggtitle("Mean LD Between QTL and Marker in Highest LD") + 
  ylab("Linkage Disequilibrum (r)") +
  xlab("Breeding Cycle") +
  scale_color_discrete(name = "Update Method",
                       labels = as.character(tp.change.factors)) +
  scale_shape_discrete(name = "Update Method",
                       labels = as.character(tp.change.factors))

ggsave(filename = file.path(figures.dir, paste("figure_", pop.makeup, "_", tp.formation, 
                                               "_mean_max_ld.jpg", sep = "")),
       height = 5, width = 6.5)



## Persistence of phase

# Calculate mean and CI for each change-cycle combination over iterations
df1 <- df %>% 
  filter(extra == "persistence") %>%
  group_by(change, cycle) %>% 
  summarize(mean = mean(value), 
            sd = sd(value),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         ci = qt(p = 1 - (0.05 / 2), df = n - 1) * se)

ggplot(data = df1, aes(x = cycle, y = mean, col = change, shape = change)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), col = "black", width = 0.10) +
  ggtitle("Persistence of LD Phase Between Training Population\nand Selection Candidates") + 
  ylab("Persistence of Phase (cor of r)") +
  xlab("Breeding Cycle") +
  ylim(c(0.0,0.65)) +
  scale_color_discrete(name = "Update Method",
                       labels = as.character(tp.change.factors)) +
  scale_shape_discrete(name = "Update Method",
                       labels = as.character(tp.change.factors))

ggsave(filename = file.path(figures.dir, paste("figure_", pop.makeup, "_", tp.formation, 
                                               "_persistence_of_phase.jpg", sep = "")),
       height = 5, width = 6.5)


```

Here is another heavy-hitter. The trend aligns really well with what we saw in the
realized prediction accuracy, with "best" maintaining a good advantage and "no change"
dropping off quickly. 


### Change in the average relationship between TP and candidates

``` {r}

# Plot
df <- lapply(X = collective.abbreviated.results, FUN = function(tpc) tpc$relationship) %>%
  bind_rows()

# Calculate mean and CI for each change-cycle combination over iterations
df1 <- df %>% 
  group_by(change, cycle) %>% 
  summarize(mean = mean(value), 
            sd = sd(value),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         ci = qt(p = 1 - (0.05 / 2), df = n - 1) * se)

ggplot(data = df1, aes(x = cycle, y = mean, col = change, shape = change)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), col = "black", width = 0.10) +
  ggtitle("Additive Genomic Relationship Between Training Population\nand Selection Candidates") + 
  ylab("Additive Genomic Relationship\n(Scaled to Base Population)") +
  xlab("Breeding Cycle") +
  ylim(c(-0.1, 1.25)) +
  scale_color_discrete(name = "Update Method",
                       labels = as.character(tp.change.factors)) +
  scale_shape_discrete(name = "Update Method",
                       labels = as.character(tp.change.factors))

ggsave(filename = file.path(figures.dir, paste("figure_", pop.makeup, "_", tp.formation, 
                                               "_average_relationship.jpg", sep = "")),
       height = 5, width = 6.5)

```

Another important figure. Adding the "best" results in a higher average relationship
between the TP and the selection candidates (compared with other methods). As we
might expect, PEVmean does slightly better CDmean. It might be good to look at the
average relationship between TP additions, especially to confirm that the algorithms
are working as expected.


### Change in the Expected Heterozygosity Among TP Additions

```{r}
# Plot
df <- lapply(X = collective.abbreviated.results, FUN = function(tpc) tpc$tp.update.exp.het) %>%
  bind_rows()

# Calculate mean and CI for each change-cycle combination over iterations
df1 <- df %>% 
  group_by(change, cycle) %>% 
  summarize(mean = mean(value), 
            sd = sd(value),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         ci = qt(p = 1 - (0.05 / 2), df = n - 1) * se)

ggplot(data = df1, aes(x = cycle, y = mean, col = change, shape = change)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), col = "black", width = 0.10) +
  ggtitle("Expected Heterozygosity of Training Population Additions") + 
  ylab("Expected Heterozygosity") +
  xlab("Breeding Cycle") +
  ylim(c(0,0.25)) +
  scale_color_discrete(name = "Update Method",
                       labels = as.character(tp.change.factors)) +
  scale_shape_discrete(name = "Update Method",
                       labels = as.character(tp.change.factors))

ggsave(filename = file.path(figures.dir, paste("figure_", pop.makeup, "_", tp.formation, 
                                               "_exp_het.jpg", sep = "")),
       height = 5, width = 6.5)

```

This is less impressive in terms of differences between the methods, but it at 
least confirms that the CDmean algorithm is working as expected. Of course, no
data is available for the "no change" method.





``` {r, eval = F, echo = F}

# Save data lists to a file
# Gather the lists
data.lists <- neyhart::find('\\.list$', class = "list")

# Create a new filename
save.file <- filename %>%
  str_replace(pattern = "simulation_results", replacement = "plot_data")

# Save the data
save(list = data.lists, file = save.file)
```