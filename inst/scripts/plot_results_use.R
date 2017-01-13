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

df1.acc <- df1 <- plotting_data$accuracy %>%
  sim.summarize() %>%
  mutate(variable = "Prediction Accuracy") %>%
  filter(heritability == "0.5") %>%
  select(-heritability)
  
# Designate labels for the individual plots within facets
n.facets <- df1 %>% 
  select(exp_name, variable) %>%
  distinct() %>% 
  nrow()

gp <- df1 %>%
  ggplot(aes(x = cycle.offset, y = mean, col = change, shape = change)) +
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), col = "black", width = 0.10) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  ylab(expression(Prediction~Accuracy~(italic(r)[MG]))) +
  xlab("Breeding Cycle") +
  scale_color_discrete(name = "Update Method",
                       labels = as.character(tp.change.factors)) +
  scale_shape_discrete(name = "Update Method",
                       labels = as.character(tp.change.factors)) +
  scale_x_continuous(breaks = seq(1, n.cycles + 1, 3)) +
  # Faceting
  facet_wrap("exp_name")

# Geom text data
gt.df <- ggplot_build(gp)$plot$data %>%
  group_by(variable, exp_name) %>%
  summarize(ymax = max(mean + ci)) %>%
  ungroup() %>%
  mutate(
    y = max(ymax) * 1.05,
    x = 1,
    label = LETTERS[seq_len(n.facets)])
         

# Modify fonts
gp1 <- gp + 
  theme(
    strip.text.x = element_text(face = "bold", size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 14),
    legend.position = "bottom") +
  geom_text(data = gt.df, aes(x = x, y = y, label = label, fontface = 2), 
            inherit.aes = FALSE, size = 5)

ggsave(plot = gp1, filename = "pred_acc.jpg", height = 6, width = 10.5)



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

# Designate labels for the individual plots within facets
n.facets <- df1 %>% 
  select(exp_name, variable) %>%
  distinct() %>% 
  nrow()

gp <- df1 %>%
  ggplot(aes(x = cycle.offset, y = mean, col = change, shape = change)) +
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), col = "black", width = 0.10) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  ylab("") +
  xlab("Breeding Cycle") +
  scale_color_discrete(name = "Update Method",
                       labels = as.character(tp.change.factors)) +
  scale_shape_discrete(name = "Update Method",
                       labels = as.character(tp.change.factors)) +
  scale_x_continuous(breaks = seq(1, n.cycles + 1, 3)) +
  # Faceting
  facet_grid(variable ~ exp_name, scale = "free_y", switch = "y")

# Geom text data
gt.df <- ggplot_build(gp)$plot$data %>%
  group_by(variable, exp_name) %>%
  summarize(ymax = max(mean + ci)) %>%
  mutate(y = max(ymax) * 1.05,
         x = 1) %>%
  ungroup() %>%
  mutate(label = LETTERS[seq_len(n.facets)])


# Modify fonts
gp1 <- gp + 
  theme(
    strip.text = element_text(face = "bold", size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 14),
    legend.position = "bottom") +
  # Add labels
  geom_text(data = gt.df, aes(x = x, y = y, label = label, fontface = 2), 
            inherit.aes = FALSE, size = 5)

ggsave(plot = gp1, filename = "gen_val_gen_var.jpg", height = 10.5, width = 10.5)






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


# For manuscript
df1 <- bind_rows(df1.persistence, df1.relatioship, df1.inbreeding, df1.qtlfreq) %>%
  mutate(variable = factor(variable, 
                           levels = c("Average Relationship\n(Scaled to Base Population)", 
                                      "Persistence of LD Phase",
                                      "Inbreeding\n(Scaled to Base Population)", 
                                      "Proportion of Fixed QTL")))

# Designate labels for the individual plots within facets
n.facets <- df1 %>% 
  select(exp_name, variable) %>%
  distinct() %>% 
  nrow()

gp <- df1 %>%
  ggplot(aes(x = cycle.offset, y = mean, col = change, shape = change)) +
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), col = "black", width = 0.10) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  ylab("") +
  xlab("Breeding Cycle") +
  scale_color_discrete(name = "Update Method",
                       labels = as.character(tp.change.factors)) +
  scale_shape_discrete(name = "Update Method",
                       labels = as.character(tp.change.factors)) +
  scale_x_continuous(breaks = seq(1, n.cycles + 1, 3)) +
  # Faceting
  facet_grid(variable ~ exp_name, scale = "free_y", switch = "y")

# Geom text data
gt.df <- ggplot_build(gp)$plot$data %>%
  group_by(variable, exp_name) %>%
  summarize(ymax = max(mean + ci)) %>%
  mutate(y = max(ymax) * 1.05,
         x = 1) %>%
  ungroup() %>%
  mutate(label = LETTERS[seq_len(n.facets)])


# Modify fonts
gp1 <- gp + 
  theme(
    strip.text = element_text(face = "bold", size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 14),
    legend.position = "bottom") +
  # Add labels
  geom_text(data = gt.df, aes(x = x, y = y, label = label, fontface = 2), 
            inherit.aes = FALSE, size = 5)

ggsave(plot = gp1, filename = "explan_vars.jpg", height = 13.5, width = 10.5)
