#' Summarize a results data frame
#' 
#' @param df A data frame of results that includes values of a variable
#' and summary metadata (i.e. the cycle, change, iteration, etc.)
#' 
#' @import dplyr
#' 
#' @export
#' 
sim.summarize <- function(df) {
  
  # Pull out colnames
  df.groups <- names(df) %>% .[. != "value"]
  
  # Make sure the variables are in correct format
  df1 <- df %>% mutate(iter = as.integer(iter),
                       cycle = as.integer(cycle),
                       change = as.character(change))
  
  # Remove duplicated change-iter-cycle pairs
  df1 <- df1 %>% 
    group_by_(.dots = df.groups) %>% 
    filter(row_number() == 1)
  
  df.groups <- names(df1) %>% .[!. %in% c("value", "iter")]
  
  df2 <- df1 %>%
    group_by_(.dots = df.groups) %>%
    summarize(mean = mean(value, na.rm = T), 
              sd = sd(value, na.rm = T),
              n = n()) %>%
    mutate(se = sd / sqrt(n),
           ci = qt(p = 1 - (0.05 / 2), df = n - 1) * se)
  
  # Add in an offset to the cycles
  df2$cycle.offset <- df2$cycle + (as.numeric(as.factor(df2$change)) - 1) * 0.1
  
  # Ungroup
  df2 %>%
    ungroup()
}

#' Plot data from a data frame summary
#' 
#' @param df.summary The summary data frame produced by the \code{sim.summarize}
#' function.
#' @param main The plot title.
#' @param ylab The y-axis title.
#' @param xlab The x-axis title. Defaults to "Breeding Cycle."
#' @param col.factors The factors used in coloring the points / lines on the plot.
#' @param text.y.scaling The scale by which to increase the maximum y value in
#' the data in order to add annotations. Defaults to 1.05 (or 5 per-cent).
#' @param print.plot Logical. Should the plot be displayed?
#' 
#' @export
#' 
sim.ggplot <- function(df.summary, main, ylab, xlab = "Breeding Cycle", 
                       col.factors, text.y.scaling = 1.05, print.plot = TRUE,
                       facet.vars = c("exp_name")) {
  
  # Designate labels for the individual plots within facets
  n.facets <- df.summary %>% 
    select_(.dots = c(facet.vars, "variable")) %>%
    distinct() %>% 
    nrow()
  
  gp <- ggplot(data = df.summary, aes(x = cycle.offset, y = mean, col = change, 
                                      shape = change)) +
    geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), col = "black", width = 0.10) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    ggtitle(main) + 
    ylab(ylab) +
    xlab(xlab) +
    # ylim(ylim) +
    # xlim(xlim) +
    scale_color_discrete(name = "Update Method",
                         labels = as.character(col.factors)) +
    scale_shape_discrete(name = "Update Method",
                         labels = as.character(col.factors)) +
    scale_x_continuous(breaks = seq(1, n.cycles + 1, 3))
  
  # Determine how to map the facets
  if (n.facets <= 2) {
    gp1 <- gp + facet_grid(as.formula(c("~", str_c(facet.vars, collapse = " + "))))
  }
  
  if (n.facets > 2) {
    gp1 <- gp + facet_grid(as.formula(c("variable ~", str_c(facet.vars, collapse = " + "))), 
                           scales = "free_y", switch = "y")
  }
  
  gp.build <- ggplot_build(gp1)
  
  # Find the y-ranges
  facet.df <- sapply(X = gp.build$panel$ranges, FUN = function(rangei) 
    rangei$y.range %>% max() ) %>%
    # convert to df
    tbl_df() %>%
    # Add 10%
    mutate(value = value * text.y.scaling)
  
  facet.df1 <- cbind(facet.df, 
                     gp.build$panel$layout %>% select_(.dots = c(facet.vars, "variable")))
  
  # Determine annotation locations
  facet.df2 <- facet.df1 %>%
    mutate(x = 1, label = LETTERS[seq_len(n.facets)]) %>%
    rename(y = value)
  
  # Modify fonts
  gp2 <- gp1 + theme(
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold")
  ) + 
    geom_text(data = facet.df2, aes(x = x, y = y, label = label, fontface = 2), 
              inherit.aes = FALSE)
    
  # Print the plot
  if (print.plot) print(gp2)
  
  # Return the plot
  return(gp2)
  
}
