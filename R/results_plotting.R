#' Summarize a results data frame
#' 
#' @param df A data frame of results that includes values of a variable
#' and summary metadata (i.e. the cycle, change, iteration, etc.)
#' 
#' @export
#' 
sim.summarize <- function(df) {
  
  # Pull out colnames
  df.groups <- names(df) %>% .[. != "value"]
  
  # Remove duplicated change-iter-cycle pairs
  df1 <- df %>% 
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
#' @param ylim The y-axis scale range.
#' @param xlim The x-axis scale range. Defaults to \code{c(1, (n.cycles + 1))}.
#' @param col.factors The factors used in coloring the points / lines on the plot.
#' @param print.plot Logical. Should the plot be displayed?
#' 
#' @export
#' 
sim.ggplot <- function(df.summary, main, ylab, xlab = "Breeding Cycle", ylim = NULL, 
                       xlim = c(1, n.cycles + 1), col.factors, print.plot = TRUE) {
  
  # Set the range in ylim
  if (is.null(ylim)) {
    ylim <- df.summary %>% 
      summarize(max = max(mean + ci, na.rm = T), min = min(mean - ci, na.rm = T)) %>% 
      summarize(max = max(max, na.rm = T), min = min(min, na.rm = T)) %>% 
      pretty() %>% 
      range()
    
  }
  
  gp <- ggplot(data = df.summary, aes(x = cycle.offset, y = mean, col = change, 
                                      shape = change)) +
    geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), col = "black", width = 0.10) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    ggtitle(main) + 
    ylab(ylab) +
    xlab(xlab) +
    ylim(ylim) +
    xlim(xlim) +
    scale_color_discrete(name = "Update Method",
                         labels = as.character(col.factors)) +
    scale_shape_discrete(name = "Update Method",
                         labels = as.character(col.factors)) +
    facet_grid(~ exp_name)
  
  # Print the plot
  if (print.plot) print(gp)
  
  # Return the plot
  return(gp)
  
}