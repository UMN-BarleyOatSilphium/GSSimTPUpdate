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