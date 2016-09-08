#' Create a crossing block from known parents
#' 
#' @description 
#' Create a crossing block based on a vector of line names for each parents 
#' (this can be the same vector) and a number of crosses, while choosing whether 
#' or not to use parents once.
#' 
#' @param parent1.lines A character vector of lines for the first parents.
#' @param parent2.lines See \code{parent1.lines}.
#' @param n.crosses The number of crosses to make.
#' @param method The mating scheme. Can be \code{"random"} (default) or 
#' \code{"all.pairwise"}. If \code{"all.pairwise"}, parents will be used
#' more than once.
#' @param use.parents.once A logical whether to only use parents once.
#' 
#' @import dplyr
#' 
#' @export
#' 
make.crossing.block <- function(parent1.lines, parent2.lines, n.crosses, 
                                method = "random", use.parents.once = FALSE) {
  
  # Find the intersection of parent1.lines and parent2.lines
  par.intersect <- intersect(parent1.lines, parent2.lines)
  n.intersect <- length(par.intersect)
  # Find the number of parent1 and parent2 lines
  n.par1 <- length(parent1.lines)
  n.par2 <- length(parent2.lines)
  
  if ( all(n.intersect == n.par1, n.intersect == n.par2) ) {
    n.possible.crosses <- (n.par1 * (n.par2 - 1)) / 2
    same.lines = T
  } else {
    n.possible.crosses <- n.par1 * n.par2
    same.lines = F
  }
  
  # Create all pairwise crosses
  sample.crosses <- expand.grid(parent1.lines, parent2.lines) %>%
    select(Var2, Var1)
  
  # Remove selfs
  sample.crosses1 <- sample.crosses %>%
    filter(apply(X = sample.crosses, MARGIN = 1, FUN = function(cross) length(unique(cross)) > 1))
  
  # If statements for methods
  if (method == "all.pairwise") {
    
    return(sample.crosses1)
    
  } else {
  
    if (n.crosses > n.possible.crosses) stop("n.crosses is more than possible.")
    
    # If parents should only be used once, sample without replacement all lines 
    # and put them into a matrix
    
    if (use.parents.once) {
      
      # Sample into a matrix
      if (same.lines) {
        
        # Create random crosses from the parent lines
        random.crosses <- sample(parent1.lines) %>%
          matrix(ncol = 2, byrow = T) %>%
          as.data.frame()
        
        # Sample n.crosses from the random crosses
        selected.crosses <- random.crosses %>%
          sample_n(size = n.crosses)

        } else {
          
        # Are there any overlapping lines in the two parent vectors?
        if (n.intersect == 0) {
          
          # If not just sample each vector into a matrix
          random.crosses <- cbind( sample(parent1.lines), sample(parent2.lines) ) %>%
            as.data.frame()
          
          # Sample n.crosses from the random crosses
          selected.crosses <- random.crosses %>%
            sample_n(size = n.crosses)
          
        } else {
          stop("Function not available.")
        }
      }
      
    } else {
      
      # Remove reciprocal
      sample.crosses2 <- sample.crosses1 %>%
        filter(!apply(X = sample.crosses1, MARGIN = 1, FUN = sort) %>% 
                 t() %>% 
                 duplicated() )
      
      # Sample n.crosses from the sample crosses
      selected.crosses <- sample.crosses2 %>%
        sample_n(size = n.crosses)
    
    }
    
    # Rename
    colnames(selected.crosses) <- c("Parent1", "Parent2")
    row.names(selected.crosses) <- NULL
    return(selected.crosses)
  }
  
} # Close the function