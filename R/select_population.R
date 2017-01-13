#' Apply selection on a population
#' 
#' @param value.mat A matrix of values (phenotypic, genotypic, or breeding) on
#' which to apply selection.
#' @param sel.intensity Either the proportion to select (i.e. from 0 to 1) or
#' the number of individuals to select.
#' @param selection The method of selection. Can be \code{"best"} (default),
#' \code{"worst"}, \code{"random"}, or \code{"tails"}.
#' @param exlcusion A vector of line names to exclude when making selection.
#' 
#' @details 
#' The \code{"best"} method selects the top \code{sel.intensity} entries based on
#' the values in the \code{value.mat}, the \code{"worst"} method selects the
#' bottom entries, the \code{"random"} method select entries randomly, and the
#' \code{"tails"} method select the top \eqn{sel.intensity/2} entries
#' and the bottom \eqn{sel.intensity/2} entries.
#' 
#' @import dplyr
#' 
#' @export
#' 
select.population <- function(value.mat, sel.intensity, selection = "best",
                              exclusion = NULL) {
  
  # Exclude based on the provided exclusion vector
  value.mat <- value.mat[!row.names(value.mat) %in% exclusion,] %>%
    as.matrix()
  
  # Number in the pheno.vec
  N.tot <- nrow(value.mat)
  
  # If the sel.intensity is between 0 and 1, find the total number to keep
  if (all(sel.intensity > 0, sel.intensity < 1)) {
    # Number to select
    N.sel <- sel.intensity * N.tot
  } else {
    N.sel <- sel.intensity
  }
  
  if (selection == "best") {
    # Sort the values
    value.mat.ordered <- value.mat[order(value.mat, decreasing = T),] %>%
      as.matrix()
    value.sel <- value.mat.ordered[seq(N.sel),]
  }
  
  # If using the "worst" method, select the bottom N.sel
  if (selection == "worst") {
    # Sort the values
    value.mat.ordered <- value.mat[order(value.mat, decreasing = F),] %>%
      as.matrix()
    value.sel <- value.mat.ordered[seq(N.sel),]
  }
  
  # If using "random", randomly select the lines
  if (selection == "random") {
    value.sel <- value.mat[sample(seq(N.tot), size = N.sel),]
  }
  
  # If using "tails", select the top and bottom
  if (selection == "tails") {
    # Top
    value.mat.ordered <- value.mat[order(value.mat, decreasing = T),] %>%
      as.matrix()
    value.sel.top <- value.mat.ordered[seq(N.sel/2),]
    
    value.mat.ordered <- value.mat[order(value.mat, decreasing = F),] %>%
      as.matrix()
    value.sel.bottom <- value.mat.ordered[seq(N.sel/2),]
    
    value.sel <- c(value.sel.top, value.sel.bottom)
    
  }
  
  # Sort the selections based on name
  value.sel <- value.sel[order(names(value.sel))]
  # Retrieve the names
  line.names <- names(value.sel)
  
  # Create output list
  output.list <- list(value.sel = as.matrix(value.sel), lines.sel = line.names)
  
  # Return the selections
  return(output.list)
} # Close the function
