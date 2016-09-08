#' Round limit
#' 
#' @description 
#' A convenience function to round a number down or up if it is outside a 
#' lower or upper limit
#' 
#' 
round.limit <- function(x, limit, option) {
  if (option == "upper") {
    if (x > limit) {
      return(limit)
    } else {return(x) }
  }
  if (option == "lower") {
    if (x < limit) {
      return(limit)
    } else {return(x) }
  }
} # Close