#' Determine if a locus is polymorphic
is.polymorphic <- function(x) {
  mean(x) %>%
    abs() != 1
} # Close the function