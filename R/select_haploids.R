#' Subset the haploid genotypes matrix
#' 
#' @param haploid.genos A 2n x m matrix of haploid alleles at m loci for n
#' individuals.
#' @param line.names Individuals to subset.
#' 
#' @import dplyr
#' @import stringr
#' 
#' @export 
#' 
select.haploids <- function(haploid.genos, line.names) {
  
  # Deal with inputs
  if (!is.list(haploid.genos)) {
    haploid.mat <- as.matrix(haploid.genos)
  } else {
    haploid.mat <- do.call("rbind", haploid.genos)
  }
  
  line.names <- as.character(line.names)
  
  # Extract the names of the lines in the haploid matrix
  haploid.lines <- row.names(haploid.mat) %>%
    # Remove gamete indicator
    str_replace(pattern = "\\.[0-9]$", replacement = "")
  
  # Detect the desired line names in the haploid line names
  haploid.subset <- haploid.mat %>%
    subset.matrix(subset = haploid.lines %in% line.names)
  
  # Return the subset
  return(haploid.subset)
  
} # Close the function
