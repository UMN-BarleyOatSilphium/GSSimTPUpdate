#' Hidden functions
#' 
#' @description 
#' Attempts to find the inverse of a matrix. If unsuccessful, find the
#' generalized inverse.
#' 
#' @param x A matrix to be inverted.
#' @param silent Logical whether to print the method of inversion used.
#' 
#' @export
#' 
#' @importFrom MASS ginv
#' 
invert.mat <- function(x, silent = FALSE) {
  
  if (try(expr = solve(x), silent = T) %>% class == "try-error") {
    x.inv <- ginv(x)
    method <- "Generalized"
  } else {
    # Otherwise use solve
    x.inv <- solve(x)
    method <- "Inverse"
  }
  
  if (!silent)  cat("\nMethod: ", method)
  return(x.inv)
  
} # Close the function

#' 
#' Determine if a locus is polymorphic
#' 
#' @import dplyr
#' 
is.polymorphic <- function(x) {
  mean(x) %>%
    abs() != 1
} # Close the function

#' 
#' Make a matrix of contrasts
#' 
#' @description 
#' Creates a matrix of contrasts between the individuals not in the training 
#' population (i.e. the candidates) and the mean of the whole population.
#' 
#' @param unphenotyped.index The index of the unphenotyped entries in the 
#' relationship matrix (A)
#' @param n.total The size of the whole population
#' 
#' @export
#' 
#' @references 
#' Rincent, R., Laloe, D., Nicolas, S., Altmann, T., Brunel, D., Revilla, P., 
#' Moreau, L. (2012). Maximizing the Reliability of Genomic Selection by 
#' Optimizing the Calibration Set of Reference Individuals: Comparison of 
#' Methods in Two Diverse Groups of Maize Inbreds (Zea mays L.). Genetics, 
#' 192(2), 715–728. http://doi.org/10.1534/genetics.112.141473
#' 
make.contrast <- function(unphenotyped.index, n.total) {
  
  # Number of unphenotyped individuals
  n.unphenotyped <- length(unphenotyped.index)
  
  # The positive value of the contrast
  pos.contrast = 1 - (1 / n.total)
  # The negative value of the contrast
  neg.contrast = -1 / n.total
  
  # Make the total constrast matrix
  contrast.mat <- matrix(data = neg.contrast, nrow = n.total, ncol = n.unphenotyped)
  
  # Fill in the positive values
  coord <- cbind(unphenotyped.index,
                 seq(n.unphenotyped))
  
  contrast.mat[coord] <- pos.contrast
  
  return(contrast.mat)
}

#' 
#' Finds the whole genome positions of QTL and markers
#' 
#' @param genome The list of hypred genomes.
#' 
#' @import dplyr
#' 
#' @export
#' 
find.pos <- function(genome) {
  
  # Return all the loci (this is easy - just the index of the number of columns)
  pos.loci <- sapply(genome, function(chr) 
    length(chr@pos.snp)) %>% 
    sum() %>%
    seq()
  
  # Number of chromosomes
  n.chr <- length(genome)
  
  # Running count of the number of loci up to the nth chromosome
  n.loci.total <- 0
  
  # Empty lists of qtl and snp positions
  pos.qtl <- list()
  
  # Iterate over chromosomes to find the position of markers and qtl
  for (i in seq(n.chr)) {
    
    # Extract the position of the qtl and add the loci total
    pos.qtl[[i]] <- genome[[i]]@pos.add.qtl$ID + n.loci.total
    
    # Total number of loci up to the chromosome
    n.loci.total = n.loci.total + genome[[i]]@num.snp.chr
    
  }
  
  # Unlist the qtl position list
  pos.qtl <- pos.qtl %>%
    unlist()
  
  # Position of markers
  pos.snp <- setdiff(pos.loci, pos.qtl)
  
  # Create list
  out.list <- list(pos.loci = pos.loci,
                   pos.qtl = pos.qtl,
                   pos.snp = pos.snp)
  
  # Return
  return(out.list)
  
} # Close the function

#' 
#' Round with a limit
#' 
#' @description 
#' A convenience function to round a number down or up if it is outside a 
#' lower or upper limit
#' 
#' @param x A vector of numbers to assess
#' @param limit The upper or lower limit for rounding
#' @param option Character vector of either \code{"upper"} or \code{"lower"}
#' corresponding to whether the \code{limit} is an upper limit or a lower limit.
#'
#' 
#' @export
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

#' 
#' M matrix
#' 
#' @description 
#' Creates an M matrix, or the orthagonal projector
#' 
#' @param n The number of rows of matrix X
#' 
#' @export
#' 
#' @references 
#' Rincent, R., Laloe, D., Nicolas, S., Altmann, T., Brunel, D., Revilla, P., 
#' Moreau, L. (2012). Maximizing the Reliability of Genomic Selection by 
#' Optimizing the Calibration Set of Reference Individuals: Comparison of 
#' Methods in Two Diverse Groups of Maize Inbreds (Zea mays L.). Genetics, 
#' 192(2), 715–728. http://doi.org/10.1534/genetics.112.141473
#' 
design.M <- function(n) {
  
  # Identity matrix of size n
  I <- diag(n) 
  # X matrix
  X <- matrix(1, n, 1)
  # Create M and return
  I - tcrossprod(X %*% solve(crossprod(X, X)), X)
  
} # Close the function

#' 
#' Parse results to data.frame
#' 
#' @description 
#' Convert a vector of values to a df during data parsing
#' 
#' @param x A vector of values
#' @param n.iters The number of simulation iterations
#' @param change The method of updating the training population
#' 
#' @import stringr
#' @import dplyr
#' @import tidyr
#' 
#' @export
#' 
nv_df <- function(x, change) {
  
  # Convert x to a tibble
  x.tbl <- tbl_df(x) %>%
    mutate(obs = names(x))
  
  # New columns names for separation
  # Extract the first name as an example
  name.eg <- names(x[1])
  # Find the number of period-separated elements
  n.cols <- str_count(string = name.eg, pattern = "\\.") + 1
  
  # Create a vector of new column names
  if (n.cols <= 3) {
    new.cols <- c("set", "rep", "cycle")
  } else {
    new.cols <- c("set", "rep", "cycle", str_c("extra", seq_len(n.cols - 3)))
  }
  
  # Create new columns based on the names of the values
  x.tbl1 <- x.tbl %>% 
    separate(obs, new.cols, sep = "\\.") %>%
    unite(iter, set, rep) %>% 
    mutate(iter = iter %>% as.factor() %>% as.numeric(),
           cycle = cycle %>% str_replace("cycle", ""))
  
  # Rearrange columns and add the change
  x.tbl2 <- x.tbl1 %>%
    mutate(change = change) %>%
    select(c((seq_len(n.cols) + 1), 1))
  
  return(x.tbl2)
  
} # Close the function

#'
#' Retrieve QTL effects from a genome
#' 
#' @param genome A list of hypred genomes.
#' 
#' @export
#' 
qtl.eff <- function(genome) {
  
  # Make sure the genome is a list
  if (!is.list(genome)) stop("genome must be a list.")
  
  # Iterate over chromosomes
  eff.list <- lapply(X = genome, FUN = function(chr) chr@add.and.dom.eff)
  
  # Separate into addive and dominance
  add.eff <- sapply(eff.list, FUN = function(i) i$add) %>%
    unlist()
  dom.eff <- sapply(eff.list, FUN = function(i) i$dom) %>%
    unlist()
  
  # Return
  list(add = add.eff, dom = dom.eff)
  
} # Close the function
    

