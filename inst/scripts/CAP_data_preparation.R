## Script to take genotype data from T3 and prepare it for use in the simulations
# To accompany the GSSimTPUpdate project

# Load the libraries
library(GSSimTPUpdate)
library(dplyr)

# Define a function to take a matrix and compare all columns and return an index of the identical columns
identical.columns <- function(input.matrix) {
  
  # Error handling
  input.matrix <- as.matrix(input.matrix)
  
  # Save the vector 1:ncol(input.matrix)
  col.ind <- 1:ncol(input.matrix)
  
  # Create a comparison list
  comparison.list <- list()
  
  # Loop over the number of columns
  for (i in col.ind) {
    
    # Extract column i
    col.i <- input.matrix[,i]
    # Compare it to all other columns
    compare.i <- apply(X = as.matrix(input.matrix[,setdiff(col.ind, i)]), MARGIN = 2, FUN = function(column) identical(col.i, column) )
    # Find the columns that match
    col.compare <- setdiff(col.ind, i)[compare.i]
    
    # Add to the list
    comparison.list[[i]] <- sort(c(i, col.compare))
  }
  
  # Find the unique groups
  comparison.list <- unique(comparison.list)
  
  # Return the comparison list
  return(comparison.list)
} # Close the function

round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}

# Load the marker data
## This marker data was downloaded from T3 with the following filters:
## Minimum MAF: 0
## Max marker missingness: 0.1
## Max entry missingness: 0.1
data("CAP.hmp")

# Extract the marker data, including marker name, chromosome, alleles, and position
marker.info <- CAP.hmp[,c(1:4)]

# Create a n x m genotype matrix for filtering
CAP.M <- t(CAP.hmp[,-c(1:4)])
# Set column names to marker names
colnames(CAP.M) <- marker.info$rs


# Extract line names
line.names <- row.names(CAP.M)
# Separate MN from ND
ND.lines <- line.names[grep(pattern = "^ND", x = line.names)]
MN.lines <- setdiff(line.names, ND.lines)

# Split the M matrix by breeding program
MN.M <- CAP.M[MN.lines,]
ND.M <- CAP.M[ND.lines,]

# Remove snps that are monomorphic in BOTH the MN marker matrix and the ND
## marker matrix
MN.M.poly <- which(abs(colMeans(MN.M, na.rm = T)) != 1)
ND.M.poly <- which(abs(colMeans(ND.M, na.rm = T)) != 1)

# Find the common polymorphic snps
M.poly <- intersect(MN.M.poly, ND.M.poly)
CAP.M <- CAP.M[,M.poly]

marker.info <- marker.info[M.poly,]
  
  
### Processing of the marker matrix for use in simulation
# Remove redundant markers
# These are characterized by having the same genotypes across all samples AND fall on the same cM position
marker.info.known <- marker.info[marker.info$chrom != "UNK",]
# Find all unique chromosome-cM positions
unique.positions <- unique(marker.info.known[,c(3,4)])
  
# Apply a function over the unique positions
non.redundant.marker.list <- apply(X = unique.positions, MARGIN = 1, FUN = function(uniq) {
  
  # For each position, find all markers with that position
  markers.at.uniq <- as.matrix(marker.info.known$chrom) %in% as.matrix(uniq[1]) & as.matrix(marker.info.known$pos) %in% as.numeric(uniq[2])
  
  # If the number of markers is 1, just return the marker
  if (sum(markers.at.uniq) == 1) {
    return(marker.info.known[markers.at.uniq,])
    
  } else { # Otherwise look more closely
    
    # Extract the marker names
    marker.names <- as.character(marker.info.known[markers.at.uniq,1])
    M.i <- CAP.M[,marker.names]
    # Determine if the genotype calls are the same between markers
    same.genos <- identical.columns(input.matrix = M.i)
    # Choose the first index from each group
    # Gather the names of the chosen markers
    non.redundant.marker.names <- colnames(M.i)[sapply(X = same.genos, FUN = function(group) group[1])]
    
    # Extract marker info for those markers
    marker.info.non.redundant <- marker.info.known[marker.info.known$rs %in% non.redundant.marker.names,]
    marker.info.non.redundant1 <- marker.info.non.redundant
    # Add or subtract a deviation for the marker position
    # A deviation of 0.01 cM will be sufficient because the minimum positive difference is 0.05
    # Determine the deviation to add/substract
    deviation <- as.vector(round2(scale(x = 1:nrow(marker.info.non.redundant), scale = F), 0)) * 10
    marker.info.non.redundant1$pos <- marker.info.non.redundant$pos + deviation
    # Use a while loop to ensure that the final genetic position is not negative
    while(any(marker.info.non.redundant1$pos < 0)) {
      deviation = deviation + 10
      marker.info.non.redundant1$pos <- marker.info.non.redundant$pos + deviation
    }
    # Return the new marker info
    return(marker.info.non.redundant1)
  } })

# Collapse the list
non.redundant.markers <- do.call("rbind", non.redundant.marker.list)
          
# Subset the marker matrix
CAP.M <- CAP.M[,as.character(non.redundant.markers$rs)]

# Change all heterozygotes to NA
CAP.M[CAP.M == 0] <- NA

# Impute marker data based on the mean across a marker, 
# then round the non-integer values to either -1 or 1
CAP.M.impute <- apply(X = CAP.M, MARGIN = 2, FUN = function(snp) {
  # Find the index of NA
  na.index <- is.na(snp)
  # Calculate the mean
  snp.mean <- mean(snp, na.rm = T)
  # Round the mean to -1 or 1
  if (snp.mean > 0) snp.mean = 1
  if (snp.mean < 0) snp.mean = -1
  # If the mean is 0, round to 1
  if (snp.mean == 0) snp.mean = 1
  # Replace the missing with the mean
  snp[na.index] <- snp.mean
  # Return
  return(snp)
})

# Remove monomorphic markers again. Those markers that had some heterozygotes
## but were otherwise monomorphic will become monomorphic after imputation
M.poly <- abs(colMeans(CAP.M.impute, na.rm = T)) != 1
CAP.M.impute <- CAP.M.impute[,M.poly]

non.redundant.markers <- non.redundant.markers[M.poly,]
  
### Processing of markers

# Convert the marker positions to Morgans (currently in 1000 * cM)\
non.redundant.markers$pos <- non.redundant.markers$pos / 1000 / 100
# Remove the "UNK" factor from the chrom levels
non.redundant.markers$chrom <- as.numeric(non.redundant.markers$chrom)
# Remove row.names
row.names(non.redundant.markers) <- NULL
  
  
CAP.M.final <- CAP.M.impute
# Since all genotypes are completely homozygous, it is easy to create the gamete matrix
CAP.haploids <- rbind(CAP.M.final, CAP.M.final)
# Change all -1 to zeros
CAP.haploids[CAP.haploids == -1] <- 0
# Rename the genotypes
row.names(CAP.haploids) <- paste(rep(row.names(CAP.M.final), 2), rep(1:2, each = nrow(CAP.M.final)), sep = ".")
# Sort on row.names
CAP.haploids <- CAP.haploids[order(row.names(CAP.haploids)),]

# The marker info comes from the non-redundant markers
CAP.markers <- non.redundant.markers

