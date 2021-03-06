% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convenience_functions.R
\name{make.contrast}
\alias{make.contrast}
\title{Make a matrix of contrasts}
\usage{
make.contrast(unphenotyped.index, n.total)
}
\arguments{
\item{unphenotyped.index}{The index of the unphenotyped entries in the 
relationship matrix (A)}

\item{n.total}{The size of the whole population}
}
\description{
Creates a matrix of contrasts between the individuals not in the training 
population (i.e. the candidates) and the mean of the whole population.
}
\references{
Rincent, R., Laloe, D., Nicolas, S., Altmann, T., Brunel, D., Revilla, P., 
Moreau, L. (2012). Maximizing the Reliability of Genomic Selection by 
Optimizing the Calibration Set of Reference Individuals: Comparison of 
Methods in Two Diverse Groups of Maize Inbreds (Zea mays L.). Genetics, 
192(2), 715–728. http://doi.org/10.1534/genetics.112.141473
}

