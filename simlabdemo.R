### Simulation Lab Demo
# 18 May, 2018

# Source functions
source("Code/hypred_simulation_FUNCTIONS.R")
source("../../../Genomic Selection/Population Structure/Code/GWAS.utils.R")

library(rrBLUP)
library(hypred)



# Load data
load("Files/Barley_CAP_simuation_starting_material.RData")


# Subset the marker data for 192 MN lines
CAP.haploid <- CAP.gametes.0.03[c(1:384),]
head(CAP.haploid)

# Look at the genetic map
marker.map <- sampled.markers.0.03
head(marker.map)

chr.len <- as.numeric(tapply(X = marker.map$pos, INDEX = marker.map$chrom, FUN = max))
chr.len

# Create a genome
hv.genome <- make.genome(n.chr = 7,
                         chr.len = chr.len,
                         n.chr.snps = 100,
                         genetic.map = marker.map$pos)

# Design the trait architecture
hv.genome <- trait.architecture(genome = hv.genome,
                                n.QTL = 50,
                                qtl.add.eff = "geometric")

# Summary
hv.genome

# Phenotype the CAP lines
values <- measure.values(genome = hv.genome, gametes = CAP.haploid, h2 = 0.7, V_e.scale = 8, n.env = 5, n.rep = 2)
phenos <- values$mean.pheno.values
hist(phenos)

# Look at the variance
values$var.components$true

# Find the highest and lowest performing lines
highest.line <- row.names(phenos)[order(phenos, decreasing = T)[1]]
lowest.line <- row.names(phenos)[order(phenos, decreasing = F)[1]]

# Make a crossing block
crossing.block <- make.crossing.block(parent1.lines = highest.line, parent2.lines = lowest.line, n.crosses = 1)

# Make a F2-derived mapping population
mapping.pop.haploid <- make.population(genome = hv.genome, 
                               named.parental.gametes = CAP.haploid,
                               crossing.block = crossing.block,
                               N = 500,
                               cycle.number = 1,
                               generations = 5, # F_x where generations = x - 1
                               pop.type = "inbred",
                               mutation.rate.snp = 7e-8,
                               mutation.rate.qtl = 7e-8)

# Phenotype the mapping population
mapping.values <- measure.values(genome = hv.genome, gametes = mapping.pop.haploid, h2 = 0.7, V_e.scale = 8, n.env = 3, n.rep = 2)
mapping.phenos <- mapping.values$mean.pheno.values
mapping.values$var.components$true

# Genotype the mapping population
mapping.genos <- genotype.markers(gametes = mapping.pop.haploid, genome = hv.genome, include.QTL = F)

# Setup the GWAS
gwas.pheno <- data.frame(gid = as.character(row.names(mapping.phenos)), pheno = mapping.phenos)
gwas.geno <- data.frame(rs = colnames(mapping.genos), 
                        chr = marker.map$chrom[marker.map$rs %in% colnames(mapping.genos)],
                        pos = (marker.map$pos[marker.map$rs %in% colnames(mapping.genos)] * 100),
                        t(mapping.genos))

colnames(gwas.geno) <- c("rs", "chr", "pos", row.names(mapping.genos))


# Run a GWAS
demo.GWAS <- GWAS(pheno = gwas.pheno, geno = gwas.geno, plot = F, P3D = T)
# Plot
manhattan(input = demo.GWAS, fdr.level = 0.10)


# Where are the QTL?
pos.and.eff <- cbind(slot(hv.genome, "pos.add.qtl")$M, slot(hv.genome, "add.and.dom.eff")$add, rep(1:7, each = 8))
pos.and.eff






