
# ---- packages ----

# ---- installation ----

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

bio.packages = c("snpStats", "SNPRelate","rtracklayer", "biomaRt")
#BiocManager::install(bio.packages)

other.packages = c("plyr", "LDheatmap", "doParallel", "ggplot2", "coin", "igraph", "devtools")
#install.packages(other.packages)

#install.packages("https://github.com/merns/postgwas/releases/download/1.11-2/postgwas_1.11-2.zip", repos=NULL)

# ---- loading ----

# Loading packages ...
lapply(bio.packages, require, character.only = TRUE)
lapply(other.packages, require, character.only = TRUE)
library(postgwas)

# Moodle for GenABEL errors
source("GWAS-extras/ztransform.R")
source("GWAS-extras/rntransform.R")
source("GWAS-extras/estlambda.R")
source("GWAS-extras/GWAA.R")


# ---- step1 ----
# Reading data into R

source("R/globals.R")

# ---- step1-a ----
library(snpStats)

# Read in PLINK files to create list
geno <- read.plink(gwas.fn$bed, gwas.fn$bim, gwas.fn$fam, na.strings = ("-9"))

# ---- step1-b ----
# Obtain the SnpMatrix object (genotypes) table from geno list
# Note: Phenotypes and covariates will be read from the clinical data file, below
genotype <- geno$genotype
print(genotype)                  # 861473 SNPs read in for 1401 subjects

# Obtain the SNP information from geno list
genoBim <- geno$map
colnames(genoBim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
print(head(genoBim))

# Remove raw files to open up memory
rm(geno)

# ---- step1-c ----
# Read in clinical file
clinical <- read.csv(clinical.fn,
                     colClasses=c("character", "factor", "factor", rep("numeric", 4)))
rownames(clinical) <- clinical$FamID
print(head(clinical))

# ---- step1-d ----
# Subset genotype for subject data
genotype <- genotype[clinical$FamID, ]
print(genotype)  # Tutorial: All 1401 subjects contain both clinical and genotype data


# ---- step2 ----
# SNP-level filtering (part 1)

# ---- step2-a ----
# Create SNP summary statistics (MAF, call rate, etc.)
snpsum.col <- col.summary(genotype)
print (head(snpsum.col))

# ---- step2-b ----
# Setting thresholds
call <- 0.95
minor <- 0.01

# Filter on MAF and call rate
use <- with(snpsum.col, (!is.na(MAF) & MAF > minor) & Call.rate >= call)
use[is.na(use)] <- FALSE # Remove NA's as well

cat (ncol(genotype) - sum(use), "SNPs will be removed due to low MAF or call rate.\n") #203287 SNPs will be removed

# Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
genotype <- genotype[,use]
snpsum.col <- snpsum.col[use,]

print (genotype) # 658186 SNPs remain


# ---- step3 ----
# sample-level filtering

# ---- step3-a ----
library (SNPRelate) # Estimating LD, relatedness, PCA
library(plyr)

# Create sample statistics (Call rate, Heterozygosity)
snpsum.row <- row.summary(genotype)

# Add the F stat (inbreeding coefficient) to snpsum.row
MAF <- snpsum.col$MAF
memory.limit(size = 99999999999)
callmatrix <- !is.na(genotype)
hetExp <- callmatrix %*% (2*MAF*(1-MAF))
hetObs <- with(snpsum.row, Heterozygosity*(ncol (genotype))*Call.rate)
snpsum.row$hetF <- 1 - (het0bs/hetExp)

head(snpsum.row)

# ---- step3-b =
# Setting thresholds

sampcall <- 0.95  # Sample call rate cut-off
hetcutoff <- 0.1  # Inbreeding coefficient cut-off

sampleuse <- with(snpsum.row, lis.na(Call.rate) & Call.rate > sampcall & abs(hetF) <= hetcutoff)
sampleuse[is.na(sampleuse)] <- FALSE # remove NA's as well
cat(nrow(genotype)-sum(sampleuse), "subjects will be removed due to low sample call rate or inbreeding coefficient.\n")
# 0 subjects removed

# Subset genotype and clinical data for subjects who pass call rate and heterozygosity crtieria
genotype <- genotype[sampleuse,]
clinical<- clinical[rownames(genotype), ]











