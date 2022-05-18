# Introduction


This web tutorial is derived from 'A guide to genome-wide association analysis and post-analytic interrogation' (Statistics in Medicine, in review). The tutorial presents fundamental concepts and specific software tools for implementing a complete genome wide association (GWA) analysis, as well as post-analytic visualization and interrogation of potentially novel findings. In this tutorial we use complete GWA data on 1401 individuals from [the PennCATH study of coronary artery disease (CAD).](http://www.ncbi.nlm.nih.gov/pubmed/21239051)

In the steps to follow we begin by demonstrating a method for downloading necessary R packages and setting global parameters as a means for saving progress while working through a GWA analysis. Next, we include quality control steps for both SNP and sample level filtering. The third section is split into principal component calculation for population stratification in statistical modeling, as well as imputation of non-typed SNPs using 1000 Genomes reference genotype data. We then demonstrate strategies to carry out the GWA analysis on the typed data using basic linear modeling functionality in R, as well as imputed data using functionality contained within the `snpStats` package. Finally, we demonstrate means for post-analytic interrogation, including methods for evaluating the performance of statistical models, as well as visualization of the global and subsetted GWAS output.

## Installing necessary packages

```r
# Run this once interactively to download and install BioConductor packages and other packages.

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

bio.packages = c("snpStats", "SNPRelate","rtracklayer", "biomaRt")
BiocManager::install(bio.packages)

other.packages = c("plyr", "LDheatmap", "doParallel", "ggplot2", "coin", "igraph", "devtools", "downloader")
install.packages(other.packages)

# GenABEL has moved to CRAN archive. The below command for local installation from CRAN archive.
install.packages("https://cran.r-project.org/src/contrib/Archive/GenABEL/GenABEL_1.7-6.tar.gz", 
type = "source", repos = NULL)

install.packages("https://github.com/merns/postgwas/releases/download/1.11-2/postgwas_1.11-2.zip", repos=NULL)

```

## Loading packages

```r
# Loading packages ...
lapply(bio.packages, require, character.only = TRUE)
lapply(other.packages, require, character.only = TRUE)
library(GenABEL)
library(postgwas)

# Moodle for GenABEL errors
source("R/ztransform.R")
source("R/rntransform.R")
source("R/estlambda.R")
source("R/GWAA.R")

```


## Configuring global parameters

Customize and Run [globals.R](R/globals.R)

```r
source("R/globals.R")

# Downloading support files
# Download and unzip data needed for this tutorial

library(downloader)

download(urlSupport, zipSupport.fn)
unzip(zipSupport.fn, exdir = data.dir)

```

