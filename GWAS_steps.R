
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
snpsum.row$hetF <- 1 - (hetObs/hetExp)

head(snpsum.row)

# ---- step3-b ----
# Setting thresholds

sampcall <- 0.95  # Sample call rate cut-off
hetcutoff <- 0.1  # Inbreeding coefficient cut-off

sampleuse <- with(snpsum.row, !is.na(Call.rate) & Call.rate > sampcall & abs(hetF) <= hetcutoff)
sampleuse[is.na(sampleuse)] <- FALSE # remove NA's as well
cat(nrow(genotype)-sum(sampleuse), "subjects will be removed due to low sample call rate or inbreeding coefficient.\n")
# 0 subjects removed

# Subset genotype and clinical data for subjects who pass call rate and heterozygosity crtieria
genotype <- genotype[sampleuse,]
clinical<- clinical[rownames(genotype), ]

# ---- step3-c ----
# Checking for Relatedness

ld.thresh <- 0.2 # LD cut-off
kin.thresh <- 0.1 # Kinship cut-off

# Create gds file, required for SNPRelate functions
snpgdsBED2GDS(gwas.fn$bed, gwas.fn$fam, gwas.fn$bim, gwas.fn$gds)
genofile <- openfn.gds(gwas.fn$gds, readonly = FALSE)

# Automatically added "-1" sample suffixes are removed
gds.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
gds.ids <- sub("-1", "", gds.ids)
add.gdsn(genofile, "sample.id", gds.ids, replace = TRUE)

# Prune SNPs for IBD analysis
set.seed(1000)

geno.sample.ids <- rownames (genotype)
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.thresh,
                          sample.id = geno.sample.ids, # Only analyze the filtered samples
                          snp.id = colnames(genotype)) # Only analyze the filtered SNPs
snpset.ibd <- unlist(snpSUB, use.names=FALSE)
cat(length(snpset.ibd), "will be used in IBD analysis\n") # expect 72812 SNPs


# ----step3-d ----
# Find IBD coefficients using Method of Moments procedure. Include pairwise kinship.
ibd <- snpgdsIBDMoM(genofile, kinship=TRUE,
                    sample.id = geno.sample.ids,
                    snp.id = snpset.ibd,
                    num.thread = 1)
ibdcoeff <- snpgdsIBDSelection(ibd) # Pairwise sample comparison
head(ibdcoeff)

# ---- step3-e ----
# Check if there are any candidates for relatedness
ibdcoeff <- ibdcoeff[ ibdcoeff$kinship >= kin.thresh, ]

# iteratively remove samples with high kinship starting with the sample with the most pairings
related.samples <- NULL
while(nrow(ibdcoeff) > 0) {
  
  # count the number of occurrences of each and take the top one
  sample.counts <- arrange(count(c(ibdcoeff$ID1, ibdcoeff$ID2)), -freq)
  rm.sample <- sample.counts[1, 'x']
  cat("Removing sample", as.character(rm.sample), 'too closely related to', sample.counts[1, 'freq'],'other samples. \n')

  # remove from ibdcoeff and add to list
  ibdcoeff <- ibdcoeff[ibdcoeff$ID1 != rm.sample & ibdcoeff$ID2 != rm.sample,]
  related.samples <- c(as.character(rm.sample), related.samples)
}

# filter genotype and clinical to include only unrelated samples
genotype <- genotype[ !(rownames(genotype) %in% related.samples), ]
clinical <- clinical[ !(clinical$FamID %in% related.samples), ]

geno.sample.ids <- rownames(genotype)

cat(length(related.samples), "similar samples removed due to correlation coefficient >=", kin.thresh,"\n")
print(genotype) # expect all 1401 subjects remain

# ---- step3-f ----
# Checking for ancestry

# Find PCA matrix
pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids, snp.id = snpset.ibd, num.thread=1)

# Create data frame of first two principal components
pctab <- data.frame(sample.id = pca$sample.id,
                    PC1 = pca$eigenvect[,1], # the first eigenvector
                    PC2 = pca$eigenvect[,2], # the second eigenvector
                    stringsAsFactors = FALSE)

# Plot the first two principal components
plot(pctab$PC1, pctab$PC2, xlab="Principal Component 1", ylab="Principal Component 2", main="Ancestry Plot")

# ---- step4 ---
# Hardy-Weinberg SNP filtering on CAD controls

hardy <- 10^-6 # HWE cut-off

CADcontrols <- clinical[clinical$CAD==0, 'FamID']
snpsum.colCont <- col.summary(genotype[CADcontrols,])
HWEuse <- with(snpsum.colCont, !is.na(z.HWE) & ( abs(z.HWE) < abs( qnorm(hardy/2) ) ) )
rm(snpsum.colCont)

HWEuse[is.na(HWEuse)] <- FALSE # Remove NA's as well
cat(ncol(genotype)-sum(HWEuse) ,"SNPs will be removed due to high HWE.\n") # 1296 SNPs removed
                         
# Subset genotype and SNP summary data for SNPs that pass HWE criteria
genotype <- genotype[,HWEuse]
                         
print (genotype) # 656890 SNPs remain



# --- steps5 ----
# Generating principal components for modeling

# ---- step5-a ---
# Set LD threshold to 0.2

ld.thresh <- 0.2

set.seed(1000)
geno.sample.ids <- rownanes(genotype)
SnpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.thresh,
                          sample.id = geno.sample.ids, # Only analyze the filtered samples
                          snp.id = colnames(genotype)) # Only analyze the filtered SNPs
snpset.pca <- unlist(snpSUB, use.names=FALSE)
cat(length(snpset.pca),"\n") # 72578 SNPs will be used in PCA analysis

pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids, snp.id = snpset.pca, num.thread=1)

# Find and record first 10 principal components
# pcs will be a N:10 matrix. Each column is a principal component.
pcs <- data.frame(FamID = pca$sample.id, pca$eigenvect[,1 : 10],
                  stringsAsFactors = FALSE)
colnames (pcs) [2:11]<-paste("pc", 1:10, sep = "")

print(head(pcs))

# ---- step6 ----
# Genotype imputation

#--- step6-a ---
# Read in 1000g data for given chromosome 16
thougeno <- read.pedfile(onethou.fn$ped, snps = onethou.fn$info, which=1)

# Obtain genotype data for given chromosome
genoMatrix <- thougeno$genotypes

# Obtain the chromosome position for each SNP
support <- thougeno$map

colnames(support) <- c("SNP", "position", "A1", "A2")
head(support)

# Imputation of non-typed 1000g SNPs
presSnps <- colnames(genotype)

# Subset for SNPs on given chromosome
presSnps <- colnames(genotype)
presDatChr <- genoBim[genoBim$SNP %in% presSnps & genoBim$chr==16, ]
targetSnps <- presDatChr$SNP

# Subset 1000g data for our SNPs
# "missing" and "present" are snpMatrix objects needed for imputation rules
is.present <- colnames(genoMatrix) %in% targetSnps

missing <- genoMatrix[,!is.present]
print (missing) # Almost 400,000 SNPs

present <- genoMatrix[,is.present]
print (present) # Our typed SNPs

# Obtain positions of SNPs to be used for imputation rules
pos.pres <- support$position[is.present]
pos.miss <- support$position[!is.present]

# Calculate and store imputation rules using snp.imputation()
rules <- snp.imputation(present, missing, pos.pres, pos.miss)

# ---- step6-b
# Remove failed imputations
rules <- rules[can.impute(rules)]
cat("Imputation rules for", length(rules), "SNPs were estimated\n") # Imputation rules for 197888 SNPs were estimated

# Quality control for imputation certainty and MAF
# Set thresholds
r2threshold <- 0.7
minor <- 0.01

# Filter on imputation certainty and MAF
rules <- rules[imputation.r2(rules) >= r2threshold]

cat(length(rules),"imputation rules remain after uncertain imputations were removed\n")
# 162565 imputation rules remain after uncertain impuations were removed

rules <- rules[imputation.maf(rules) >= minor]
cat(length(rules) ,"imputation rules remain after MAF filtering\n") # 162565 imputation rules remain after MAF filtering

# Obtain posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]

imputed <- impute.snps(rules, target, as.numeric=FALSE)
print(imputed) # 162565 SNPs were imputed

# ---- step7 ----
# Association analysis of typed SNPs

# ---- step7-a ----
# Merge clinical data and principal components to create phenotype table
phenoSub <- merge(clinical,pcs) # data.frame => [ FamID CAD sex age hdl pci pc2 ... pc10 ]

# We will do a rank-based inverse normal transformation of hdl
phenoSub$phenotype <- rntransform(phenoSub$hdl, family="gaussian")

# Show that the assumptions of normality met after transformation
par(mfrow=c(1,2))
hist(phenoSub$hdl, main="Histogram of HDL", xlab="HDL")
hist(phenoSub$phenotype, main="Histogram of Tranformed HDL", xlab="Transformed HDL")

# Remove unnecessary columns from table
phenoSub$hdl <- NULL
phenoSub$ldl <- NULL
phenoSub$tg <- NULL
phenoSub$CAD <- NULL

# Rename columns to match names necessary for GWAS() function
phenoSub <- rename(phenoSub, replace=c(FamID="id"))

# Include only subjects with hdl data
phenoSub <- phenoSub[!is.na(phenoSub$phenotype), ]
# 1309 subjects included with phenotype data

print (head (phenoSub))

# ---- step7-b ----
# Run GWAS analysis (using parallel processing)

# Note: This function writes a file, but does not produce an R object
start <- Sys.time()
GWAA(genodata=genotype, phenodata=phenoSub, filename=gwaa.fname)
end <- Sys.time()
print(end-start)

# ---- step8 ----
# Association analysis of imputed SNPs

# ---- step8-a ----
# Carry out association testing for imputed SNPs using snp.rhs.tests()
rownames (phenoSub) <- phenoSub$id

imp <- snp.rhs.tests(phenotype ~ sex + age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                     family = "Gaussian", data = phenoSub, snp.data = target, rules = rules)

# Obtain p values for imputed SNPs by calling methods on the returned GlmTests object.
results <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results <- results[!is.na(results$p.value),]

# Write a file containing the results
write.csv(results, impute.out.fname, row.names=FALSE)

# Merge imputation testing results with support to obtain coordinates
imputeOut<-merge(results, support[, c("SNP", "position")])
imputeOut$chr <- 16

imputeOut$type <- "imputed"

# Find the -log_10 of the p-values
imputeOut$Neg_logP <- -log10(imputeOut$p.value)

# Order by p-value
imputeOut <- arrange(imputeOut, p.value)
print(head(imputeOut))

# ---- step8-b ----

# Read in file containing protein coding genes coords
genes <- read.csv(protein.coding.coords.fname, stringsAsFactors = FALSE)

# Subset for CETP SNPs
impCETP <- map2gene("CETP", coords = genes, SNPs = imputeOut)

# Filter only the imputed CETP SNP genotypes
impCETPgeno <- imputed[, impCETP$SNP ]


# ---- step9 ----
# Data Integration

# ---- step9-a ----
# Read in GWAS output that was produced by GWAA function
GWASout <- read.table(gwaa.fname, header=TRUE, colClasses=c("character", rep("numeric",4)))

# Find the -log_10 of the p-values
GWASout$Neg_logP <- -log10(GWASout$p.value)

# Merge output with genoBim by SNP name to add position and chromosome number
GWASout <- merge(GWASout, genoBim[,c("SNP", "chr", "position")])

# Order SNPs by significance
GWASout <- arrange(GWASout, -Neg_logP)
print(head(GWASout))

# ---- step9-b ----
# Combine typed and imputed
GWASout$type <- "typed"

GWAScomb<-rbind.fill(GWASout, imputeOut)
head(GWAScomb)
tail(GWAScomb)

# Subset for CETP SNPs
typCETP <- map2gene("CETP", coords = genes, SNPs = GWASout)

# Combine CETP SNPs from imputed and typed analysis
CETP <- rbind.fill(typCETP, impCETP) [,c("SNP","p.value","Neg_logP","chr", "position", "type", "gene")]
print (CETP)

# ---- step10 ----
# Visualizing and QC of GWA findings

par(mfrow=c(1,1))

# ---- step10-a ----
# Create Manhattan Plot
GWAS_Manhattan(GWAScomb)

# ---- step10-b ----
# Rerun the GWAS using unadjusted model
phenoSub2 <- phenoSub[,c("id","phenotype")] # remove all extra factors, leave only phenotype

GWAA(genodata=genotype, phenodata=phenoSub2, filename=gwaa.unadj.fname)
GWASoutUnadj <- read.table(gwaa.unadj.fname, header=TRUE, colClasses=c("character", rep("numeric",4)))

# Create QQ plots for adjusted and unadjusted model outputs
par(mfrow=c(1,2))
lambdaUnadj <- estlambda(GWASoutUnadj$t.value^2, plot=TRUE, method="median")
lambdaAdj <- estlambda(GWASout$t.value^2, plot=TRUE, method="median")
cat(sprintf("Unadjusted lambda: %s\nAdjusted lambda: %s\n", lambdaUnadj$estimate, lambdaAdj$estimate))
                           
# ---- step10-c ----
# Calculate standardized lambda
lambdaUnadj_1000 <- 1+(lambdaUnadj$estimate-1)/nrow(phenoSub) * 1000
lambdaAdj_1000 <- 1+(lambdaAdj$estimate-1) /nrow(phenoSub) * 1000
cat(sprintf("Standardized unadjusted lambda: %s\nStandardized adjusted lambda: %s\n", lambdaUnadj_1000, lambdaAdj_1000))


# ---- step10-d ----

# Add "rs247617" to CETP
CETP <- rbind.fill(GWASout [GWASout$SNP == "rs247617",], CETP)

# Combine genotypes and imputed genotypes for CETP region
subgen <- cbind(genotype[,colnames(genotype) %in% CETP$SNP], impCETPgeno) # CETP subsets from typed and imputed SNPs
                
# Subset SNPs for only certain genotypes
certain <- apply(as(subgen, 'numeric'), 2, function(x) { all(x %in% c(0,1,2,NA)) })
subgen <- subgen[,certain]

# Subset and order CETP SNPs by position
CETP <- CETP[CETP$SNP %in% colnames(subgen),]
CETP <- arrange(CETP, position)
subgen <- subgen[, order(match(colnames (subgen) ,CETP$SNP)) ]

# Create LDheatmap
ld <- ld(subgen, subgen, stats="R.squared") # Find LD map of CETP SNPs
                                    
ll <- LDheatmap(ld, CETP$position, flip=TRUE, name="myLDgrob", title=NULL)

# Add genes, recombination
llplusgenes <- LDheatmap.addGenes(ll, chr = "chr16", genome = "hg19", genesLocation = 0.01)

# Add plot of -log(p)
plot.new()
llQplot2 <- LDheatmap.addGrob(llplusgenes, rectGrob(gp = gpar(col = "white")),height = .34)
pushViewport(viewport(x = 0.483, y= 0.76, width = .91 ,height = .4))

grid.draw(ggplotGrob({
  qplot(position, Neg_logP, data = CETP, xlab="", ylab = "Negative Log P-value", xlim = range(CETP$position),
        asp = 1/10, color = factor(type), colour=c("#000000", "#D55E00")) +
    theme(axis.text.x = element_blank(),
          axis.title.y = element_text(size = rel(0.75)), legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    scale_color_manual(values = c("red", "black"))
}))

# ---- step10-e ----
# Create regional association plot

library(postgwas)

# Create data.frame of most significant SNP only
snps<-data.frame(SNP=c("rs1532625"))

# Change column names necessary to run regionalplot function
GWAScomb <- rename(GWAScomb, replace=c(p.value="P", chr="CHR", position="BP"))

# Edit biomartConfigs so regionalplot function
# pulls from human genome build 37/hg19

myconfig <- biomartConfigs$hsapiens
myconfig$hsapiens$gene$host <- "grch37.ensembl.org"
myconfig$hsapiens$gene$mart <- "ENSEMBL_MART_ENSEMBL"
myconfig$hsapiens$snpéhost <- "grch37.ensembl.org"
myconfig$hsapiens$snpémart <- "ENSEMBL_MART_SNP"

# Run regionalplot using HAPMAP data (pop = CEU)
regionalplot(snps, GWAScomb, biomart.config = myconfig, window.size = 400000, draw.snpname = data.frame (
  snps = c("rs1532625", "rs247617"),
  text = c("rs1532625", "rs247617"),
  angle = c(20, 160),
  length = c(1, 1),
  cex = c(0.8)
),
ld.options = list(
  gts.source = 2,
  max.snps.per.window = 2000,
  rsquare.min = 0.8,
  show.rsquare.text = FALSE
),
out.format = list(file = "png", panels.per.page = 4))




