### Map QTLs 1 of 3
pop <- "NBH"
source("/home/jmiller1/Sex_Map/code/control_file.R")
## For plotting marker_dens <- list()

# Table of Chroms with sig QTLs test.QTLs <- read.table(file.path(basedir,
# 'rQTL/metadata/QTLs.txt'), sep = '\t', header = T)

## Get chrom number vector test.QTLs$chrm.n <- gsub('chr', '', test.QTLs$chrom)

# print(paste(pop, X, sep = ' '))

## read in the QTL cross
cross.18 <- read.cross.jm(file = file.path(indir, "NBH_CHR5_filt_conv.csvr"), format = "csvr", 
  geno = c(1:3), estimate.map = FALSE)

### Pull names from plinkfile
path <- file.path(plinkdir, "NBH_CHR5_filt_conv.ped")
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
cross.18$pheno$ID <- paste(popname, indname, sep = "_")

## Remove problematic individuals (found by kinship analysis)
popdir <- file.path("/home/jmiller1/QTL_Map_Raw/popgen", "rQTL", pop, "REMAPS")
con <- file(file.path(popdir, "kinship.keep.ind.txt"), open = "r")
keepers <- readLines(con)
close(con)

print("Dropping kinship outliers")
cross.18 <- subset(cross.18, ind = cross.18$pheno$ID %in% keepers)
cross.18 <- subset(cross.18, ind = !is.na(cross.18$pheno$Phen))


sex <- read.table("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/data/sex.txt")
rownames(sex) <- sex$ID
cross.18$pheno$sex <- sex[cross.18$pheno$ID, 2]

cross.a <- subset(cross.18, ind = cross.18$pheno$ID[which(cross.18$pheno$sex == 1)])
cross.b <- subset(cross.18, ind = cross.18$pheno$ID[which(cross.18$pheno$sex == 0)])

### Table before missing filter
gt.a <- geno.table(cross.a)
gt.b <- geno.table(cross.b)

pval.a <- log10(gt.a$P.value)
pval.b <- log10(gt.b$P.value)

#### Filter Conservative
drop.missing <- function(cross, M) {
  gt <- geno.table(cross)
  todrop <- rownames(gt[which(gt$missing > M), ])
  paste(length(todrop), "markers dropped")
  cross <- drop.markers(cross, unlist(todrop))
  return(cross)
}
## missing more than 35/46 filtered missing more than 33/43 filtered
cross.a <- drop.missing(cross.a, 11)
cross.b <- drop.missing(cross.b, 10)
gt.missing.a <- geno.table(cross.a)
gt.missing.b <- geno.table(cross.b)

gt.cross.par <- NA
cutoff <- 0.1

cross.a <- drop.markers(cross.a, rownames(gt.missing.a[gt.missing.a$P.value < cutoff, 
  ]))
gt.pval.a <- geno.table(cross.a)
sum(nmar(cross.a))

cross.b <- drop.markers(cross.b, rownames(gt.missing.b[gt.missing.b$P.value < cutoff, 
  ]))
gt.pval.b <- geno.table(cross.b)
sum(nmar(cross.b))

print("forming initial linkage groups to fix phase...")
grpRf <- 0.3
grpLod <- 4  ## Standard LG form LOD
finLod <- 5  ## Higher final NBH LOD
finRf <- 0.15

save.image("/home/jmiller1/Sex_Map/datasetsexmap.rsave")
cross.test.a <- formLinkageGroups(cross.a, max.rf = grpRf, min.lod = grpLod, reorgMarkers = TRUE)
### Flip phase and retest
cross.test.a <- switchAlleles(cross.test.a, markernames(cross.test.a, chr = chrnames(cross.test.a)[1]))
cross.test.a <- formLinkageGroups(cross.test.a, max.rf = grpRf, min.lod = grpLod, 
  reorgMarkers = TRUE)
cross.test.a <- switchAlleles(cross.test.a, markernames(cross.test.a, chr = chrnames(cross.test.a)[1]))
cross.test.a <- formLinkageGroups(cross.test.a, max.rf = grpRf, min.lod = grpLod, 
  reorgMarkers = TRUE)
cross.test.a <- subset(cross.test.a, chr = which.max(nmar(cross.test.a)))
cross.test.a <- orderMarkers(cross.test.a, chr = 1, window = 5, use.ripple = T, error.prob = 0.05, 
  map.function = "kosambi", sex.sp = F, maxit = 1000, tol = 0.01)
cross.test.a <- removeDoubleXO(cross.test.a, verbose = T)
POS.map.a <- est.map(cross.test.a, error.prob = 0.01, map.function = "kosambi", chr = 1, 
  maxit = 1000)
cross.test.a <- replace.map(cross.test.a, POS.map.a)
save.image("/home/jmiller1/Sex_Map/datasetsexmap.rsave")

cross.test.b <- formLinkageGroups(cross.b, max.rf = grpRf, min.lod = grpLod, reorgMarkers = TRUE)
### Flip phase a couple times to see if additional markers fall into the largest LG
### and retest
cross.test.b <- switchAlleles(cross.test.b, markernames(cross.test.b, chr = chrnames(cross.test.b)[1]))
cross.test.b <- formLinkageGroups(cross.test.b, max.rf = grpRf, min.lod = grpLod, 
  reorgMarkers = TRUE)
cross.test.b <- switchAlleles(cross.test.b, markernames(cross.test.b, chr = chrnames(cross.test.b)[1]))
cross.test.b <- formLinkageGroups(cross.test.b, max.rf = grpRf, min.lod = grpLod, 
  reorgMarkers = TRUE)
cross.test.b <- subset(cross.test.b, ind = cross.test.b$pheno$ID[which(cross.test.b$pheno$sex == 
  0)])
cross.test.b <- subset(cross.test.b, chr = which.max(nmar(cross.test.b)))
cross.test.b <- orderMarkers(cross.test.b, chr = 1, window = 5, use.ripple = T, error.prob = 0.05, 
  map.function = "kosambi", sex.sp = F, maxit = 1000, tol = 0.01)
cross.test.b <- removeDoubleXO(cross.test.b, verbose = T)
POS.map.b <- est.map(cross.test.b, error.prob = 0.01, map.function = "kosambi", chr = 1, 
  maxit = 1000)
cross.test.b <- replace.map(cross.test.b, POS.map.b)
save.image("/home/jmiller1/Sex_Map/datasetsexmap.rsave")
