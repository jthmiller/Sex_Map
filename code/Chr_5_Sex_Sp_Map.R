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

cross.a5 <- subset(cross.18, ind = cross.18$pheno$ID[which(cross.18$pheno$sex == 1)])
cross.b5 <- subset(cross.18, ind = cross.18$pheno$ID[which(cross.18$pheno$sex == 0)])



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
save.image("/home/jmiller1/Sex_Map/dataset/sexmap.rsave")

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
###save.image("/home/jmiller1/Sex_Map/dataset/sexmap.rsave")







load("/home/jmiller1/Sex_Map/dataset/sexmap.rsave")
chr <- 5
##### MAP only chr 5 scaffolds ######
chr5 <- read.table("/home/jmiller1/Sex_Map/metadata/chrm5_only.txt", header=F)
chr5 <- as.character(unlist(chr5$V1))

### those that have markers in the map
mapping <- unique(gsub(":.*",'',markernames(cross.18)))
chr5.mapped <- chr5[chr5 %in% mapping]

###cross.18.mapped.scafs <- subset(cross.18,chr)

chr5.mapped.only <- drop.markers(cross.18,markernames(cross.18)[!gsub(":.*",'',markernames(cross.18)) %in% chr5.mapped])

### Table before missing filter
gt <- geno.table(chr5.mapped.only)
pval.a <- log10(gt$P.value)

#### Filter Conservative
drop.missing <- function(cross, M) {
  gt <- geno.table(cross)
  todrop <- rownames(gt[which(gt$missing > M), ])
  paste(length(todrop), "markers dropped")
  cross <- drop.markers(cross, unlist(todrop))
  return(cross)
}
chr5.mapped.only.a <- subset(chr5.mapped.only, ind = chr5.mapped.only$pheno$ID[which(chr5.mapped.only$pheno$sex == 1)])
chr5.mapped.only.b <- subset(chr5.mapped.only, ind = chr5.mapped.only$pheno$ID[which(chr5.mapped.only$pheno$sex == 0)])

## missing more than 35/46 filtered missing more than 33/43 filtered
chr5.mapped.only.a  <- drop.missing(chr5.mapped.only.a , 5)
chr5.mapped.only.b  <- drop.missing(chr5.mapped.only.b , 5)
gt.missing.a <- geno.table(chr5.mapped.only.a)
gt.missing.b <- geno.table(chr5.mapped.only.b)
sum(nmar(chr5.mapped.only.b))
sum(nmar(chr5.mapped.only.a))

#########
### No segregation distortion filter
#########

grpRf <- 0.3
grpLod <- 3  ## Standard LG form LOD
finLod <- 4  ## Higher final NBH LOD
finRf <- 0.05



cross.test.b <- formLinkageGroups(chr5.mapped.only.b, max.rf = grpRf, min.lod = grpLod, reorgMarkers = TRUE)
### Flip phase a couple times to see if additional markers fall into the largest LG
### and retest
cross.test.b <- switchAlleles(cross.test.b, markernames(cross.test.b, chr = chrnames(cross.test.b)[1]))
cross.test.b <- formLinkageGroups(cross.test.b, max.rf = grpRf, min.lod = grpLod,
  reorgMarkers = TRUE)
cross.test.b <- switchAlleles(cross.test.b, markernames(cross.test.b, chr = chrnames(cross.test.b)[1]))

cross.test.b <- formLinkageGroups(cross.test.b, max.rf = finRf, min.lod = finLod,
  reorgMarkers = TRUE)


cross.test.bx1 <- subset(cross.test.b, chr = 1) ## All BB scaffolds
cross.test.bx2 <- subset(cross.test.b, chr = 2) ## All AB or BB
cross.test.bx3 <- subset(cross.test.b, chr = 3) ## All AB or BB
cross.test.bx4 <- subset(cross.test.b, chr = 4) ## 1:4:3 AA:AB:BB



cross.test.b <- orderMarkers(cross.test.b, chr = 1, window = 5, use.ripple = T, error.prob = 0.05,
  map.function = "kosambi", sex.sp = F, maxit = 1000, tol = 0.01)
cross.test.b <- removeDoubleXO(cross.test.b, verbose = T)
POS.map.b <- est.map(cross.test.b, error.prob = 0.01, map.function = "kosambi", chr = 1,
  maxit = 1000)
cross.test.b <- replace.map(cross.test.b, POS.map.b)



##### Seems like male
cross.test.a <- formLinkageGroups(chr5.mapped.only.a, max.rf = grpRf, min.lod = grpLod, reorgMarkers = TRUE)
### Flip phase and retest
cross.test.a <- switchAlleles(cross.test.a, markernames(cross.test.a, chr = chrnames(cross.test.a)[1]))
cross.test.a <- formLinkageGroups(cross.test.a, max.rf = grpRf, min.lod = grpLod,
  reorgMarkers = TRUE)
cross.test.a <- switchAlleles(cross.test.a, markernames(cross.test.a, chr = chrnames(cross.test.a)[1]))
cross.test.a <- formLinkageGroups(cross.test.a, max.rf = grpRf, min.lod = grpLod,
  reorgMarkers = TRUE)


cross.test.aX <- subset(cross.test.a, chr = 1) ### Normally seg 1:2:1
cross.test.ax <- orderMarkers(cross.test.a, chr = 1, window = 5, use.ripple = T, error.prob = 0.05,
  map.function = "kosambi", sex.sp = F, maxit = 1000, tol = 0.01)
cross.test.ax <- removeDoubleXO(cross.test.a, verbose = T)
POS.map.ax <- est.map(cross.test.a, error.prob = 0.01, map.function = "kosambi", chr = 1,
  maxit = 1000)
cross.test.ax <- replace.map(cross.test.a, POS.map.a)

cross.test.ay <- subset(cross.test.a, chr = 2) ### ALl hets
cross.test.ay <- orderMarkers(cross.test.ay, chr = 1, window = 5, use.ripple = T, error.prob = 0.05,
  map.function = "kosambi", sex.sp = F, maxit = 1000, tol = 0.01)
cross.test.ay <- removeDoubleXO(cross.test.ay, verbose = T)
POS.map.ay <- est.map(cross.test.ay, error.prob = 0.01, map.function = "kosambi", chr = 1,
  maxit = 1000)
cross.test.ay <- replace.map(cross.test.ay, POS.map.a)




cross.test.bx1 <- subset(cross.test.b, chr = 1) ## All BB scaffolds
cross.test.bx2 <- subset(cross.test.b, chr = 2) ## All AB or BB ()
cross.test.bx3 <- subset(cross.test.b, chr = 3) ## All AB or BB
cross.test.bx4 <- subset(cross.test.b, chr = 4) ## 1:4:3 AA:AB:BB (might be partial recombine?)


cross.test.a1 <- subset(cross.test.a, chr = 1) ### Normally seg 1:2:1 (snps fixed in each parent)
cross.test.a2 <- subset(cross.test.a, chr = 2) ### All hets

scaftell <- function(Z){
  sort(table(gsub(":.*",'',markernames(Z))))
}


scaftell(cross.test.a1) ## AA AB BB (1:2:1) NW_012224575.1 (X)
scaftell(cross.test.a2) ## AB (0:1:0) NW_012224575.1 (Z)


scaftell(cross.test.bx1)### BB (0:0:1) NW_012224575.1 (Z)
scaftell(cross.test.bx2)## AB BB (0:1:1) NW_012224575.1 (X)



cross.test.a1 <- orderMarkers(cross.test.a1, chr = 1, window = 5, use.ripple = T, error.prob = 0.05,
  map.function = "kosambi", sex.sp = F, maxit = 1000, tol = 0.01)
cross.test.a1 <- removeDoubleXO(cross.test.a1, verbose = T)
POS.map.a1 <- est.map(cross.test.a1, error.prob = 0.01, map.function = "kosambi", chr = 1,
  maxit = 1000)
cross.test.a1 <- replace.map(cross.test.a1, POS.map.a)



save.image("/home/jmiller1/Sex_Map/dataset/sexmap.mapped_scafs_only.rsave")
