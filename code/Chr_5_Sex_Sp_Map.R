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

## Remove problematic individuals (found by kinship analysis)
popdir <- file.path("/home/jmiller1/QTL_Map_Raw/popgen", "rQTL", pop, "REMAPS")
con <- file(file.path(popdir, "kinship.keep.ind.txt"), open = "r")
keepers <- readLines(con)
close(con)

print("Dropping kinship outliers")
cross.18 <- subset(cross.18, ind = cross.18$pheno$ID %in% keepers)

sex <- read.table("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/data/sex.txt",stringsAsFactors=F)
rownames(sex) <- sex$ID
sex.vec <- sex[as.character(cross.18$pheno$ID), 'sex']
cross.18$pheno$sex <- sex.vec

a.parent <- geno.table(subset(cross.18, ind = cross.18$pheno$ID=='NBH_NBH1M'))
b.parent <- geno.table(subset(cross.18, ind = cross.18$pheno$ID=='NBH_NBH1F'))
parents <- geno.table(subset(cross.18, ind = (cross.18$pheno$ID=='NBH_NBH1F' | cross.18$pheno$ID=='NBH_NBH1M')))
parents.ab <- subset(cross.18, ind = (cross.18$pheno$ID=='NBH_NBH1F' | cross.18$pheno$ID=='NBH_NBH1M'))

## drop invar
drops.invar <- rownames(parents)[which(parents$AA==2 | parents$BB==2 | parents$AB>0 | parents$missing==2)]
cross.18 <- drop.markers(cross.18, drops.invar)
parents.ab <- drop.markers(parents.ab, drops.invar)

##switch phase relative to parents
toswitch <- c(rownames(a.parent)[which(a.parent$BB==1)],rownames(b.parent)[which(a.parent$AA==1)])
cross.18 <- switchAlleles(cross.18,toswitch)
parents.ab  <- switchAlleles(parents.ab,toswitch)

### those conf in pars
par.table <- geno.table(parents.ab)
parents.fixed <- rownames(par.table)[which(par.table$AA==1 & par.table$BB==1)]
###geno.table(cross.18)[parents.fixed,2:8]

## a is male and 1
## b is female and 0

cross.a.all <- subset(cross.18, ind = cross.18$pheno$ID[which(cross.18$pheno$sex == 1)])
#cross.a.parent <- subset(cross.a.all, ind = is.na(cross.a.all$pheno$Pheno))

cross.b.all <- subset(cross.18, ind = cross.18$pheno$ID[which(cross.18$pheno$sex == 0)])
#cross.b.parent <- subset(cross.b.all, ind = is.na(cross.b.all$pheno$Pheno))

cross.a <- subset(cross.a.all, ind = !is.na(cross.a.all$pheno$Pheno))
cross.b <- subset(cross.b.all, ind = !is.na(cross.b.all$pheno$Pheno))

### Table before missing filter 2.943293e-02
gt.a <- geno.table(cross.a)
gt.b <- geno.table(cross.b)

#### Filter Conservative
drop.missing <- function(cross, M) {
  gt <- geno.table(cross)
  todrop <- rownames(gt[which(gt$missing > M), ])
  paste(length(todrop), "markers dropped")
  cross <- drop.markers(cross, unlist(todrop))
  return(cross)
}
## missing more than 35/46 filtered missing more than 33/43 filtered
cross.a <- drop.missing(cross.a, 8)
cross.b <- drop.missing(cross.b, 6)
gt.missing.a <- geno.table(cross.a)
gt.missing.b <- geno.table(cross.b)

## drop down to fixed in pars
#cross.a.fix <- drop.markers(cross.a,markernames(cross.a)[!which(markernames(cross.a)%in% parents.fixed)])
###head(geno.table(cross.a.fix),100)
#cross.b.fix <- drop.markers(cross.b,markernames(cross.b)[!which(markernames(cross.b)%in% parents.fixed)])

print("forming initial linkage groups to fix phase...")
grpRf <- 0.2
grpLod <- 4  ## Standard LG form LOD
finLod <- 5  ## Higher final NBH LOD
finRf <- 0.15

cross.test.a <- formLinkageGroups(cross.a, max.rf = grpRf, min.lod = grpLod, reorgMarkers = TRUE)
cross.test.b <- formLinkageGroups(cross.b, max.rf = grpRf, min.lod = grpLod, reorgMarkers = TRUE)

cross.test.a.linked <- subset(cross.test.a, chr=c(1,2,3))
cross.test.b.linked <- subset(cross.test.b, chr=c(1,2,3))

marks <- c(markernames(cross.test.a.linked),markernames(cross.test.b.linked))

cross.18 <- drop.markers(cross.18, markernames(cross.18)[!markernames(cross.18)%in% marks])

cross.18 <- formLinkageGroups(cross.18, max.rf = grpRf, min.lod = grpLod, reorgMarkers = TRUE)













cross.18



marks <- c(markernames(cross.test.a),markernames(cross.test.b))

cross.a.parent <- drop.markers(cross.a.parent, markernames(cross.a.parent)[!markernames(cross.a.parent)%in% marks])
cross.b.parent <- drop.markers(cross.b.parent, markernames(cross.b.parent)[!markernames(cross.b.parent)%in% marks])


cross.a.parent




geno.table(cross.test.a.linked)
geno.table(cross.test.b.linked)




save.image("/home/jmiller1/Sex_Map/dataset/parent_filtered_sexmap.rsave")

cross.test.a <- formLinkageGroups(cross.a, max.rf = finRf, min.lod = grpLod, reorgMarkers = TRUE)
cross.test.b <- formLinkageGroups(cross.b, max.rf = finRf, min.lod = grpLod, reorgMarkers = TRUE)



cross.test.a <- switchAlleles(cross.test.a, checkAlleles(cross.test.a, threshold=4)[,1])

a[]

##smaller group
cross.test.a.linked <- subset(cross.test.a, chr=1)
cross.test.b.linked <- subset(cross.test.b, chr=1)

cross.test.a <- formLinkageGroups(cross.test.a.linked, max.rf = 0.1, min.lod = grpLod, reorgMarkers = TRUE)
cross.test.b <- formLinkageGroups(cross.test.b.linked, max.rf = 0.1, min.lod = grpLod, reorgMarkers = TRUE)

save.image("/home/jmiller1/Sex_Map/dataset/parent_filtered_sexmap.rsave")


cross.test.a.linked <- subset(cross.test.a, chr=c(1,2,3))
cross.test.b.linked <- subset(cross.test.b, chr=c(1,2,3))

marks <- c(markernames(cross.test.a.linked),markernames(cross.test.b.linked))
marks.a <- markernames(cross.a.parent)
marks.b <- markernames(cross.b.parent)

cross.a.parent <- drop.markers(cross.a.parent,marks.a[!marks.a %in% marks])
cross.b.parent <- drop.markers(cross.b.parent,marks.b[!marks.b %in% marks])

geno.table(cross.a.parent)


cross.a.parent <- subset(cross.a.parent, ind = is.na(cross.a.all$pheno$Pheno))
cross.a.parent <- subset(cross.b.parent, ind = is.na(cross.a.all$pheno$Pheno))

get.scafs <- function(X){
  unique(gsub(":.*",'',markernames(cross.test.b)))
}
a <- get.scafs(cross.test.a)

comp.scafs <- function(X,Y){
  sum(get.scafs(X) %in% get.scafs(Y))
}











cross.test.a.linked <- subset(cross.test.a, chr=c(1,2))
cross.test.b.linked <- subset(cross.test.b, chr=c(1,2))




cross.test.a.linked <- geno.table(subset(cross.test.a, chr=1)

### Flip phase and retest





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



for pval

##gt.cross.par <- NA
##cutoff <- 0.01



head(geno.table(cross.b.fix),100)
NW_012225098.1


NW_012225098.1:36365  NW_012225098.1       5 19  8 11      0      0
NW_012225098.1:133639 NW_012225098.1       2 20  7 14      0      0
NW_012225098.1:133860 NW_012225098.1       5 19  8 11      0      0
NW_012225098.1:133870 NW_012225098.1       5 19  8 11


X <- rels(cross.18)
diag(X) <- 1
cols <- unlist(sapply(strsplit(rownames(X), "_"), "[[", 1))
labs <- unlist(sapply(strsplit(rownames(X), "_"), "[[", 2))

sx <- cross.18$pheno$sex
names(sx) <- cross.18$pheno$ID
sx[names(x)]

fit <- cmdscale(as.dist(1 - X), eig = TRUE, k = 2)
x <- fit$points[, 1]
y <- fit$points[, 2]



pdf("/home/jmiller1/public_html/chr5.pdf", width = 20, height = 20)
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Metric MDS for Chromosome 5 (sex chromosome, unfiltered snps)",
  type = "n")
text(x, y, labels = labs, col = brewer.pal(11, "Spectral")[as.factor(sx[names(x)])], cex = 2)
dev.off()


###gt.a[rownames(parents)[parents$AB==1],]
#pval.a <- log10(gt.a$P.value)
#pval.b <- log10(gt.b$P.value)
## do not use
##cross.a <- drop.markers(cross.a, rownames(gt.missing.a[gt.missing.a$P.value < cutoff,
##  ]))
##gt.pval.a <- geno.table(cross.a)
##sum(nmar(cross.a))

#cross.b <- drop.markers(cross.b, rownames(gt.missing.b[gt.missing.b$P.value < cutoff,
##  ]))
##gt.pval.b <- geno.table(cross.b)
##sum(nmar(cross.b))
