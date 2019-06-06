### Map QTLs 1 of 3
pop <- "NBH"
source("/home/jmiller1/Sex_Map/code/control_file.R")

## Due to dropout, ~6/92 genotypes are wrong, error rate of 6.5%

## For plotting marker_dens <- list()

# Table of Chroms with sig QTLs test.QTLs <- read.table(file.path(basedir,
# 'rQTL/metadata/QTLs.txt'), sep = '\t', header = T)

## Get chrom number vector test.QTLs$chrm.n <- gsub('chr', '', test.QTLs$chrom)

# print(paste(pop, X, sep = ' '))
## read in the QTL cross
cross.18 <- read.cross.jm(file = file.path(indir, "NBH_CHR5_filt_conv.csvr"), format = "csvr",
  geno = c(1:3), estimate.map = FALSE)

### Pull names from plinkfile
#path <- file.path(plinkdir, "NBH_CHR5_filt_conv.ped")
#popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
#indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
#newname <- paste(popname, indname, sep = "_")
#cross.18$pheno$ID <- paste(popname, indname, sep = "_")

## Remove problematic individuals (found by kinship analysis)
popdir <- file.path("/home/jmiller1/QTL_Map_Raw/popgen", "rQTL", pop, "REMAPS")
con <- file(file.path(popdir, "kinship.keep.ind.txt"), open = "r")
keepers <- readLines(con)
close(con)

print("Dropping kinship outliers")
cross.18 <- subset(cross.18, ind = cross.18$pheno$ID %in% keepers)


#### Add sex of individuals from kinship analysis
sex <- read.table("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/data/sex.txt",stringsAsFactors=F)
rownames(sex) <- sex$ID
sex.vec <- sex[as.character(cross.18$pheno$ID), 'sex']
cross.18$pheno$sex <- sex.vec

#### SUBSET to two different cross objects by sex
cross.male <- subset(cross.18, ind = cross.18$pheno$ID[which(cross.18$pheno$sex == 1)])
cross.female <- subset(cross.18, ind = cross.18$pheno$ID[which(cross.18$pheno$sex == 0)])
parent.both <- subset(cross.18, ind = cross.18$pheno$ID[is.na(cross.18$pheno$Pheno)])

### Female can be mapped by regular filtering 1:2:1
parent.female <- subset(cross.female, ind = cross.female$pheno$ID[is.na(cross.male$pheno$Pheno)])
parent.male <- subset(cross.male, ind = cross.male$pheno$ID[is.na(cross.male$pheno$Pheno)])

### Drop parents from offspring cross
cross.male <- subset(cross.male, ind = !is.na(cross.male$pheno$Phen))
cross.female <- subset(cross.female, ind = !is.na(cross.female$pheno$Phen))

## Extra subsetted frame (after mapping loci for non-missing data, determine what is happening
## with dropped markers
cross.a5 <- subset(cross.18, ind = cross.18$pheno$ID[which(cross.18$pheno$sex == 1)])
cross.b5 <- subset(cross.18, ind = cross.18$pheno$ID[which(cross.18$pheno$sex == 0)])

### Table before missing filter
gt.a <- geno.table(cross.male)
gt.b <- geno.table(cross.female)

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
## ~75% genotypes per locus

cross.male <- drop.missing(cross.male, 5)
cross.female <- drop.missing(cross.female, 5)
gt.missing.a <- geno.table(cross.male)
gt.missing.b <- geno.table(cross.female)
sum(nmar(cross.male))
sum(nmar(cross.female))

################################################################################
################################################################################
################################################################################




################################################################################
### CHR 5 ONLY #################################################################
################################################################################

### Only CHR5
chr5 <- read.table("/home/jmiller1/Sex_Map/metadata/chrm5_only.txt", header=F)
chr5 <- as.character(unlist(chr5$V1))

### those that have markers in the map
mapping <- unique(gsub(":.*",'',markernames(cross.female)))
chr5.mapped <- chr5[chr5 %in% mapping]

###cross.18.mapped.scafs <- subset(cross.18,chr)
################################################################################
### CHR 5 ONLY #################################################################
################################################################################

cross.male.RF08.LOD5 <- formLinkageGroups(cross.male, max.rf = 0.1, min.lod = 5,reorgMarkers = TRUE)

save.image('NOAHS')

gt.male.RF08.LOD5 <- geno.table(cross.male.RF08.LOD5)

save.image('NOAHS')

index <- gsub(":.*",'',rownames(gt.male.RF08.LOD5)) %in% chr5

save.image('NOAHS')


## 1 has most markers
## 23 has many AA/AB markers

## 26 is the deleted regions
## 25 is also the deleted regions

#####################################################################
male.deleted <- subset(cross.male.RF08.LOD5, chr = c(25,26))

male.AB <- subset(cross.male.RF08.LOD5, chr = c(1,23))

male.sex.linked <- subset(cross.male.RF08.LOD5, chr = c(1,23,25,26))
#####################################################################


#####################################################################
male.sex.linked <- switchAlleles(male.sex.linked,markernames(male.sex.linked, chr=c(23,25))

male.sex.linked <- formLinkageGroups(male.sex.linked, max.rf = 0.08, min.lod = 5,reorgMarkers = TRUE)

male.sex.linked <- subset(male.sex.linked, chr = keep)

male.sex.linked <- orderMarkers(male.sex.linked, window = 5, use.ripple = T, error.prob = 0.15,
  map.function = "kosambi", sex.sp = F, maxit = 1000, tol = 0.01)

POS.map.us <- est.map(male.sex.linked, error.prob = 0.15, map.function = "kosambi", chr = 1,maxit = 1000)

male.sex.linked <- replace.map(male.sex.linked, POS.map.us)
#####################################################################


#####################################################################
male.AB <- switchAlleles(male.AB,markernames(male.AB, chr=23))

male.AB <- formLinkageGroups(male.AB, max.rf = 0.05, min.lod = 5,reorgMarkers = TRUE)

gt.male.AB <- geno.table(male.AB)

index <- gsub(":.*",'',rownames(gt.male.AB)) %in% chr5
#####################################################################


##tail(sort(table(gt.male.AB[index,'chr'])),8)
##tail(sort(table(gt.male.AB[,'chr'])),8)

###AAxAB
gt.male.AB[gt.male.AB$chr=='3',]
gt.male.AB[gt.male.AB$chr=='7',]
gt.male.AB[gt.male.AB$chr=='34',]
### Switch AB
gt.male.AB[gt.male.AB$chr=='2',]
gt.male.AB[gt.male.AB$chr=='4',]
gt.male.AB[gt.male.AB$chr=='16',]
gt.male.AB[gt.male.AB$chr=='37',]
gt.male.AB[gt.male.AB$chr=='62',]

### All AB
gt.male.AB[gt.male.AB$chr=='1',]

#male.AB.sw <- switchAlleles(male.AB,markernames(male.AB, chr=c(2,4,16,37,62)))
##male.AB.sw <- subset(male.AB, chr = c(3,7,2,4,16,34,37,62))


gt.male.AB <- geno.table(male.AB)
chrs.map <- names(table(gt.male.AB[index,'chr'])[table(gt.male.AB[index,'chr']) > 0])
male.AB.sw <- subset(male.AB, chr = chrs.map[-1]) ### ~ 500 scaffolds
gt.male.AB.sw  <- geno.table(male.AB.sw)
gt.male.AB.sw [gt.male.AB.sw$chr=='4',]


male.AB.sw.lg <- switchAlleles(male.AB.sw,markernames(male.AB.sw, chr=c(2,4)))
test <- formLinkageGroups(male.AB.sw.lg , max.rf = 0.1, min.lod = 5,reorgMarkers = TRUE)
gt.test <- geno.table(test)
gt.test[gt.test$chr=='3',]
test <- switchAlleles(test,markernames(test, chr=3))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='4',]
test <- switchAlleles(test,markernames(test, chr=4))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='5',]
test <- switchAlleles(test,markernames(test, chr=5))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='7',]
test <- drop.markers(test,markernames(test, chr=7))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='8',]
test <- switchAlleles(test,markernames(test, chr=8))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='9',]
test <- switchAlleles(test,markernames(test, chr=9))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='12',]
test <- switchAlleles(test,markernames(test, chr=12))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='13',]
test <- switchAlleles(test,markernames(test, chr=13))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='14',]
test <- switchAlleles(test,markernames(test, chr=14))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='15',]
test <- drop.markers(test,markernames(test, chr=15))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='16',]
test <- drop.markers(test,markernames(test, chr=16))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='17',]
test <- drop.markers(test,markernames(test, chr=17))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='18',]
test <- drop.markers(test,markernames(test, chr=18))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='19',]
test <- drop.markers(test,markernames(test, chr=19))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='21',]
test <- switchAlleles(test,markernames(test, chr=21))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='22',]
test <- drop.markers(test,markernames(test, chr=22))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='23',]
test <- drop.markers(test,markernames(test, chr=23))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='25',]
test <- switchAlleles(test,markernames(test, chr=25))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='26',]
test <- switchAlleles(test,markernames(test, chr=26))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='27',]
test <- drop.markers(test,markernames(test, chr=27))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='28',]
test <- drop.markers(test,markernames(test, chr=28))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='29',]
test <- drop.markers(test,markernames(test, chr=29))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='30',]
test <- switchAlleles(test,markernames(test, chr=30))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='31',]
test <- switchAlleles(test,markernames(test, chr=31))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='32',]
test <- drop.markers(test,markernames(test, chr=32))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='32',]
test <- switchAlleles(test,markernames(test, chr=33))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='34',]
test <- drop.markers(test,markernames(test, chr=34))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='35',]
test <- drop.markers(test,markernames(test, chr=35))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='36',]
test <- switchAlleles(test,markernames(test, chr=36))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='37',]
test <- switchAlleles(test,markernames(test, chr=37))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='37',]
test <- drop.markers(test,markernames(test, chr=38))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='39',]
test <- drop.markers(test,markernames(test, chr=39))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='40',]
test <- switchAlleles(test,markernames(test, chr=40))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='42',]
test <- drop.markers(test,markernames(test, chr=42))
gt.test <- geno.table(test)
gt.test[gt.test$chr=='43',]

#sapply(as.character(c(75:84)),function(x){
# print(x)
#print(gt.test[gt.test$chr==x,])
#})


test <- drop.markers(test,markernames(test, chr=c(43,45,46,47,48,49,50,52,53,54,55,56,58,59,60,62,64,65,66,67,68,69,70,71,72,73,75,76,77,78,79,80,82)))
test <- switchAlleles(test,markernames(test, chr=c(57,61,63,74,81)))
gt.test <- geno.table(test)

test2 <- subset(test, chr = names(table(gt.test$chr)[table(gt.test$chr) > 1]))
test2.relaxed <- formLinkageGroups(test2 , max.rf = 0.3, min.lod = 3,reorgMarkers = TRUE)

save.image('NOAHS')

gt.test2 <- geno.table(test2.relaxed)
gt.test2[gt.test2$chr=='1',]

test2.relaxed <- orderMarkers(test2.relaxed, window = 5, use.ripple = T, error.prob = 0.2,
  map.function = "kosambi", sex.sp = F, maxit = 1000, tol = 0.01)

POS.map.us <- est.map(test2.relaxed, error.prob = 0.2, map.function = "kosambi", chr = 1,maxit = 1000)

test2.relaxed <- replace.map(test2.relaxed, POS.map.us)

save.image('NOAHS')

write.cross(test2.relaxed, "/home/jmiller1/Sex_Map/dataset/test2.relaxed",format="tidy")


##write.table(cbind(sort(table(gsub(':.*','',rownames(gt.test2))))), 'test2.relaxed')



#test2.hi <- formLinkageGroups(test2 , max.rf = 0.1, min.lod = 5,reorgMarkers = TRUE)
#test.relaxed <- formLinkageGroups(test , max.rf = 0.3, min.lod = 3,reorgMarkers = TRUE)


male.AB.bk <- subset(male.AB, chr = c(3,7,2,4,16,34,37,62))
male.allhets <- subset(male.AB, chr = 1)




male.AB.sw.lg <- formLinkageGroups(male.AB.sw , max.rf = 0.3, min.lod = 3,reorgMarkers = TRUE)

"
gt.male.sw.lg <- geno.table(male.AB.sw.lg)

male.AB.sw.lg <- switchAlleles(male.AB.sw.lg,markernames(male.AB.sw.lg, chr=1))

male.AB.sw.lg <- formLinkageGroups(male.AB.sw.lg, max.rf = 0.3, min.lod = 3,reorgMarkers = TRUE)

male.AB.sw.lg <- orderMarkers(male.AB.sw.lg, window = 5, use.ripple = T, error.prob = 0.2,
  map.function = "kosambi", sex.sp = F, maxit = 1000, tol = 0.01)

POS.map.us <- est.map(male.AB.sw.lg, error.prob = 0.15, map.function = "kosambi", chr = 1,maxit = 1000)

male.AB.sw.lg <- replace.map(male.AB.sw.lg, POS.map.us)

save.image('NOAHS')



write.cross(male.AB.sw.lg, "/home/jmiller1/Sex_Map/dataset/map.chr5.more.scaffs",format="tidy")

write.csv(gt.missing.a,"/home/jmiller1/Sex_Map/dataset/male.gtable")

####### Simple map

yourdir <- "/home/jmiller1/Sex_Map/dataset/"

map <- read.csv(file.path(yourdir,'map.chr5.more.scaffs_map.csv'),stringsAsFactors=F)

gts <- read.csv(file.path(yourdir,'male.gtable'),stringsAsFactors=F,row.names=1)

## Scaffold_ID scaffold_position LG genetic_position

Scaffold_ID <- gsub(':.*','',map$X)

scaffold_position <- as.numeric(gsub('.*:','',map$X))

map <- data.frame(Scaffold_ID,scaffold_position,chr=5,pos=as.numeric(map$pos))

av.pos <- sapply(unique(map$Scaffold_ID),
            function(X){mean(map[map$Scaffold_ID==X,4])})

orient <- sapply(unique(map$Scaffold_ID),
            function(X){ cor( map[map$Scaffold_ID==X,4] ,map[map$Scaffold_ID==X,2])})

dir <- ifelse(orient > 0,"+","-")

simp.map <- data.frame(Scaffold_ID=unique(map$Scaffold_ID),av.pos=as.numeric(av.pos),orient=as.numeric(orient),dir,row.names=unique(map$Scaffold_ID))

simp.map <- simp.map[order(simp.map$av.pos),]


### Total scaffs in MALE

male <- table(gts$chr)
male.mapped <- table(map$Scaffold_ID)
m <- names(male.mapped)

evid <- data.frame(all.loci=as.numeric(male[m]), num.mapped=as.numeric(male.mapped[m]),row.names=m)
evid$proportion_scaff_linked <- as.numeric(evid$num.mapped/evid$all.loci)

simp.map$proportion_scaff_linked <- evid[rownames(simp.map),3]

write.csv(simp.map,file.path(yourdir,'simp.map'))












































map <- read.csv('/home/jmiller1/genomes_jm/mapped/scripts/chr5_map.csv')



'screen here

screen -dr  7971.pts-97.farm

male.AB <- subset(male.AB, chr =    )

male.sex.linked <- orderMarkers(male.sex.linked, window = 5, use.ripple = T, error.prob = 0.15,
  map.function = "kosambi", sex.sp = F, maxit = 1000, tol = 0.01)

POS.map.us <- est.map(male.sex.linked, error.prob = 0.15, map.function = "kosambi", chr = 1,maxit = 1000)

male.sex.linked <- replace.map(male.sex.linked, POS.map.us)

















gt.male.RF08.LOD5[gt.male.RF08.LOD5$chr==1,  ]

keep <- unique(gt.male.RF08.LOD5[index, 'chr' ])

gt.male.RF08.LOD5 <- subset(gt.male.RF08.LOD5, chr = keep)

gt.male.RF08.LOD5 <- orderMarkers(gt.male.RF08.LOD5, window = 5, use.ripple = T, error.prob = 0.15,
  map.function = "kosambi", sex.sp = F, maxit = 1000, tol = 0.01)

POS.map.us <- est.map(gt.male.RF08.LOD5, error.prob = 0.15, map.function = "kosambi", chr = 1,maxit = 1000)

gt.male.RF08.LOD5 <- replace.map(gt.male.RF08.LOD5, POS.map.us)

"screen is here

################################################################################
################################################################################














##parent cross types
parent.female.5 <- drop.markers(parent.female,markernames(parent.female)[!gsub(":.*",'',markernames(parent.female)) %in% chr5.mapped])
parent.male.5 <- drop.markers(parent.male,markernames(parent.male)[!gsub(":.*",'',markernames(parent.male)) %in% chr5.mapped])
parent.both.5 <- drop.markers(parent.both,markernames(parent.both)[!gsub(":.*",'',markernames(parent.both)) %in% chr5.mapped])

gt.parents.both <- geno.table(parent.both)

gt.parents.mal <- geno.table(parent.male.5)
gt.parents.fem <- geno.table(parent.female.5)

gt.parents.mal.not.het <- rownames(gt.parents.mal)[!gt.parents.mal$AB==1]
gt.parents.mal.hom <- rownames(gt.parents.mal)[!gt.parents.mal$AB==1]

gt.parents.fem.marks <- rownames(gt.parents.fem)[gt.parents.fem$AB==1]

marks.female.mappable <- gt.parents.fem.marks[gt.parents.fem.marks %in% gt.parents.mal.not.het]


### Regions with 1:0:1 seg in males

## Look for high link in males with the 1:0:1 region
cross.male.5 <- drop.markers(cross.male,markernames(cross.male)[!gsub(":.*",'',markernames(cross.male)) %in% chr5.mapped])
cross.female.5 <- drop.markers(cross.female,markernames(cross.female)[!gsub(":.*",'',markernames(cross.female)) %in% chr5.mapped])

switch male GP BB to AA
mal.BB <- rownames(gt.parents.mal)[gt.parents.mal$BB==1]
cross.test.male <- switchAlleles(cross.male.5,mal.BB)
cross.test.male <- formLinkageGroups(cross.test.male, max.rf = 0.1, min.lod = 3,reorgMarkers = TRUE)

cross.test.5.unswitched <- subset(cross.test.male, chr = c(2,3)) ####  AA:2.7  AB:94.5  BB:2.8 (522 markers)

cross.test.5.unswitched <- orderMarkers(cross.test.5.unswitched, window = 5, use.ripple = T, error.prob = 0.15,
  map.function = "kosambi", sex.sp = F, maxit = 1000, tol = 0.01)
POS.map.us <- est.map(cross.test.5.unswitched, error.prob = 0.15, map.function = "kosambi", chr = 1,maxit = 1000)
cross.males.1_0_1 <- replace.map(cross.test.5.unswitched, POS.map.us)

write.cross(cross.males.1_0_1, "/home/jmiller1/Sex_Map/dataset/map.chr5",format="tidy")











###### Try to switch phase
cross.test.male.switch <- switchAlleles(cross.test.male,markernames(cross.test.male,chr=2))
cross.test.male.switch <- formLinkageGroups(cross.test.male.switch, max.rf = 0.1, min.lod = 6,reorgMarkers = TRUE)
cross.test.male.switch <- subset(cross.test.male.switch, chr = c(2)) #### AA:48.2  AB:49.2  BB:2.5 (122 markers)

cross.males.1_0_1 <- orderMarkers(cross.test.male.switch, chr = 1, window = 5, use.ripple = T, error.prob = 0.05,
  map.function = "kosambi", sex.sp = F, maxit = 1000, tol = 0.01)

POS.map <- est.map(cross.test.male.switch, error.prob = 0.01, map.function = "kosambi", chr = 1,maxit = 1000)
cross.test.male.switch <- replace.map(cross.test.male.switch, POS.map)

##save.image('NOAHS')
##load('NOAHS')




##### FEMALE
cross.test.female <- formLinkageGroups(cross.female.5, max.rf = 0.05, min.lod = 3,reorgMarkers = TRUE)


male.5.8 <-
head(geno.table(cross.male.5,chr=1)[,1:5],200)



##save.image('NOAHS')
##load('NOAHS')






### Low reco rate
cross.test.male <- formLinkageGroups(cross.male.5, max.rf = 0.05, min.lod = 3,reorgMarkers = TRUE)

cross.males.1_0_1 <- subset(cross.test.male, chr = c(17,20,22,23,24,41))
cross.males.1_0_1 <- formLinkageGroups(cross.males.1_0_1, max.rf = 0.25, min.lod = 2,reorgMarkers = TRUE)

switching <- markernames(cross.males.1_0_1, chr = chrnames(cross.males.1_0_1)[1])

### Switch in testcross
cross.males.1_0_1 <- switchAlleles(cross.males.1_0_1, switching)
### Must also switch alleles in fem parent, and others to compare
cross.female.5 <- switchAlleles(cross.female.5, switching)
cross.male.5 <- switchAlleles(cross.male.5, switching)

parent.male.5 <- switchAlleles(parent.male.5, switching)
parent.female.5 <- switchAlleles(parent.female.5, switching)

gt.parents.mal <- geno.table(parent.male.5)
gt.parents.fem <- geno.table(parent.female.5)
gt.female.5 <- geno.table(cross.female.5)
gt.male.5 <- geno.table(cross.male.5)

## Link mcr1 genes and map region to begin building the X in males
cross.males.1_0_1 <- formLinkageGroups(cross.males.1_0_1, max.rf = 0.1, min.lod = 3,reorgMarkers = TRUE)

cross.males.1_0_1 <- orderMarkers(cross.males.1_0_1, chr = 1, window = 5, use.ripple = T, error.prob = 0.05,
  map.function = "kosambi", sex.sp = F, maxit = 1000, tol = 0.01)

POS.map.males.1_0_1 <- est.map(cross.males.1_0_1, error.prob = 0.01, map.function = "kosambi", chr = 1,maxit = 1000)
cross.males.1_0_1 <- replace.map(cross.males.1_0_1, POS.map.males.1_0_1)


### Which markers are in linkage with the odd segregation ratios in males
mcr1 <- c('NW_012224869.1:3233519','NW_012224869.1:3233364','NW_012224869.1:3233460','NW_012224869.1:2090106','NW_012224869.1:2996642')
noahs <- c('NW_012224575.1','NW_012224869.1')

## NW_012224575.1:1809000 1816000
## NW_012224575.1:1880000 1909000
## NW_012224869.1:227000 232000

gt.male.5[mcr1,]
gt.female.5[mcr1,]
gt.parents.mal[mcr1,]
gt.parents.fem[mcr1,]

mapping <-
gt.male.5[gsub(":.*",'',rownames(gt.male.5)) %in% noahs,]


##try to map in males
cross.males.5.map <- formLinkageGroups(cross.male.5, max.rf = 0.1, min.lod = 5,reorgMarkers = TRUE)
gt.cross.male.map <- geno.table(cross.males.5.map)
geno.table(cross.males.5.map,chr=7)[mcr1,]
### chr7 contains mc1r regions


cross.males.5.mc1r <- subset(cross.males.5.map, chr = 7)
geno.table(cross.males.5.map,chr=7)




cross.males.5.relax <- formLinkageGroups(cross.male.5, max.rf = 0.13, min.lod = 3,reorgMarkers = TRUE)
geno.table(cross.males.5.relax)[mcr1,] ## on chr 1
cross.males.5.relax <- subset(cross.males.5.relax, chr = 2) ## subset to map

cross.males.5.relax <- orderMarkers(cross.males.5.relax, chr = 1, window = 5, use.ripple = T, error.prob = 0.05,
  map.function = "kosambi", sex.sp = F, maxit = 1000, tol = 0.01)
POS.map.males.relax <- est.map(cross.males.5.relax, error.prob = 0.01, map.function = "kosambi", chr = 1,maxit = 1000)
cross.males.5.relax <- replace.map(cross.males.5.relax, POS.map.males.relax)




cross.males.5.relax <- switchAlleles(cross.males.5.relax, markernames(cross.males.5.relax,chr=6))
cross.males.5.relax <- formLinkageGroups(cross.males.5.relax, max.rf = 0.15, min.lod = 3,reorgMarkers = TRUE)
geno.table(cross.males.5.relax)[mcr1,]

## Map chr1 with mcr1 genes on it
cross.males.5.relax <- orderMarkers(cross.males.5.relax, chr = 1, window = 5, use.ripple = T, error.prob = 0.05,
  map.function = "kosambi", sex.sp = F, maxit = 1000, tol = 0.01)



POS.map.males.relax <- est.map(cross.males.5.relax, error.prob = 0.01, map.function = "kosambi", chr = 1,maxit = 1000)
cross.males.5.relax <- replace.map(cross.males.5.relax, POS.map.males.relax)


##### Switch markers at the begining to see if there are other scaffolds that make the transistion make sense
cross.male.all.switch <- switchAlleles(cross.male, switching)
cross.male.all.switch <- formLinkageGroups(cross.male.all.switch, max.rf = 0.1, min.lod = 3,reorgMarkers = TRUE)
geno.table(cross.males.5.map)[mcr1,]










## Female recombination (no seg distortion among females.)
hets.female.none.male <- drop.markers(cross.female,markernames(cross.female)[!markernames(cross.female) %in% marks.female.mappable])
nothets.male.none.male <- drop.markers(cross.male,markernames(cross.male)[!markernames(cross.male) %in% marks.female.mappable])

geno.table(nothets.male.none.male)

geno.table(hets.female.none.male)









In male offspring, X will be from mother meiosis.
In female offspring, one X is low recomb. from male
the other, is recombining.

option 1, Use high linkage, to pull out the


###before filters
save.image("/home/jmiller1/Sex_Map/dataset/sexmap_scratch.rsave")


A A    B B
A A    A B
A B    A B
B B    A B

 M  F
 AA BB
 BB AA


### cross.female <- drop.markers(cross.female, rownames(gt.missing.b[gt.missing.b$P.value < cutoff,
##  ]))
## gt.pval.b <- geno.table(cross.female)
## sum(nmar(cross.female)) ### 2024 in females





gt.cross.par <- NA

## Drop by segregation distortion. Skip (instead by parent genotypes)
##cutoff <- 0.1
##cross.male <- drop.markers(cross.male, rownames(gt.missing.a[gt.missing.a$P.value < cutoff,
##  ]))
## gt.pval.a <- geno.table(cross.male)
### sum(nmar(cross.male)) ### 2378 in males

## Get parents and switch to correct phase







print("forming initial linkage groups to fix phase...")
grpRf <- 0.3
grpLod <- 4  ## Standard LG form LOD
finLod <- 5  ## Higher final NBH LOD
finRf <- 0.15

save.image("/home/jmiller1/Sex_Map/datasetsexmap_scratch.rsave")
cross.test.male <- formLinkageGroups(cross.male, max.rf = grpRf, min.lod = grpLod, reorgMarkers = TRUE)
### Flip phase and retest
cross.test.male <- switchAlleles(cross.test.male, markernames(cross.test.male, chr = chrnames(cross.test.male)[1]))
cross.test.male <- formLinkageGroups(cross.test.male, max.rf = grpRf, min.lod = grpLod,
  reorgMarkers = TRUE)
cross.test.male <- switchAlleles(cross.test.male, markernames(cross.test.male, chr = chrnames(cross.test.male)[1]))
cross.test.male <- formLinkageGroups(cross.test.male, max.rf = grpRf, min.lod = grpLod,
  reorgMarkers = TRUE)
cross.test.male <- subset(cross.test.male, chr = which.max(nmar(cross.test.male)))
cross.test.male <- orderMarkers(cross.test.male, chr = 1, window = 5, use.ripple = T, error.prob = 0.05,
  map.function = "kosambi", sex.sp = F, maxit = 1000, tol = 0.01)
cross.test.male <- removeDoubleXO(cross.test.male, verbose = T)
POS.map.a <- est.map(cross.test.male, error.prob = 0.01, map.function = "kosambi", chr = 1,
  maxit = 1000)
cross.test.male <- replace.map(cross.test.a, POS.map.a)
save.image("/home/jmiller1/Sex_Map/dataset/sexmap.rsave")

cross.test.female <- formLinkageGroups(cross.female, max.rf = grpRf, min.lod = grpLod, reorgMarkers = TRUE)
### Flip phase a couple times to see if additional markers fall into the largest LG
### and retest
cross.test.female <- switchAlleles(cross.test.female, markernames(cross.test.female, chr = chrnames(cross.test.female)[1]))
cross.test.female <- formLinkageGroups(cross.test.female, max.rf = grpRf, min.lod = grpLod,
  reorgMarkers = TRUE)
cross.test.female <- switchAlleles(cross.test.female, markernames(cross.test.female, chr = chrnames(cross.test.female)[1]))
cross.test.female <- formLinkageGroups(cross.test.female, max.rf = grpRf, min.lod = grpLod,
  reorgMarkers = TRUE)
cross.test.female <- subset(cross.test.female, ind = cross.test.female$pheno$ID[which(cross.test.female$pheno$sex ==
  0)])
cross.test.female <- subset(cross.test.female, chr = which.max(nmar(cross.test.female)))
cross.test.female <- orderMarkers(cross.test.female, chr = 1, window = 5, use.ripple = T, error.prob = 0.05,
  map.function = "kosambi", sex.sp = F, maxit = 1000, tol = 0.01)
cross.test.female <- removeDoubleXO(cross.test.female, verbose = T)
POS.map.b <- est.map(cross.test.female, error.prob = 0.01, map.function = "kosambi", chr = 1,
  maxit = 1000)
cross.test.female <- replace.map(cross.test.female, POS.map.b)
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



cross.test.female <- formLinkageGroups(chr5.mapped.only.b, max.rf = grpRf, min.lod = grpLod, reorgMarkers = TRUE)
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





map <- read.csv('/home/jmiller1/genomes_jm/mapped/scripts/chr5_map.csv')

av.pos <- sapply(unique(map$Scaffold_ID),function(X){mean(map[map$Scaffold_ID==X,4])})
orient <- sapply(unique(map$Scaffold_ID),function(X){ cor( map[map$Scaffold_ID==X,4] ,map[map$Scaffold_ID==X,2] )})
dir <- ifelse(orient > 0,"+","-")

cbind(unique(as.character(map$Scaffold_ID))[order(av.pos)],dir)
