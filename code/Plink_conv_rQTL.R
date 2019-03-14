#!/bin/R
library(qtl)
source("/home/jmiller1/Sex_Map/PLINK2RQTL.f2.R")
setwd("/home/jmiller1/Sex_Map")
dir <- "/home/jmiller1/Sex_Map"

name <- "NBH_CHR5_filt_conv"

ped <- paste(dir, name, ".ped", sep = "")
map <- paste(dir, name, ".map", sep = "")
PLINKtoCSVR(ped = ped, map = map, out = paste(dir, name, ".parents.csvr", sep = ""))
