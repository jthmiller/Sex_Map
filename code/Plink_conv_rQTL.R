#!/bin/R
library(qtl)
source("/home/jmiller1/Sex_Map/code/PLINK2RQTL.f2.R")
setwd("/home/jmiller1/Sex_Map")

indir <- "/home/jmiller1/Sex_Map/plink_files"
outdir <- "/home/jmiller1/Sex_Map/dataset"
name <- "NBH_CHR5_filt_conv"

ped <- paste(indir, name, ".ped", sep = "")
map <- paste(indir, name, ".map", sep = "")
PLINKtoCSVR(ped = ped, map = map, out = paste(outdir, name, ".parents.csvr", sep = ""))
