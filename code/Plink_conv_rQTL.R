#!/bin/R
library(qtl)
source("/home/jmiller1/Sex_Map/PLINK2RQTL.f2.R")
setwd("/home/jmiller1/Sex_Map")
dir <- "/home/jmiller1/Sex_Map"

name <- ""

ped <- paste(dir, name, ".filt.pk.recode.ped", sep = "")
map <- paste(dir, name, ".filt.pk.recode.map", sep = "")
PLINKtoCSVR(ped = ped, map = map, out = paste(dir, name, ".parents.csvr", sep = ""))
