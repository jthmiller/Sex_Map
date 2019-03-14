#!/bin/R

## Directories
basedir <- "/home/jmiller1/Sex_Map/"
indir <- "/home/jmiller1/Sex_Map/dataset/"

## Funtions for processing rQTL map data
source(file.path(basedir, "code/source_file.R"))
# source(file.path(basedir, 'rQTL/scripts/QTL_remap/QTL/model_source_file.R'))

## Libraries
flib <- "/share/apps/rmodules"
fpacks <- c("devtools", "httr", "RColorBrewer", "qtl")
lapply(fpacks, require, character.only = TRUE, lib.loc = flib)

mylib <- "/home/jmiller1/R/x86_64-pc-linux-gnu-library/3.5"
mpacks <- c("qtl", "foreach", "qtl2", "qtlTools", "doParallel", "plyr")
lapply(mpacks, require, character.only = TRUE, lib.loc = mylib)

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if (trace) 
      cat(nm, ":")
    source(file.path(path, nm), ...)
    if (trace) 
      cat("\n")
  }
}

sourceDir("doParallel/R")

### Phenotype translation
trsl.bin <- c(0, 0, 0, 1, 1, 1)
names(trsl.bin) <- as.character(0:5)

## Parameters for rQTL for population specific datasets (NBH markers require at
## least 70% genotypes )
if (pop == "NBH") {
  ns <- "North"
  confirmed = T
  reorder.marks <- F
  mapped.only = TRUE
  grpLod <- 10  ## Standard LG form LOD
  finLod <- 12  ## Higher final NBH LOD
  grpRf <- 0.3
  finRf <- 0.15
  cutoff <- 1e-05
  miss <- 12
  miss1 <- 10
  miss2 <- 8
  droppo <- 3
} else if (pop == "ELR") {
  ns <- "South"
  confirmed = T
  reorder.marks <- F
  mapped.only = TRUE
  missing <- 0.9
  grpLod <- 10  ## Standard LG form LOD
  finLod <- 12  ## Higher final ELR LOD
  grpRf <- 0.3
  finRf <- 0.15
  cutoff <- 1e-04  ## Higher, need more power to detect seg distortion
  miss <- 10
  miss1 <- 10
  miss2 <- 10
  droppo <- 15
} else if (pop == "NEW") {
  ns <- "North"
  confirmed = T
  reorder.marks <- F
  mapped.only = TRUE
  inds <- c(NA)  # determined to be dropped low cov
  missing <- 0.8
  grpLod <- 10  ## Standard LG form LOD
  finLod <- 12  ## Higher final ELR LOD
  grpRf <- 0.3
  finRf <- 0.15
  cutoff <- 1e-05
  miss <- 15
  miss1 <- 10
  miss2 <- 8
  droppo <- 2
} else if (pop == "BRP") {
  ns <- "North"
  confirmed = T
  reorder.marks <- F
  mapped.only = TRUE
  inds <- c(NA)  # determined to be dropped low cov
  missing <- 0.8
  grpLod <- 10  ## Standard LG form LOD
  finLod <- 12  ## Higher final ELR LOD
  grpRf <- 0.3
  finRf <- 0.15
  cutoff <- 1e-05
  miss <- 15
  miss1 <- 10
  miss2 <- 8
  droppo <- 2
}
