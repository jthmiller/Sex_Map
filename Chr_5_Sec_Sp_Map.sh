#!/bin/bash -l
#SBATCH -o /home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/array_error_out/scan/out_%x_%A_%a.txt
#SBATCH -e /home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/array_error_out/scan/err_%x_%A_%a.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=QTL.Remap
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4000
#SBATCH -p high

####QTLs are on chrm '1 2 8 13 18 24'

scriptdir='/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap'

pop=$1

Rscript $scriptdir/QTL/18.8.22.qtl_SCAN.R --vanilla $pop $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
