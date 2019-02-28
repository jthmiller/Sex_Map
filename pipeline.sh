### set base dir
basedir=/home/jmiller1/sex_chrom

## align to ncbi genome
sbatch -p low -t 48:00:00 --export=basedir="/home/jmiller1/Sex_Map" ${basedir}/align.sh

## calculate depth to find rad sites
sbatch -t 48:00:00 --mem=60000 depth.sh

sbatch  -t 48:00:00 --mem=60000 callgt.sh
