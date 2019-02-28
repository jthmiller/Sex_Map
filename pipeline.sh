### set base dir
basedir=/home/jmiller1/Sex_Map

## align to ncbi genome
sbatch -p low -t 48:00:00 --export=basedir="/home/jmiller1/Sex_Map" ${basedir}/align.sh

## calculate depth to find rad sites
sbatch -t 48:00:00 --mem=64000 --export=basedir="/home/jmiller1/Sex_Map" ${basedir}/depth.sh

## call genotypes
sbatch  -t 48:00:00 --mem=64000 --export=basedir="/home/jmiller1/Sex_Map" ${basedir}/callgt.sh
