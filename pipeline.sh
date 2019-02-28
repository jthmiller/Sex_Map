### set base dir
basedir=/home/jmiller1/Sex_Map

## align to ncbi genome
sbatch -p low -t 48:00:00 --export=basedir="/home/jmiller1/Sex_Map" ${basedir}/align.sh

## calculate depth to find rad sites
sbatch -t 48:00:00 --mem=64000 --export=basedir="/home/jmiller1/Sex_Map" ${basedir}/depth.sh

## call genotypes
sbatch  -t 48:00:00 --mem=64000 --export=basedir="/home/jmiller1/Sex_Map" ${basedir}/callgt.sh


9180990_39       low array_jo jmiller1  R     1:00:56      1 6   2G     c9-95
9180990_40       low array_jo jmiller1  R     1:00:56      1 6   2G     c9-95
9180990_41       low array_jo jmiller1  R     1:00:56      1 6   2G     c9-95

9180990_42
