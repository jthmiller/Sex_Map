#!/bin/bash -l

my_bamtools=/home/jmiller1/bin/bamtools-master/bin/bamtools
my_samtools=/home/jmiller1/bin/samtools-1.3/samtools
bamlist=/home/jmiller1/QTL_Map_Raw/align/bam.list.txt
BT=/home/nreid/bin/bedtools/bin/bedtools2

###ls ${basedir}/alignments/*.bam | xargs -n1 -P5 $my_samtools index
ls ${basedir}/alignments/*.bam > ${basedir}/metadata/bamlist.txt

### Previously determined that depth should be above 100 for a decent gt, and less than 150000
### find all bases that meet criteria

$my_bamtools merge -list ${basedir}/metadata/bamlist.txt |
$my_bamtools filter -in stdin -mapQuality ">30" -isProperPair true |
$my_samtools depth -m 50 -d 15000 /dev/stdin |
gzip -c > ${basedir}/metadata/depth.gz

zcat ${basedir}/metadata/depth.gz | awk '$3>384' |
awk '{OFS="\t"}{print $1,$2-1,$2,$3}' |
$BT merge -i stdin -d 10 -c 4 -o mean,max |
awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$3-$2}' |
awk '$6>$min' | awk '$5<$max' > ${basedir}/metadata/radsites.bed


### try all at once
#$my_bamtools merge -list ${basedir}/metadata/bamlist.txt |
#$my_bamtools filter -in stdin -mapQuality ">30" -isProperPair true |
#$my_samtools depth -m 50 -d 15000 /dev/stdin |
#awk '$3>384' |
#awk '{OFS="\t"}{print $1,$2-1,$2,$3}' |
#$BT merge -i stdin -d 10 -c 4 -o mean,max |
#awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$3-$2}' |
#awk '$6>$min' |
#awk '$5<$max' > ${basedir}/metadata/radsites.bed
