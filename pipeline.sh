### set base dir
basedir=/home/jmiller1/Sex_Map

## align to ncbi genome
sbatch -p low -t 48:00:00 --export=basedir="/home/jmiller1/Sex_Map" ${basedir}/align.sh

## calculate depth to find rad sites
sbatch -p high -t 48:00:00 --mem=64000 --export=basedir="/home/jmiller1/Sex_Map" ${basedir}/depth.sh

## call genotypes
sbatch -t 48:00:00 --export=basedir="/home/jmiller1/Sex_Map" ${basedir}/callgt.sh

REMOVE EMPTY VCFs

## get all chr5 and unmapped contigs
grep '^chain' ${basedir}/metadata/meiotic_map.chainfile | cut -f3,8 | less -S | sort | uniq > ${basedir}/metadata/Scaf_chrm.txt
grep 'chr5' ${basedir}/metadata/Scaf_chrm.txt | cut -f1 > ${basedir}/metadata/Scaf_chrm5_only.txt
grep -v 'chr' ${basedir}/metadata/Scaf_chrm.txt | cut -f1 >> ${basedir}/metadata/Scaf_chrm5_only.txt
find "${basedir}/genotypes" -type f -size +10k -name '*.vcf' |
grep -f ${basedir}/metadata/Scaf_chrm5_only.txt > ${basedir}/metadata/Scafolds_to_map.txt

while read F
do
	dest=$(basename $F)
	bgzip -@ 6 -f ${F} -c > ${basedir}/chr5_bcfs/"${dest%\.vcf}".bcf
	bcftools index -f ${basedir}/chr5_bcfs/"${dest%\.vcf}".bcf --threads 6
done < "${basedir}/metadata/Scafolds_to_map.txt"


find ${basedir}/chr5_bcfs/ -type f -size +10k -name '*1.bcf' > ${basedir}/metadata/filelist.txt
bcftools concat -D -a --file-list ${basedir}/metadata/filelist.txt -O b -o ${basedir}/chr5_bcfs/NBH_CHR5.bcf
bcftools sort ${basedir}/chr5_bcfs/NBH_CHR5.bcf -O b -o ${basedir}/chr5_bcfs/NBH_CHR5_sorted.bcf
bcftools index ${basedir}/chr5_bcfs/NBH_CHR5_sorted.bcf
bcftools stats --depth 10,5000,10 ${basedir}/chr5_bcfs/NBH_CHR5_sorted.bcf > NBH_CHR5_unfilt.depthstats

bcftools view -m2 -M2 -v snps -q 0.1 -Q 0.9 -e 'INFO/AC<10 & FMT/DP<10 & FMT/DP>1000 & INFO/MQM<30 & INFO/NS<10' \
	${basedir}/chr5_bcfs/NBH_CHR5_sorted.bcf -O v -o ${basedir}/chr5_bcfs/NBH_CHR5_filt.vcf
bcftools stats --depth 1,1000,5 ${basedir}/chr5_bcfs/NBH_CHR5_filt.vcf > NBH_CHR5_filt.depthstats
