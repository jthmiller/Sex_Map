### set base dir
basedir=/home/jmiller1/Sex_Map

## align to ncbi genome
sbatch -p low -t 48:00:00 --export=basedir="/home/jmiller1/Sex_Map" ${basedir}/code/align.sh

## calculate depth to find rad sites
sbatch -p high -t 48:00:00 --mem=64000 --export=basedir="/home/jmiller1/Sex_Map" ${basedir}/code/depth.sh

## call genotypes
sbatch -t 48:00:00 --export=basedir="/home/jmiller1/Sex_Map" ${basedir}/code/callgt.sh

REMOVE EMPTY VCFs

## get all chr5 and unmapped contigs
grep '^chain' ${basedir}/metadata/meiotic_map.chainfile | cut -f3,8 | less -S | sort | uniq > ${basedir}/metadata/Scaf_chrm.txt
grep 'chr5' ${basedir}/metadata/Scaf_chrm.txt | cut -f1 > ${basedir}/metadata/chrm5_only.txt

grep 'chr5' ${basedir}/metadata/Scaf_chrm.txt | cut -f1 > ${basedir}/metadata/Scaf_chrm5_only.txt
grep -v 'chr' ${basedir}/metadata/Scaf_chrm.txt | cut -f1 >> ${basedir}/metadata/Scaf_chrm5_only.txt
find "${basedir}/genotypes" -type f -size +10k -name '*.vcf' |
grep -f ${basedir}/metadata/Scaf_chrm5_only.txt > ${basedir}/metadata/Scafolds_to_map.txt

## bgzip non-empty contig vcfs and index
while read F
do
	dest=$(basename $F)
	bgzip -@ 6 -f ${F} -c > ${basedir}/chr5_bcfs/"${dest%\.vcf}".bcf
	bcftools index -f ${basedir}/chr5_bcfs/"${dest%\.vcf}".bcf --threads 6
done < "${basedir}/metadata/Scafolds_to_map.txt"

### Calculate some depth statistics to determine upper coverage cutoff
find ${basedir}/chr5_bcfs/ -type f -size +10k -name '*1.bcf' > ${basedir}/metadata/filelist.txt
bcftools concat -D -a --file-list ${basedir}/metadata/filelist.txt -O b -o ${basedir}/chr5_bcfs/NBH_CHR5.bcf
bcftools sort ${basedir}/chr5_bcfs/NBH_CHR5.bcf -O b -o ${basedir}/chr5_bcfs/NBH_CHR5_sorted.bcf
bcftools index ${basedir}/chr5_bcfs/NBH_CHR5_sorted.bcf
bcftools stats --depth 10,5000,10 ${basedir}/chr5_bcfs/NBH_CHR5_sorted.bcf > NBH_CHR5_unfilt.depthstats

### Filter bcfs and write a vcf for plink, and convert to csv for rQTL
bcftools view -m2 -M2 -v snps -q 0.1 -Q 0.9 -e 'INFO/AC<10 & FMT/DP<10 & FMT/DP>1000 & INFO/MQM<30 & INFO/NS<10' \
	${basedir}/chr5_bcfs/NBH_CHR5_sorted.bcf -O v -o ${basedir}/chr5_bcfs/NBH_CHR5_filt.vcf
bcftools stats --depth 1,1000,5 ${basedir}/chr5_bcfs/NBH_CHR5_filt.vcf > NBH_CHR5_filt.depthstats

### Plink conversion PED and MAP files for rQTL
plink=~/bin/plink
##init_flagset='--make-bed  --allow-extra-chr --autosome-num 24 --allow-no-sex --family'
pop='NBH'

### Update barcode id to samp id
$plink \
	--vcf ${basedir}/chr5_bcfs/NBH_CHR5_filt.vcf \
	--out ${basedir}/plink_files/NBH_CHR5_filt \
	--update-ids ${basedir}/plink_files/SOMM.txt \
	--allow-no-sex \
	--allow-extra-chr \
	--set-missing-var-ids @:# \
	--make-bed

### relate samp id to phenotypes
$plink \
	--bfile ${basedir}/plink_files/NBH_CHR5_filt \
	--out ${basedir}/plink_files/NBH_CHR5_pheno \
	--pheno ${basedir}/plink_files/SOMM.FAM.2.txt \
	--allow-no-sex \
	--allow-extra-chr \
	--all-pheno \
	--family \
	--keep-cluster-names ${pop} \
	--make-bed

### recode to new map and ped for conversion
$plink \
	--bfile ${basedir}/plink_files/NBH_CHR5_pheno \
	--out ${basedir}/plink_files/NBH_CHR5_filt_conv \
	--pheno ${basedir}/plink_files/SOMM.FAM.2.txt \
	--allow-no-sex \
	--allow-extra-chr \
	--all-pheno \
	--family \
	--keep-cluster-names ${pop} \
	--biallelic-only strict \
	--snps-only just-acgt \
	--nonfounders \
	--recode


Rscript ${basedir}/code/Plink_conv_rQTL.R
Rscript ${basedir}/code/Chr_5_Sex_Sp_Map.R



### PARENTS #################################################
## align to ncbi genome
sbatch -p low -t 48:00:00 --export=basedir="/home/jmiller1/Sex_Map" ${basedir}/code/align_parents.sh
## calculate depth to find rad sites
sbatch -p high -t 48:00:00 --mem=64000 --export=basedir="/home/jmiller1/Sex_Map" ${basedir}/code/depth.sh
## call genotypes
sbatch -t 48:00:00 --export=basedir="/home/jmiller1/Sex_Map" ${basedir}/code/callgt.sh
