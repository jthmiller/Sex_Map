#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o /home/jmiller1/QTL_Map_Raw/QTL_scripts/out_er/array_job_out_%A_%a.txt
#SBATCH -e /home/jmiller1/QTL_Map_Raw/QTL_scripts/out_er/array_job_err_%A_%a.txt
#SBATCH --array=1-8577
#SBATCH -p med
#SBATCH --mem=3000
#SBATCH -t 48:00:00
###### number of nodes
###### number of processors
#SBATCH --cpus-per-task=3

refdir=/home/jmiller1/bin/code/NCBI_Genome
bwagenind=${refdir}/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fna
all_scaf=${refdir}/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.genomefile
indir=/home/jmiller1/Sex_Map/alignments

my_freebayes=/home/jmiller1/bin/freebayes/bin/freebayes
my_bedtools=/home/jmiller1/bin/bedtools2/bin/bedtools
my_bamtools=/home/jmiller1/bin/bamtools-master/bin/bamtools

###SOMM087_ACAAGCTA        NBH     NBH1M
###SOMM087_AACCGAGA        NBH     NBH1F

##scaf to align
scaf=$(sed -n "$SLURM_ARRAY_TASK_ID p" ${refdir}/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.scaffolds)
endpos=$(expr $(grep -P "$scaf\t" ${refdir}/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fna.fai | cut -f 2) - 1)
region=$scaf:1..$endpos
echo $region

outfile=$scaf.vcf.gz
echo $outfile

bam_dir=${refdir}/alignments
vcf_out=${basedir}/genotypes/parents
bed_regions=${basedir}/metadata/radsites.bed
##bam_list=${basedir}/metadata/bamlist.txt
bam_list=$(cat <(find $indir -name 'SOMM087_ACAAGCTA.bam') <(find $indir -name 'SOMM087_AACCGAGA.bam'))

$my_bamtools merge -list $bam_list -region $region| \
$my_bamtools filter -in "${indir}"/SOMM087_ACAAGCTA.bam -mapQuality '>30' -isProperPair true | \
	$my_bedtools intersect -sorted -a stdin -b $bed_regions -g $all_scaf | \
	$my_freebayes -f $bwagenind --use-best-n-alleles 4 --pooled-discrete --stdin |
bgzip > $vcf_out/$outfile
