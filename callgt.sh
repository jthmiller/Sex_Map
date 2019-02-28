#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o /home/jmiller1/QTL_Map_Raw/QTL_scripts/out_er/array_job_out_%A_%a.txt
#SBATCH -e /home/jmiller1/QTL_Map_Raw/QTL_scripts/out_er/array_job_err_%A_%a.txt
#SBATCH --array=1-8577
#SBATCH -p med
#SBATCH --mem=30000
#SBATCH -t 48:00:00
###### number of nodes
###### number of processors
#SBATCH --cpus-per-task=6

refdir=/home/jmiller1/bin/code/NCBI_Genome
bwagenind=${refdir}/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fna

my_freebayes=/home/jmiller1/bin/freebayes/bin/freebayes
my_bedtools=/home/jmiller1/bin/bedtools2/bin/bedtools
my_bamtools=/home/jmiller1/bin/bamtools-master/bin/bamtools

##scaf to align
scaf=$(sed -n "$SLURM_ARRAY_TASK_ID p" ${refdir}/unique.scaf.lg.txt)

endpos=$(expr $(grep -P "$scaf\t" ${refdir}/unsplit_merge.fasta.fai | cut -f 2) - 1)

region=$scaf:1..$endpos
echo $region

outfile=$scaf.vcf
echo $outfile

bam_dir=/home/jmiller1/QTL_Map_Raw/align
vcf_out=/home/jmiller1/QTL_Map_Raw/vcf/freebayes.array
bed_regions=~/QTL_Map_Raw/align/radsites.sorted.bed
bam_list=/home/jmiller1/QTL_Map_Raw/align/bamlist.txt
merged_bams=/home/jmiller1/QTL_Map_Raw/align/SOMM0_ALL.bam
all_scaf=/home/jmiller1/bin/code/ALLMAPS_OUT/unsplit_merge.fasta.genomefile.txt

$my_bamtools merge -list $bam_list -region $region| \
	$my_bamtools filter -in stdin -mapQuality '>30' -isProperPair true | \
	$my_bedtools intersect -sorted -a stdin -b $bed_regions -g $all_scaf | \
	$my_freebayes -f $bwagenind --use-best-n-alleles 4 --pooled-discrete --stdin > $vcf_out/$outfile
