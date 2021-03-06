#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o /home/jmiller1/Sex_Map/alignments/er_out/array_job_out_%A_%a.txt
#SBATCH -e /home/jmiller1/Sex_Map/alignments/er_out/array_job_err_%A_%a.txt
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=2gb
#SBATCH -p lo
#SBATCH --array=1-2

my_samtools=/home/jmiller1/bin/samtools-1.3/samtools
my_bwa=/home/jmiller1/bin/bwa-0.7.12/bwa
my_samblstr=/home/jmiller1/bin/samblaster-master/samblaster

## ncbi genome
refdir=/home/jmiller1/bin/code/NCBI_Genome
bwagenind=${refdir}/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fna

## bamdir
outdir=${basedir}/alignments

### Parents in BRP lib  fastq
indir=/home/jmiller1/QTL_Map_Raw/demulti_out/BP/

###SOMM087 ACAAGCTA        NBH     NBH1M
###SOMM087 AACCGAGA        NBH     NBH1F
###SOMM088 AATCCGTC        BLI     BI1124M
###SOMM088 CACCTTAC        NEW     NEW911M


p1=$(find $indir -name "*RA*ACAAGCTA*")
p2=$(find $indir -name "*RA*AACCGAGA*")
q1=$(printf "%s\n" "$p1" "$p2")

fq1=$(echo "$q1" | sed -n $(echo $SLURM_ARRAY_TASK_ID)p)
fq2=$(echo $fq1 | sed 's/RA/RB/')

lib='SOMM087'
root=$(echo $fq1 | sed 's/.*\///' | sed 's/_R._/_/' | cut -c 1-8,11-18)
rg=$(echo \@RG\\tID:$root\\tPL:Illumina\\tPU:x\\tLB:$lib\\tSM:$root)
tempsort=$root.temp
outfile=$outdir/$root.bam

echo $root
echo $fq1
echo $fq2
echo $rg
echo $tempsort
echo $outfile

$my_bwa mem $bwagenind -t 6 -R $rg $fq1 $fq2 | $my_samblstr | $my_samtools view -S -h -u - | $my_samtools sort - -T /scratch/$tempsort -O bam -o $outfile

$my_samtools index $outfile
