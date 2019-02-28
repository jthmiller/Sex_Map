#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o /home/jmiller1/Sex_Map/alignments/er_out/array_job_out_%A_%a.txt
#SBATCH -e /home/jmiller1/Sex_Map/alignments/er_out/array_job_err_%A_%a.txt
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=2gb
#SBATCH -p lo
#SBATCH --array=1-96

my_samtools=/home/jmiller1/bin/samtools-1.3/samtools
my_bwa=/home/jmiller1/bin/bwa-0.7.12/bwa
my_samblstr=/home/jmiller1/bin/samblaster-master/samblaster

### NBH fastq
indir=/home/jmiller1/QTL_Map_Raw/demulti_out/NBH/

## ncbi genome
refdir=/home/jmiller1/bin/code/NCBI_Genome
bwagenind=${refdir}/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fna

## bamdir
outdir=${basedir}/alignments

lib=$(echo $fq1 | sed 's/_/\t/'|cut -f 1)
fq1=$(find $indir -name "*RA*fastq.gz" | sed -n $(echo $SLURM_ARRAY_TASK_ID)p)
fq2=$(echo $fq1 | sed 's/RA/RB/')
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
