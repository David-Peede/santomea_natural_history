#!/bin/bash
#
#SBATCH -J fastq_2_bam # Job name.
#SBATCH -n 12 # Number of cores.
#SBATCH -N 1 # Ensure that all cores are on one machine.
#SBATCH -t 7-0 # Runtime in D-HH:MM (or use minutes).
#SBATCH --mem 96g # Memory in GB.
#SBATCH -o fastq_2_bam-%A.out # File for STDOUT (with jobid = %j).
#SBATCH -e fastq_2_bam-%A.err # File for STDERR (with jobid = %j).
#SBATCH --mail-type=ALL # Type of email notification: BEGIN,END,FAIL,ALL.
#SBATCH --mail-user= # Email where notifications will be sent.


module add trim_galore
module add bwa/0.7.15
module add samtools/1.4
module add picard/2.2.4
module add gatk/3.7


FQDIR= # Path to your fastq directory
MYDIR= # Path to your desired output directory.
REF= # Path to the reference fasta.


for SAMPLE in # List of sample prefixes.
do
        echo ${SAMPLE}
trim_galore -q 20 --paired --illumina \
$FQDIR/${SAMPLE}_R1_001.fastq.gz \
$FQDIR/${SAMPLE}_R2_001.fastq.gz \
-o $MYDIR/

bwa mem -M -t 6 -k 20 -R '@RG\tID:'${SAMPLE}':'${SAMPLE}'\tSM:'${SAMPLE}'\tLB:ga-'${SAMPLE}'\tPL:Illumina' \
$REF \
$MYDIR/${SAMPLE}_R1_001_val_1.fq.gz \
$MYDIR/${SAMPLE}_R2_001_val_2.fq.gz | samtools sort -@ 6 \
-o $MYDIR/${SAMPLE}_R1R2.bam

samtools view -bF 4 -@ 6 -q 1 \
$MYDIR/${SAMPLE}_R1R2.bam > $MYDIR/${SAMPLE}_sorted.bam 

java -jar /nas/longleaf/apps/picard/2.2.4/picard-tools-2.2.4/picard.jar MarkDuplicates \
INPUT=$MYDIR/${SAMPLE}_sorted.bam \
OUTPUT=$MYDIR/${SAMPLE}_dedup.bam \
METRICS_FILE=$MYDIR/${SAMPLE}_dedup_metrics.txt 

java -jar /nas/longleaf/apps/picard/2.2.4/picard-tools-2.2.4/picard.jar BuildBamIndex \
INPUT=$MYDIR/${SAMPLE}_dedup.bam

java -Xmx16g -jar /nas/longleaf/apps/gatk/3.7/GenomeAnalysisTK.jar \
-nt 6 \
-T RealignerTargetCreator \
-R $REF \
-I $MYDIR/${SAMPLE}_dedup.bam \
-o $MYDIR/${SAMPLE}_realignertargetcreator.intervals \

java -Xmx64g -Djava.io.tmpdir=/proj/matutelb/projects/drosophila/zaprionus/temp_trash \
-jar /nas/longleaf/apps/gatk/3.7/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $REF \
-I $MYDIR/${SAMPLE}_dedup.bam \
-targetIntervals $MYDIR/${SAMPLE}_realignertargetcreator.intervals \
-o $MYDIR/${SAMPLE}_realigned.bam \
-LOD 2.0
done;