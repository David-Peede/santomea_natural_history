#!/bin/bash
for SAMP in [sample_prefix_list]; do
sbatch -J ${SAMP}_mpileup -N 1 -n 12 -t 1-0 --mem=2G -o ${SAMP}_mpileup-%A.out -e ${SAMP}_mpileup-%A.err --mail-type=ALL --mail-user=EMAIL@live.unc.edu --wrap="module add samtools; bcftools mpileup --threads 12 -f REF.fna ${SAMP}_realigned.bam | bcftools call -m -f GQ -Oz -o ${SAMP}_unfiltered.vcf.gz"
done