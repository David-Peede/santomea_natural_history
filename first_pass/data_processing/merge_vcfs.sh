#!/bin/bash
for CHR in 2L 2R 3L 3R 4 X; do for SPECIES in san tei yak; do
sbatch -J ${SPECIES}_${CHR}_merge -N 1 -n 6 -t 1-0 --mem=2G -o  ${SPECIES}_${CHR}_merge-%A.out -e  ${SPECIES}_${CHR}_merge-%A.err --mail-type=ALL --mail-user=EMAIL@live.unc.edu --wrap="module add samtools; bcftools merge -l all_${SPECIES}_${CHR}_vcf_list.txt -Oz -o all_${SPECIES}_${CHR}_unfiltered.vcf.gz"
done; done