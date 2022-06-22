#!/bin/bash
for CHR in 2L 2R 3L 3R 4 X; do for SPECIES in san tei yak; do
sbatch -J ${SPECIES}_${CHR}_filter -N 1 -n 6 -t 1-0 --mem=2G -o ${SPECIES}_${CHR}_filter-%A.out -e ${SPECIES}_${CHR}_filter-%A.err --mail-type=ALL --mail-user=dpeede@live.unc.edu --wrap="module add samtools; bcftools view -M2 -e 'INFO/INDEL="INDEL" && F_MISSING>0.1' -Oz -o all_${SPECIES}_${CHR}_filtered.vcf.gz all_${SPECIES}_${CHR}_unfiltered.vcf.gz"
done; done