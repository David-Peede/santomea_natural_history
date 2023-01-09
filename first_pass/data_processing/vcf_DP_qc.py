import allel
import numpy as np
import sys

# sys.argv[1] = sample prefix
# sys.argv[2] = reference genome
# sys.argv[3] = chromosome

unfiltered_vcf = str(sys.argv[1])+'_unfiltered_m'+str(sys.argv[2])+'_chr'+str(sys.argv[3])+'.vcf.gz'
callset = allel.read_vcf(unfiltered_vcf, fields=['DP'])

DP = callset['variants/DP']

lower_bound = np.percentile(DP, 2.5)
upper_bound = np.percentile(DP, 97.5)

filtered_vcf = str(sys.argv[1])+'_filtered_DP_m'+str(sys.argv[2])+'_chr'+str(sys.argv[3])+'.vcf.gz'
print("bcftools view -i 'DP>="+str(lower_bound.astype(int))+" && DP<="+str(upper_bound.astype(int))+" && QUAL>=20 && MQ>=30' -Oz -o "+filtered_vcf+" "+unfiltered_vcf)