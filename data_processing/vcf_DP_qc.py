import allel
import numpy as np
import sys

VCF = str(sys.argv[1])
callset = allel.read_vcf(VCF, fields=['DP'])

DP = callset['variants/DP']

lower_bound = np.percentile(DP, 2.5)
upper_bound = np.percentile(DP, 97.5)

print("bcftools view -i 'DP>="+lower_bound+" && DP<="+upper_bound+" && MQ>=30' -Oz -o "+VCF)