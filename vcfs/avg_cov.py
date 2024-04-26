# Import packages
import allel
import numpy as np
import sys

# sys.argv[1] = sample prefix
# sys.argv[2] = species

# Load the data.
filtered_vcf = './'+str(sys.argv[2])+'/'+str(sys.argv[1])+'_dp_mq_gq_filtered.vcf.gz'
callset = allel.read_vcf(filtered_vcf, fields=['DP'])

# Extract the read depth information.
dp = callset['variants/DP']

# Calculate the average.
avg_cov = np.mean(dp)

# Write the information to an output file.
outfile = sys.stdout
outfile.write('\t'.join([str(sys.argv[1]), str(avg_cov)])+'\n')