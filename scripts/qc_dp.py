# Import packages
import allel
import numpy as np
import sys

# sys.argv[1] = sample prefix
# sys.argv[2] = species

# Load the data.
unfiltered_vcf = './'+str(sys.argv[2])+'/'+str(sys.argv[1])+'_unfiltered.vcf.gz'
callset = allel.read_vcf(unfiltered_vcf, fields=['DP'])

# Extract the read depth information.
dp = callset['variants/DP']

# Calculate lower and upper bounds.
lower_bound = np.percentile(DP, 2.5)
upper_bound = np.percentile(DP, 97.5)

# Write the information to an output file.
outfile = sys.stdout
outfile.write('\t'.join([str(sys.argv[1]), str(lower_bound.astype(int)), str(upper_bound.astype(int))])+'\n')