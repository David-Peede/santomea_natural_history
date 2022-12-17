# Import packages.
import allel
import numcodecs
import numpy as np
import sys
import zarr

# sys.argv[1] = refrence species
# sys.argv[1] = chromosome

# Define file paths.
vcf_path = './{0}/san_yak_tei_filtered_chr{1}.vcf.gz'.format(str(sys.argv[1]), str(sys.argv[2]))
zarr_path = '../zarrs/{0}/chr{1}.zarr'.format(str(sys.argv[1]), str(sys.argv[2]))

# Convert the vcf file to a zarr array.
allel.vcf_to_zarr(
    vcf_path, zarr_path, group=str(sys.argv[2]),
    fields=['samples', 'DP', 'POS', 'CHROM', 'GT'], log=sys.stdout, overwrite=True,
)