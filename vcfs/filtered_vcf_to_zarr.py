# Import packages.
import allel
import numcodecs
import numpy as np
import sys
import zarr

# sys.argv[1] = refrence
# sys.argv[2] = chromosome
# sys.argv[3] = called genotype threshold

# Define file paths.
vcf_path = f'./{sys.argv[1]}/san_yak_tei_mel_{sys.argv[3]}_filtered_chr{sys.argv[2]}.vcf.gz'
zarr_path = f'../zarrs/{sys.argv[3]}/{sys.argv[1]}/chr{sys.argv[2]}.zarr'

# Convert the vcf file to a zarr array.
allel.vcf_to_zarr(
    vcf_path, zarr_path, group=str(sys.argv[2]),
    fields=['samples', 'DP', 'POS', 'CHROM', 'GT'], log=sys.stdout, overwrite=True,
)
