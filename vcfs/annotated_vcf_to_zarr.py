# Import packages.
import allel
import numcodecs
import numpy as np
import sys
import zarr

# sys.argv[1] = chromosome
# sys.argv[2] = called genotype threshold
# sys.argv[3] = annotation type

# Define file paths.
vcf_path = f'./mel_annotations/san_yak_tei_mel_{sys.argv[2]}_filtered_{sys.argv[3]}_chr{sys.argv[1]}.vcf.gz'
zarr_path = f'../zarrs/{sys.argv[2]}/mel/{sys.argv[3]}/chr{sys.argv[1]}.zarr'

# Convert the vcf file to a zarr array.
allel.vcf_to_zarr(
    vcf_path, zarr_path, group=str(sys.argv[1]),
    fields=['samples', 'DP', 'POS', 'CHROM', 'GT'], log=sys.stdout, overwrite=True,
)
