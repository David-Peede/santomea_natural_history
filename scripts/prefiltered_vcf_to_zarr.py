# Import packages.
import allel
import numcodecs
import numpy as np
import sys
import zarr

# sys.argv[1] = refrence species

# Define file paths.
vcf_path = './{0}/san_yak_tei_prefiltered.vcf.gz'.format(str(sys.argv[1]))
zarr_path = '../zarrs/{0}/san_yak_tei_prefiltered.zarr'.format(str(sys.argv[1]))

# Convert the vcf file to a zarr array.
allel.vcf_to_zarr(vcf_path, zarr_path, fields=['samples', 'DP', 'POS', 'CHROM', 'GT'], log=sys.stdout, overwrite=True)