# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

# sys.argv[1] = reference
# sys.argv[2] = missing data threshold


# Define a function to load genotyope and positions arrays.
def load_callset_pos(prefix, miss, chrom):
    # Intialize the file path.
    path = f'../zarrs/{miss}/{prefix}/chr{chrom}.zarr'
    # Load the vcf file.
    callset = zarr.open_group(path, mode='r')
    # Extract the genotypes.
    geno = callset[f'{chrom}/calldata/GT']
    # Load the positions.
    pos = allel.SortedIndex(callset[f'{chrom}/variants/POS'])
    return geno, pos

# Define a function to calculate nucleotide diversity.
def pixy(gt, pop_idx):
    # Extract the alternative allele count matrix.
    aac = gt.take(pop_idx, axis=1).count_alleles()
    # If the alternative allele count matrix is mono-allelic.
    if aac.shape[1] == 1:
        # Then there are no pairwise differences.
        return 0
    # Else, the site is bi-allelic.
    else:
        # Determine the number of called chromosomes for each site.
        called_chromosomes = np.nansum(aac, axis=1)
        # Create a mask where there are no called genotypes.
        mask = called_chromosomes == 0
        # Determine the allele counts of the derived/alternative allele.
        derived_allele_count = aac[:, 1]
        # Determine the allele counts of the ancestral/reference allele.
        ancestral_allele_count = called_chromosomes - derived_allele_count
        # Determine the number of comparisons per site.
        nC2 = np.array([((n * (n - 1)) / 2) for n in called_chromosomes])
        # Calculate the numerator.
        numerator = np.nansum((derived_allele_count[~mask] * ancestral_allele_count[~mask]))
        # Calculate the denominator.
        denominator = np.nansum(nC2[~mask])
        # Calculate pixy.
        pi_pixy = numerator / denominator
        return pi_pixy

# Read in meta data has a pandas dataframe.
meta_df = pd.read_csv('./dros_meta_data.txt', sep='\t')
# Intialize a dictionary.
idx_dicc = {}
# For all populations.
for pop in ['san', 'yak_symp', 'yak_allo', 'tei', 'mel']:
    # Fill the dictionary.
    idx_dicc[pop] = meta_df[meta_df['population'] == pop].index.values
# Fill the dictionary for the last population.
idx_dicc['yak'] = np.concatenate((idx_dicc['yak_allo'], idx_dicc['yak_symp']))

# Extract the reference assembly.
ref = str(sys.argv[1])
miss = str(sys.argv[2])

# Read in the window dataframe.
window_df = pd.read_csv(f'./tables/{ref}_{miss}_5kb_windows_qced.txt', sep='\t')
# Extract window information.
chroms = window_df.chr.values
starts = window_df.start.values
ends = window_df.end.values
# Intialize a results matrix to store the results.
results_mat = np.empty((window_df.shape[0], 6))

# For every non-overlapping window.
for idx in range(window_df.shape[0]):
    # Extract the window information.
    chrom = chroms[idx]
    start = starts[idx]
    end = ends[idx]
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos(ref, miss, chrom)
    # Identify the window to extract.
    wind_loc = all_pos.locate_range(start, end)
    # Compute pi.
    san_pi = pixy(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_idx=idx_dicc['san'],
    )
    symp_pi = pixy(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_idx=idx_dicc['yak_symp'],
    )
    allo_pi = pixy(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_idx=idx_dicc['yak_allo'],
    )
    yak_pi = pixy(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_idx=idx_dicc['yak'],
    )
    tei_pi = pixy(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_idx=idx_dicc['tei'],
    )
    mel_pi = pixy(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_idx=idx_dicc['mel'],
    )
    # Append the results.
    results_mat[idx, :] = np.array([san_pi, symp_pi, allo_pi, yak_pi, tei_pi, mel_pi])
        
# Export the the results matrix.
np.savetxt(
    f'./results/{ref}/{miss}/pi_5kb_windows.txt.gz',
    results_mat, fmt='%1.15f', delimiter='\t', newline='\n',
)