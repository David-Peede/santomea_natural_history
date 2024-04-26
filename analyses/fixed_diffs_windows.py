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

# Define a function to calculate alternative allele frequencies.
def calc_alt_freqs(gt):
    # If there are no altenative alleles...
    if (gt.count_alleles().shape[1] == 1):
        # Calculate alternative allele frequencies.
        alt_freqs = gt.count_alleles().to_frequencies()[:, 0] - 1
    # Else...
    else:
        # Calculate alternative allele frequencies.
        alt_freqs = gt.count_alleles().to_frequencies()[:, 1]
    return alt_freqs

# Define a function to count the number of fixed differences between two populations.
def count_fixed_diffs(gt, pop_a, pop_b):
    # Compute the allele frequencies.
    pop_a_freqs = calc_alt_freqs(gt=gt.take(pop_a, axis=1))
    pop_b_freqs = calc_alt_freqs(gt=gt.take(pop_b, axis=1))
    # Count the number of fixed differences.
    n_fixed_diffs = (np.abs(pop_a_freqs - pop_b_freqs) == 1).sum()
    return n_fixed_diffs

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
results_mat = np.empty((window_df.shape[0], 13))

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
    # Count the number of fixed differences.
    san_mel = count_fixed_diffs(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_a=idx_dicc['san'], pop_b=idx_dicc['mel'],
    )
    san_tei = count_fixed_diffs(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_a=idx_dicc['san'], pop_b=idx_dicc['tei'],
    )
    san_yak = count_fixed_diffs(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_a=idx_dicc['san'], pop_b=idx_dicc['yak'],
    )
    san_symp = count_fixed_diffs(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_a=idx_dicc['san'], pop_b=idx_dicc['yak_symp'],
    )
    san_allo = count_fixed_diffs(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_a=idx_dicc['san'], pop_b=idx_dicc['yak_allo'],
    )
    tei_mel = count_fixed_diffs(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_a=idx_dicc['tei'], pop_b=idx_dicc['mel'],
    )
    tei_yak = count_fixed_diffs(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak'],
    )
    tei_symp = count_fixed_diffs(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak_symp'],
    )
    tei_allo = count_fixed_diffs(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak_allo'],
    )
    yak_mel = count_fixed_diffs(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_a=idx_dicc['yak'], pop_b=idx_dicc['mel'],
    )
    symp_allo = count_fixed_diffs(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_a=idx_dicc['yak_symp'], pop_b=idx_dicc['yak_allo'],
    )
    symp_mel = count_fixed_diffs(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_a=idx_dicc['yak_symp'], pop_b=idx_dicc['mel'],
    )
    allo_mel = count_fixed_diffs(
        gt=allel.GenotypeArray(callset[wind_loc]),
        pop_a=idx_dicc['yak_allo'], pop_b=idx_dicc['mel'],
    )
    # Append the results.
    results_mat[idx, :] = np.array([
        san_mel, san_tei, san_yak, san_symp, san_allo,
        tei_mel, tei_yak, tei_symp, tei_allo,
        yak_mel, symp_allo, symp_mel, allo_mel,
    ])
        
# Export the the results matrix.
np.savetxt(
    f'./results/{ref}/{miss}/fixed_diffs_5kb_windows.txt.gz',
    results_mat, fmt='%1.15f', delimiter='\t', newline='\n',
)