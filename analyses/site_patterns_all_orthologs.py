# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

# sys.argv[1] = reference
# sys.argv[2] = missing data threshold
# sys.argv[3] = p1 population
# sys.argv[4] = p2 population
# sys.argv[5] = p3 population
# sys.argv[6] = p4 population

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

# Define a site pattern function.
def site_patterns(p1, p2, p3, p4):
    # Calculate site pattern counts.
    abba = np.nansum((1 - p1) * (p2) * (p3) * (1 - p4))
    baba = np.nansum((p1) * (1 - p2) * (p3) * (1 - p4))
    bbaa = np.nansum((p1) * (p2) * (1 - p3) * (1 - p4))
    baaa = np.nansum((p1) * (1 - p2) * (1 - p3) * (1 - p4))
    abaa = np.nansum((1 - p1) * (p2) * (1 - p3) * (1 - p4))
    aaba = np.nansum((1 - p1) * (1 - p2) * (p3) * (1 - p4))
    return abba, baba, bbaa, baaa, abaa, aaba

# Define a function to calculate site patterns.
def dros_site_patterns(
    gt,
    p1_idx, p2_idx, p3_idx, p4_idx,
):
    # Determine the indicies where each population has called genotypes.
    p1_mask = (gt.take(p1_idx, axis=1).is_called() == True).any(axis=1)
    p2_mask = (gt.take(p2_idx, axis=1).is_called() == True).any(axis=1)
    p3_mask = (gt.take(p3_idx, axis=1).is_called() == True).any(axis=1)
    p4_mask = (gt.take(p4_idx, axis=1).is_called() == True).any(axis=1)
    # Determine the indicied where all populations have called genotypes.
    called_mask = (p1_mask & p2_mask & p3_mask & p4_mask)
    # If there are no sites called between all three samples...
    if (called_mask.sum() == 0):
        # Set the results to np.nan since we don't have any sites to perform computations on.
        results = np.full(6, np.nan)
    # Else...
    else:
        # Determine the indicies where we have varibale sites.
        var_mask = gt.compress(called_mask, axis=0).count_alleles().is_variant()
        # If there are no variable sites...
        if (var_mask.sum() == 0):
            # Set the results to 0 since we are iterating over QC'ed regions.
            results = np.zeros(6)
        # Else...
        else:
            # Calculate the alternative allele frequencies.
            p1_alt_freqs = calc_alt_freqs(gt.take(p1_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p2_alt_freqs = calc_alt_freqs(gt.take(p2_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p3_alt_freqs = calc_alt_freqs(gt.take(p3_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p4_alt_freqs = calc_alt_freqs(gt.take(p4_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            # Polarize the allele frequencies based on the most common allele in the outgroup.
            p1_der_freqs = np.where(p4_alt_freqs > 0.5, np.abs(p1_alt_freqs - 1), p1_alt_freqs)
            p2_der_freqs = np.where(p4_alt_freqs > 0.5, np.abs(p2_alt_freqs - 1), p2_alt_freqs)
            p3_der_freqs = np.where(p4_alt_freqs > 0.5, np.abs(p3_alt_freqs - 1), p3_alt_freqs)
            p4_der_freqs = np.where(p4_alt_freqs > 0.5, np.abs(p4_alt_freqs - 1), p4_alt_freqs)
            # Calculate the site pattern counts.
            abba, baba, bbaa, baaa, abaa, aaba = site_patterns(
                p1_der_freqs, p2_der_freqs, p3_der_freqs, p4_der_freqs,
            )
            # Intialize a results array.
            results = np.array([abba, baba, bbaa, baaa, abaa, aaba])
    return results

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

# Extract the reference assembly and focal populations.
ref = str(sys.argv[1])
miss = str(sys.argv[2])
p1 = str(sys.argv[3])
p2 = str(sys.argv[4])
p3 = str(sys.argv[5])
p4 = str(sys.argv[6])

# Read in the ortholog dataframe.
ortholog_df = pd.read_csv(f'./tables/{ref}_{miss}_ortholog_coords.txt', sep='\t')
# Extract ortholog information.
chroms = ortholog_df.chr.values
starts = ortholog_df.start.values
ends = ortholog_df.end.values
# Intialize a results matrix to store the results.
results_mat = np.empty((ortholog_df.shape[0], 6))

# For every ortholog window.
for idx in range(ortholog_df.shape[0]):
    # Extract the ortholog information.
    chrom = chroms[idx]
    start = starts[idx]
    end = ends[idx]
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos(ref, miss, chrom)
    # Determine the positions for this gene.
    gene_idx = np.where(((start <= all_pos) & (all_pos <= end)))[0]
    # If there are no sites in this gene.
    if gene_idx.size == 0:
        # Append the results.
        results_mat[idx, :] = np.full(6, np.nan)
    # Else there are sites to do computations on.
    else:
        # Identify the window to extract.
        wind_loc = all_pos.locate_range(start, end)
        # Calculate site patterns.
        sps = dros_site_patterns(
            gt=allel.GenotypeArray(callset[wind_loc]),
            p1_idx=idx_dicc[p1], p2_idx=idx_dicc[p2],
            p3_idx=idx_dicc[p3], p4_idx=idx_dicc[p4],
        )
        # Append the results.
        results_mat[idx, :] = sps
        
# Export the the results matrix.
np.savetxt(
    f'./results/{ref}/{miss}/p1_{p1}_p2_{p2}_p3_{p3}_p4_{p4}_sps_orthologs.txt.gz',
    results_mat, fmt='%1.15f', delimiter='\t', newline='\n',
)