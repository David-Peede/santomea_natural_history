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

# Read in the ortholog dataframe.
ortholog_df = pd.read_csv(f'./tables/{ref}_{miss}_ortholog_coords.txt', sep='\t')
# Extract ortholog information.
chroms = ortholog_df.chr.values
starts = ortholog_df.start.values
ends = ortholog_df.end.values
# Intialize a results matrix to store the results.
results_mat = np.empty((ortholog_df.shape[0], 6))

# For every non-overlapping window.
for idx in range(ortholog_df.shape[0]):
    # Extract the window information.
    chrom = chroms[idx]
    start = starts[idx]
    end = ends[idx]
    print(chrom, start, end, ortholog_df.gene.values[idx])
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
        # Identify the gene window to extract.
        gene_loc = all_pos.locate_range(start, end)
        # Compute tajimas d.
        san_d = allel.tajima_d(allel.GenotypeArray(callset[gene_loc]).take(idx_dicc['san'], axis=1).count_alleles())
        symp_d = allel.tajima_d(allel.GenotypeArray(callset[gene_loc]).take(idx_dicc['yak_symp'], axis=1).count_alleles())
        allo_d = allel.tajima_d(allel.GenotypeArray(callset[gene_loc]).take(idx_dicc['yak_allo'], axis=1).count_alleles())
        yak_d = allel.tajima_d(allel.GenotypeArray(callset[gene_loc]).take(idx_dicc['yak'], axis=1).count_alleles())
        tei_d = allel.tajima_d(allel.GenotypeArray(callset[gene_loc]).take(idx_dicc['tei'], axis=1).count_alleles())
        mel_d = allel.tajima_d(allel.GenotypeArray(callset[gene_loc]).take(idx_dicc['mel'], axis=1).count_alleles())
        # Append the results.
        results_mat[idx, :] = np.array([san_d, symp_d, allo_d, yak_d, tei_d, mel_d])
        
# Export the the results matrix.
np.savetxt(
    f'./results/{ref}/{miss}/taj_d_orthologs.txt.gz',
    results_mat, fmt='%1.15f', delimiter='\t', newline='\n',
)