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

# Extract the ref assembly.
ref = str(sys.argv[1])
miss = str(sys.argv[2])

# Load the ortholog meta data.
ortholog_coords_df = pd.read_csv(f'./tables/{ref}_ortholog_coords.txt', sep='\t')

# Intialize a dictionary to store coordinates.
ortho_coord_dicc = {}
# For all chromosomes.
for chrom in ortholog_coords_df['chr'].unique():
    # Intialize the dictionary.
    ortho_coord_dicc[chrom] = {}
    # Subset the dataframe.
    chrom_df = ortholog_coords_df[ortholog_coords_df['chr'] == chrom]
    # Fill the dictionary.
    ortho_coord_dicc[chrom]['start'] = chrom_df['start'].values
    ortho_coord_dicc[chrom]['end'] = chrom_df['end'].values

# Intialize lists.
n_tot = []
n_seg = []

# For every chromosome.
for chrom in ortholog_coords_df['chr'].unique():
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos(ref, miss, chrom)
    # Extract the starts, ends, and all positions.
    starts = ortho_coord_dicc[chrom]['start']
    ends = ortho_coord_dicc[chrom]['end']
    # For all start and end positions.
    for idx in range(starts.size):
        # Extract the start and end positions.
        start = starts[idx]
        end = ends[idx]
        # Determine the positions for this gene.
        gene_idx = np.where(((start <= all_pos) & (all_pos <= end)))[0]
        # If there are sites in this gene.
        if gene_idx.size > 0:
            # Identify the window to extract.
            wind_loc = all_pos.locate_range(start, end)
            # Update the lists.
            n_tot.append(gene_idx.size)
            n_seg.append(allel.GenotypeArray(callset[wind_loc]).count_alleles().count_segregating())
        # Else, there are no sites that passed qc.
        else:
            # Update the lists.
            n_tot.append(0)
            n_seg.append(0)
            
# Update the orthlog data frame.
ortholog_coords_df['tot_sites'] = n_tot
ortholog_coords_df['seg_sites'] = n_seg

# Export the annotated ortholog informtaion.
ortholog_coords_df.to_csv(f'./tables/{ref}_{miss}_ortholog_coords.txt', sep='\t', index=False)