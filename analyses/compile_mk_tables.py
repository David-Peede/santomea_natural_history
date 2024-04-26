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
def load_callset_pos(prefix, miss, chrom, ann):
    # Intialize the file path.
    path = f'../zarrs/{miss}/{prefix}/{ann}/chr{chrom}.zarr'
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

# Load the ortholog meta data.
ortholog_coords_df = pd.read_csv(f'./tables/{ref}_{miss}_ortholog_coords.txt', sep='\t')

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

# Intialize a mk dictionary.
mk_dicc = {}
# For every population.
for pop in idx_dicc:
    # Intialize the subdictionary.
    mk_dicc[pop] = {
        'ds': [], 'dn': [],
        'ps': [], 'pn': [],
    }
    
# For every chromosome.
for chrom in ortholog_coords_df['chr'].unique():
    # Extract the genotype callsets and positions.
    s_callset, s_pos = load_callset_pos(ref, miss, chrom, 'synonymous_variant')
    n_callset, n_pos = load_callset_pos(ref, miss, chrom, 'missense_variant')
    # Extract the starts, ends, and all positions.
    starts = ortho_coord_dicc[chrom]['start']
    ends = ortho_coord_dicc[chrom]['end']
    # For all start and end positions.
    for idx in range(starts.size):
        # Extract the start and end positions.
        start = starts[idx]
        end = ends[idx]
        # Determine the positions for this gene.
        s_idx = np.where(((start <= s_pos) & (s_pos <= end)))[0]
        n_idx = np.where(((start <= n_pos) & (n_pos <= end)))[0]
        # If there are syn-muts in this gene.
        if s_idx.size > 0:
            # Identify the window to extract.
            s_loc = s_pos.locate_range(start, end)
            # For every population.
            for pop in mk_dicc:
                # Count the number of syn-muts.
                mk_dicc[pop]['ds'].append(allel.GenotypeArray(s_callset[s_loc]).take(idx_dicc[pop], axis=1).count_alleles().count_non_segregating())
                mk_dicc[pop]['ps'].append(allel.GenotypeArray(s_callset[s_loc]).take(idx_dicc[pop], axis=1).count_alleles().count_segregating())
        # Else, there are no syn-muts that passed qc.
        else:
            # For every population.
            for pop in mk_dicc:
                # Update the lists.
                mk_dicc[pop]['ds'].append(0)
                mk_dicc[pop]['ps'].append(0)
        # If there are non-syn-muts in this gene.
        if n_idx.size > 0:
            # Identify the window to extract.
            n_loc = n_pos.locate_range(start, end)
            # For every population.
            for pop in mk_dicc:
                # Count the number of non-syn-muts.
                mk_dicc[pop]['dn'].append(allel.GenotypeArray(n_callset[n_loc]).take(idx_dicc[pop], axis=1).count_alleles().count_non_segregating())
                mk_dicc[pop]['pn'].append(allel.GenotypeArray(n_callset[n_loc]).take(idx_dicc[pop], axis=1).count_alleles().count_segregating())
        # Else, there are no non-syn-sites that passed qc.
        else:
            # For every population.
            for pop in mk_dicc:
                # Update the lists.
                mk_dicc[pop]['dn'].append(0)
                mk_dicc[pop]['pn'].append(0)

# For every population.
for pop in mk_dicc:
    for mut in mk_dicc[pop]:          
        # Update the orthlog data frame.
        ortholog_coords_df[f'{pop}_{mut}'] = mk_dicc[pop][mut]

# Export the annotated ortholog informtaion.
ortholog_coords_df.to_csv(f'./tables/{ref}_{miss}_ortholog_mk.txt', sep='\t', index=False)