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

# Read in the ortholog dataframe.
ortholog_df = pd.read_csv(f'./tables/{ref}_{miss}_ortholog_coords.txt', sep='\t')
# Extract ortholog information.
chroms = ortholog_df.chr.values
starts = ortholog_df.start.values
ends = ortholog_df.end.values
# Intialize a results matricies to store the results.
all_mat = np.empty((ortholog_df.shape[0], 6))

# For every ortholog window.
for idx in range(ortholog_df.shape[0]):
    # Extract the ortholog information.
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
        all_mat[idx, :] = np.full(6, np.nan)
    # Else there are sites to do computations on.
    else:
        # Identify the gene window to extract.
        gene_loc = all_pos.locate_range(start, end)
        # Compute pi.
        all_san_pi = pixy(
            gt=allel.GenotypeArray(callset[gene_loc]),
            pop_idx=idx_dicc['san'],
        )
        all_symp_pi = pixy(
            gt=allel.GenotypeArray(callset[gene_loc]),
            pop_idx=idx_dicc['yak_symp'],
        )
        all_allo_pi = pixy(
            gt=allel.GenotypeArray(callset[gene_loc]),
            pop_idx=idx_dicc['yak_allo'],
        )
        all_yak_pi = pixy(
            gt=allel.GenotypeArray(callset[gene_loc]),
            pop_idx=idx_dicc['yak'],
        )
        all_tei_pi = pixy(
            gt=allel.GenotypeArray(callset[gene_loc]),
            pop_idx=idx_dicc['tei'],
        )
        all_mel_pi = pixy(
            gt=allel.GenotypeArray(callset[gene_loc]),
            pop_idx=idx_dicc['mel'],
        )
        # Append the results.
        all_mat[idx, :] = np.array([
            all_san_pi, all_symp_pi, all_allo_pi,
            all_yak_pi, all_tei_pi, all_mel_pi,
        ])

# Export the the results matricies.
np.savetxt(
    f'./results/{ref}/{miss}/pi_all_orthologs.txt.gz',
    all_mat, fmt='%1.15f', delimiter='\t', newline='\n',
)