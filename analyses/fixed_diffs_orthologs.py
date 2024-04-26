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

# Define a function to load genotyope and positions arrays for coding mutations.
def load_coding_mut_callset_pos(prefix, miss, chrom, ann):
    # Intialize the file path.
    path = f'../zarrs/{miss}/{prefix}/{ann}/chr{chrom}.zarr'
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

# Read in the ortholog dataframe.
ortholog_df = pd.read_csv(f'./tables/{ref}_{miss}_ortholog_coords.txt', sep='\t')
# Extract ortholog information.
chroms = ortholog_df.chr.values
starts = ortholog_df.start.values
ends = ortholog_df.end.values
# Intialize a results matricies to store the results.
all_mat = np.empty((ortholog_df.shape[0], 11))
int_mat = np.empty((ortholog_df.shape[0], 11))
syn_mat = np.empty((ortholog_df.shape[0], 11))
mis_mat = np.empty((ortholog_df.shape[0], 11))

# For every ortholog window.
for idx in range(ortholog_df.shape[0]):
    # Extract the ortholog information.
    chrom = chroms[idx]
    start = starts[idx]
    end = ends[idx]
    # Extract the genotype callset and positions.
    callset, all_pos = load_callset_pos(ref, miss, chrom)
    _, syn_sites = load_coding_mut_callset_pos(ref, miss, chrom, 'synonymous_variant')
    _, mis_sites = load_coding_mut_callset_pos(ref, miss, chrom, 'missense_variant')
    # Determine the positions for this gene.
    gene_idx = np.where(((start <= all_pos) & (all_pos <= end)))[0]
    # If there are no sites in this gene.
    if gene_idx.size == 0:
        # Append the results.
        all_mat[idx, :] = np.full(11, np.nan)
        int_mat[idx, :] = np.full(11, np.nan)
        syn_mat[idx, :] = np.full(11, np.nan)
        mis_mat[idx, :] = np.full(11, np.nan)
    # Else there are sites to do computations on.
    else:
        # Identify the gene window to extract.
        gene_loc = all_pos.locate_range(start, end)
        # Count the number of fixed differences.
        all_san_mel = count_fixed_diffs(
            gt=allel.GenotypeArray(callset[gene_loc]),
            pop_a=idx_dicc['san'], pop_b=idx_dicc['mel'],
        )
        all_san_tei = count_fixed_diffs(
            gt=allel.GenotypeArray(callset[gene_loc]),
            pop_a=idx_dicc['san'], pop_b=idx_dicc['tei'],
        )
        all_san_yak = count_fixed_diffs(
            gt=allel.GenotypeArray(callset[gene_loc]),
            pop_a=idx_dicc['san'], pop_b=idx_dicc['yak'],
        )
        all_san_symp = count_fixed_diffs(
            gt=allel.GenotypeArray(callset[gene_loc]),
            pop_a=idx_dicc['san'], pop_b=idx_dicc['yak_symp'],
        )
        all_san_allo = count_fixed_diffs(
            gt=allel.GenotypeArray(callset[gene_loc]),
            pop_a=idx_dicc['san'], pop_b=idx_dicc['yak_allo'],
        )
        all_tei_mel = count_fixed_diffs(
            gt=allel.GenotypeArray(callset[gene_loc]),
            pop_a=idx_dicc['tei'], pop_b=idx_dicc['mel'],
        )
        all_tei_yak = count_fixed_diffs(
            gt=allel.GenotypeArray(callset[gene_loc]),
            pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak'],
        )
        all_tei_symp = count_fixed_diffs(
            gt=allel.GenotypeArray(callset[gene_loc]),
            pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak_symp'],
        )
        all_tei_allo = count_fixed_diffs(
            gt=allel.GenotypeArray(callset[gene_loc]),
            pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak_allo'],
        )
        all_yak_mel = count_fixed_diffs(
            gt=allel.GenotypeArray(callset[gene_loc]),
            pop_a=idx_dicc['yak'], pop_b=idx_dicc['mel'],
        )
        all_symp_allo = count_fixed_diffs(
            gt=allel.GenotypeArray(callset[gene_loc]),
            pop_a=idx_dicc['yak_symp'], pop_b=idx_dicc['yak_allo'],
        )
        # Append the results.
        all_mat[idx, :] = np.array([
            all_san_mel, all_san_tei, all_san_yak, all_san_symp, all_san_allo,
            all_tei_mel, all_tei_yak, all_tei_symp, all_tei_allo,
            all_yak_mel, all_symp_allo,
        ])
        # Locate the coding mutation indicies.
        syn_idx = np.where(((start <= syn_sites) & (syn_sites <= end)))[0]
        mis_idx = np.where(((start <= mis_sites) & (mis_sites <= end)))[0]
        # If there are no coding mutations.
        if syn_idx.size == 0 and mis_idx.size == 0:
            # Append the results.
            int_mat[idx, :] = np.array([
                all_san_mel, all_san_tei, all_san_yak, all_san_symp, all_san_allo,
                all_tei_mel, all_tei_yak, all_tei_symp, all_tei_allo,
                all_yak_mel, all_symp_allo,
            ])
            syn_mat[idx, :] = np.full(11, np.nan)
            mis_mat[idx, :] = np.full(11, np.nan)
        # Else-if there are only synonmous coding mutations.
        elif syn_idx.size != 0 and mis_idx.size == 0:
            # Extract the position of the coding mutations.
            gene_pos = all_pos[gene_idx]
            syn_pos = syn_sites[syn_idx]
            # Generate mutation masks.
            syn_mask = np.isin(gene_pos, syn_pos)
            int_mask = ~syn_mask
            # Count the number of fixed differences.
            syn_san_mel = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['mel'],
            )
            syn_san_tei = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['tei'],
            )
            syn_san_yak = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['yak'],
            )
            syn_san_symp = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['yak_symp'],
            )
            syn_san_allo = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['yak_allo'],
            )
            syn_tei_mel = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['tei'], pop_b=idx_dicc['mel'],
            )
            syn_tei_yak = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak'],
            )
            syn_tei_symp = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak_symp'],
            )
            syn_tei_allo = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak_allo'],
            )
            syn_yak_mel = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['yak'], pop_b=idx_dicc['mel'],
            )
            syn_symp_allo = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['yak_symp'], pop_b=idx_dicc['yak_allo'],
            )
            # Append the results.
            syn_mat[idx, :] = np.array([
                syn_san_mel, syn_san_tei, syn_san_yak, syn_san_symp, syn_san_allo,
                syn_tei_mel, syn_tei_yak, syn_tei_symp, syn_tei_allo,
                syn_yak_mel, syn_symp_allo,
            ])
            mis_mat[idx, :] = np.full(11, np.nan)
            # If there are no intronic mutations.
            if int_mask.sum() == 0:
                # Append the results.
                int_mat[idx, :] = np.full(11, np.nan)
            # Else there are intronic mutations.
            else:
                # Count the number of fixed differences.
                int_san_mel = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['san'], pop_b=idx_dicc['mel'],
                )
                int_san_tei = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['san'], pop_b=idx_dicc['tei'],
                )
                int_san_yak = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['san'], pop_b=idx_dicc['yak'],
                )
                int_san_symp = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['san'], pop_b=idx_dicc['yak_symp'],
                )
                int_san_allo = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['san'], pop_b=idx_dicc['yak_allo'],
                )
                int_tei_mel = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['tei'], pop_b=idx_dicc['mel'],
                )
                int_tei_yak = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak'],
                )
                int_tei_symp = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak_symp'],
                )
                int_tei_allo = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak_allo'],
                )
                int_yak_mel = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['yak'], pop_b=idx_dicc['mel'],
                )
                int_symp_allo = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['yak_symp'], pop_b=idx_dicc['yak_allo'],
                )
                # Append the results.
                int_mat[idx, :] = np.array([
                    int_san_mel, int_san_tei, int_san_yak, int_san_symp, int_san_allo,
                    int_tei_mel, int_tei_yak, int_tei_symp, int_tei_allo,
                    int_yak_mel, int_symp_allo,
                ])
        # Else-if there are only missense coding mutations.
        elif syn_idx.size == 0 and mis_idx.size != 0:
            # Extract the position of the coding mutations.
            gene_pos = all_pos[gene_idx]
            mis_pos = mis_sites[mis_idx]
            # Generate mutation masks.
            mis_mask = np.isin(gene_pos, mis_pos)
            int_mask = ~mis_mask
            # Count the number of fixed differences.
            mis_san_mel = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['mel'],
            )
            mis_san_tei = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['tei'],
            )
            mis_san_yak = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['yak'],
            )
            mis_san_symp = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['yak_symp'],
            )
            mis_san_allo = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['yak_allo'],
            )
            mis_tei_mel = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['tei'], pop_b=idx_dicc['mel'],
            )
            mis_tei_yak = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak'],
            )
            mis_tei_symp = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak_symp'],
            )
            mis_tei_allo = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak_allo'],
            )
            mis_yak_mel = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['yak'], pop_b=idx_dicc['mel'],
            )
            mis_symp_allo = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['yak_symp'], pop_b=idx_dicc['yak_allo'],
            )
            # Append the results.
            mis_mat[idx, :] = np.array([
                mis_san_mel, mis_san_tei, mis_san_yak, mis_san_symp, mis_san_allo,
                mis_tei_mel, mis_tei_yak, mis_tei_symp, mis_tei_allo,
                mis_yak_mel, mis_symp_allo,
            ])
            syn_mat[idx, :] = np.full(11, np.nan)
            # If there are no intronic mutations.
            if int_mask.sum() == 0:
                # Append the results.
                int_mat[idx, :] = np.full(11, np.nan)
            # Else there are intronic mutations.
            else:
                # Count the number of fixed differences.
                int_san_mel = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['san'], pop_b=idx_dicc['mel'],
                )
                int_san_tei = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['san'], pop_b=idx_dicc['tei'],
                )
                int_san_yak = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['san'], pop_b=idx_dicc['yak'],
                )
                int_san_symp = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['san'], pop_b=idx_dicc['yak_symp'],
                )
                int_san_allo = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['san'], pop_b=idx_dicc['yak_allo'],
                )
                int_tei_mel = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['tei'], pop_b=idx_dicc['mel'],
                )
                int_tei_yak = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak'],
                )
                int_tei_symp = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak_symp'],
                )
                int_tei_allo = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak_allo'],
                )
                int_yak_mel = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['yak'], pop_b=idx_dicc['mel'],
                )
                int_symp_allo = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['yak_symp'], pop_b=idx_dicc['yak_allo'],
                )
                # Append the results.
                int_mat[idx, :] = np.array([
                    int_san_mel, int_san_tei, int_san_yak, int_san_symp, int_san_allo,
                    int_tei_mel, int_tei_yak, int_tei_symp, int_tei_allo,
                    int_yak_mel, int_symp_allo,
                ])
        # Else, there ar both types of coding mutations.
        else:
            # Extract the position of the coding mutations.
            gene_pos = all_pos[gene_idx]
            syn_pos = syn_sites[syn_idx]
            mis_pos = mis_sites[mis_idx]
            # Generate mutation masks.
            syn_mask = np.isin(gene_pos, syn_pos)
            mis_mask = np.isin(gene_pos, mis_pos)
            int_mask = ~(syn_mask | mis_mask)
            # Count the number of fixed differences.
            syn_san_mel = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['mel'],
            )
            syn_san_tei = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['tei'],
            )
            syn_san_yak = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['yak'],
            )
            syn_san_symp = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['yak_symp'],
            )
            syn_san_allo = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['yak_allo'],
            )
            syn_tei_mel = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['tei'], pop_b=idx_dicc['mel'],
            )
            syn_tei_yak = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak'],
            )
            syn_tei_symp = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak_symp'],
            )
            syn_tei_allo = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak_allo'],
            )
            syn_yak_mel = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['yak'], pop_b=idx_dicc['mel'],
            )
            syn_symp_allo = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['yak_symp'], pop_b=idx_dicc['yak_allo'],
            )
            mis_san_mel = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['mel'],
            )
            mis_san_tei = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['tei'],
            )
            mis_san_yak = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['yak'],
            )
            mis_san_symp = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['yak_symp'],
            )
            mis_san_allo = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['san'], pop_b=idx_dicc['yak_allo'],
            )
            mis_tei_mel = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['tei'], pop_b=idx_dicc['mel'],
            )
            mis_tei_yak = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak'],
            )
            mis_tei_symp = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak_symp'],
            )
            mis_tei_allo = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak_allo'],
            )
            mis_yak_mel = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['yak'], pop_b=idx_dicc['mel'],
            )
            mis_symp_allo = count_fixed_diffs(
                gt=allel.GenotypeArray(callset[gene_loc]),
                pop_a=idx_dicc['yak_symp'], pop_b=idx_dicc['yak_allo'],
            )
            # Append the results.
            syn_mat[idx, :] = np.array([
                syn_san_mel, syn_san_tei, syn_san_yak, syn_san_symp, syn_san_allo,
                syn_tei_mel, syn_tei_yak, syn_tei_symp, syn_tei_allo,
                syn_yak_mel, syn_symp_allo,
            ])
            mis_mat[idx, :] = np.array([
                mis_san_mel, mis_san_tei, mis_san_yak, mis_san_symp, mis_san_allo,
                mis_tei_mel, mis_tei_yak, mis_tei_symp, mis_tei_allo,
                mis_yak_mel, mis_symp_allo,
            ])
            # If there are no intronic mutations.
            if int_mask.sum() == 0:
                # Append the results.
                int_mat[idx, :] = np.full(11, np.nan)
            # Else there are intronic mutations.
            else:
                # Count the number of fixed differences.
                int_san_mel = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['san'], pop_b=idx_dicc['mel'],
                )
                int_san_tei = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['san'], pop_b=idx_dicc['tei'],
                )
                int_san_yak = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['san'], pop_b=idx_dicc['yak'],
                )
                int_san_symp = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['san'], pop_b=idx_dicc['yak_symp'],
                )
                int_san_allo = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['san'], pop_b=idx_dicc['yak_allo'],
                )
                int_tei_mel = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['tei'], pop_b=idx_dicc['mel'],
                )
                int_tei_yak = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak'],
                )
                int_tei_symp = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak_symp'],
                )
                int_tei_allo = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['tei'], pop_b=idx_dicc['yak_allo'],
                )
                int_yak_mel = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['yak'], pop_b=idx_dicc['mel'],
                )
                int_symp_allo = count_fixed_diffs(
                    gt=allel.GenotypeArray(callset[gene_loc]),
                    pop_a=idx_dicc['yak_symp'], pop_b=idx_dicc['yak_allo'],
                )
                # Append the results.
                int_mat[idx, :] = np.array([
                    int_san_mel, int_san_tei, int_san_yak, int_san_symp, int_san_allo,
                    int_tei_mel, int_tei_yak, int_tei_symp, int_tei_allo,
                    int_yak_mel, int_symp_allo,
                ])

# Export the the results matricies.
np.savetxt(
    f'./results/{ref}/{miss}/fixed_diffs_all_orthologs.txt.gz',
    all_mat, fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/fixed_diffs_int_orthologs.txt.gz',
    int_mat, fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/fixed_diffs_syn_orthologs.txt.gz',
    syn_mat, fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/fixed_diffs_mis_orthologs.txt.gz',
    mis_mat, fmt='%1.15f', delimiter='\t', newline='\n',
)