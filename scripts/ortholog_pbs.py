# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

# sys.argv[1] = refrence

# Ignore divide by 0/np.nan error and encode as np.nan's.
np.seterr(divide='ignore', invalid='ignore')

# Define a function to load genotyope and positions arrays.
def load_callset_pos(prefix, chrom):
    # Intialize the file path.
    path = f'../zarrs/{prefix}/chr{chrom}.zarr'
    # Load the vcf file.
    callset = zarr.open_group(path, mode='r')
    # Extract the genotypes.
    geno = callset[f'{chrom}/calldata/GT']
    # Load the positions.
    pos = allel.SortedIndex(callset[f'{chrom}/variants/POS'])
    return geno, pos

# Define a function to calculate PBS.
def calc_pbs_n1(gt, pop_a, pop_b, pop_c):
    # Determine allele counts.
    a_ac = gt.take(pop_a, axis=1).count_alleles()
    b_ac = gt.take(pop_b, axis=1).count_alleles()
    c_ac = gt.take(pop_c, axis=1).count_alleles()
    # Calculate the numerator and denominator for Hudson's Fst estimator.
    a_b_num, a_b_den = allel.hudson_fst(a_ac, b_ac)
    a_c_num, a_c_den = allel.hudson_fst(a_ac, c_ac)
    c_b_num, c_b_den = allel.hudson_fst(c_ac, b_ac)
    # Calculate Fst.
    a_b_raw_fst = a_b_num / a_b_den
    a_c_raw_fst = a_c_num / a_c_den
    c_b_raw_fst = c_b_num / c_b_den
    # Correct for Fst values of 1 that will lead to inf.
    a_b_fst_raw = np.where(a_b_raw_fst == 1, 0.99999, a_b_raw_fst)
    a_c_fst_raw = np.where(a_c_raw_fst == 1, 0.99999, a_c_raw_fst)
    c_b_fst_raw = np.where(c_b_raw_fst == 1, 0.99999, c_b_raw_fst)
    # Correct for negative Fst values.
    a_b_fst = np.where(a_b_fst_raw < 0, 0, a_b_fst_raw)
    a_c_fst = np.where(a_c_fst_raw < 0, 0, a_c_fst_raw)
    c_b_fst = np.where(c_b_fst_raw < 0, 0, c_b_fst_raw)
    # Calculate divergence.
    a_b_t = -np.log(1 - a_b_fst)
    a_c_t = -np.log(1 - a_c_fst)
    c_b_t = -np.log(1 - c_b_fst)
    # Calculate the raw PBS.
    raw_pbs = (a_b_t + a_c_t - c_b_t) / 2
    # Calculate the normalization factor.
    norm = 1 + (a_b_t + a_c_t + c_b_t) / 2
    # Calculate PBS.
    pbs = raw_pbs / norm
    # If pbs is an empty array...
    if pbs.size == 0:
        # Return the results.
        return pbs, np.nan
    # Else...
    else:
        return pbs, np.nanmean(pbs)

# Define a function to calculate pbs for each orthologous gene.
def pbs_genes(refrence):
    # Read in meta data has a pandas dataframe.
    meta_df = pd.read_csv(
        './matute_meta_data.txt',
        sep='\t',
        names=['strain', 'population'],
    )
    # Intialize a dictionary.
    idx_dicc = {}
    # For all populations.
    for pop in ['san', 'yak_symp', 'yak_allo', 'tei']:
        # Fill the dictionary.
        idx_dicc[pop] = meta_df[meta_df['population'] == pop].index.values
    # Fill the dictionary for the last population.
    idx_dicc['yak'] = np.concatenate((idx_dicc['yak_allo'], idx_dicc['yak_symp']))
    # Load the ortholog data frame.
    ortho_df = pd.read_csv(
        f'../annotations/{refrence}_ortholog_qc_passed_coords.txt', sep='\t',
        names=['gene_id', 'chr', 'start', 'end', 'length', 'tot_sites', 'seg_sites'],
    )
    # Extract the chromosomes, start, and end positions.
    chromosomes = ortho_df['chr'].values
    starts = ortho_df['start'].values
    ends = ortho_df['end'].values
    # Intialize a results matrix to store the results.
    results_mat = np.empty((ortho_df.shape[0], 5))
    # For every gene.
    for gene in range(ortho_df.shape[0]):
        # Grab the chromsome, start, and end positions.
        chromosome = chromosomes[gene]
        start = starts[gene]
        end = ends[gene]
        # Extract the genotype callset and positions.
        callset, all_pos = load_callset_pos(refrence, chromosome)
        # Locate the gene.
        gene_loc = all_pos.locate_range(start, end)
       # Calculate pbs.
        _, san_yak_tei = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[gene_loc]), pop_a=idx_dicc['san'],
            pop_b=idx_dicc['yak'], pop_c=idx_dicc['tei'],
        )
        _, san_yak_symp_yak_allo = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[gene_loc]), pop_a=idx_dicc['san'],
            pop_b=idx_dicc['yak_symp'], pop_c=idx_dicc['yak_allo'],
        )
        _, san_yak_symp_tei = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[gene_loc]), pop_a=idx_dicc['san'],
            pop_b=idx_dicc['yak_symp'], pop_c=idx_dicc['tei'],
        )
        _, san_yak_allo_tei = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[gene_loc]), pop_a=idx_dicc['san'],
            pop_b=idx_dicc['yak_allo'], pop_c=idx_dicc['tei'],
        )
        _, yak_symp_yak_allo_tei = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[gene_loc]), pop_a=idx_dicc['yak_symp'],
            pop_b=idx_dicc['yak_allo'], pop_c=idx_dicc['tei'],
        )
        # Append the result matrix.
        results_mat[gene, :] = np.array([
            san_yak_tei,
            san_yak_symp_yak_allo, san_yak_symp_tei,
            san_yak_allo_tei, yak_symp_yak_allo_tei,
        ])
    # Export the the results matrix.
    np.savetxt(
        f'./results/{refrence}/ortholog_pbs.csv',
        results_mat, fmt='%1.15f', delimiter=',', newline='\n',
    )
    return

# Calculate pbs.
pbs_genes(refrence=str(sys.argv[1]))
