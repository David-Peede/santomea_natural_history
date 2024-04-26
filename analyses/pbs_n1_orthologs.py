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
    a_b_fst = np.nansum(a_b_num) / np.nansum(a_b_den)
    a_c_fst = np.nansum(a_c_num) / np.nansum(a_c_den)
    c_b_fst = np.nansum(c_b_num) / np.nansum(c_b_den)
    # Correct for Fst values of 1 that will lead to inf.
    if a_b_fst == 1:
        a_b_fst = 0.99999
    if a_c_fst == 1:
        a_c_fst = 0.99999
    if c_b_fst == 1:
        c_b_fst = 0.99999
    # Calculate divergence.
    a_b_t = -np.log(1 - a_b_fst)
    a_c_t = -np.log(1 - a_c_fst)
    c_b_t = -np.log(1 - c_b_fst)
    # Calculate the PBS per focal population.
    pbs_a = (a_b_t + a_c_t - c_b_t) / 2
    pbs_b = (a_b_t + c_b_t - a_c_t) / 2
    pbs_c = (a_c_t + c_b_t - a_b_t) / 2
    # Calculate the normalization factor.
    norm = 1 + (pbs_a + pbs_b + pbs_c)
    # Calculate PBS.
    pbs = pbs_a / norm
    return pbs

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
results_mat = np.empty((ortholog_df.shape[0], 17))

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
        results_mat[idx, :] = np.full(17, np.nan)
    # Else there are sites to do computations on.
    else:
        # Identify the window to extract.
        wind_loc = all_pos.locate_range(start, end)
        # Calculate pbs.
        san_tei_mel = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['san'],
            pop_b=idx_dicc['tei'], pop_c=idx_dicc['mel'],
        )
        san_yak_mel = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['san'],
            pop_b=idx_dicc['yak'], pop_c=idx_dicc['mel'],
        )
        san_allo_mel = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['san'],
            pop_b=idx_dicc['yak_allo'], pop_c=idx_dicc['mel'],
        )
        san_symp_mel = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['san'],
            pop_b=idx_dicc['yak_symp'], pop_c=idx_dicc['mel'],
        )
        san_yak_tei = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['san'],
            pop_b=idx_dicc['yak'], pop_c=idx_dicc['tei'],
        )
        san_allo_tei = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['san'],
            pop_b=idx_dicc['yak_allo'], pop_c=idx_dicc['tei'],
        )
        san_symp_tei = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['san'],
            pop_b=idx_dicc['yak_symp'], pop_c=idx_dicc['tei'],
        )
        san_symp_allo = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['san'],
            pop_b=idx_dicc['yak_symp'], pop_c=idx_dicc['yak_allo'],
        )
        tei_san_mel = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['tei'],
            pop_b=idx_dicc['san'], pop_c=idx_dicc['mel'],
        )
        tei_yak_mel = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['tei'],
            pop_b=idx_dicc['yak'], pop_c=idx_dicc['mel'],
        )
        tei_yak_san = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['tei'],
            pop_b=idx_dicc['yak'], pop_c=idx_dicc['san'],
        )
        symp_allo_mel = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['yak_symp'],
            pop_b=idx_dicc['yak_allo'], pop_c=idx_dicc['mel'],
        )
        symp_allo_tei = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['yak_symp'],
            pop_b=idx_dicc['yak_allo'], pop_c=idx_dicc['tei'],
        )
        symp_allo_san = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['yak_symp'],
            pop_b=idx_dicc['yak_allo'], pop_c=idx_dicc['san'],
        )
        yak_san_tei = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['yak'],
            pop_b=idx_dicc['san'], pop_c=idx_dicc['tei'],
        )
        yak_tei_mel = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['yak'],
            pop_b=idx_dicc['tei'], pop_c=idx_dicc['mel'],
        )
        yak_san_mel = calc_pbs_n1(
            gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['yak'],
            pop_b=idx_dicc['san'], pop_c=idx_dicc['mel'],
        )
        # Append the results.
        results_mat[idx, :] = np.array([
            san_tei_mel, san_yak_mel, san_allo_mel, san_symp_mel,
            san_yak_tei, san_allo_tei, san_symp_tei, san_symp_allo,
            tei_san_mel, tei_yak_mel, tei_yak_san,
            symp_allo_mel, symp_allo_tei, symp_allo_san,
            yak_san_tei, yak_tei_mel, yak_san_mel,
        ])
        
# Export the the results matrix.
np.savetxt(
    f'./results/{ref}/{miss}/pbs_n1_orthologs.txt.gz',
    results_mat, fmt='%1.15f', delimiter='\t', newline='\n',
)