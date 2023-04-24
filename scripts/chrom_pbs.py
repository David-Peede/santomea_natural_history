# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

# sys.argv[1] = refrence
# sys.argv[1] = chromosome

# Ignore divide by 0/np.nan error and encode as np.nan's.
np.seterr(divide='ignore', invalid='ignore')

# Define a function to load genotyope and positions arrays.
def load_gt_pos(prefix, chrom):
    # Intialize the file path.
    path = f'../zarrs/{prefix}/chr{chrom}.zarr'
    # Load the vcf file.
    callset = zarr.open_group(path, mode='r')
    # Convert the genotype callset to an array.
    gt = allel.GenotypeArray(callset[f'{chrom}/calldata/GT'])
    # Load the positions.
    pos = allel.SortedIndex(callset[f'{chrom}/variants/POS'])
    return gt, pos

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

# Define a function to calculate pbs per snp per chromosome.
def pbs_genome(refrence, chromosome):
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
    # Extract the genotype callset and positions.
    gt, _ = load_gt_pos(refrence, chromosome)
    # Calculate pbs.
    san_yak_tei, _ = calc_pbs_n1(
        gt=gt, pop_a=idx_dicc['san'],
        pop_b=idx_dicc['yak'], pop_c=idx_dicc['tei'],
    )
    san_yak_symp_yak_allo, _ = calc_pbs_n1(
        gt=gt, pop_a=idx_dicc['san'],
        pop_b=idx_dicc['yak_symp'], pop_c=idx_dicc['yak_allo'],
    )
    san_yak_symp_tei, _ = calc_pbs_n1(
        gt=gt, pop_a=idx_dicc['san'],
        pop_b=idx_dicc['yak_symp'], pop_c=idx_dicc['tei'],
    )
    san_yak_allo_tei, _ = calc_pbs_n1(
        gt=gt, pop_a=idx_dicc['san'],
        pop_b=idx_dicc['yak_allo'], pop_c=idx_dicc['tei'],
    )
    yak_symp_yak_allo_tei, _ = calc_pbs_n1(
        gt=gt, pop_a=idx_dicc['yak_symp'],
        pop_b=idx_dicc['yak_allo'], pop_c=idx_dicc['tei'],
    )
    # Export the the results matrix.
    np.savetxt(
        f'./results/{refrence}/san_yak_tei_chr{chromosome}_pbs.csv',
        [san_yak_tei], fmt='%1.15f', delimiter=',', newline='\n',
    )
    np.savetxt(
        f'./results/{refrence}/san_yak_symp_yak_allo_chr{chromosome}_pbs.csv',
        [san_yak_symp_yak_allo], fmt='%1.15f', delimiter=',', newline='\n',
    )
    np.savetxt(
        f'./results/{refrence}/san_yak_symp_tei_chr{chromosome}_pbs.csv',
        [san_yak_symp_tei], fmt='%1.15f', delimiter=',', newline='\n',
    )
    np.savetxt(
        f'./results/{refrence}/san_yak_allo_tei_chr{chromosome}_pbs.csv',
        [san_yak_allo_tei], fmt='%1.15f', delimiter=',', newline='\n',
    )
    np.savetxt(
        f'./results/{refrence}/yak_symp_yak_allo_tei_chr{chromosome}_pbs.csv',
        [yak_symp_yak_allo_tei], fmt='%1.15f', delimiter=',', newline='\n',
    )
    return


# Calculate pbs.
pbs_genome(refrence=str(sys.argv[1]), chromosome=str(sys.argv[2]))
