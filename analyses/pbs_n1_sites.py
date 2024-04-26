# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

# sys.argv[1] = ref
# sys.argv[2] = missing data threshold
# sys.argv[3] = chromosome


# Define a function to load genotyope and positions arrays.
def load_gt_pos(prefix, miss, chrom):
    # Intialize the file path.
    path = f'../zarrs/{miss}/{prefix}/chr{chrom}.zarr'
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
    a_b_fst = np.where(a_b_raw_fst == 1, 0.99999, a_b_raw_fst)
    a_c_fst = np.where(a_c_raw_fst == 1, 0.99999, a_c_raw_fst)
    c_b_fst = np.where(c_b_raw_fst == 1, 0.99999, c_b_raw_fst)
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
    # If pbs is an empty array...
    if pbs.size == 0:
        # Return the results.
        return pbs, np.nan
    # Else...
    else:
        return pbs, np.nanmean(pbs)

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
chrom = str(sys.argv[3])

# Extract the genotype callset and positions.
gt, _ = load_gt_pos(ref, miss, chrom)

# Calculate pbs.
san_tei_mel = calc_pbs_n1(
    gt=gt, pop_a=idx_dicc['san'],
    pop_b=idx_dicc['tei'], pop_c=idx_dicc['mel'],
)
san_yak_mel = calc_pbs_n1(
    gt=gt, pop_a=idx_dicc['san'],
    pop_b=idx_dicc['yak'], pop_c=idx_dicc['mel'],
)
san_allo_mel = calc_pbs_n1(
    gt=gt, pop_a=idx_dicc['san'],
    pop_b=idx_dicc['yak_allo'], pop_c=idx_dicc['mel'],
)
san_symp_mel = calc_pbs_n1(
    gt=gt, pop_a=idx_dicc['san'],
    pop_b=idx_dicc['yak_symp'], pop_c=idx_dicc['mel'],
)
san_yak_tei = calc_pbs_n1(
    gt=gt, pop_a=idx_dicc['san'],
    pop_b=idx_dicc['yak'], pop_c=idx_dicc['tei'],
)
san_allo_tei = calc_pbs_n1(
    gt=gt, pop_a=idx_dicc['san'],
    pop_b=idx_dicc['yak_allo'], pop_c=idx_dicc['tei'],
)
san_symp_tei = calc_pbs_n1(
    gt=gt, pop_a=idx_dicc['san'],
    pop_b=idx_dicc['yak_symp'], pop_c=idx_dicc['tei'],
)
san_symp_allo = calc_pbs_n1(
    gt=gt, pop_a=idx_dicc['san'],
    pop_b=idx_dicc['yak_symp'], pop_c=idx_dicc['yak_allo'],
)
tei_san_mel = calc_pbs_n1(
    gt=gt, pop_a=idx_dicc['tei'],
    pop_b=idx_dicc['san'], pop_c=idx_dicc['mel'],
)
tei_yak_mel = calc_pbs_n1(
    gt=gt, pop_a=idx_dicc['tei'],
    pop_b=idx_dicc['yak'], pop_c=idx_dicc['mel'],
)
tei_yak_san = calc_pbs_n1(
    gt=gt, pop_a=idx_dicc['tei'],
    pop_b=idx_dicc['yak'], pop_c=idx_dicc['san'],
)
symp_allo_mel = calc_pbs_n1(
    gt=gt, pop_a=idx_dicc['yak_symp'],
    pop_b=idx_dicc['yak_allo'], pop_c=idx_dicc['mel'],
)
symp_allo_tei = calc_pbs_n1(
    gt=gt, pop_a=idx_dicc['yak_symp'],
    pop_b=idx_dicc['yak_allo'], pop_c=idx_dicc['tei'],
)
symp_allo_san = calc_pbs_n1(
    gt=gt, pop_a=idx_dicc['yak_symp'],
    pop_b=idx_dicc['yak_allo'], pop_c=idx_dicc['san'],
)
yak_san_tei = calc_pbs_n1(
    gt=gt, pop_a=idx_dicc['yak'],
    pop_b=idx_dicc['san'], pop_c=idx_dicc['tei'],
)
yak_tei_mel = calc_pbs_n1(
    gt=gt, pop_a=idx_dicc['yak'],
    pop_b=idx_dicc['tei'], pop_c=idx_dicc['mel'],
)
yak_san_mel = calc_pbs_n1(
    gt=gt, pop_a=idx_dicc['yak'],
    pop_b=idx_dicc['san'], pop_c=idx_dicc['mel'],
)

# Export the the results matrix.
np.savetxt(
    f'./results/{ref}/{miss}/san_tei_mel_chr{chrom}_pbs_n1.txt.gz',
    [san_tei_mel], fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/san_yak_mel_chr{chrom}_pbs_n1.txt.gz',
    [san_yak_mel], fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/san_allo_mel_chr{chrom}_pbs_n1.txt.gz',
    [san_allo_mel], fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/san_symp_mel_chr{chrom}_pbs_n1.txt.gz',
    [san_symp_mel], fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/san_yak_tei_chr{chrom}_pbs_n1.txt.gz',
    [san_yak_tei], fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/san_allo_tei_chr{chrom}_pbs_n1.txt.gz',
    [san_allo_tei], fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/san_symp_tei_chr{chrom}_pbs_n1.txt.gz',
    [san_symp_tei], fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/san_symp_allo_chr{chrom}_pbs_n1.txt.gz',
    [san_symp_allo], fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/tei_san_mel_chr{chrom}_pbs_n1.txt.gz',
    [tei_san_mel], fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/tei_yak_mel_chr{chrom}_pbs_n1.txt.gz',
    [tei_yak_mel], fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/tei_yak_san_chr{chrom}_pbs_n1.txt.gz',
    [tei_yak_san], fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/symp_allo_mel_chr{chrom}_pbs_n1.txt.gz',
    [symp_allo_mel], fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/symp_allo_tei_chr{chrom}_pbs_n1.txt.gz',
    [symp_allo_tei], fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/symp_allo_san_chr{chrom}_pbs_n1.txt.gz',
    [symp_allo_san], fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/yak_san_tei_chr{chrom}_pbs_n1.txt.gz',
    [yak_san_tei], fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/yak_tei_mel_chr{chrom}_pbs_n1.txt.gz',
    [yak_tei_mel], fmt='%1.15f', delimiter='\t', newline='\n',
)
np.savetxt(
    f'./results/{ref}/{miss}/yak_san_mel_chr{chrom}_pbs_n1.txt.gz',
    [yak_san_mel], fmt='%1.15f', delimiter='\t', newline='\n',
)