# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr

# sys.argv[1] = reference


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
def calc_pbs(gt, pop_a, pop_b, pop_c):
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
    # Correct for negative Fst values.
    a_b_fst = np.where(a_b_raw_fst < 0, 0, a_b_raw_fst)
    a_c_fst = np.where(a_c_raw_fst < 0, 0, a_c_raw_fst)
    c_b_fst = np.where(c_b_raw_fst < 0, 0, c_b_raw_fst)
    # Calculate PBS.
    pbs = (
        ((np.log(1.0 - a_b_fst) * -1.0) +\
         (np.log(1.0 - a_c_fst) * -1.0) -\
         (np.log(1.0 - c_b_fst) * -1.0)) / float(2)
    )
    return pbs, np.nanmean(pbs)

# Define a function to perform bootstrapping.
def pbs_bootstraps(reference):
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
    # Intialize a dictionary of chromosome lengths.
    if reference == 'san':
        chromosome_dicc = {
            '2L': 27915280, '2R': 22646313,
            '3L': 25621868, '3R': 29038104,
            '4': 1429326, 'X': 22631256,
        }
    elif reference == 'yak':
        chromosome_dicc = {
            '2L': 31052931, '2R': 23815334,
            '3L': 25180761, '3R': 30730773,
            '4': 1429802, 'X': 24674056,
        }
    # Intialize a counter. 
    idx = 0
    # Intialize a results matrix to store the results.
    results_mat = np.empty((10_000, 7))
    # For one billion tries....
    for _ in range(1_000_000_000):
        # Randomly select a chromosome.
        chromosome = np.random.choice(list(chromosome_dicc.keys()))
        # Randomly generate a start and end position.
        start = np.random.randint((chromosome_dicc[chromosome] - (5_000 - 1)))
        end = start + 5_000
        # Extract the genotype callset and positions.
        callset, all_pos = load_callset_pos(reference, chromosome)
        # Determine the positions in this window.
        wind_idx = np.where(((start <= all_pos) & (all_pos <= end)))[0]
        # If there are no sites in this window...
        if (wind_idx.size == 0):
            continue
        # Else...
        else:
            # Identify the window to extract.
            wind_loc = all_pos.locate_range(start, end)
            # Calculate pbs.
            _, san_yak_symp_yak_allo = calc_pbs(
                gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['san'],
                pop_b=idx_dicc['yak_symp'], pop_c=idx_dicc['yak_allo'],
            )
            _, san_yak_symp_tei = calc_pbs(
                gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['san'],
                pop_b=idx_dicc['yak_symp'], pop_c=idx_dicc['tei'],
            )
            _, san_yak_allo_tei = calc_pbs(
                gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['san'],
                pop_b=idx_dicc['yak_allo'], pop_c=idx_dicc['tei'],
            )
            _, yak_symp_yak_allo_tei = calc_pbs(
                gt=allel.GenotypeArray(callset[wind_loc]), pop_a=idx_dicc['yak_symp'],
                pop_b=idx_dicc['yak_allo'], pop_c=idx_dicc['tei'],
            )
            # Append the results.
            results_mat[idx, :] = np.array([
                chromosome, start, end,
                san_yak_symp_yak_allo, san_yak_symp_tei,
                san_yak_allo_tei, yak_symp_yak_allo_tei,
            ])
            # Move the index counter forward.
            idx += 1
            # If we just finished the 10000th window.
            if (idx == 10_000):
                break
            # Else...
            else:
                # Continue until you have bootstrapped 10000 windows.
                continue
    # Export the the results matrix.
    np.savetxt(
        f'./results/{refrence}/bootstrap_pbs.csv',
        results_mat, fmt='%s', delimiter=',', newline='\n',
    )
    return


# Conduct bootstrapping.
pbs_bootstraps(reference=str(sys.argv[1]))