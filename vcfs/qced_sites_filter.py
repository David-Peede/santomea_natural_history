# Import packages.
import gzip
import pandas as pd
import sys


# sys.argv[1] = reference
# sys.argv[2] = chromosome
# sys.argv[3] = missingness threshold



# Define a function to filter a vcf file given a list of known qc'ed positions.
def qced_sites_filter(vcf, ref, chrom, miss):
    # Read in the qc'ed sites information.
    qc_df = pd.read_csv(
        f'./{ref}_{miss}_qc_passed_sites_chr{chrom}.txt',
        sep='\t', names=['chrom', 'pos'],
    )
    # Intialize a hash table with all the qc'ed sites.
    qced_sites = set(qc_df.pos.values)
    # Intialize the filtered vcf file to be outputted to stdout.
    new_vcf = sys.stdout
    # Open the original vcf file.
    with gzip.open(vcf, 'rt') as data:
        # Iterate through every line in the original vcf file.
        for line in data:
            # If the line is a meta info line...
            if line.startswith('#'):
                # Write it to the new vcf.
                new_vcf.write(line)
            # Else...
            else:
                # Grab the position.
                pos = int(line.split()[1])
                # If the current position is a qc'ed site.
                if pos in qced_sites:
                    # Write the this position's entry to the new vcf file.
                    new_vcf.write(line)
                    # Remove the current position from the hash table.
                    qced_sites.remove(pos)
    return

# Perform the final filtering step.
qced_sites_filter(
    vcf=f'./{sys.argv[1]}/san_yak_tei_mel_multi_filtered_chr{sys.argv[2]}.vcf.gz', 
    ref=str(sys.argv[1]),
    chrom=str(sys.argv[2]),
    miss=str(sys.argv[3]),
)