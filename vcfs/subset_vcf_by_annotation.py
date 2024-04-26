# Import packages.
import gzip
import sys

### sys.argv[1] = vcf file ###
### sys.argv[2] = chromosome ###
### sys.argv[3] = annotation type ###


# Define function to filter the joint vcf file by a specified annotation.
def vcf_annotation_subset(vcf, chrom, ann):
    """
    ###########################################################################
    INPUT
        vcf: A gzipped annotated VCF file.
        ann: Annotation type.
    ---------------------------------------------------------------------------
    OUTPUT: Filtered VCF file to standard out.
    ###########################################################################
    """
    # Intialize the filtered vcf file to be outputted to stdout.
    new_vcf = sys.stdout
    # Intialize a report file.
    report = open(f'./mel_annotations/mel_75_{ann}_chr{chrom}_report.txt', 'w')
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
                # Split the line by tabs.
                spline = line.split()
                # Grab the position.
                pos = spline[1]
                # Grab the info field.
                info = spline[7]
                # If the info field contains the annotation.
                if ann in info:
                    # Write it to the new vcf.
                    new_vcf.write(line)
                    # Write to the report file.
                    report.write('\t'.join([chrom, pos, info])+'\n')
    # Close the report file.
    report.close()
    return

# Subset the VCF file by annotation type.
vcf_annotation_subset(vcf=str(sys.argv[1]), chrom=str(sys.argv[2]), ann=str(sys.argv[3]))