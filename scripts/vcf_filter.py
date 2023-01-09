# Import packages.
import gzip
import sys

# sys.argv[1] = sample prefix
# sys.argv[2] = species

def vcf_filter(vcf, refernce_genome, prefix, header_1, header_2):
    """
    vcf: VCF filtered for DP and MQ
    --------------------------------------------------------------------------
    output: A VCF filtered for SNPs with a GQ of 40 or higher
    """
    # Intilize the output.
    outfile = sys.stdout
    # Intialize a naming dictionary for the chromosomes.
    if (refernce_genome == 'san'):
        chrom_dicc = {
            'NC_053016.2': '2L', 'NC_053017.2': '2R',
            'NC_053018.2': '3L', 'NC_053019.2': '3R',
            'NC_053020.2': '4', 'NC_053021.2': 'X',
        }
        # Write the header.
        outfile.write(header_1+'\n')
    elif (refernce_genome == 'yak'):
        chrom_dicc = {
            'NC_052527.2': '2L', 'NC_052528.2': '2R',
            'NC_052529.2': '3L', 'NC_052530.2': '3R',
            'NC_052531.2': '4', 'NC_052526.2': 'X',
        }
        # Write the header.
        outfile.write(header_2+'\n')
    # Open the gzipped vcf.
    with gzip.open(vcf, 'rt') as data:
        # For every line in the file.
        for line in data:
            # If the line is a meta data line.
            if line.startswith('##'):
                # Continue to the next line.
                continue
            # Else-if the line is is the header line.
            elif line.startswith('#'):
                # Split the line by tabs.
                spline = line.split()
                # Change the sample name.
                spline[-1] = prefix
                # Write it to output.
                outfile.write('\t'.join(spline)+'\n')
            else:
                # Split the line by tabs.
                spline = line.split()
                # Grab specific column info.
                chrom = spline[0]
                ref = spline[3]
                alt = spline[4]
                info = spline[7]
                format_field = spline[8]
                geno = spline[9]
                # If the site is multi-allelic or and indel.
                if (len(ref) + len(alt)) > 2:
                    # Continue to the next line.
                    continue
                # Else-if the site was called an indel.
                elif 'INDEL' in format_field:
                    # Continue to the next line.
                    continue
                # Else-if the site is not for one of our chromosomes of interest.
                elif chrom not in chrom_dicc:
                    # Continue to the next line.
                    continue
                # Else.
                else:
                    # Change the chromosome naming such that it doesn't throw errors.
                    spline[0] = chrom_dicc[chrom]
                    # If the formar field only has GT.
                    if len(format_field) == 2:
                        # Write it to output.
                        outfile.write('\t'.join(spline)+'\n')
                    # Else.
                    else:
                        # Extract the format information.
                        gt, pl, gq = geno.split(':')
                        # If the genotype quality is less than 40.
                        if (int(gq) < 40):
                            # Continue to the next line
                            continue
                        # Else.
                        else:
                            # Write it to output.
                            outfile.write('\t'.join(spline)+'\n')
    return

# Intialize new headers
san_header = '''##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##bcftoolsVersion=1.13+htslib-1.13
##bcftoolsCommand=mpileup --threads 4 -f ../assemblies/Dsan.genomic.fna ./san/Zimbabwe_sorted_chrom_subset.bam
##reference=file://../assemblies/Dsan.genomic.fna
##contig=<ID=2L,length=27915280>
##contig=<ID=2R,length=22646313>
##contig=<ID=3L,length=25621868>
##contig=<ID=3R,length=29038104>
##contig=<ID=4,length=1429326>
##contig=<ID=x,length=22631256>
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of raw reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of raw reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
##INFO=<ID=RPBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Read Position Bias (closer to 0 is better)">
##INFO=<ID=MQBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Mapping Quality Bias (closer to 0 is better)">
##INFO=<ID=BQBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Base Quality Bias (closer to 0 is better)">
##INFO=<ID=MQSBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Mapping Quality vs Strand Bias (closer to 0 is better)">
##INFO=<ID=SCBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Soft-Clip Length Bias (closer to 0 is better)">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Phred-scaled Genotype Quality">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">'''

yak_header = '''##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##bcftoolsVersion=1.13+htslib-1.13
##bcftoolsCommand=mpileup --threads 4 -f ../assemblies/Dyak.genomic.fna ./yak/Zimbabwe_sorted_chrom_subset.bam
##reference=file://../assemblies/Dyak.genomic.fna
##contig=<ID=X,length=24674056>
##contig=<ID=2L,length=31052931>
##contig=<ID=2R,length=23815334>
##contig=<ID=3L,length=25180761>
##contig=<ID=3R,length=30730773>
##contig=<ID=4,length=1429802>
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of raw reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of raw reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
##INFO=<ID=RPBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Read Position Bias (closer to 0 is better)">
##INFO=<ID=MQBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Mapping Quality Bias (closer to 0 is better)">
##INFO=<ID=BQBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Base Quality Bias (closer to 0 is better)">
##INFO=<ID=MQSBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Mapping Quality vs Strand Bias (closer to 0 is better)">
##INFO=<ID=SCBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Soft-Clip Length Bias (closer to 0 is better)">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Phred-scaled Genotype Quality">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">'''

# Compile the vcf file.
vcf_file = '{0}_dp_mq_filtered.vcf.gz'.format(str(sys.argv[1]))

# Filter the vcf.
vcf_filter(vcf_file, str(sys.argv[2]), str(sys.argv[1]), san_header, yak_header)