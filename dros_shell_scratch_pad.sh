### PHYLOGEN ###

# Download files from NCBI REFSEQ.
sbatch -J yak_faa -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/yak_faa-%A.out -e logs/yak_faa-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_yakuba/latest_assembly_versions/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1_translated_cds.faa.gz"
sbatch -J yak_fna -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/yak_fna-%A.out -e logs/yak_fna-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_yakuba/latest_assembly_versions/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1_cds_from_genomic.fna.gz"
sbatch -J yak_gff -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/yak_gff-%A.out -e logs/yak_gff-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_yakuba/latest_assembly_versions/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1_genomic.gff.gz"
sbatch -J san_faa -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/san_faa-%A.out -e logs/san_faa-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_santomea/latest_assembly_versions/GCF_016746245.2_Prin_Dsan_1.1/GCF_016746245.2_Prin_Dsan_1.1_translated_cds.faa.gz"
sbatch -J san_fna -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/san_fna-%A.out -e logs/san_fna-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_santomea/latest_assembly_versions/GCF_016746245.2_Prin_Dsan_1.1/GCF_016746245.2_Prin_Dsan_1.1_cds_from_genomic.fna.gz"
sbatch -J san_gff -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/san_gff-%A.out -e logs/san_gff-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_santomea/latest_assembly_versions/GCF_016746245.2_Prin_Dsan_1.1/GCF_016746245.2_Prin_Dsan_1.1_genomic.gff.gz"
sbatch -J tei_faa -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/tei_faa-%A.out -e logs/tei_faa-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_teissieri/latest_assembly_versions/GCF_016746235.2_Prin_Dtei_1.1/GCF_016746235.2_Prin_Dtei_1.1_translated_cds.faa.gz"
sbatch -J tei_fna -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/tei_fna-%A.out -e logs/tei_fna-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_teissieri/latest_assembly_versions/GCF_016746235.2_Prin_Dtei_1.1/GCF_016746235.2_Prin_Dtei_1.1_cds_from_genomic.fna.gz"
sbatch -J tei_gff -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/tei_gff-%A.out -e logs/tei_gff-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_teissieri/latest_assembly_versions/GCF_016746235.2_Prin_Dtei_1.1/GCF_016746235.2_Prin_Dtei_1.1_genomic.gff.gz"
sbatch -J mel_faa -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/mel_faa-%A.out -e logs/mel_faa-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_melanogaster/latest_assembly_versions/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_translated_cds.faa.gz"
sbatch -J mel_fna -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/mel_fna-%A.out -e logs/mel_fna-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_melanogaster/latest_assembly_versions/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_cds_from_genomic.fna.gz"
sbatch -J mel_gff -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/mel_gff-%A.out -e logs/mel_gff-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_melanogaster/latest_assembly_versions/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.gz"


# Rename the files.
cat species_prefixes.txt | while read prefix species ; do gunzip -c ${prefix}translated_cds.faa.gz > ${species}.AA.fa ; done
cat species_prefixes.txt | while read prefix species ; do gunzip -c ${prefix}cds_from_genomic.fna.gz > ${species}.CDS.fa ; done
cat species_prefixes.txt | while read prefix species ; do mv ${prefix}genomic.gff.gz ${species}.genomic.gff.gz ; done


# Make a list of protein lengths.
for dros in Dmel Dtei Dyak Dsan; do
seqkit fx2tab -n -l ${dros}.AA.fa | awk 'BEGIN{FS="\t"}{gsub("^.*protein_id=", "", $1); gsub("].*$", "", $1); print}' | sort > ${dros}.AA_len.tsv
done


# Make a file linking gene name in the first column and the protein name in the second column.
for dros in Dmel Dtei Dyak Dsan; do
grep "^>" ${dros}.AA.fa | awk 'BEGIN{FS="\t"}{gsub("^.*GeneID:", "", $1); gsub("].*protein_id=", "\t", $1); gsub("].*$", "", $1); print}' | sort -k 2,2 | uniq > ${dros}.gene_AA.tsv
done


# Determine the longest protein for each gene.
for dros in Dmel Dtei Dyak Dsan; do
join -1 2 ${dros}.gene_AA.tsv -2 1 ${dros}.AA_len.tsv | awk 'len[$2]<$3 {len[$2]=$3; id[$2]=$1} END {for (i in id) {print id[i]}}' > ${dros}.longest_AA.txt
done


# Subset out the longest protein per gene.
for dros in Dmel Dtei Dyak Dsan; do
sbatch -J ${dros}_longest_aa -N 1 -n 1 -t 3:00:00 --mem=2G -o logs/${dros}_longest_aa-%A.out -e logs/${dros}_longest_aa-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="modlue load seqkit; seqkit grep -r -n -f ${dros}.longest_AA.txt ${dros}.AA.fa > ${dros}.longest.AA.fa"
sbatch -J ${dros}_longest_cds -N 1 -n 1 -t 3:00:00 --mem=2G -o logs/${dros}_longest_cds-%A.out -e logs/${dros}_longest_cds-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="modlue load seqkit; seqkit grep -r -n -f ${dros}.longest_AA.txt ${dros}.CDS.fa > ${dros}.longest.CDS.fa"
done


# Clean up the headers.
for dros in Dmel Dtei Dyak Dsan; do
cp ${dros}.longest.CDS.fa ${dros}.longest.simpleheader.CDS.fa
cp ${dros}.longest.AA.fa ${dros}.longest.simpleheader.AA.fa
sed -i "s/>.*gene=/>${dros}_/g" ${dros}.longest.simpleheader.CDS.fa
sed -i "s/>.*gene=/>${dros}_/g" ${dros}.longest.simpleheader.AA.fa
sed -i 's/].*$//g' ${dros}.longest.simpleheader.CDS.fa
sed -i 's/].*$//g' ${dros}.longest.simpleheader.AA.fa
done


# Copy the proteomes to the new directory.
for dros in Dmel Dtei Dyak Dsan; do
cp ${dros}.longest.simpleheader.AA.fa ./dros_proteomes
done


# Run OrthoFinder.
sbatch -J ortho_finder_run -N 1 -n 15 -t 1-0 --mem=300G -o logs/ortho_finder_run-%A.out -e logs/ortho_finder_run-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate py2; python ./OrthoFinder_source/orthofinder.py -t 24 -a 6 -f ./dros_proteomes"


# Move all the single copy ortholog protein fastas into the alignment directory and name them .faa to specify amino acid fasta.
cat ./dros_proteomes/OrthoFinder/Results_Oct31/Orthogroups/Orthogroups_SingleCopyOrthologues.txt | while read Orthogroup ; do cp ./dros_proteomes/OrthoFinder/Results_Oct31/Orthogroup_Sequences/${Orthogroup}.fa ./alignments/${Orthogroup}.faa ; done


# Extract the single copy orthologues from the cds fasta files.
cat ./dros_proteomes/OrthoFinder/Results_Oct31/Orthogroups/Orthogroups_SingleCopyOrthologues.txt | while read Orthogroup ; do cat ./species_prefixes.txt | while read prefix species ; do grep ${species} ./alignments/${Orthogroup}.faa | sed 's/>//g' | seqkit grep -n -f - ${species}.longest.simpleheader.CDS.fa >> ./alignments/${Orthogroup}.fna ; done ; done


# Align all of the amino acids.
while IFS= read -r line; do
sbatch -J ${line}_align -N 1 -n 1 -t 1:00:00 --mem=5G -o logs/${line}_align-%A.out -e logs/${line}_align-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load clustal_omega/1.2.4; clustalo -i ./alignments/${line}.faa -o ./alignments/${line}.aln.faa"
done < ./SingleCopyOrthologues00

while IFS= read -r line; do
sbatch -J ${line}_align -N 1 -n 1 -t 1:00:00 --mem=5G -o logs/${line}_align-%A.out -e logs/${line}_align-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load clustal_omega/1.2.4; clustalo -i ./alignments/${line}.faa -o ./alignments/${line}.aln.faa"
done < ./SingleCopyOrthologues01

while IFS= read -r line; do
sbatch -J ${line}_align -N 1 -n 1 -t 1:00:00 --mem=5G -o logs/${line}_align-%A.out -e logs/${line}_align-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load clustal_omega/1.2.4; clustalo -i ./alignments/${line}.faa -o ./alignments/${line}.aln.faa"
done < ./SingleCopyOrthologues02

while IFS= read -r line; do
sbatch -J ${line}_align -N 1 -n 1 -t 1:00:00 --mem=5G -o logs/${line}_align-%A.out -e logs/${line}_align-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load clustal_omega/1.2.4; clustalo -i ./alignments/${line}.faa -o ./alignments/${line}.aln.faa"
done < ./SingleCopyOrthologues03

while IFS= read -r line; do
sbatch -J ${line}_align -N 1 -n 1 -t 1:00:00 --mem=5G -o logs/${line}_align-%A.out -e logs/${line}_align-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load clustal_omega/1.2.4; clustalo -i ./alignments/${line}.faa -o ./alignments/${line}.aln.faa"
done < ./SingleCopyOrthologues04

while IFS= read -r line; do
sbatch -J ${line}_align -N 1 -n 1 -t 1:00:00 --mem=5G -o logs/${line}_align-%A.out -e logs/${line}_align-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load clustal_omega/1.2.4; clustalo -i ./alignments/${line}.faa -o ./alignments/${line}.aln.faa"
done < ./SingleCopyOrthologues05

while IFS= read -r line; do
sbatch -J ${line}_align -N 1 -n 1 -t 1:00:00 --mem=5G -o logs/${line}_align-%A.out -e logs/${line}_align-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load clustal_omega/1.2.4; clustalo -i ./alignments/${line}.faa -o ./alignments/${line}.aln.faa"
done < ./SingleCopyOrthologues06

while IFS= read -r line; do
sbatch -J ${line}_align -N 1 -n 1 -t 1:00:00 --mem=5G -o logs/${line}_align-%A.out -e logs/${line}_align-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load clustal_omega/1.2.4; clustalo -i ./alignments/${line}.faa -o ./alignments/${line}.aln.faa"
done < ./SingleCopyOrthologues07

while IFS= read -r line; do
sbatch -J ${line}_align -N 1 -n 1 -t 1:00:00 --mem=5G -o logs/${line}_align-%A.out -e logs/${line}_align-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load clustal_omega/1.2.4; clustalo -i ./alignments/${line}.faa -o ./alignments/${line}.aln.faa"
done < ./SingleCopyOrthologues08

while IFS= read -r line; do
sbatch -J ${line}_align -N 1 -n 1 -t 1:00:00 --mem=5G -o logs/${line}_align-%A.out -e logs/${line}_align-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load clustal_omega/1.2.4; clustalo -i ./alignments/${line}.faa -o ./alignments/${line}.aln.faa"
done < ./SingleCopyOrthologues09

while IFS= read -r line; do
sbatch -J ${line}_align -N 1 -n 1 -t 1:00:00 --mem=5G -o logs/${line}_align-%A.out -e logs/${line}_align-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load clustal_omega/1.2.4; clustalo -i ./alignments/${line}.faa -o ./alignments/${line}.aln.faa"
done < ./SingleCopyOrthologues10

while IFS= read -r line; do
sbatch -J ${line}_align -N 1 -n 1 -t 1:00:00 --mem=5G -o logs/${line}_align-%A.out -e logs/${line}_align-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load clustal_omega/1.2.4; clustalo -i ./alignments/${line}.faa -o ./alignments/${line}.aln.faa"
done < ./SingleCopyOrthologues11

sbatch -J OG0002317_align_extra_mem -N 1 -n 1 -t 1:00:00 --mem=25G -o logs/OG0002317_align_extra_mem-%A.out -e logs/OG0002317_align_extra_mem-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load clustal_omega/1.2.4; clustalo -i ./alignments/OG0002317.faa -o ./alignments/OG0002317.aln.faa"
sbatch -J OG0011257_align_extra_mem -N 1 -n 1 -t 1:00:00 --mem=25G -o logs/OG0011257_align_extra_mem-%A.out -e logs/OG0011257_align_extra_mem-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load clustal_omega/1.2.4; clustalo -i ./alignments/OG0011257.faa -o ./alignments/OG0011257.aln.faa"


# Check to see if all the alignments are done.
while IFS= read -r line; do
if [ -f ./alignments/${line}.aln.faa ]
then
	echo ok
else
	echo ${line}
fi
done < ./dros_proteomes/OrthoFinder/Results_Oct31/Orthogroups/Orthogroups_SingleCopyOrthologues.txt


# Create codon alignments for PAML.
while IFS= read -r line; do
./pal2nal.v14/pal2nal.pl ./alignments/${line}.aln.faa ./alignments/${line}.fna -output paml -nogap 1> ./alignments/${line}.pal2nal 2> ./pal2nal_logs/${line}.pal2nal.log
done < ./dros_proteomes/OrthoFinder/Results_Oct31/Orthogroups/Orthogroups_SingleCopyOrthologues.txt

while IFS= read -r line; do
./pal2nal.v14/pal2nal.pl ./alignments/${line}.aln.faa ./alignments/${line}.fna -output paml -nogap 1> ./alignments/${line}.pal2nal 2> ./pal2nal_logs/${line}.pal2nal.log
done < qc_failed_pal2nal


# Make a list of alignments with errors.
cat ./dros_proteomes/OrthoFinder/Results_Oct31/Orthogroups/Orthogroups_SingleCopyOrthologues.txt | while read Orthogroup ; do grep "ERROR" ./pal2nal_log/${Orthogroup}.pal2nal.log | sed 's/.pal2nal.*//g' >> failed_alignments ; done

while IFS= read -r line; do
if grep -q "ERROR" ./pal2nal_logs/${line}.pal2nal.log
then
	echo ${line} >> qc_failed_pal2nal
fi
done < ./dros_proteomes/OrthoFinder/Results_Oct31/Orthogroups/Orthogroups_SingleCopyOrthologues.txt

cat ./qc_failed_pal2nal | while read Orthogroup ; do cat ./species_prefixes.txt | while read prefix species ; do grep ${species} ./alignments/${Orthogroup}.faa | sed 's/>//g' | seqkit grep -n -f - ${species}.longest.simpleheader.CDS.fa >> ./alignments/${Orthogroup}.fna ; done ; done

NOTE: I needed to replace paranthese with _ ie sed -i 's/(/_/g; s/)/_/g' Dmel.longest.simpleheader.CDS.fa


# Run PAML on all orthologues.
while IFS= read -r line; do
prefix=${line}
infile=../alignments/${line}.pal2nal
outfile=${line}.out
awk 'NR==1{print $1,$2,"'''${infile}'''"}NR==2{print $1,$2,"'''${outfile}'''"}NR>2{print $0}' yn00_pre.ctl > yn00.ctl;
../paml-4.10.5/bin/yn00 yn00.ctl;
done < Orthogroups_SingleCopyOrthologues.txt


# Generate trees for CODEML.
while IFS= read -r line; do
species=$(grep '>' ./alignments/${line}.faa | sed 's/>//g')
mel=$(echo ${species} | awk '{ print $1 }')
san=$(echo ${species} | awk '{ print $2 }')
tei=$(echo ${species} | awk '{ print $3 }')
yak=$(echo ${species} | awk '{ print $4 }')
echo "(${mel},(${tei},(${yak},${san})));" > ./alignments/${line}.tree
done < Orthogroups_SingleCopyOrthologues.txt


# Run CODEML on all orthologues.
while IFS= read -r line; do
echo "seqfile = ../alignments/${line}.pal2nal   * sequence data filename
outfile = ../${line}.out   * main result file name
treefile = ../alignments/${line}.tree
noisy = 0      * 0,1,2,3,9: how much rubbish on the screen
verbose = 0      * 1:detailed output
runmode = -2     * -2:pairwise
seqtype = 1      * 1:codons
CodonFreq = 2      * 0:equal, 1:F1X4, 2:F3X4, 3:F61
model = 1      *
NSsites = 0      *
icode = 0      * 0:universal code
fix_kappa = 1      * 1:kappa fixed, 0:kappa to be estimated
kappa = 1      * initial or fixed kappa
fix_omega = 0      * 1:omega fixed, 0:omega to be estimated
omega = 0.5    * initial omega value" > codeml.ctl
../paml-4.10.5/bin/codeml codeml.ctl
done < Orthogroups_SingleCopyOrthologues.txt

# Create a mapping file.
while IFS= read -r line; do
gene=$(grep '>Dmel_' ./alignments/${line}.faa | sed 's/>Dmel_//g')
echo -e "${line}\t${gene}\n ">> ortholog_to_gene.txt
done < Orthogroups_SingleCopyOrthologues.txt


# Rename and move all of the ortholologous alignment files.
cat obp_orthologs.txt | while read ortholog gene ; do cp ./alignments/${ortholog}.aln.faa ./test_alignments/${gene}.aln.faa ; done
cat obp_orthologs.txt | while read ortholog gene ; do cp ./alignments/${ortholog}.fna ./test_alignments/${gene}.fna ; done
cat obp_orthologs.txt | while read ortholog gene ; do cp ./alignments/${ortholog}.pal2nal ./test_alignments/${gene}.pal2nal ; done
cat five_random_orthologs.txt | while read ortholog gene ; do cp ./alignments/${ortholog}.aln.faa ./test_alignments/${gene}.aln.faa ; done
cat five_random_orthologs.txt | while read ortholog gene ; do cp ./alignments/${ortholog}.fna ./test_alignments/${gene}.fna ; done
cat five_random_orthologs.txt | while read ortholog gene ; do cp ./alignments/${ortholog}.pal2nal ./test_alignments/${gene}.pal2nal ; done


# Parse all of the paml outputs.
paml_file = './{0}.out'.format(str(sys.argv[1]))
# Open the paml output.
with open(paml_file, 'r') as data:
    # Loop through every line.
    for line in data:
        # If the line starts with Dmel_.
        if line.startswith('Dmel_'):
            # Print the gene.
            gene = line[5:].strip()
        else:
            # Split the line by tabs.
            spline = line.split()
            # If the line has information.
            if len(spline) > 1:
                # Grab the first two elements.
                idx_1 = spline[0]
                idx_2 = spline[1]
                # If it is Dsan v Dmel.
                if (idx_1 == '2') & (idx_2 == '1'):
                    # Record the omega value.
                    omega_san_mel = spline[6].strip()
                # Else if it is Dyak v Dmel.
                elif (idx_1 == '4') & (idx_2 == '1'):
                    # Record the omega value.
                    omega_yak_mel = spline[6].strip()
                # Else if it is Dsan v Dtei.
                elif (idx_1 == '3') & (idx_2 == '2'):
                    # Record the omega value.
                    omega_san_tei = spline[6].strip()
                # Else if it is Dsan v Dtei.
                elif (idx_1 == '4') & (idx_2 == '3'):
                    # Record the omega value.
                    omega_yak_tei = spline[6].strip()
                # Else...
                else:
                    continue
            # Else...
            else:
                continue
# Construct an output list.
output_list = [gene, omega_san_mel, omega_yak_mel, omega_san_tei, omega_yak_tei]
# Open the output file.
out_file = open('./{0}_omega.txt'.format(gene), 'w')
# Write the output file.
out_file.write('\t'.join(output_list)+'\n')
# Close the output file.
out_file.close()
# Print the output.
print(gene)


# Create condon alignments in fasta format.
while IFS= read -r line; do
./pal2nal.v14/pal2nal.pl ./alignments/${line}.aln.faa ./alignments/${line}.fna -output fasta -nogap 1> ./alignments/${line}.pal2nal.fa 2> ./pal2nal_logs/${line}.pal2nal.fa.log
done < ./dros_proteomes/OrthoFinder/Results_Oct31/Orthogroups/Orthogroups_SingleCopyOrthologues.txt


# Make a list of alignments with errors.
while IFS= read -r line; do
if grep -q "ERROR" ./pal2nal_logs/${line}.pal2nal.fa.log
then
        echo ${line} >> qc_failed_pal2nal_fa
fi
done < ./dros_proteomes/OrthoFinder/Results_Oct31/Orthogroups/Orthogroups_SingleCopyOrthologues.txt

# Rename and alignments for download.
cat cleaned_ortholog_to_gene.txt | while read ortholog gene ; do cp ./alignments/${ortholog}.aln.faa ./cleaned_alns/${gene}.aln.faa ; cp ./alignments/${ortholog}.fna ./cleaned_alns/${gene}.fna ; cp ./alignments/${ortholog}.pal2nal ./cleaned_alns/${gene}.pal2nal; cp ./alignments/${ortholog}.pal2nal.fa ./cleaned_alns/${gene}.pal2nal.fa; done


# Clean all of the files.
cat cleaned_ortholog_to_gene.txt | while read ortholog gene ; do sed -i 's/[*_].*//g' ./cleaned_alns/${gene}.aln.faa ; sed -i 's/[*_].*//g' ./cleaned_alns/${gene}.fna ; sed -i 's/[*_].*//g' ./cleaned_alns/${gene}.pal2nal; sed -i 's/[*_].*//g' ./cleaned_alns/${gene}.pal2nal.fa; done

# Create a mapping file for ortholog coordinates.
while IFS= read -r line; do
mel=$(grep '>Dmel' ./alignments/${line}.faa | sed 's/>Dmel_//g')
san=$(grep '>Dsan' ./alignments/${line}.faa | sed 's/>Dsan_//g')
yak=$(grep '>Dyak' ./alignments/${line}.faa | sed 's/>Dyak_//g')
echo -e "${line}\t${mel}\t${san}\t${yak}">> orthogroup_mel_san_yak_mapping.txt
done < Orthogroups_SingleCopyOrthologues.txt



### POPGEN ###

# Rename fastqs by strain.
cat matute_pe_fq.txt | while read sra strain dros ; do mv ${sra}_1.fastq.gz ${strain}_1.fastq.gz ; done
cat matute_pe_fq.txt | while read sra strain dros ; do mv ${sra}_2.fastq.gz ${strain}_2.fastq.gz ; done
cat matute_se_fq.txt | while read sra strain dros ; do mv ${sra}.fastq.gz ${strain}.fastq.gz ; done
### rename_fastqs.sh ###


# Download the reference genomes.
sbatch -J san_fna -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/san_fna-%A.out -e logs/san_fna-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_santomea/latest_assembly_versions/GCF_016746245.2_Prin_Dsan_1.1/GCF_016746245.2_Prin_Dsan_1.1_genomic.fna.gz"
sbatch -J yak_fna -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/yak_fna-%A.out -e logs/yak_fna-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_yakuba/latest_assembly_versions/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1_genomic.fna.gz"
sbatch -J tei_fna -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/tei_fna-%A.out -e logs/tei_fna-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_teissieri/latest_assembly_versions/GCF_016746235.2_Prin_Dtei_1.1/GCF_016746235.2_Prin_Dtei_1.1_genomic.fna.gz"
sbatch -J mel_fna -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/mel_fna-%A.out -e logs/mel_fna-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_melanogaster/latest_assembly_versions/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz"


# Index the reference genomes for bowtie2.
for dros in Dmel Dtei Dyak Dsan; do
sbatch -J ${dros}_fna_index -N 1 -n 1 -t 12:00:00 --mem=5G -o logs/${dros}_fna_index-%A.out -e logs/${dros}_fna_index-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load bowtie2/2.3.0; bowtie2-build ${dros}.genomic.fna ${dros}.genomic"
done

# Index the reference genomes for samtools.
for dros in Dmel Dtei Dyak Dsan; do
sbatch -J ${dros}_fna_index -N 1 -n 1 -t 12:00:00 --mem=5G -o logs/${dros}_fna_index-%A.out -e logs/${dros}_fna_index-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load samtools; samtools faidx ${dros}.genomic.fna"
done


# Trim all paired end reads.
cat matute_pe_fq.txt | while read sra strain dros; do 
sbatch -J ${strain}_pe_trim -N 1 -n 1 -t 3:00:00 --mem=2G -o logs/${strain}_pe_trim-%A.out -e logs/${strain}_pe_trim-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load trimgalore fastqc cutadapt; trim_galore -q 20 --paired --illumina ${strain}_1.fastq.gz ${strain}_2.fastq.gz"
done


# Trim all single end reads.
cat matute_se_fq.txt | while read sra strain dros; do 
sbatch -J ${strain}_se_trim -N 1 -n 1 -t 3:00:00 --mem=2G -o logs/${strain}_se_trim-%A.out -e logs/${strain}_se_trim-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load trimgalore fastqc cutadapt; trim_galore -q 20 ${strain}.fastq.gz"
done


# Align paired-end reads.
cat matute_pe_fq.txt | while read sra strain dros; do 
sbatch -J ${strain}_yak_align_pe -N 1 -n 2 -t 12:00:00 --mem=2G -o logs/${strain}_yak_align_pe-%A.out -e logs/${strain}_yak_align_pe-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load bowtie2; bowtie2 --no-unal -p 4 --rg-id ${strain} -x ../assemblies/Dyak.genomic -1 ${strain}_1_val_1.fq.gz -2 ${strain}_2_val_2.fq.gz -S ../alignments/yak/${strain}.sam"
sbatch -J ${strain}_san_align_pe -N 1 -n 2 -t 12:00:00 --mem=2G -o logs/${strain}_san_align_pe-%A.out -e logs/${strain}_san_align_pe-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load bowtie2; bowtie2 --no-unal -p 4 --rg-id ${strain} -x ../assemblies/Dsan.genomic -1 ${strain}_1_val_1.fq.gz -2 ${strain}_2_val_2.fq.gz -S ../alignments/san/${strain}.sam"
sbatch -J ${strain}_tei_align_pe -N 1 -n 2 -t 12:00:00 --mem=2G -o logs/${strain}_tei_align_pe-%A.out -e logs/${strain}_tei_align_pe-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load bowtie2; bowtie2 --no-unal -p 4 --rg-id ${strain} -x ../assemblies/Dtei.genomic -1 ${strain}_1_val_1.fq.gz -2 ${strain}_2_val_2.fq.gz -S ../alignments/tei/${strain}.sam"
sbatch -J ${strain}_mel_align_pe -N 1 -n 2 -t 12:00:00 --mem=2G -o logs/${strain}_mel_align_pe-%A.out -e logs/${strain}_mel_align_pe-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load bowtie2; bowtie2 --no-unal -p 4 --rg-id ${strain} -x ../assemblies/Dmel.genomic -1 ${strain}_1_val_1.fq.gz -2 ${strain}_2_val_2.fq.gz -S ../alignments/mel/${strain}.sam"
done


# Align single-end reads.
cat matute_se_fq.txt | while read sra strain dros; do 
sbatch -J ${strain}_yak_align_se -N 1 -n 2 -t 12:00:00 --mem=2G -o logs/${strain}_yak_align_se-%A.out -e logs/${strain}_yak_align_se-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load bowtie2; bowtie2 --no-unal -p 4 --rg-id ${strain} -x ../assemblies/Dyak.genomic -U ${strain}_trimmed.fq.gz -S ../alignments/yak/${strain}.sam"
sbatch -J ${strain}_san_align_se -N 1 -n 2 -t 12:00:00 --mem=2G -o logs/${strain}_san_align_se-%A.out -e logs/${strain}_san_align_se-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load bowtie2; bowtie2 --no-unal -p 4 --rg-id ${strain} -x ../assemblies/Dsan.genomic -U ${strain}_trimmed.fq.gz -S ../alignments/san/${strain}.sam"
sbatch -J ${strain}_tei_align_se -N 1 -n 2 -t 12:00:00 --mem=2G -o logs/${strain}_tei_align_se-%A.out -e logs/${strain}_tei_align_se-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load bowtie2; bowtie2 --no-unal -p 4 --rg-id ${strain} -x ../assemblies/Dtei.genomic -U ${strain}_trimmed.fq.gz -S ../alignments/tei/${strain}.sam"
sbatch -J ${strain}_mel_align_se -N 1 -n 2 -t 12:00:00 --mem=2G -o logs/${strain}_mel_align_se-%A.out -e logs/${strain}_mel_align_se-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load bowtie2; bowtie2 --no-unal -p 4 --rg-id ${strain} -x ../assemblies/Dmel.genomic -U ${strain}_trimmed.fq.gz -S ../alignments/mel/${strain}.sam"
done


# Using 4 threads (-@ 4) filter out unmapped reads (-F 4) and reads with a MAPQ less than 30 (-q 30) then convert to a BAM. Next sort the BAM by read name (-n).
cat matute_samps.txt | while read strain; do for dros in yak san tei mel; do
sbatch -J ${strain}_${dros}_sam_2_bam -N 1 -n 2 -t 3:00:00 --mem=5G --account=ccmb-condo -o logs/${strain}_${dros}_sam_2_bam-%A.out -e logs/${strain}_${dros}_sam_2_bam-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load samtools; samtools view -@ 4 -b -q 30 ./${dros}/${strain}.sam | samtools sort -@ 4 -n -O BAM > ./${dros}/${strain}_read_name_sorted.bam"
done; done


# Add the proper tag to mark duplicates.
cat matute_samps.txt | while read strain; do for dros in yak san tei mel; do
sbatch -J ${strain}_${dros}_fixmate -N 1 -n 2 -t 3:00:00 --mem=5G --account=ccmb-condo -o logs/${strain}_${dros}_fixmate-%A.out -e logs/${strain}_${dros}_fixmate-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load samtools; samtools fixmate -@ 4 -m -O BAM ./${dros}/${strain}_read_name_sorted.bam ./${dros}/${strain}_fixmate.bam"
done; done


# Sort the BAM by coordinates.
cat matute_samps.txt | while read strain; do for dros in yak san tei mel; do
sbatch -J ${strain}_${dros}_sort_coords -N 1 -n 2 -t 3:00:00 --mem=5G --account=ccmb-condo -o logs/${strain}_${dros}_sort_coords-%A.out -e logs/${strain}_${dros}_sort_coords-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load samtools; samtools sort -@ 4 -O BAM -o ./${dros}/${strain}_coord_sorted.bam ./${dros}/${strain}_fixmate.bam"
done; done


# Remove duplicate records.
cat matute_samps.txt | while read strain; do for dros in yak san tei mel; do
sbatch -J ${strain}_${dros}_rmv_dups -N 1 -n 2 -t 3:00:00 --mem=5G --account=ccmb-condo -o logs/${strain}_${dros}_rmv_dups-%A.out -e logs/${strain}_${dros}_rmv_dups-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load samtools; samtools markdup -@ 4 -r -s -O BAM ./${dros}/${strain}_coord_sorted.bam ./${dros}/${strain}_coord_sorted_duprmv.bam"
done; done


# Index the bam file.
cat matute_samps.txt | while read strain; do for dros in yak san tei mel; do
sbatch -J ${strain}_${dros}_bam_idx -N 1 -n 2 -t 3:00:00 --mem=5G --account=ccmb-condo -o logs/${strain}_${dros}_bam_idx-%A.out -e logs/${strain}_${dros}_bam_idx-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load samtools; samtools index -b -@ 4 ./${dros}/${strain}_coord_sorted_duprmv.bam"
sbatch -J ${strain}_${dros}_bam_idx -N 1 -n 2 -t 3:00:00 --mem=5G --account=ccmb-condo -o logs/${strain}_${dros}_bam_idx-%A.out -e logs/${strain}_${dros}_bam_idx-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load samtools; samtools index -b -@ 4 ./${dros}/${strain}_sorted.bam"
done; done


# Subset the four major Dmel chromosomes.
cat matute_samps.txt | while read strain; do for dros in mel; do
sbatch -J ${strain}_${dros}_chrom-subset -N 1 -n 2 -t 3:00:00 --mem=5G --account=ccmb-condo -o logs/${strain}_${dros}_chrom-subset-%A.out -e logs/${strain}_${dros}_chrom-subset-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load samtools; samtools view -@ 4 -b -o ./${dros}/${strain}_coord_sorted_duprmv_chrom_subset.bam ./${dros}/${strain}_coord_sorted_duprmv.bam NC_004354.4 NT_033779.5 NT_033778.4 NT_037436.4 NT_033777.3 NC_004353.4"
sbatch -J ${strain}_${dros}_chrom-subset -N 1 -n 2 -t 3:00:00 --mem=5G --account=ccmb-condo -o logs/${strain}_${dros}_chrom-subset-%A.out -e logs/${strain}_${dros}_chrom-subset-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load samtools; samtools view -@ 4 -b -o ./${dros}/${strain}_sorted_chrom_subset.bam ./${dros}/${strain}_sorted.bam NC_004354.4 NT_033779.5 NT_033778.4 NT_037436.4 NT_033777.3 NC_004353.4"
done; done


# Subset the four major Dtei chromosomes.
cat matute_samps.txt | while read strain; do for dros in tei; do
sbatch -J ${strain}_${dros}_chrom-subset -N 1 -n 2 -t 3:00:00 --mem=5G --account=ccmb-condo -o logs/${strain}_${dros}_chrom-subset-%A.out -e logs/${strain}_${dros}_chrom-subset-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load samtools; samtools view -@ 4 -b -o ./${dros}/${strain}_coord_sorted_duprmv_chrom_subset.bam ./${dros}/${strain}_coord_sorted_duprmv.bam NC_053029.1 NC_053030.1 NC_053031.1 NC_053032.1 NC_053033.1 NC_053034.1"
sbatch -J ${strain}_${dros}_chrom-subset -N 1 -n 2 -t 3:00:00 --mem=5G --account=ccmb-condo -o logs/${strain}_${dros}_chrom-subset-%A.out -e logs/${strain}_${dros}_chrom-subset-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load samtools; samtools view -@ 4 -b -o ./${dros}/${strain}_sorted_chrom_subset.bam ./${dros}/${strain}_sorted.bam NC_053029.1 NC_053030.1 NC_053031.1 NC_053032.1 NC_053033.1 NC_053034.1"
done; done


# Subset the four major Dyak chromosomes.
cat matute_samps.txt | while read strain; do for dros in yak; do
sbatch -J ${strain}_${dros}_chrom-subset -N 1 -n 2 -t 3:00:00 --mem=5G --account=ccmb-condo -o logs/${strain}_${dros}_chrom-subset-%A.out -e logs/${strain}_${dros}_chrom-subset-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load samtools; samtools view -@ 4 -b -o ./${dros}/${strain}_coord_sorted_duprmv_chrom_subset.bam ./${dros}/${strain}_coord_sorted_duprmv.bam NC_052526.2 NC_052527.2 NC_052528.2 NC_052529.2 NC_052530.2 NC_052531.2"
sbatch -J ${strain}_${dros}_chrom-subset -N 1 -n 2 -t 3:00:00 --mem=5G --account=ccmb-condo -o logs/${strain}_${dros}_chrom-subset-%A.out -e logs/${strain}_${dros}_chrom-subset-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load samtools; samtools view -@ 4 -b -o ./${dros}/${strain}_sorted_chrom_subset.bam ./${dros}/${strain}_sorted.bam NC_052526.2 NC_052527.2 NC_052528.2 NC_052529.2 NC_052530.2 NC_052531.2"
done; done


# Subset the four major Dsan chromosomes.
cat matute_samps.txt | while read strain; do for dros in san; do
sbatch -J ${strain}_${dros}_chrom-subset -N 1 -n 2 -t 3:00:00 --mem=5G --account=ccmb-condo -o logs/${strain}_${dros}_chrom-subset-%A.out -e logs/${strain}_${dros}_chrom-subset-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load samtools; samtools view -@ 4 -b -o ./${dros}/${strain}_coord_sorted_duprmv_chrom_subset.bam ./${dros}/${strain}_coord_sorted_duprmv.bam NC_053016.2 NC_053017.2 NC_053018.2 NC_053019.2 NC_053020.2 NC_053021.2"
sbatch -J ${strain}_${dros}_chrom-subset -N 1 -n 2 -t 3:00:00 --mem=5G --account=ccmb-condo -o logs/${strain}_${dros}_chrom-subset-%A.out -e logs/${strain}_${dros}_chrom-subset-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load samtools; samtools view -@ 4 -b -o ./${dros}/${strain}_sorted_chrom_subset.bam ./${dros}/${strain}_sorted.bam NC_053016.2 NC_053017.2 NC_053018.2 NC_053019.2 NC_053020.2 NC_053021.2"
done; done


# Calculate genome-wide coverage.
cat matute_samps.txt | while read strain; do for dros in yak san tei mel; do
sbatch -J ${strain}_${dros}_cov_calc -N 1 -n 1 -t 12:00:00 --mem=2G --account=ccmb-condo -o logs/${strain}_${dros}_cov_calc-%A.out -e logs/${strain}_${dros}_cov_calc-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load samtools; samtools depth -a ./${dros}/${strain}_coord_sorted_duprmv_chrom_subset.bam | awk '{sum+=\$3; sumsq+=\$3*\$3} END { print \"Average = \",sum/NR; print \"Stdev = \",sqrt(sumsq/NR - (sum/NR)**2)}' > ./coverage/${strain}_${dros}_coord_sorted_duprmv_chrom_subset.txt"
sbatch -J ${strain}_${dros}_cov_calc -N 1 -n 1 -t 12:00:00 --mem=2G --account=ccmb-condo -o logs/${strain}_${dros}_cov_calc-%A.out -e logs/${strain}_${dros}_cov_calc-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load samtools; samtools depth -a ./${dros}/${strain}_sorted_chrom_subset.bam | awk '{sum+=\$3; sumsq+=\$3*\$3} END { print \"Average = \",sum/NR; print \"Stdev = \",sqrt(sumsq/NR - (sum/NR)**2)}' > ./coverage/${strain}_${dros}_sorted_chrom_subset.txt"
done; done


# Call variants using bcftools.
cat matute_samps.txt | while read strain; do for dros in yak san tei mel; do
sbatch -J ${strain}_${dros}_bam_2_vcf -N 1 -n 2 -t 3-0 --mem=4G -o logs/${strain}_${dros}_bam_2_vcf-%A.out -e logs/${strain}_${dros}_bam_2_vcf-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load bcftools/1.13 gsl/2.5 gcc/8.3 perl/5.24.1; bcftools mpileup --threads 4 -f ../assemblies/D${dros}.genomic.fna ./${dros}/${strain}_sorted_chrom_subset.bam | bcftools call -m -f GQ -Oz -o ../vcfs/${dros}/${strain}_unfiltered.vcf.gz"
done; done


# Determine the read depth cut-offs for every sample.
cat matute_samps.txt | while read strain; do for dros in yak san; do
sbatch -J ${strain}_${dros}_dp_qc -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/${strain}_${dros}_dp_qc-%A.out -e logs/${strain}_${dros}_dp_qc-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python qc_dp.py ${strain} ${dros} > ./${dros}/${strain}_dp_cutoffs.txt"
done; done


# Filter each sample for read depth, quality, and mapping quality.
cat san_dp_cutoffs.txt | while read strain dp_min dp_max; do
sbatch -J ${strain}_san_dp_filter -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/${strain}_san_dp_filter-%A.out -e logs/${strain}_san_dp_filter-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load bcftools/1.13 gsl/2.5 gcc/8.3 perl/5.24.1; bcftools view -i 'DP>${dp_min} && DP<${dp_max} && QUAL>20 && MQ>30' -Oz -o ./san/${strain}_dp_mq_filtered.vcf.gz ./san/${strain}_unfiltered.vcf.gz"
done


# Filter each sample for read depth, quality, and mapping quality.
cat yak_dp_cutoffs.txt | while read strain dp_min dp_max; do
sbatch -J ${strain}_yak_dp_filter -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/${strain}_yak_dp_filter-%A.out -e logs/${strain}_yak_dp_filter-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load bcftools/1.13 gsl/2.5 gcc/8.3 perl/5.24.1; bcftools view -i 'DP>${dp_min} && DP<${dp_max} && QUAL>20 && MQ>30' -Oz -o ./yak/${strain}_dp_mq_filtered.vcf.gz ./yak/${strain}_unfiltered.vcf.gz"
done


# Filter out INDELs and sites that do not have a GQ of 40 or higher.
cat matute_samps.txt | while read strain; do for dros in yak san; do
sbatch -J ${strain}_${dros}_vcf_filter -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/${strain}_${dros}_vcf_filter-%A.out -e logs/${strain}_${dros}_vcf_filter-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load tabix python/3.7.4; cd ./${dros}; python3 ../vcf_filter.py ${strain} ${dros} | bgzip > ${strain}_final_filtered.vcf.gz"
done; done


# Index the filtered vcf files.
cat matute_samps.txt | while read strain; do for dros in yak san; do
sbatch -J ${strain}_${dros}_vcf_index -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/${strain}_${dros}_vcf_index-%A.out -e logs/${strain}_${dros}_vcf_index-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load tabix; tabix -p vcf ./${dros}/${strain}_final_filtered.vcf.gz"
done; done


# Generate Dsan sample lists.
cat san_samps.txt | while read strain meta; do for dros in yak san; do
printf "${strain}_final_filtered.vcf.gz\n" >> ./${dros}/san_vcf_list.txt
done; done

# Generate Dyak sympatric sample lists.
cat yak_symp_samps.txt | while read strain meta; do for dros in yak san; do
printf "${strain}_final_filtered.vcf.gz\n" >> ./${dros}/yak_symp_vcf_list.txt
done; done

# Generate Dyak allopatric sample lists.
cat yak_allo_samps.txt | while read strain meta; do for dros in yak san; do
printf "${strain}_final_filtered.vcf.gz\n" >> ./${dros}/yak_allo_vcf_list.txt
done; done

# Generate Dtei sample lists.
cat tei_samps.txt | while read strain meta; do for dros in yak san; do
printf "${strain}_final_filtered.vcf.gz\n" >> ./${dros}/tei_vcf_list.txt
done; done


# Merge the individually filtered vcf files.
for pop in san yak_symp yak_allo tei; do for dros in san yak; do
sbatch -J ${pop}_${dros}_vcf_merge -N 1 -n 1 -t 3:00:00 --mem=2G --account=ccmb-condo -o logs/${pop}_${dros}_vcf_merge-%A.out -e logs/${pop}_${dros}_vcf_merge-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load bcftools/1.13 gsl/2.5 gcc/8.3 perl/5.24.1; cd ./${dros}; bcftools merge -m all -l ${pop}_vcf_list.txt -Oz -o ${pop}_unfiltered.vcf.gz"
done; done


# Filter out sites that don't have 90% genotype information.
for pop in san yak_symp yak_allo tei; do for dros in san yak; do
sbatch -J ${pop}_${dros}_vcf_filter -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/${pop}_${dros}_vcf_filter-%A.out -e logs/${pop}_${dros}_vcf_filter-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load bcftools/1.13 gsl/2.5 gcc/8.3 perl/5.24.1; bcftools view -M2 -e 'F_MISSING>0.1' -Oz -o ./${dros}/${pop}_filtered.vcf.gz ./${dros}/${pop}_unfiltered.vcf.gz"
done; done


# Index the filtered merged vcf files.
for pop in san yak_symp yak_allo tei; do for dros in san yak; do
sbatch -J ${pop}_${dros}_index -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/${pop}_${dros}_index-%A.out -e logs/${pop}_${dros}_index-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load tabix; tabix -p vcf ./${dros}/${pop}_filtered.vcf.gz"
done; done


# Merge the population filtered vcf files.
for dros in san yak; do
sbatch -J ${dros}_vcf_merge -N 1 -n 1 -t 12:00:00 --mem=2G --account=ccmb-condo -o logs/${dros}_vcf_merge-%A.out -e logs/${dros}_vcf_merge-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load bcftools/1.13 gsl/2.5 gcc/8.3 perl/5.24.1; cd ./${dros}; bcftools merge -m all -Oz -o san_yak_tei_unfiltered.vcf.gz san_filtered.vcf.gz yak_symp_filtered.vcf.gz yak_allo_filtered.vcf.gz tei_filtered.vcf.gz"
done


# Index the merged vcf files.
for dros in san yak; do
sbatch -J ${dros}_merge_index -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/${dros}_merge_index-%A.out -e logs/${dros}_merge_index-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load tabix; tabix -p vcf ./${dros}/san_yak_tei_unfiltered.vcf.gz"
done


# Filter out sites that don't have 90% genotype information.
for dros in san yak; do
sbatch -J ${dros}_vcf_filter -N 1 -n 1 -t 3:00:00 --mem=2G --account=ccmb-condo -o logs/${dros}_vcf_filter-%A.out -e logs/${dros}_vcf_filter-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load bcftools/1.13 gsl/2.5 gcc/8.3 perl/5.24.1; cd ./${dros}; bcftools view -M2 -e 'F_MISSING>0.1' -Oz -o san_yak_tei_prefiltered.vcf.gz san_yak_tei_unfiltered.vcf.gz"
done


# Index the prefiltered merged vcf files.
for dros in san yak; do
sbatch -J ${dros}_prefilter_index -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/${dros}_prefilter_index-%A.out -e logs/${dros}_prefilter_index-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load tabix; tabix -p vcf ./${dros}/san_yak_tei_prefiltered.vcf.gz"
done


# Convert the prefiltered data to zarr arrays.
for dros in san yak; do
sbatch -J ${dros}_prefiltered_zarr -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/${dros}_prefiltered_zarr-%A.out -e logs/${dros}_prefiltered_zarr-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python prefiltered_vcf_to_zarr.py ${dros}"
done


# Filter out sites that don't have 90% genotype information across all populations by chromosomes.
for dros in san yak; do for chrom in 2L 2R 3L 3R 4 X; do
sbatch -J ${dros}_${chrom}_subset -N 1 -n 1 -t 1-0 --mem=2G --account=ccmb-condo -o logs/${dros}_${chrom}_subset-%A.out -e logs/${dros}_${chrom}_subset-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load bcftools/1.13 gsl/2.5 gcc/8.3 perl/5.24.1; bcftools view -R ./${dros}_qc_passed_sites_chr${chrom}.txt -Oz -o ./${dros}/san_yak_tei_filtered_chr${chrom}.vcf.gz ./${dros}/san_yak_tei_prefiltered.vcf.gz"
done; done


# Index the final vcf files.
for dros in san yak; do for chrom in 2L 2R 3L 3R 4 X; do
sbatch -J ${dros}_${chrom}_index -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/${dros}_${chrom}_index-%A.out -e logs/${dros}_${chrom}_index-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load tabix; tabix -p vcf ./${dros}/san_yak_tei_filtered_chr${chrom}.vcf.gz"
done; done


# Convert the filtered data to zarr arrays.
for dros in san yak; do for chrom in 2L 2R 3L 3R 4 X; do
sbatch -J ${dros}_${chrom}_filtered_zarr -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/${dros}_${chrom}_filtered_zarr-%A.out -e logs/${dros}_${chrom}_filtered_zarr-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python filtered_vcf_to_zarr.py ${dros} ${chrom}"
done; done










### MISC META DATA ###
>NC_004354.4 Drosophila melanogaster chromosome X
>NT_033779.5 Drosophila melanogaster chromosome 2L
>NT_033778.4 Drosophila melanogaster chromosome 2R
>NT_037436.4 Drosophila melanogaster chromosome 3L
>NT_033777.3 Drosophila melanogaster chromosome 3R
>NC_004353.4 Drosophila melanogaster chromosome 4
>NC_024512.1 Drosophila melanogaster chromosome Y


>NC_053029.1 Drosophila teissieri strain GT53w chromosome 2L, Prin_Dtei_1.1, whole genome shotgun sequence
>NC_053030.1 Drosophila teissieri strain GT53w chromosome 2R, Prin_Dtei_1.1, whole genome shotgun sequence
>NC_053031.1 Drosophila teissieri strain GT53w chromosome 3L, Prin_Dtei_1.1, whole genome shotgun sequence
>NC_053032.1 Drosophila teissieri strain GT53w chromosome 3R, Prin_Dtei_1.1, whole genome shotgun sequence
>NC_053033.1 Drosophila teissieri strain GT53w chromosome 4, Prin_Dtei_1.1, whole genome shotgun sequence
>NC_053034.1 Drosophila teissieri strain GT53w chromosome X, Prin_Dtei_1.1, whole genome shotgun sequence


>NC_052526.2 Drosophila yakuba strain Tai18E2 chromosome X, Prin_Dyak_Tai18E2_2.1, whole genome shotgun sequence
>NC_052527.2 Drosophila yakuba strain Tai18E2 chromosome 2L, Prin_Dyak_Tai18E2_2.1, whole genome shotgun sequence
>NC_052528.2 Drosophila yakuba strain Tai18E2 chromosome 2R, Prin_Dyak_Tai18E2_2.1, whole genome shotgun sequence
>NC_052529.2 Drosophila yakuba strain Tai18E2 chromosome 3L, Prin_Dyak_Tai18E2_2.1, whole genome shotgun sequence
>NC_052530.2 Drosophila yakuba strain Tai18E2 chromosome 3R, Prin_Dyak_Tai18E2_2.1, whole genome shotgun sequence
>NC_052531.2 Drosophila yakuba strain Tai18E2 chromosome 4, Prin_Dyak_Tai18E2_2.1, whole genome shotgun sequenc


>NC_053016.2 Drosophila santomea strain STO CAGO 1482 chromosome 2L, Prin_Dsan_1.1, whole genome shotgun sequence
>NC_053017.2 Drosophila santomea strain STO CAGO 1482 chromosome 2R, Prin_Dsan_1.1, whole genome shotgun sequence
>NC_053018.2 Drosophila santomea strain STO CAGO 1482 chromosome 3L, Prin_Dsan_1.1, whole genome shotgun sequence
>NC_053019.2 Drosophila santomea strain STO CAGO 1482 chromosome 3R, Prin_Dsan_1.1, whole genome shotgun sequence
>NC_053020.2 Drosophila santomea strain STO CAGO 1482 chromosome 4, Prin_Dsan_1.1, whole genome shotgun sequence
>NC_053021.2 Drosophila santomea strain STO CAGO 1482 chromosome X, Prin_Dsan_1.1, whole genome shotgun sequence





BS14	san
C550_39	san
C650_14	san
CAR1600	san
Qiuja630.39	san
Quija37	san
Rain42	san
sanC1350.14	san
sanCAR1490.5	san
sanCOST1250.5	san
sanCOST1270.6	san
sanOBAT1200.13	san
sanOBAT1200.5	san
sanRain39	san
sanSTO7	san
sanThena5	san
san_Field3	san
1_19	yak_symp
1_5	yak_symp
1_6	yak_symp
1_7	yak_symp
2_11	yak_symp
2_14	yak_symp
2_6	yak_symp
2_8	yak_symp
3_16	yak_symp
3_2	yak_symp
3_23	yak_symp
4_21	yak_symp
BAR_1000_2	yak_symp
Bosu_1235_14	yak_symp
Cascade_SN6_1	yak_symp
COST_1235_2	yak_symp
COST_1235_3	yak_symp
Montecafe_17_17	yak_symp
OBAT_1200_5	yak_symp
SA_3	yak_symp
SN7	yak_symp
SN_Cascade_22	yak_symp
BIOKO_NE_4_6	yak_allo
Cascade_19_16	yak_allo
Cascade_21	yak_allo
Anton_1_Principe	yak_allo
Anton_2_Principe	yak_allo
Abidjan_12	yak_allo
Tai_18	yak_allo
Airport_16_5	yak_allo
Cascade_18	yak_allo
SanTome_city_14_26	yak_allo
SJ14	yak_allo
SJ4	yak_allo
SJ7	yak_allo
SJ_1	yak_allo
Balancha_1	tei
cascade_2_1	tei
cascade_2_2	tei
cascade_2_4	tei
cascade_4_1	tei
cascade_4_2	tei
cascade_4_3	tei
House_Bioko	tei
Bata2	tei
Bata8	tei
La_Lope_Gabon	tei
Selinda	tei
Zimbabwe	tei


BS14	san_st
C550_39	san_st
C650_14	san_st
CAR1600	san_st
Qiuja630.39	san_st
Quija37	san_st
Rain42	san_st
sanC1350.14	san_st
sanCAR1490.5	san_st
sanCOST1250.5	san_st
sanCOST1270.6	san_st
sanOBAT1200.13	san_st
sanOBAT1200.5	san_st
sanRain39	san_st
sanSTO7	san_st
sanThena5	san_st
san_Field3	san_st

Balancha_1	tei_bi
cascade_2_1	tei_bi
cascade_2_2	tei_bi
cascade_2_4	tei_bi
cascade_4_1	tei_bi
cascade_4_2	tei_bi
cascade_4_3	tei_bi
House_Bioko	tei_bi
Bata2	tei_eg
Bata8	tei_eg
La_Lope_Gabon	tei_ga
Selinda	tei_zi
Zimbabwe	tei_zi

BIOKO_NE_4_6	yak_bi
Cascade_19_16	yak_bi
Cascade_21	yak_bi
Anton_1_Principe	yak_pr
Anton_2_Principe	yak_pr
Abidjan_12	yak_ic
Tai_18	yak_ic
Airport_16_5	yak_ll
Cascade_18	yak_ll
SanTome_city_14_26	yak_ll
SJ14	yak_ll
SJ4	yak_ll
SJ7	yak_ll
SJ_1	yak_ll

1_19	yak_hz
1_5	yak_hz
1_6	yak_hz
1_7	yak_hz
2_11	yak_hz
2_14	yak_hz
2_6	yak_hz
2_8	yak_hz
3_16	yak_hz
3_2	yak_hz
3_23	yak_hz
4_21	yak_hz
BAR_1000_2	yak_hz
Bosu_1235_14	yak_hz
Cascade_SN6_1	yak_hz
COST_1235_2	yak_hz
COST_1235_3	yak_hz
Montecafe_17_17	yak_hz
OBAT_1200_5	yak_hz
SA_3	yak_hz
SN7	yak_hz
SN_Cascade_22	yak_hz




























SRR5860641	BS14	san
SRR5860642	C550_39	san
SRR5860635	C650_14	san
SRR5860637	CAR1600	san
SRR5860636	Rain42	san
SRR5860638	san_Field3	san
SRR5860615	Balancha_1	tei
SRR5860623	cascade_2_1	tei
SRR5860616	cascade_2_2	tei
SRR5860622	cascade_2_4	tei
SRR5860618	cascade_4_1	tei
SRR5860572	cascade_4_2	tei
SRR5860617	cascade_4_3	tei
SRR5860621	House_Bioko	tei
SRR5860571	La_Lope_Gabon	tei
SRR5860620	Selinda	tei
SRR5860619	Zimbabwe	tei
SRR5860605	Qiuja630.39	san
SRR5860610	Quija37	san
SRR5860606	sanC1350.14	san
SRR5860607	sanCAR1490.5	san
SRR5860624	sanCOST1250.5	san
SRR5860609	sanCOST1270.6	san
SRR5860612	sanOBAT1200.13	san
SRR5860614	sanOBAT1200.5	san
SRR5860613	sanRain39	san
SRR5860611	sanSTO7	san
SRR5860608	sanThena5	san
SRR5860576	Bata2	tei
SRR5860577	Bata8	tei
SRR5860595	BIOKO_NE_4_6	yak
SRR5860604	Cascade_19_16	yak
SRR5860578	Cascade_21	yak
SRR5860601	1_19	yak
SRR5860600	1_5	yak
SRR5860598	1_6	yak
SRR5860649	1_7	yak
SRR5860596	2_11	yak
SRR5860599	2_14	yak
SRR5860593	2_6	yak
SRR5860603	2_8	yak
SRR5860602	3_16	yak
SRR5860590	3_2	yak
SRR5860655	3_23	yak
SRR5860647	4_21	yak
SRR5860579	BAR_1000_2	yak
SRR5860646	Bosu_1235_14	yak
SRR5860588	Cascade_SN6_1	yak
SRR5860650	COST_1235_2	yak
SRR5860597	COST_1235_3	yak
SRR5860587	Montecafe_17_17	yak
SRR5860589	OBAT_1200_5	yak
SRR5860652	SA_3	yak
SRR5860580	SN7	yak
SRR5860581	SN_Cascade_22	yak
SRR5860591	Airport_16_5	yak
SRR5860592	Cascade_18	yak
SRR5860574	SanTome_city_14_26	yak
SRR5860653	SJ14	yak
SRR5860654	SJ4	yak
SRR5860594	SJ7	yak
SRR5860575	SJ_1	yak
SRR5860651	Anton_1_Principe	yak
SRR5860585	Anton_2_Principe	yak
SRR5860586	Abidjan_12	yak
SRR5860648	Tai_18	yak



OG0001737	Obp8a
OG0002724	Obp19a
OG0002725	Obp19b
OG0002726	Obp19c
OG0002727	Obp19d
OG0003026	Obp22a
OG0003634	Obp28a
OG0005211	Obp44a
OG0005484	Obp46a
OG0005546	Obp47a
OG0005589	Obp47b
OG0005784	Obp49a
OG0005972	Obp50a
OG0005973	Obp50b
OG0005974	Obp50c
OG0005975	Obp50d
OG0005977	Obp50e
OG0006048	Obp51a
OG0006619	Obp56a
OG0006620	Obp56b
OG0006621	Obp56c
OG0006622	Obp56d
OG0006624	Obp56i
OG0006667	Obp57c
OG0006668	Obp57b
OG0006669	Obp57a
OG0006679	Obp57e
OG0006680	Obp57d
OG0006981	Obp58b
OG0006982	Obp58c
OG0006983	Obp58d
OG0008680	Obp69a
OG0009180	Obp73a
OG0010005	Obp83a
OG0010006	Obp83b
OG0010012	Obp83cd
OG0010013	Obp83ef
OG0010014	Obp83g
OG0010125	Obp84a
OG0010242	Obp85a
OG0011570	Obp93a
OG0012468	Obp99a
OG0012474	Obp99c
OG0012476	Obp99d
OG0012477	Obp99b



OG0002689       zld
OG0002690       pico
OG0002691       Nup205
OG0002692       THADA
OG0002693       Hers



#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/001/SRR5860571/SRR5860571_1.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/001/SRR5860571/SRR5860571_2.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/002/SRR5860572/SRR5860572_1.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/002/SRR5860572/SRR5860572_2.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/004/SRR5860574/SRR5860574.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/005/SRR5860575/SRR5860575.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/006/SRR5860576/SRR5860576.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/007/SRR5860577/SRR5860577.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/008/SRR5860578/SRR5860578.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/009/SRR5860579/SRR5860579.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/000/SRR5860580/SRR5860580.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/001/SRR5860581/SRR5860581.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/005/SRR5860585/SRR5860585.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/006/SRR5860586/SRR5860586.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/007/SRR5860587/SRR5860587.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/008/SRR5860588/SRR5860588.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/009/SRR5860589/SRR5860589.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/000/SRR5860590/SRR5860590.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/001/SRR5860591/SRR5860591.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/002/SRR5860592/SRR5860592.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/003/SRR5860593/SRR5860593.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/004/SRR5860594/SRR5860594.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/005/SRR5860595/SRR5860595.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/006/SRR5860596/SRR5860596.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/007/SRR5860597/SRR5860597.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/008/SRR5860598/SRR5860598.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/009/SRR5860599/SRR5860599.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/000/SRR5860600/SRR5860600.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/001/SRR5860601/SRR5860601.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/002/SRR5860602/SRR5860602.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/003/SRR5860603/SRR5860603.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/004/SRR5860604/SRR5860604.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/005/SRR5860605/SRR5860605.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/006/SRR5860606/SRR5860606.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/007/SRR5860607/SRR5860607.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/008/SRR5860608/SRR5860608.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/009/SRR5860609/SRR5860609.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/000/SRR5860610/SRR5860610.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/001/SRR5860611/SRR5860611.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/002/SRR5860612/SRR5860612.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/003/SRR5860613/SRR5860613.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/004/SRR5860614/SRR5860614.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/005/SRR5860615/SRR5860615_1.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/005/SRR5860615/SRR5860615_2.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/006/SRR5860616/SRR5860616_1.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/006/SRR5860616/SRR5860616_2.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/007/SRR5860617/SRR5860617_1.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/007/SRR5860617/SRR5860617_2.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/008/SRR5860618/SRR5860618_1.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/008/SRR5860618/SRR5860618_2.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/009/SRR5860619/SRR5860619_1.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/009/SRR5860619/SRR5860619_2.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/000/SRR5860620/SRR5860620_1.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/000/SRR5860620/SRR5860620_2.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/001/SRR5860621/SRR5860621_1.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/001/SRR5860621/SRR5860621_2.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/002/SRR5860622/SRR5860622_1.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/002/SRR5860622/SRR5860622_2.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/003/SRR5860623/SRR5860623_1.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/003/SRR5860623/SRR5860623_2.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/004/SRR5860624/SRR5860624.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/005/SRR5860635/SRR5860635_1.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/005/SRR5860635/SRR5860635_2.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/006/SRR5860636/SRR5860636_1.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/006/SRR5860636/SRR5860636_2.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/007/SRR5860637/SRR5860637_1.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/007/SRR5860637/SRR5860637_2.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/008/SRR5860638/SRR5860638_1.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/008/SRR5860638/SRR5860638_2.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/001/SRR5860641/SRR5860641_1.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/001/SRR5860641/SRR5860641_2.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/002/SRR5860642/SRR5860642_1.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/002/SRR5860642/SRR5860642_2.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/006/SRR5860646/SRR5860646.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/007/SRR5860647/SRR5860647.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/008/SRR5860648/SRR5860648.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/009/SRR5860649/SRR5860649.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/000/SRR5860650/SRR5860650.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/001/SRR5860651/SRR5860651.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/002/SRR5860652/SRR5860652.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/003/SRR5860653/SRR5860653.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/004/SRR5860654/SRR5860654.fastq.gz
#ftp.sra.ebi.ac.uk/vol1/fastq/SRR586/005/SRR5860655/SRR5860655.fastq.gz



SRR5860641	BS14	san
SRR5860642	C550_39	san
SRR5860635	C650_14	san
SRR5860637	CAR1600	san
SRR5860636	Rain42	san
SRR5860638	san_Field3	san
SRR5860615	Balancha_1	tei
SRR5860623	cascade_2_1	tei
SRR5860616	cascade_2_2	tei
SRR5860622	cascade_2_4	tei
SRR5860618	cascade_4_1	tei
SRR5860572	cascade_4_2	tei
SRR5860617	cascade_4_3	tei
SRR5860621	House_Bioko	tei
SRR5860571	La_Lope_Gabon	tei
SRR5860620	Selinda	tei
SRR5860619	Zimbabwe	tei



SRR5860605	Qiuja630.39	san
SRR5860610	Quija37	san
SRR5860606	sanC1350.14	san
SRR5860607	sanCAR1490.5	san
SRR5860624	sanCOST1250.5	san
SRR5860609	sanCOST1270.6	san
SRR5860612	sanOBAT1200.13	san
SRR5860614	sanOBAT1200.5	san
SRR5860613	sanRain39	san
SRR5860611	sanSTO7	san
SRR5860608	sanThena5	san
SRR5860576	Bata2	tei
SRR5860577	Bata8	tei
SRR5860595	BIOKO_NE_4_6	yak
SRR5860604	Cascade_19_16	yak
SRR5860578	Cascade_21	yak
SRR5860601	1_19	yak
SRR5860600	1_5	yak
SRR5860598	1_6	yak
SRR5860649	1_7	yak
SRR5860596	2_11	yak
SRR5860599	2_14	yak
SRR5860593	2_6	yak
SRR5860603	2_8	yak
SRR5860602	3_16	yak
SRR5860590	3_2	yak
SRR5860655	3_23	yak
SRR5860647	4_21	yak
SRR5860579	BAR_1000_2	yak
SRR5860646	Bosu_1235_14	yak
SRR5860588	Cascade_SN6_1	yak
SRR5860650	COST_1235_2	yak
SRR5860597	COST_1235_3	yak
SRR5860587	Montecafe_17_17	yak
SRR5860589	OBAT_1200_5	yak
SRR5860652	SA_3	yak
SRR5860580	SN7	yak
SRR5860581	SN_Cascade_22	yak
SRR5860591	Airport_16_5	yak
SRR5860592	Cascade_18	yak
SRR5860574	SanTome_city_14_26	yak
SRR5860653	SJ14	yak
SRR5860654	SJ4	yak
SRR5860594	SJ7	yak
SRR5860575	SJ_1	yak
SRR5860651	Anton_1_Principe	yak
SRR5860585	Anton_2_Principe	yak
SRR5860586	Abidjan_12	yak
SRR5860648	Tai_18	yak