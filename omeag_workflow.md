# Workflow



## 1. Download Data

All data is downloaded from NCBI RefSeq.

```bash
# Download .faa, .fna, and gff files.
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
```



## 2. Filter to only keep the longest isoform of each protein.

```bash
# Make a list of protein lengths.
for dros in Dmel Dtei Dyak Dsan; do
seqkit fx2tab -n -l ${dros}.AA.fa | awk 'BEGIN{FS="\t"}{gsub("^.*protein_id=", "", $1); gsub("].*$", "", $1); print}' | sort > ${dros}.AA_len.tsv
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
```



## 3. Rename the fasta headers.

```bash
# Clean up the headers.
for dros in Dmel Dtei Dyak Dsan; do
cp ${dros}.longest.CDS.fa ${dros}.longest.simpleheader.CDS.fa
cp ${dros}.longest.AA.fa ${dros}.longest.simpleheader.AA.fa
sed -i "s/>.*gene=/>${dros}_/g" ${dros}.longest.simpleheader.CDS.fa
sed -i "s/>.*gene=/>${dros}_/g" ${dros}.longest.simpleheader.AA.fa
sed -i 's/].*$//g' ${dros}.longest.simpleheader.CDS.fa
sed -i 's/].*$//g' ${dros}.longest.simpleheader.AA.fa
done
```



## 4. Find 1:1 orthologs.

```bash
# Run OrthoFinder.
sbatch -J ortho_finder_run -N 1 -n 15 -t 1-0 --mem=300G -o logs/ortho_finder_run-%A.out -e logs/ortho_finder_run-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate py2; python ./OrthoFinder_source/orthofinder.py -t 24 -a 6 -f ./dros_proteomes"
```



## 5. Align 1:1 orthologs.

```bash
# Extract the single copy orthologues from the cds fasta files.
cat ./dros_proteomes/OrthoFinder/Results_Oct31/Orthogroups/Orthogroups_SingleCopyOrthologues.txt | while read Orthogroup ; do cat ./species_prefixes.txt | while read prefix species ; do grep ${species} ./alignments/${Orthogroup}.faa | sed 's/>//g' | seqkit grep -n -f - ${species}.longest.simpleheader.CDS.fa >> ./alignments/${Orthogroup}.fna ; done ; done
# Align the orthologs.
while IFS= read -r line; do
sbatch -J ${line}_align -N 1 -n 1 -t 1:00:00 --mem=5G -o logs/${line}_align-%A.out -e logs/${line}_align-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load clustal_omega/1.2.4; clustalo -i ./alignments/${line}.faa -o ./alignments/${line}.aln.faa"
done < ./dros_proteomes/OrthoFinder/Results_Oct31/Orthogroups/Orthogroups_SingleCopyOrthologues.txt 
```



## 6. Create codon-based alignments.

```bash
# Create codon alignments for PAML.
while IFS= read -r line; do
./pal2nal.v14/pal2nal.pl ./alignments/${line}.aln.faa ./alignments/${line}.fna -output paml -nogap 1> ./alignments/${line}.pal2nal 2> ./pal2nal_logs/${line}.pal2nal.log
done < ./dros_proteomes/OrthoFinder/Results_Oct31/Orthogroups/Orthogroups_SingleCopyOrthologues.txt
```



## 7. Run yn00.

```bash
### yn00_pre.ctl ###
seqfile = anyfile  * sequence data file name
outfile = yn * main result file
verbose = 0 * 1: detailed output (list sequences), 0: concise output
icode = 0 * 0:universal code; 1:mammalian mt; 2-10:see below
weighting = 0 * weighting pathways between codons (0/1)?
commonf3x4 = 0 * use one set of codon freqs for all pairs (0/1)? 
```

```bash
# Run PAML on all orthologues.
while IFS= read -r line; do
prefix=${line}
infile=../alignments/${line}.pal2nal
outfile=${line}.out
awk 'NR==1{print $1,$2,"'''${infile}'''"}NR==2{print $1,$2,"'''${outfile}'''"}NR>2{print $0}' yn00_pre.ctl > yn00.ctl;
../paml-4.10.5/bin/yn00 yn00.ctl;
done < Orthogroups_SingleCopyOrthologues.txt
```



## 8. Run codeml.

```bash
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
```



