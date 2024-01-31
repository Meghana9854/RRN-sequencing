# Using QIIME2 to taxonomically classify the OTU generated  # Performed only on PacBio data

module load qiime2/2023.2
source activate qiime2-2023.2 

## 11.1. Importing the OTU biom file from step 10 (from 1.vsearch.sh) into QIIME2 ##

#(i) for option i: GTDB-based chimeras removed #
mkdir ./9.tax_assign/qiime_input/GTDB
for i in ./8.otu_binning/biom_files_formatted/GTDB/*.biomv210.biom; 
do
  qiime tools import \
     --input-path "$i" \
     --type 'FeatureTable[Frequency]' \
     --input-format BIOMV210Format \
     --output-path ./9.tax_assign/qiime_input/GTDB/"$(basename "$i" .biomv210.biom)".qza
done    

#(ii) for option ii: RefSeq-based chimeras removed #
mkdir ./9.tax_assign/qiime_input/RefSeq
for i in ./8.otu_binning/biom_files_formatted/RefSeq/*.biomv210.biom; 
do
   qiime tools import \
     --input-path "$i" \
     --type 'FeatureTable[Frequency]' \
     --input-format BIOMV210Format \
     --output-path ./9.tax_assign/qiime_input/RefSeq/"$(basename "$i" .biomv210.biom)".qza
done   


## 11.2. Importing the fasta files from step 9 into QIIME2 ## 
mkdir ./8.otu_binning/fasta_files

#(i) for option i: GTDB-based chimeras removed #
mkdir ./8.otu_binning/fasta_files/GTDB
find ./8.otu_binning/GTDB/ -name "*.fasta" -type f -exec cp {} ./8.otu_binning/fasta_files/GTDB \;
for i in ./8.otu_binning/fasta_files/GTDB/*fasta; 
do
qiime tools import \
    --input-path "$i" \
    --output-path ./9.tax_assign/qiime_input/GTDB/"$(basename "$i" .fasta)".qza \
    --type 'FeatureData[Sequence]' \
    --input-format 'MixedCaseDNAFASTAFormat'
done    

#(ii) for option ii: RefSeq-based chimeras removed #
mkdir ./8.otu_binning/fasta_files/RefSeq
find ./8.otu_binning/RefSeq/ -name "*.fasta" -type f -exec cp {} ./8.otu_binning/fasta_files/RefSeq \;
for i in ./8.otu_binning/fasta_files/RefSeq/*fasta; 
do
qiime tools import \
    --input-path "$i" \
    --output-path ./9.tax_assign/qiime_input/RefSeq/"$(basename "$i" .fasta)".qza \
    --type 'FeatureData[Sequence]' \
    --input-format 'MixedCaseDNAFASTAFormat'
done 

#(ii) for option iii: rrnDB-based chimeras removed #
mkdir ./8.otu_binning/fasta_files/rrnDB
find ./8.otu_binning/RefSeq/ -name "*.fasta" -type f -exec cp {} ./8.otu_binning/fasta_files/RefSeq \;
for i in ./8.otu_binning/fasta_files/RefSeq/*fasta; 
do
qiime tools import \
    --input-path "$i" \
    --output-path ./9.tax_assign/qiime_input/RefSeq/"$(basename "$i" .fasta)".qza \
    --type 'FeatureData[Sequence]' \
    --input-format 'MixedCaseDNAFASTAFormat'
done 

