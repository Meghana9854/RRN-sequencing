## Assigning taxonomy to the OTU tables using QIIME2 with NB trained classifiers for SILVA and Greengenes2 ##

## importing the biom and fasta files from Step 10 (script vsearch.sh) into QIIME2 ##
# first fomatting the .biom files #
mkdir ./11.OTU_binning/biom_formatted_files/
find ./11.OTU_binning/ -name "*.biom" -type f -exec cp {} ./11.OTU_binning/biom_formatted_files \;
module load biom/2.7.1
for i in ./11.OTU_binning/biom_formatted_files/*.biom
do 
biom convert -i "$i" -o "${i%.*}.biomv210.biom" --table-type="OTU table" --to-hdf5
done
rm ./11.OTU_binning/biom_formatted_files/*_pc.biom
module unload biom/2.7.1  

## Load QIIME2 ##
module load qiime2/2023.2
source activate qiime2-2023.2

## Step 1.1: importing the biom files into QIIME2 ##
mkdir ./12.tax_assign/qiime_inputs
for i in ./11.OTU_binning/biom_formatted_files/*.biomv210.biom
do
   qiime tools import \
     --input-path "$i" \
     --type 'FeatureTable[Frequency]' \
     --input-format BIOMV210Format \
     --output-path ./12.tax_assign/qiime_inputs/"$(basename "$i" .biomv210.biom)".qza
done       

## Step 1.2: importing the fasta files into QIIME2 ##
for i in ./11.OTU_binning/fasta_files/*.fasta
do
qiime tools import \
    --input-path "$i" \
    --output-path ./12.tax_assign/qiime_inputs/"$(basename "$i" .fasta)".qza \
    --type 'FeatureData[Sequence]' \
    --input-format 'MixedCaseDNAFASTAFormat'
done

## Assigning taxonomy with QIIME2 using NB trained classfiers for SILVA and Greengenes2 ##
## SILVA NB classifier for the v3-v4 region was generated as per https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494
## Greengenes2 NB classifier for the v3-v4 region was generated as per https://github.com/biocore/greengenes2/blob/main/scripts/classifier.sh#L31-L41, and backbone files for reference reads (2022.10.backbone.full-length.fna.qza) and taxonomy (2022.10.backbone.tax.qza) were taken from http://ftp.microbio.me/greengenes_release/2022.10/ 

# general over view of the NB training scripts is as below
mkdir ./nb_training
mkdir ./nb_training/output_files
# provide the path to the refernce reads (in QIIME2 .qza format)
ref_reads = ./nb_training/database_seqs.fna.qza   #e.g. for the Greengnes2 database file path would be ./nb_training/2022.10.backbone.full-length.fna.qza
# provide the path to the taxonomy file (in QIIME2 .qza format)
ref_tax = ./nb_training/database_tax.qza        #e.g. for the Greengnes2 database file path would be ./nb_training/2022.10.backbone.tax.qza

# extracting the v3-v4 reads from the reference reads #
qiime feature-classifier extract-reads \
    --i-sequences $ref_reads \
    --p-f-primer CCTACGGGNGGCWGCAG --p-r-primer GACTACHVGGGTATCTAATCC \
    --p-read-orientation both \
    --o-reads ./nb_training/output_files/reads_v3_v4.qza \
    --p-n-jobs 17

# NB training the classifier #
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads ./nb_training/output_files/reads_v3_v4.qza \
    --i-reference-taxonomy $ref_tax \
    --o-classifier ./nb_training/output_files/classifier_v3_v4.qza

conda deactivate
