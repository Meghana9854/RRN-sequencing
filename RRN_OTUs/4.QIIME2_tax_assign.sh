## Assigning taxonomy using QIIME2 classifiers ## # Performed only on PacBio data
## classifiers used: QIIME2 BLAST (QB), NB trained classifiers (QNBT), and VSEARCH exact match + sklearn (QVPSK)

## I. OPTION 1: using QIIME2 BLAST (QB) ##
mkdir ./9.tax_assign/qiime_output/
mkdir ./9.tax_assign/qiime_output/blast

# depending on the database that is to be used provide the following: 
# path to the location of the input fasta file 
input = ./9.tax_assign/qiime_input/database_name  #e.g. file path for GTDB-based chimera removed sequences would be ./9.tax_assign/qiime_input/GTDB (from step 11.2 in 2.QIIME2_import.sh)

# path to reference reads in qza format (needs to be imported into QIIME2 prior to this step)
ref_reads = ./rrn_database/database_seqs.qza       #e.g. file path for GTDB reference reads would be ./FANGORN/GTDB_seqs.qza  

# path to reference taxonomoy in qza format (needs to be imported into QIIME2 prior to this step)
ref_tax = ./rrn_databases/for_qiime/database_taxonomy.qza #e.g. file path for GTDB reference taxonomy would be ./rrn_databases/for_qiime/GTDB_taxonomy.qza

for i in $input/*_otus.qza
do
qiime feature-classifier classify-consensus-blast \
  --i-query "$i" \
  --i-reference-reads $ref_reads \
  --i-reference-taxonomy $ref_tax \
  --o-classification ./9.tax_assign/qiime_output/blast/"$(basename "$i" _otus.qza)"_BLAST.qza \
  --o-search-results ./9.tax_assign/qiime_output/blast/"$(basename "$i" _otus.qza)"_out_BLAST.qza
done

## II. OPTION 2: using QIIME2 NB_training (QNBT) ##
mkdir ./9.tax_assign/qiime_output/nb_trained

# individually performed for each RRN primer pair and for each database, for this provide the following:
# path to the location of the input fasta file 
input = ./9.tax_assign/qiime_input/database_name  #e.g. file path for rrnDB-based chimera removed sequences would be ./9.tax_assign/qiime_input/rrnDB (from step 11.2 in 2.QIIME2_import.sh)

# path to the NB trained classifier (details provided in 3.QIIME2_NB_classifiers)
class = ./nb_classifier/database_name_primer_pair.qza #e.g. file path for rrnDB pair A (519F-2428R) primer pair ./nb_classifier/rrnDB_pairA.qza

## 1. pair A: 519F-2428R ##
for i in $input/519F_2428R_*_otus.qza
do
qiime feature-classifier classify-sklearn \
  --i-classifier $class \
  --i-reads "$i" \
  --o-classification ./9.tax_assign/qiime_output/nb_trained/"$(basename "$i" _otus.qza)"_NBt.qza
done

## 2. pair B: 27Fv2-2241R ##
for i in $input/27Fv2_2241R_*_otus.qza
do
qiime feature-classifier classify-sklearn \
  --i-classifier $class \
  --i-reads "$i" \
  --o-classification ./9.tax_assign/qiime_output/nb_trained/"$(basename "$i" _otus.qza)"_NBt.qza
done

## 3. pair C: 27Fv2-2428R ##
for i in $input/27Fv2_2428R_*_otus.qza
do
qiime feature-classifier classify-sklearn \
  --i-classifier $input \
  --i-reads "$i" \
  --o-classification ./9.tax_assign/qiime_output/nb_trained/"$(basename "$i" _otus.qza)"_NBt.qza
done

## 4. pair D: 519F-2241R ##
for i in $input/519F_2241R_*_otus.qza
do
qiime feature-classifier classify-sklearn \
  --i-classifier $class \
  --i-reads "$i" \
  --o-classification ./9.tax_assign/qiime_output/nb_trained/"$(basename "$i" _otus.qza)"_NBt.qza
done

## III. OPTION 3: using QIIME2 vsearch exact match + sklearn (QVPSK) ##
mkdir ./9.tax_assign/qiime_output/vsearch_plus_sklearn

# individually performed for each RRN primer pair and for each database, for this provide the following:
# path to the location of the input fasta file 
input = ./9.tax_assign/qiime_input/database_name  #e.g. file path for rrnDB-based chimera removed sequences would be ./9.tax_assign/qiime_input/rrnDB (from step 11.2 in 2.QIIME2_import.sh)

# path to reference reads in qza format (needs to be imported into QIIME2 prior to this step)
ref_reads = ./rrn_database/database_seqs.qza       #e.g. file path for rrnDB reference reads would be ./rrn_databases/GTDB_seqs.qza  

# path to reference taxonomoy in qza format (needs to be imported into QIIME2 prior to this step)
ref_tax = ./rrn_databases/for_qiime/database_taxonomy.qza #e.g. file path for rrnDB reference taxonomy would be ./rrn_databases/for_qiime/rrnDB_taxonomy.qza

# path to the NB trained classifier (details provided in 3.QIIME2_NB_classifiers)
class = ./nb_classifier/database_name_primer_pair.qza #e.g. file path for rrnDB pair A (519F-2428R) primer pair ./nb_classifier/rrnDB_pairA.qza

## 1. pair A: 519F-2428R ##
for i in $input/519F_2428R_*_otus.qza
do
qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query "$i" \
  --i-reference-reads $ref_reads \
  --i-reference-taxonomy $ref_tax \
  --i-classifier $class \
  --o-classification ./9.tax_assign/qiime_output/vsearch_plus_sklearn/"$(basename "$i" .qza)"_vpsk.qza
done

## 2. pair B: 27Fv2-2241R ##
for i in $input/27Fv2_2241R_*_otus.qza
do
qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query "$i" \
  --i-reference-reads $ref_reads \
  --i-reference-taxonomy $ref_tax \
  --i-classifier $class \
  --o-classification ./9.tax_assign/qiime_output/vsearch_plus_sklearn/"$(basename "$i" .qza)"_vpsk.qza
done

## 3. pair C: 27Fv2-2428R ##
for i in $input/27Fv2_2428R_rrnDBv2_*_otus.qza
do
qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query "$i" \
  --i-reference-reads $ref_reads \
  --i-reference-taxonomy $ref_tax \
  --i-classifier $class \
  --o-classification ./9.tax_assign/qiime_output/vsearch_plus_sklearn/"$(basename "$i" .qza)"_vpsk.qza
done

## 4. pair D: 519F-2241R ##
for i in $input/519F_2241R_*_otus.qza
do
qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query "$i" \
  --i-reference-reads $ref_reads \
  --i-reference-taxonomy $ref_tax \
  --i-classifier $class \
  --o-classification ./9.tax_assign/qiime_output/vsearch_plus_sklearn/"$(basename "$i" .qza)"_vpsk.qza
done

conda deactivate
