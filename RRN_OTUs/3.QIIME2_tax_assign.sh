## Assigning taxonomy using QIIME2 classifiers ## # Performed only on PacBio data
## classifiers used: QIIME2 BLAST (QB), NB trained classifiers (QNBT), and VSEARCH exact match + sklearn (QVPSK)

## option 1: using QIIME2 BLAST ##
#mkdir ./9.tax_assign/qiime_output/rrn_DBv2/
#mkdir ./9.tax_assign/qiime_output/rrn_DBv2/blast

#for i in ./9.tax_assign/qiime_input/rrn_DBv2/*_otus.qza
#do
#qiime feature-classifier classify-consensus-blast \
  #--i-query "$i" \
  #--i-reference-reads /data/Food/analysis/R1150_biotransformation/other_rrn_databases/for_qiime/rrn_DBv2_seqs.qza \
  #--i-reference-taxonomy /data/Food/analysis/R1150_biotransformation/other_rrn_databases/for_qiime/rrn_DBv2_taxonomy_qiime_formatted.qza \
  #--o-classification ./9.tax_assign/qiime_output/rrn_DBv2/blast/"$(basename "$i" _otus.qza)"_rrnDBv2_BLAST.qza \
  #--o-search-results ./9.tax_assign/qiime_output/rrn_DBv2/blast/"$(basename "$i" _otus.qza)"_out_rrnDBv2_BLAST.qza
#done

## option 2: using QIIME2 NB_training ##
#mkdir ./9.tax_assign/qiime_output/rrn_DBv2/nb_trained

## 1. pair A: 519F-2428R ##
for i in ./9.tax_assign/qiime_input/rrn_DBv2/519F_2428R_rrnDBv2_*_otus.qza
do
qiime feature-classifier classify-sklearn \
  --i-classifier /data/Food/analysis/R1150_biotransformation/other_rrn_databases/for_qiime/NB_classifers/classifier_rrn_DBv2_5L_qiime_formatted.qza \
  --i-reads "$i" \
  --o-classification ./9.tax_assign/qiime_output/rrn_DBv2/nb_trained/"$(basename "$i" _otus.qza)"_rrnDBv2_NBt.qza
done

## 2. PairB 27Fv2-2241R ##
for i in ./9.tax_assign/qiime_input/rrn_DBv2/27Fv2_2241R_rrnDBv2_*_otus.qza
do
qiime feature-classifier classify-sklearn \
  --i-classifier /data/Food/analysis/R1150_biotransformation/other_rrn_databases/for_qiime/NB_classifers/classifier_rrn_DBv2_2S_qiime_formatted.qza \
  --i-reads "$i" \
  --o-classification ./9.tax_assign/qiime_output/rrn_DBv2/nb_trained/"$(basename "$i" _otus.qza)"_rrnDBv2_NBt.qza
done

## 3. PairC 27Fv2-2428R ##
for i in ./9.tax_assign/qiime_input/rrn_DBv2/27Fv2_2428R_rrnDBv2_*_otus.qza
do
qiime feature-classifier classify-sklearn \
  --i-classifier /data/Food/analysis/R1150_biotransformation/other_rrn_databases/for_qiime/NB_classifers/classifier_rrn_DBv2_2L_qiime_formatted.qza \
  --i-reads "$i" \
  --o-classification ./9.tax_assign/qiime_output/rrn_DBv2/nb_trained/"$(basename "$i" _otus.qza)"_rrnDBv2_NBt.qza
done

## 4. pair D: 519F-2241R ##
for i in ./9.tax_assign/qiime_input/rrn_DBv2/519F_2241R_rrnDBv2_*_otus.qza
do
qiime feature-classifier classify-sklearn \
  --i-classifier /data/Food/analysis/R1150_biotransformation/other_rrn_databases/for_qiime/NB_classifers/classifier_rrn_DBv2_5S_qiime_formatted.qza \
  --i-reads "$i" \
  --o-classification ./9.tax_assign/qiime_output/rrn_DBv2/nb_trained/"$(basename "$i" _otus.qza)"_rrnDBv2_NBt.qza
done

##(ii) using the vsearch exact match + sklearn using FANGORN database ## 
mkdir ./9.tax_assign/qiime_output/rrn_DBv2/vsearch_plus_sklearn

rrn_db=/data/Food/analysis/R1150_biotransformation/other_rrn_databases/for_qiime/rrn_DBv2_seqs.qza
rrn_tax=/data/Food/analysis/R1150_biotransformation/other_rrn_databases/for_qiime/rrn_DBv2_taxonomy_qiime_formatted.qza

## 1. pair A: 519F-2428R ##
for i in ./9.tax_assign/qiime_input/rrn_DBv2/519F_2428R_rrnDBv2_*_otus.qza
do
qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query "$i" \
  --i-reference-reads $rrn_db \
  --i-reference-taxonomy $rrn_tax \
  --i-classifier /data/Food/analysis/R1150_biotransformation/other_rrn_databases/for_qiime/NB_classifers/classifier_rrn_DBv2_5L_qiime_formatted.qza \
  --o-classification ./9.tax_assign/qiime_output/rrn_DBv2/vsearch_plus_sklearn/"$(basename "$i" .qza)"_rrnDBv2_vpsk.qza
done

## 2. PairB 27Fv2-2241R ##
for i in ./9.tax_assign/qiime_input/rrn_DBv2/27Fv2_2241R_rrnDBv2_*_otus.qza
do
qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query "$i" \
  --i-reference-reads $rrn_db \
  --i-reference-taxonomy $rrn_tax \
  --i-classifier /data/Food/analysis/R1150_biotransformation/other_rrn_databases/for_qiime/NB_classifers/classifier_rrn_DBv2_2S_qiime_formatted.qza \
  --o-classification ./9.tax_assign/qiime_output/rrn_DBv2/vsearch_plus_sklearn/"$(basename "$i" .qza)"_rrnDBv2_vpsk.qza
done

## 3. PairC 27Fv2-2428R ##
for i in ./9.tax_assign/qiime_input/rrn_DBv2/27Fv2_2428R_rrnDBv2_*_otus.qza
do
qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query "$i" \
  --i-reference-reads $rrn_db \
  --i-reference-taxonomy $rrn_tax \
  --i-classifier /data/Food/analysis/R1150_biotransformation/other_rrn_databases/for_qiime/NB_classifers/classifier_rrn_DBv2_2L_qiime_formatted.qza \
  --o-classification ./9.tax_assign/qiime_output/rrn_DBv2/vsearch_plus_sklearn/"$(basename "$i" .qza)"_rrnDBv2_vpsk.qza
done

## 4. pair D: 519F-2241R ##
for i in ./9.tax_assign/qiime_input/rrn_DBv2/519F_2241R_rrnDBv2_*_otus.qza
do
qiime feature-classifier classify-hybrid-vsearch-sklearn \
  --i-query "$i" \
  --i-reference-reads $rrn_db \
  --i-reference-taxonomy $rrn_tax \
  --i-classifier /data/Food/analysis/R1150_biotransformation/other_rrn_databases/for_qiime/NB_classifers/classifier_rrn_DBv2_5S_qiime_formatted.qza \
  --o-classification ./9.tax_assign/qiime_output/rrn_DBv2/vsearch_plus_sklearn/"$(basename "$i" .qza)"_rrnDBv2_vpsk.qza
done

conda deactivate
