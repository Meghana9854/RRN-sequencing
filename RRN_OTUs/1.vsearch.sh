## Generating OTUs from RRN amplicons using vsearch   # Performed only on PacBio data
## as per the vserach github page (https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline)

### for this script to work you need to have your sample names in the working directory with the unique part of all samples names (in samplenames.txt)
### and a file for the primer pairs (in primerpairs.txt)

# give the location of where vsearch is stored so it can be called by $vsearch later throughout the script #
vsearch = ./vsearch-2.22.1/bin/vsearch
# this script runs from the 2.new_vsearch directory #
mkdir ./2.new_vsearch
cd ./2.new_vsearch/

mkdir 1.length_filtered
mkdir 2.per_sample_derep_fulllength
mkdir ./2.per_sample_derep_fulllength/uc_files

for i in $(cat samplenames.txt)
do
## 1. keeping reads only between 3500 to 5000 bp ##
   $vsearch --fastq_filter ./1.data/"$i"_qls.fastq.gz \
       --fastq_qmax 93 \
       --fastq_minlen 3500 \
       --fastq_maxlen 5000 \
       --fastaout ./1.length_filtered/"$i"_pb_lenfiltered.fasta\
       --fasta_width 0   
## 2. derepelicating reads with each sample ##
   $vsearch --derep_fulllength ./1.length_filtered/"$i"_pb_lenfiltered.fasta \
       --strand both \
       --output ./2.per_sample_derep_fulllength/"$i"_pb_derep.fasta \
       --sizeout \
       --uc ./2.per_sample_derep_fulllength/uc_files/"$i"_pb_derep.uc \
       --relabel "$i". \
       --fasta_width 0         
done       

## using nanoplot to look at the number of sequences and distribution of the read lengths ##
mkdir ./1.length_filtered/nanoplots
cat ./1.length_filtered/*_pb_lenfiltered.fasta >> ./1.length_filtered/nanoplots/onefile.fasta
mkdir ./2.per_sample_derep_fulllength/nanoplots
cat ./2.per_sample_derep_fulllength/*_pb_derep.fasta >> ./2.per_sample_derep_fulllength/nanoplots/onefile.fasta
## running nanoplot on the cat files ##
module load nanoplot/1.28.2
NanoPlot --fasta ./1.length_filtered/nanoplots/onefile.fasta --huge -o ./1.length_filtered/nanoplots
NanoPlot --fasta ./2.per_sample_derep_fulllength/nanoplots/onefile.fasta --huge -o ./2.per_sample_derep_fulllength/nanoplots
module unload nanoplot/1.28.2

## 3. merging all the derep.fasta files from each sample into one fasta file, but keeping the RRN primer pairs separate to generate OTU tables for each pair (needed for downstream analyis) ##

mkdir ./2.per_sample_derep_fulllength/per_primer/
## for the 27Fv2_2241R (2S) primer pair ##
cat ./2.per_sample_derep_fulllength/*_27Fv2_2241R_pb_derep.fasta >> ./2.per_sample_derep_fulllength/per_primer/27Fv2_2241R_derep.fasta
## for the 27Fv2_2428R (2L) primer pair ##
cat ./2.per_sample_derep_fulllength/*_27Fv2_2428R_pb_derep.fasta >> ./2.per_sample_derep_fulllength/per_primer/27Fv2_2428R_derep.fasta
## for the 519F_2241R (5S) primer pair ##
cat ./2.per_sample_derep_fulllength/*_519F_2241R_pb_derep.fasta >> ./2.per_sample_derep_fulllength/per_primer/519F_2241R_derep.fasta
## for the 519F_2428R (5L) primer pair ##
cat ./2.per_sample_derep_fulllength/*_519F_2428R_pb_derep.fasta >> ./2.per_sample_derep_fulllength/per_primer/519F_2428R_derep.fasta

mkdir ./3.all_derep_fulllength_pb
mkdir ./3.all_derep_fulllength_pb/uc_files
mkdir ./4.precluster_before_chimera_detection
mkdir ./4.precluster_before_chimera_detection/uc_file
mkdir ./5.de_novo_chimera_removal
mkdir ./6.ref_based_chimera_removal
mkdir ./6.ref_based_chimera_removal/GTDB_NR
mkdir ./6.ref_based_chimera_removal/RefSeq_NR
mkdir ./6.ref_based_chimera_removal/rrnDB

for i in $(cat primerpairs.txt)
do
## 4. dereplicating across all the samples using the combined fasta files for each primer pair ## 
   $vsearch --derep_fulllength ./2.per_sample_derep_fulllength/per_primer/"$i"_derep.fasta \
      --minuniquesize 2 \
      --sizein \
      --sizeout \
      --fasta_width 0 \
      --uc ./3.all_derep_fulllength_pb/uc_files/"$i"_all_derep.uc \
      --output ./3.all_derep_fulllength_pb/"$i"_all_derep.fasta 

## 5. precluster at 98% before chimera detection ##
   $vsearch --cluster_size ./3.all_derep_fulllength_pb/"$i"_all_derep.fasta \
     --threads $10 \
     --id 0.98 \
     --strand both \
     --sizein \
     --sizeout \
     --fasta_width 0 \
     --uc ./4.precluster_before_chimera_detection/uc_file/"$i"_all_preclustered.uc \
     --centroids ./4.precluster_before_chimera_detection/"$i"_all_preclustered.fasta

## 6. de novo chimera detection ##
    $vsearch --uchime_denovo ./4.precluster_before_chimera_detection/"$i"_all_preclustered.fasta \
     --sizein \
     --sizeout \
     --fasta_width 0 \
     --nonchimeras ./5.de_novo_chimera_removal/"$i"_all_denovo_nonchimeras.fasta

## 7. reference based chimera removal 
## option i: with FANGORN GTDB NR nrRep as the reference database
   $vsearch --uchime_ref ./5.de_novo_chimera_removal/"$i"_all_denovo_nonchimeras.fasta \
    --db ./Fangorn_db/GTDB_NR_207/nrRep.fasta \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --nonchimeras ./6.ref_based_chimera_removal/GTDB_NR/"$i"_all_GTDB_ref_nonchimeras.fasta    

## option ii: with FANGORN RefSeq NR nrRep as the reference database
   $vsearch --uchime_ref ./5.de_novo_chimera_removal/"$i"_all_denovo_nonchimeras.fasta \
    --db ./Fangorn_db/RefSeq_NR_207/nrRep.fasta \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --nonchimeras ./6.ref_based_chimera_removal/RefSeq_NR/"$i"_all_RefSeq_ref_nonchimeras.fasta
    
## option iii: with rrnDBv2 as the reference database
   $vsearch --uchime_ref ./5.de_novo_chimera_removal/"$i"_all_denovo_nonchimeras.fasta \
    --db ./rrn_databases/rrn_DBv2.fasta \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --nonchimeras ./6.ref_based_chimera_removal/rrnDB/"$i"_all_rrnDB_ref_nonchimeras.fasta

done     

## using nanoplot to look at the number of sequences and distribution of the read lengths after chimera removal (de novo and ref-based) ##
mkdir ./5.de_novo_chimera_removal/nanoplots
mkdir ./6.ref_based_chimera_removal/GTDB_nanoplots
mkdir ./6.ref_based_chimera_removal/RefSeq_nanoplots
mkdir ./6.ref_based_chimera_removal/rrnDB_nanoplots
module load nanoplot/1.28.2
for i in $(cat primerpairs.txt)
do
NanoPlot --fasta ./5.de_novo_chimera_removal/"$i"_all_denovo_nonchimeras.fasta --huge -o ./5.de_novo_chimera_removal/nanoplots/"$i"/
NanoPlot --fasta  ./6.ref_based_chimera_removal/GTDB_NR/"$i"_all_GTDB_ref_nonchimeras.fasta --huge -o ./6.ref_based_chimera_removal/GTDB_nanoplots/"$i"/
NanoPlot --fasta ./6.ref_based_chimera_removal/RefSeq_NR/"$i"_all_RefSeq_ref_nonchimeras.fasta --huge -o ./6.ref_based_chimera_removal/RefSeq_nanoplots/"$i"/
NanoPlot --fasta ./6.ref_based_chimera_removal/rrnDB/"$i"_all_rrnDB_ref_nonchimeras.fasta --huge -o ./6.ref_based_chimera_removal/rrnDB_nanoplots/"$i"/
done
module unload nanoplot/1.28.2

## 8. now that we have the unique sequences for all the reads across samples, we're putting those back against the individual reads ##
mkdir ./7.map_pl_output
mkdir ./7.map_pl_output/intermediate_files

# NOTE: for the following part of the script to run the pearl script map.pl is needed to be present in the working directory # 
# which was obatined from the vserach github page (https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline) #

## (i) for option i: GTDB-based chimera removal ##
## for the 27Fv2_2241R (2S) primer pair 
echo Extract all non-chimeric, non-singleton sequences, dereplicate
perl ./map.pl ./3.all_derep_fulllength_pb/27Fv2_2241R_all_derep.fasta ./4.precluster_before_chimera_detection/uc_file/27Fv2_2241R_all_preclustered.uc ./6.ref_based_chimera_removal/GTDB_NR/27Fv2_2241R_all_GTDB_ref_nonchimeras.fasta > ./7.map_pl_output/intermediate_files/27Fv2_2241R_GTDB_all_nonchimeras_derep.fasta
echo Extract all non-chimeric, non-singleton sequences in each sample
perl ./map.pl ./2.per_sample_derep_fulllength/per_primer/27Fv2_2241R_derep.fasta ./3.all_derep_fulllength_pb/uc_files/27Fv2_2241R_all_derep.uc ./7.map_pl_output/intermediate_files/27Fv2_2241R_GTDB_all_nonchimeras_derep.fasta > ./7.map_pl_output/27Fv2_2241R_GTDB_all_nonchimeras.fasta

## for the 27Fv2_2428R (2L) primer pair 
echo Extract all non-chimeric, non-singleton sequences, dereplicate
perl ./map.pl ./3.all_derep_fulllength_pb/27Fv2_2428R_all_derep.fasta ./4.precluster_before_chimera_detection/uc_file/27Fv2_2428R_all_preclustered.uc ./6.ref_based_chimera_removal/GTDB_NR/27Fv2_2428R_all_GTDB_ref_nonchimeras.fasta > ./7.map_pl_output/intermediate_files/27Fv2_2428R_GTDB_all_nonchimeras_derep.fasta
echo Extract all non-chimeric, non-singleton sequences in each sample
perl ./map.pl ./2.per_sample_derep_fulllength/per_primer/27Fv2_2428R_derep.fasta ./3.all_derep_fulllength_pb/uc_files/27Fv2_2428R_all_derep.uc ./7.map_pl_output/intermediate_files/27Fv2_2428R_GTDB_all_nonchimeras_derep.fasta > ./7.map_pl_output/27Fv2_2428R_GTDB_all_nonchimeras.fasta

## for the 519F_2241R (5S) primer pair 
echo Extract all non-chimeric, non-singleton sequences, dereplicate
perl ./map.pl ./3.all_derep_fulllength_pb/519F_2241R_all_derep.fasta ./4.precluster_before_chimera_detection/uc_file/519F_2241R_all_preclustered.uc ./6.ref_based_chimera_removal/GTDB_NR/519F_2241R_all_GTDB_ref_nonchimeras.fasta > ./7.map_pl_output/intermediate_files/519F_2241R_GTDB_all_nonchimeras_derep.fasta
echo Extract all non-chimeric, non-singleton sequences in each sample
perl ./map.pl ./2.per_sample_derep_fulllength/per_primer/519F_2241R_derep.fasta ./3.all_derep_fulllength_pb/uc_files/519F_2241R_all_derep.uc ./7.map_pl_output/intermediate_files/519F_2241R_GTDB_all_nonchimeras_derep.fasta > ./7.map_pl_output/519F_2241R_GTDB_all_nonchimeras.fasta

## for the 519F_2428R (5L) primer pair 
echo Extract all non-chimeric, non-singleton sequences, dereplicate
perl ./map.pl ./3.all_derep_fulllength_pb/519F_2428R_all_derep.fasta ./4.precluster_before_chimera_detection/uc_file/519F_2428R_all_preclustered.uc ./6.ref_based_chimera_removal/GTDB_NR/519F_2428R_all_GTDB_ref_nonchimeras.fasta > ./7.map_pl_output/intermediate_files/519F_2428R_GTDB_all_nonchimeras_derep.fasta
echo Extract all non-chimeric, non-singleton sequences in each sample
perl ./map.pl ./2.per_sample_derep_fulllength/per_primer/519F_2428R_derep.fasta ./3.all_derep_fulllength_pb/uc_files/519F_2428R_all_derep.uc ./7.map_pl_output/intermediate_files/519F_2428R_GTDB_all_nonchimeras_derep.fasta > ./7.map_pl_output/519F_2428R_GTDB_all_nonchimeras.fasta


## (ii) for option ii: RefSeq-based chimera removal ##
## for the 27Fv2_2241R (2S) primer pair 
echo Extract all non-chimeric, non-singleton sequences, dereplicate
perl ./map.pl ./3.all_derep_fulllength_pb/27Fv2_2241R_all_derep.fasta ./4.precluster_before_chimera_detection/uc_file/27Fv2_2241R_all_preclustered.uc ./6.ref_based_chimera_removal/RefSeq_NR/27Fv2_2241R_all_RefSeq_ref_nonchimeras.fasta > ./7.map_pl_output/intermediate_files/27Fv2_2241R_RefSeq_all_nonchimeras_derep.fasta
echo Extract all non-chimeric, non-singleton sequences in each sample
perl ./map.pl ./2.per_sample_derep_fulllength/per_primer/27Fv2_2241R_derep.fasta ./3.all_derep_fulllength_pb/uc_files/27Fv2_2241R_all_derep.uc ./7.map_pl_output/intermediate_files/27Fv2_2241R_RefSeq_all_nonchimeras_derep.fasta > ./7.map_pl_output/27Fv2_2241R_RefSeq_all_nonchimeras.fasta

## for the 27Fv2_2428R (2L) primer pair 
echo Extract all non-chimeric, non-singleton sequences, dereplicate
perl ./map.pl ./3.all_derep_fulllength_pb/27Fv2_2428R_all_derep.fasta ./4.precluster_before_chimera_detection/uc_file/27Fv2_2428R_all_preclustered.uc ./6.ref_based_chimera_removal/RefSeq_NR/27Fv2_2428R_all_RefSeq_ref_nonchimeras.fasta > ./7.map_pl_output/intermediate_files/27Fv2_2428R_RefSeq_all_nonchimeras_derep.fasta
echo Extract all non-chimeric, non-singleton sequences in each sample
perl ./map.pl ./2.per_sample_derep_fulllength/per_primer/27Fv2_2428R_derep.fasta ./3.all_derep_fulllength_pb/uc_files/27Fv2_2428R_all_derep.uc ./7.map_pl_output/intermediate_files/27Fv2_2428R_RefSeq_all_nonchimeras_derep.fasta > ./7.map_pl_output/27Fv2_2428R_RefSeq_all_nonchimeras.fasta

## for the 519F_2241R (5S) primer pair 
echo Extract all non-chimeric, non-singleton sequences, dereplicate
perl ./map.pl ./3.all_derep_fulllength_pb/519F_2241R_all_derep.fasta ./4.precluster_before_chimera_detection/uc_file/519F_2241R_all_preclustered.uc ./6.ref_based_chimera_removal/RefSeq_NR/519F_2241R_all_RefSeq_ref_nonchimeras.fasta > ./7.map_pl_output/intermediate_files/519F_2241R_RefSeq_all_nonchimeras_derep.fasta
echo Extract all non-chimeric, non-singleton sequences in each sample
perl ./map.pl ./2.per_sample_derep_fulllength/per_primer/519F_2241R_derep.fasta ./3.all_derep_fulllength_pb/uc_files/519F_2241R_all_derep.uc ./7.map_pl_output/intermediate_files/519F_2241R_RefSeq_all_nonchimeras_derep.fasta > ./7.map_pl_output/519F_2241R_RefSeq_all_nonchimeras.fasta

## for the 519F_2428R (5L) primer pair 
echo Extract all non-chimeric, non-singleton sequences, dereplicate
perl ./map.pl ./3.all_derep_fulllength_pb/519F_2428R_all_derep.fasta ./4.precluster_before_chimera_detection/uc_file/519F_2428R_all_preclustered.uc ./6.ref_based_chimera_removal/RefSeq_NR/519F_2428R_all_RefSeq_ref_nonchimeras.fasta > ./7.map_pl_output/intermediate_files/519F_2428R_RefSeq_all_nonchimeras_derep.fasta
echo Extract all non-chimeric, non-singleton sequences in each sample
perl ./map.pl ./2.per_sample_derep_fulllength/per_primer/519F_2428R_derep.fasta ./3.all_derep_fulllength_pb/uc_files/519F_2428R_all_derep.uc ./7.map_pl_output/intermediate_files/519F_2428R_RefSeq_all_nonchimeras_derep.fasta > ./7.map_pl_output/519F_2428R_RefSeq_all_nonchimeras.fasta


## (ii) for option 3: rrnDB-based chimera removal ##
## for the 27Fv2_2241R (2S) primer pair 
echo Extract all non-chimeric, non-singleton sequences, dereplicate
perl ./map.pl ./3.all_derep_fulllength_pb/27Fv2_2241R_all_derep.fasta ./4.precluster_before_chimera_detection/uc_file/27Fv2_2241R_all_preclustered.uc ./6.ref_based_chimera_removal/rrnDB/27Fv2_2241R_all_rrnDB_ref_nonchimeras.fasta > ./7.map_pl_output/intermediate_files/27Fv2_2241R_rrnDB_all_nonchimeras_derep.fasta
echo Extract all non-chimeric, non-singleton sequences in each sample
perl ./map.pl ./2.per_sample_derep_fulllength/per_primer/27Fv2_2241R_derep.fasta ./3.all_derep_fulllength_pb/uc_files/27Fv2_2241R_all_derep.uc ./7.map_pl_output/intermediate_files/27Fv2_2241R_rrnDB_all_nonchimeras_derep.fasta > ./7.map_pl_output/27Fv2_2241R_rrnDB_all_nonchimeras.fasta

## for the 27Fv2_2428R (2L) primer pair 
echo Extract all non-chimeric, non-singleton sequences, dereplicate
perl ./map.pl ./3.all_derep_fulllength_pb/27Fv2_2428R_all_derep.fasta ./4.precluster_before_chimera_detection/uc_file/27Fv2_2428R_all_preclustered.uc ./6.ref_based_chimera_removal/rrnDB/27Fv2_2428R_all_rrnDB_ref_nonchimeras.fasta > ./7.map_pl_output/intermediate_files/27Fv2_2428R_rrnDB_all_nonchimeras_derep.fasta
echo Extract all non-chimeric, non-singleton sequences in each sample
perl ./map.pl ./2.per_sample_derep_fulllength/per_primer/27Fv2_2428R_derep.fasta ./3.all_derep_fulllength_pb/uc_files/27Fv2_2428R_all_derep.uc ./7.map_pl_output/intermediate_files/27Fv2_2428R_rrnDB_all_nonchimeras_derep.fasta > ./7.map_pl_output/27Fv2_2428R_rrnDB_all_nonchimeras.fasta

## for the 519F_2241R (5S) primer pair 
echo Extract all non-chimeric, non-singleton sequences, dereplicate
perl ./map.pl ./3.all_derep_fulllength_pb/519F_2241R_all_derep.fasta ./4.precluster_before_chimera_detection/uc_file/519F_2241R_all_preclustered.uc ./6.ref_based_chimera_removal/rrnDB/519F_2241R_all_rrnDB_ref_nonchimeras.fasta > ./7.map_pl_output/intermediate_files/519F_2241R_rrnDB_all_nonchimeras_derep.fasta
echo Extract all non-chimeric, non-singleton sequences in each sample
perl ./map.pl ./2.per_sample_derep_fulllength/per_primer/519F_2241R_derep.fasta ./3.all_derep_fulllength_pb/uc_files/519F_2241R_all_derep.uc ./7.map_pl_output/intermediate_files/519F_2241R_rrnDB_all_nonchimeras_derep.fasta > ./7.map_pl_output/519F_2241R_rrnDB_all_nonchimeras.fasta

## for the 519F_2428R (5L) primer pair 
echo Extract all non-chimeric, non-singleton sequences, dereplicate
perl ./map.pl ./3.all_derep_fulllength_pb/519F_2428R_all_derep.fasta ./4.precluster_before_chimera_detection/uc_file/519F_2428R_all_preclustered.uc ./6.ref_based_chimera_removal/rrnDB/519F_2428R_all_rrnDB_ref_nonchimeras.fasta > ./7.map_pl_output/intermediate_files/519F_2428R_rrnDB_all_nonchimeras_derep.fasta
echo Extract all non-chimeric, non-singleton sequences in each sample
perl ./map.pl ./2.per_sample_derep_fulllength/per_primer/519F_2428R_derep.fasta ./3.all_derep_fulllength_pb/uc_files/519F_2428R_all_derep.uc ./7.map_pl_output/intermediate_files/519F_2428R_rrnDB_all_nonchimeras_derep.fasta > ./7.map_pl_output/519F_2428R_rrnDB_all_nonchimeras.fasta

## 9. Generating the OTU table through clustering, checked over a range of clustering threholds (97%, 98%, 99% , 99.5% nad 99.9%) ##
mkdir ./8.otu_binning

## (i) for option i: GTDB-based chimeras removed ##
mkdir ./8.otu_binning/GTDB
mkdir ./8.otu_binning/GTDB/99.9pc
for i in $(cat primerpairs.txt)
do
## with clustering threshold at 99.9% ##  
   $vsearch --cluster_size ./7.map_pl_output/"$i"_GTDB_all_nonchimeras.fasta \
    --threads $10 \
    --id 0.999 \  #0.97 for 97%, 0.98 for 98%, 0.99 for 99%, 0.995 for 99.5%, and 0.999 for 99.9% clustering was applied#
    --strand both \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc ./8.otu_binning/GTDB/99.9pc/"$i"_GTDB_999pc_clustered.uc \
    --centroids ./8.otu_binning/GTDB/99.9pc/"$i"_GTDB_999pc_otus.fasta \
    --otutabout ./8.otu_binning/GTDB/99.9pc/"$i"_GTDB_999pc_otutab.txt \
    --biomout ./8.otu_binning/GTDB/99.9pc/"$i"_GTDB_999pc_otutab.biom       
done 

## (ii) for option ii: RefSeq-based chimeras removed ##
mkdir ./8.otu_binning/RefSeq
mkdir ./8.otu_binning/RefSeq/99.9pc   
for i in $(cat primerpairs.txt)
do
## with clustering threshold at 99.9% ##  
   $vsearch --cluster_size ./7.map_pl_output/"$i"_RefSeq_all_nonchimeras.fasta \
    --threads $10 \
    --id 0.999 \   #0.97 for 97%, 0.98 for 98%, 0.99 for 99%, 0.995 for 99.5%, and 0.999 for 99.9% clustering was applied#
    --strand both \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc ./8.otu_binning/RefSeq/99.9pc/"$i"_RefSeq_999pc_clustered.uc \
    --centroids ./8.otu_binning/RefSeq/99.9pc/"$i"_RefSeq_999pc_otus.fasta \
    --otutabout ./8.otu_binning/RefSeq/99.9pc/"$i"_RefSeq_999pc_otutab.txt \
    --biomout ./8.otu_binning/RefSeq/99.9pc/"$i"_RefSeq_999pc_otutab.biom    
done 

## (iii) for option iii: rrnDB-based chimeras removed ##
mkdir ./8.otu_binning/rrnDB
mkdir ./8.otu_binning/rrnDB/99.9pc   
for i in $(cat primerpairs.txt)
do
## with clustering threshold at 99.9% ##  
   $vsearch --cluster_size ./7.map_pl_output/"$i"_rrnDB_all_nonchimeras.fasta \
    --threads $10 \
    --id 0.999 \   #0.97 for 97%, 0.98 for 98%, 0.99 for 99%, 0.995 for 99.5%, and 0.999 for 99.9% clustering was applied#
    --strand both \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc ./8.otu_binning/rrnDB/99.9pc/"$i"_rrnDB_999pc_clustered.uc \
    --centroids ./8.otu_binning/rrnDB/99.9pc/"$i"_rrnDB_999pc_otus.fasta \
    --otutabout ./8.otu_binning/rrnDB/99.9pc/"$i"_rrnDB_999pc_otutab.txt \
    --biomout ./8.otu_binning/rrnDB/99.9pc/"$i"_rrnDb_999pc_otutab.biom    
done 

## 10. Formatting vsearch's biome file (OTU table) for input in QIIME2 ##
mkdir ./8.otu_binning/biom_files_formatted

## (i) for option i: GTDB-based chimeras removed ##
mkdir ./8.otu_binning/biom_files_formatted/GTDB
## copying the biom files to another directory before converting ##
find ./8.otu_binning/GTDB/ -name "*.biom" -type f -exec cp {} ./8.otu_binning/biom_files_formatted/GTDB \;
module load biom/2.7.1
for i in ./8.otu_binning/biom_files_formatted/GTDB/*.biom
do
  biom convert -i "$i" -o "${i%.*}.biomv210.biom" --table-type="OTU table" --to-hdf5
done
rm ./8.otu_binning/biom_files_formatted/GTDB/*otutab.biom

## (ii) for option ii: RefSeq-based chimeras removed ##
mkdir ./8.otu_binning/biom_files_formatted/RefSeq
## copying the biom files to another directory before converting ##
find ./8.otu_binning/RefSeq/ -name "*.biom" -type f -exec cp {} ./8.otu_binning/biom_files_formatted/RefSeq \;
for i in ./8.otu_binning/biom_files_formatted/RefSeq/*.biom
do
  biom convert -i "$i" -o "${i%.*}.biomv210.biom" --table-type="OTU table" --to-hdf5
done
rm ./8.otu_binning/biom_files_formatted/RefSeq/*otutab.biom

## (iii) for option iii: rrnDB-based chimeras removed ##
mkdir ./8.otu_binning/biom_files_formatted/rrn_DB
## copying the biom files to another directory before converting ##
find ./8.otu_binning/rrn_DBv2/ -name "*.biom" -type f -exec cp {} ./8.otu_binning/biom_files_formatted/rrn_DB \;
for i in ./8.otu_binning/biom_files_formatted/rrn_DBv2/*.biom
do
  biom convert -i "$i" -o "${i%.*}.biomv210.biom" --table-type="OTU table" --to-hdf5
done
rm ./8.otu_binning/biom_files_formatted/rrn_DBv2/*otutab.biom
module unload biom/2.7.1
