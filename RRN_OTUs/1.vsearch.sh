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

#mkdir ./3.all_derep_fulllength_pb
#mkdir ./3.all_derep_fulllength_pb/uc_files
#mkdir ./4.precluster_before_chimera_detection
#mkdir ./4.precluster_before_chimera_detection/uc_file
#mkdir ./5.de_novo_chimera_removal
#mkdir ./6.ref_based_chimera_removal
#mkdir ./6.ref_based_chimera_removal/GTDB_NR
#mkdir ./6.ref_based_chimera_removal/RefSeq_NR
#mkdir ./6.ref_based_chimera_removal/rrnDB

#for i in $(cat primerpairs.txt)
#do
## 4. dereplicating across all the samples i.e. from the combined fasta files for each primer pair ## 
   #$vsearch --derep_fulllength ./2.per_sample_derep_fulllength/per_primer/"$i"_derep.fasta \
      #--minuniquesize 2 \
      #--sizein \
      #--sizeout \
      #--fasta_width 0 \
      #--uc ./3.all_derep_fulllength_pb/uc_files/"$i"_all_derep.uc \
      #--output ./3.all_derep_fulllength_pb/"$i"_all_derep.fasta 

## 5. precluster at 98% before chimera detection ##
   #$vsearch --cluster_size ./3.all_derep_fulllength_pb/"$i"_all_derep.fasta \
     #--threads $10 \
     #--id 0.98 \
     #--strand both \
     #--sizein \
     #--sizeout \
     #--fasta_width 0 \
     #--uc ./4.precluster_before_chimera_detection/uc_file/"$i"_all_preclustered.uc \
     #--centroids ./4.precluster_before_chimera_detection/"$i"_all_preclustered.fasta

## 6. de novo chimera detection ##
    #$vsearch --uchime_denovo ./4.precluster_before_chimera_detection/"$i"_all_preclustered.fasta \
     #--sizein \
     #--sizeout \
     #--fasta_width 0 \
     #--nonchimeras ./5.de_novo_chimera_removal/"$i"_all_denovo_nonchimeras.fasta

## 7. reference based chimera removal ## option i: with FANGORN GTDB NR nrRep as the reference database downloaded from https://melbourne.figshare.com/articles/dataset/Fangorn_rrn_Database/20086916?file=37794402##
   #$vsearch --uchime_ref ./5.de_novo_chimera_removal/"$i"_all_denovo_nonchimeras.fasta \
    #--db /data/Food/analysis/R1150_biotransformation/Fangorn_db/GTDB_NR_207/nrRep.fasta \
    #--sizein \
    #--sizeout \
    #--fasta_width 0 \
    #--nonchimeras ./6.ref_based_chimera_removal/GTDB_NR/"$i"_all_GTDB_ref_nonchimeras.fasta    

## 7. reference based chimera removal ## option ii: with FANGORN RefSeq NR nrRep as the reference database downloaded from https://melbourne.figshare.com/articles/dataset/Fangorn_rrn_Database/20086916?file=37794402##
   #$vsearch --uchime_ref ./5.de_novo_chimera_removal/"$i"_all_denovo_nonchimeras.fasta \
    #--db /data/Food/analysis/R1150_biotransformation/Fangorn_db/RefSeq_NR_207/nrRep.fasta \
    #--sizein \
    #--sizeout \
    #--fasta_width 0 \
    #--nonchimeras ./6.ref_based_chimera_removal/RefSeq_NR/"$i"_all_RefSeq_ref_nonchimeras.fasta
    
## 7. reference based chimera removal ## option iii: with rrnDBv2 as the reference database
   #$vsearch --uchime_ref ./5.de_novo_chimera_removal/"$i"_all_denovo_nonchimeras.fasta \
    #--db /data/Food/analysis/R1150_biotransformation/other_rrn_databases/rrn_DBv2.fasta \
    #--sizein \
    #--sizeout \
    #--fasta_width 0 \
    #--nonchimeras ./6.ref_based_chimera_removal/rrnDB/"$i"_all_rrnDB_ref_nonchimeras.fasta

#done     

## nanoplot to look at the no.of sequneces and distribution of the read lengths from chimera removal (de novo and ref-based) ##
#mkdir ./5.de_novo_chimera_removal/nanoplots
#mkdir ./6.ref_based_chimera_removal/GTDB_nanoplots
#mkdir ./6.ref_based_chimera_removal/RefSeq_nanoplots

#module load nanoplot/1.28.2
#for i in $(cat primerpairs.txt)
#do
#NanoPlot --fasta ./5.de_novo_chimera_removal/"$i"_all_denovo_nonchimeras.fasta --huge -o ./5.de_novo_chimera_removal/nanoplots/"$i"/
#NanoPlot --fasta  ./6.ref_based_chimera_removal/GTDB_NR/"$i"_all_GTDB_ref_nonchimeras.fasta --huge -o ./6.ref_based_chimera_removal/GTDB_nanoplots/"$i"/
#NanoPlot --fasta ./6.ref_based_chimera_removal/RefSeq_NR/"$i"_all_RefSeq_ref_nonchimeras.fasta --huge -o ./6.ref_based_chimera_removal/RefSeq_nanoplots/"$i"/
#done

#module unload nanoplot/1.28.2

## 8. now that we have the unique sequences for all the reads across samples, we're putting those back against the individual reads ##
