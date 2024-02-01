## Generating OTUs for the v3-v4 data ## Performed only on Illumina data ##

## first removing the adpters and primer sequences from the reads using cutadapt ##
## using paired end mode with cutadapt ## for this script to run a file containing sample names (samplenmaes.txt) and diretory with raw data (1.data) need to be created in the working directory
mkdir ./2.primers_removed 
module load cutadapt/2.6 
adapter_fw='TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
adapter_rev='GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
PrimerF_GI='CCTACGGGNGGCWGCAG'
PrimerR_GI='GACTACHVGGGTATCTAATCC'
for i in $(cat samplesnames.txt)
do
cutadapt -a $adapter_fw -A $adapter_rev -b $PrimerF_GI -B $PrimerR_GI \
-o ./2.primers_removed/"$i"_trimmed.fastq -p ./2.primers_removed/"$i"_trimmed.fastq ./1.reads/"$i".fastq ./1.reads/"$i".fastq \
-e 0.1 \
--discard-untrimmed \
--cores 4
done
module unload cutadapt/2.6 

######## checking the quality of trimmed reads using FastQC #########
## for all reads ##
mkdir ./2.primers_removed/fastqc
cat ./2.primers_removed/*.fastq > ./2.primers_removed/fastqc/onefile_trimmed.fastq
module load fastqc/0.11.8
fastqc ./2.primers_removed/fastqc/onefile_trimmed.fastq -o ./2.primers_removed/fastqc/
## for forward and revrse reads separately ##
# forward reads only
cat ./2.primers_removed/*trimmed_R1.fastq > ./2.primers_removed/fastqc/allforward_trimmed.fastq
fastqc ./2.primers_removed/fastqc/allforward_trimmed.fastq -o ./2.primers_removed/fastqc/
# reverse raeds only 
cat ./2.primers_removed/*trimmed_R2.fastq > ./2.primers_removed/fastqc/allreverse_trimmed.fastq
fastqc ./2.primers_removed/fastqc/allreverse_trimmed.fastq -o ./2.primers_removed/fastqc/
module unload fastqc/0.11.8
#######################################################################

############## VSEARCH PIPELINE ##############

## merging the forward and reverse reads ##
mkdir 3.merged
# provide the path to the folder containing the adpater and primer-removed reads 
data = 2.primers_removed
# provide the path to the vsearch folder
vsearch = ./vsearch-2.22.1/bin/vsearch
# provide the maximum no of mismatches allowed
maxdiffs="15"
# provide the minimum alignmnet length
overlap="130"
#*****************************************************************************************
# Step 1: merge data
for file1 in ${raw_data}/*R1_trimmed.fastq
    do

        echo ""
        echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        echo Merging paired reads
        echo forward reads are:
        echo $(basename ${file1})
        echo reverse reads are:
        echo $(basename ${file1} R1_trimmed.fastq)R2_trimmed.fastq

    $vsearch -fastq_mergepairs ${file1} -reverse "${raw_data}/$(basename -s R1_trimmed.fastq ${file1})R2_trimmed.fastq" -fastqout "3.merged/$(basename "$file1")" -fastq_maxdiffs ${maxdiffs} -fastq_minovlen ${overlap} --fastq_eeout
done

######## checking the quality of merged reads on FastQC #########
mkdir ./3.merged/fastqc
cat ./3.merged/*.fastq > ./3.merged/fastqc/onefile_merged.fastq
module load fastqc/0.11.8
fastqc ./3.merged/fastqc/onefile_merged.fastq -o ./3.merged/fastqc/
module unload fastqc/0.11.8

#*****************************************************************************************
mkdir ./3.merged/stats
mkdir ./4.qual_filtered
mkdir ./4.qual_filtered/qual_filt
mkdir ./5.per_sample_derep
mkdir ./5.per_sample_derep/uc_files

for i in $(cat samplenames.txt)
do
## Step 2: Calculate quality statistics ## 
$vsearch --fastq_eestats ./3.merged/"$i"_R1_trimmed.fastq \
         --output ./3.merged/stats/"$i".stats

## Step 3: Quality filtering ## fastq_maxlen 446 as that is the max length we have based on FastQC ##    
$vsearch --fastq_filter ./3.merged/"$i"_R1_trimmed.fastq \
         --fastq_maxee 1 \
         --fastq_minlen 250 \
         --fastq_maxlen 446 \
         #--fastq_maxns 0 \
         #--fastaout ./4.qual_filtered/qual_filt/"$i"_qual_filtered.fasta \
         #--fasta_width 0
         
## Step 4: Dereplicate at sample level and relabel with sample_n ##
$vsearch --derep_fulllength ./4.qual_filtered/qual_filt/"$i"_qual_filtered.fasta \
         --strand plus \
         --output ./5.per_sample_derep/"$i"_ps_derep.fasta \
         --sizeout \
         --uc ./5.per_sample_derep/uc_files/"$i"_derep.uc \
         --relabel "$i" \
         --fasta_width 0         
done

## see the no. of reads at each step ##
echo From Step 3: Sum of sequences in each sample from fastq_maxee: $(cat ./4.qual_filtered/qual_filt/*.fasta | grep -c "^>")
echo From Step 4: Sum of unique sequences in each sample: $(cat ./5.per_sample_derep/*.fasta | grep -c "^>")

## Step 5: Merge all samples ## 
cat ./5.per_sample_derep/*_ps_derep.fasta > ./5.per_sample_derep/all.fasta

## Step 6: Dereplicate across samples and remove singletons ##
mkdir ./6.all_derep
$vsearch --derep_fulllength ./5.per_sample_derep/all.fasta \
         --minuniquesize 2 \
         --sizein \
         --sizeout \
         --fasta_width 0 \
         --uc ./6.all_derep/all.derep.uc \
         --output ./6.all_derep/all.derep.fasta

## Step 7: Precluster at 98% before chimera detection ##
mkdir ./7.precluster
$vsearch --cluster_size ./6.all_derep/all.derep.fasta \
         --id 0.98 \
         --strand plus \
         --sizein \
         --sizeout \
         --fasta_width 0 \
         --uc ./7.precluster/all.preclustered.uc \
         --centroids ./7.precluster/all.preclustered.fasta
         
## Step 8: De novo chimera detection ##
mkdir ./8.de_novo_chimera_rem
$vsearch --uchime_denovo ./7.precluster/all.preclustered.fasta \
         --sizein \
         --sizeout \
         --nonchimeras ./8.de_novo_chimera_rem/all_denovo_nonchimeras.fasta \
         --fasta_width 0 
         
## Step 9: Reference_based chimera detection ##
mkdir ./9.ref_based_chimera_rem
$vsearch --uchime_ref ./8.de_novo_chimera_rem/all_denovo_nonchimeras.fasta \
         --db ./reference_databases/SILVA_138.1_SSURef_NR99_tax_silva.fasta \
         --sizein \
         --sizeout \
         --fasta_width 0 \
         --nonchimeras ./9.ref_based_chimera_rem/all_ref_SILVA_nonchimeras.fasta 
                
## see the no. of reads at each step ##      
#echo From Step 6: All unique non-singleton sequences: $(grep -c "^>" ./6.all_derep/all.derep.fasta)
#echo From Step 7: Unique sequences after preclustering: $(grep -c "^>" ./7.precluster/all.preclustered.fasta) 
#echo From Step 8: Unique sequences after de novo chimera detection: $(grep -c "^>" ./8.de_novo_chimera_rem/all.denovo.nonchimeras.fasta)
#echo From Step 9: Unique sequences after reference-based chimera detection using SILVA database: $(grep -c "^>" ./9.ref_based_chimera_rem/all_ref_SILVA_nonchimeras.fasta)

## Step 10: Now that we have the unique sequences for all the reads across samples, we're putting those back against the individual reads ##
## perl script (map.pl) was obtained from https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline ##
mkdir ./10.mapping
echo Extract all non-chimeric, non-singleton sequences, dereplicated
perl ./map.pl ./6.all_derep/all.derep.fasta ./7.precluster/all.preclustered.uc ./9.ref_based_chimera_rem/all_ref_SILVA_nonchimeras.fasta > ./10.mapping/all_nonchimeras_SILVA_derep.fasta
echo Extract all non-chimeric, non-singleton sequences in each sample
perl ./map.pl ./5.per_sample_derep/all.fasta ./6.all_derep/all.derep.uc ./10.mapping/all_nonchimeras_SILVA_derep.fasta > ./10.mapping/all_nonchimeras_SILVA.fasta

## see the no. of reads at each step ##   
echo From Step 10: Unique non-chimeric, non-singleton sequences - SILVA: $(grep -c "^>" ./10.mapping/all_nonchimeras_SILVA_derep.fasta)
echo From Step 10: Sum of unique non-chimeric, non-singleton sequences in each sample - SILVA: $(grep -c "^>" ./10.mapping/all_nonchimeras_SILVA.fasta)

## Step 11: OTU clustering, generate OTU table ##
mkdir ./11.OTU_binning
mkdir ./11.OTU_binning/SILVA
mkdir ./11.OTU_binning/SILVA/97_pc
# OTU clustering tried with 97%, 98% and 99% ##
$vsearch --cluster_size ./10.mapping/all_nonchimeras_SILVA.fasta \
    --id 0.97 \  #use 0.97 for 97%, 0.98 for 98%, and 0.99 for 99%
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc ./11.OTU_binning/SILVA/97_pc/all_clustered_SILVA_97_pc.uc \
    --centroids ./11.OTU_binning/SILVA/97_pc/all_otus_SILVA_97_pc.fasta \
    --biomout ./11.OTU_binning/SILVA/97_pc/all_otutab_SILVA_97_pc.biom

## see the no. of OTUs generated ## 
echo From Step11: Number of OTUs - 97pc - SILVA: $(grep -c "^>" ./11.OTU_binning/SILVA/97_pc/all_otus_SILVA_97_pc.fasta)
