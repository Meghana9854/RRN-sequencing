## Generating OTUs for the v3-v4 data ## Performed only on Illumina data ##

############## CUTADAPT FOR ADAPTER AND PRIMER TRIMMING ##############
##### Removing primers from the raw data first based on how quality is better when primers are removed before merging (ref: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3187-5)## 
## using paired end mode with cutadapt ## reverse complements of For and Rev primer are not given because they are automatically checked for in paired end mode with cutadapt when two input files are given ##

#cd /data/Food/analysis/R1150_biotransformation/RRN_mc_Illumina/1.i.subsampling/

#mkdir ./2.primers_removed 

#module load cutadapt/2.6 

#adapter_fw='TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
#adapter_rev='GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
#PrimerF_GI='CCTACGGGNGGCWGCAG'
#PrimerR_GI='GACTACHVGGGTATCTAATCC'

#for i in $(cat subsamplesnames.txt)
#do
#cutadapt -a $adapter_fw -A $adapter_rev -b $PrimerF_GI -B $PrimerR_GI \
#-o ./2.primers_removed/"$i"_trimmed.fastq -p ./2.primers_removed/"$i"_trimmed.fastq ./1.reads/"$i".fastq ./1.reads/"$i".fastq \
#-e 0.1 \
#--discard-untrimmed \
#--cores 4
#done

## NOTE: cutadapt -a is for 3'linked adapters/sequences, -b is for 5' or 3' adapters/sequences, and -g is for 5' linked adapters/sequences 

#module unload cutadapt/2.6 

######## checking the trimmed reads on fastqc #########
## for all reads ##
#mkdir ./2.primers_removed/fastqc
#cat ./2.primers_removed/*.fastq > ./2.primers_removed/fastqc/onefile_trimmed.fastq

#module load fastqc/0.11.8
#fastqc ./2.primers_removed/fastqc/onefile_trimmed.fastq -o ./2.primers_removed/fastqc/

## for forward and revrse reads separately ##
#cat ./2.primers_removed/*trimmed_R1.fastq > ./2.primers_removed/fastqc/allforward_trimmed.fastq
#fastqc ./2.primers_removed/fastqc/allforward_trimmed.fastq -o ./2.primers_removed/fastqc/

#cat ./2.primers_removed/*trimmed_R2.fastq > ./2.primers_removed/fastqc/allreverse_trimmed.fastq
#fastqc ./2.primers_removed/fastqc/allreverse_trimmed.fastq -o ./2.primers_removed/fastqc/

#module unload fastqc/0.11.8
