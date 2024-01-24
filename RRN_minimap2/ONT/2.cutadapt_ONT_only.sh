# Removing the RRN primer sequneces using cutadapt
module load cutadapt/2.6

ls -d *.fastq | awk -F '.fastq' '{print "cutadapt -g CAGCMGCCGCGGTAA...CGTCGTGAGACAGKTYGG -g CCRAMCTGTCTCACGACG...TTACCGCGGCKGCTG \
-o ./final_trimmed_reads/"$0" /data/Food/analysis/R1150_biotransformation/RRN_mc_ONT/final_trimmed_reads/519F_2428R/"$0" -e0.25 --discard-untrimmed"}' > intermediate_script.sh

sh intermediate_script.sh
module unload cutadapt/2.6

# Checking the quality of reads (per primer pair) post demultiplexing and cutadapt using Nanoplot 
module load nanoplot/1.28.2

cat *27Fv2_2241R.fastq > 2S.fastq
NanoPlot --fastq 2S.fastq -o ./nanoplot3/2S

cat *27Fv2_2428R.fastq > 2L.fastq
NanoPlot --fastq 2L.fastq -o ./nanoplot3/2L

cat *519F_2241R.fastq > 5S.fastq
NanoPlot --fastq 5S.fastq -o ./nanoplot3/5S

cat *519F_2428R.fastq > 5L.fastq
NanoPlot --fastq 5L.fastq -o ./nanoplot3/5L

module unload nanoplot/1.28.2
