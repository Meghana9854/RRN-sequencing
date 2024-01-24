#Removing the RRN primer sequneces using cutadapt
module load cutadapt/2.6

ls -d *.fastq | awk -F '.fastq' '{print "cutadapt -g CAGCMGCCGCGGTAA...CGTCGTGAGACAGKTYGG -g CCRAMCTGTCTCACGACG...TTACCGCGGCKGCTG \
-o ./final_trimmed_reads/"$0" /data/Food/analysis/R1150_biotransformation/RRN_mc_ONT/final_trimmed_reads/519F_2428R/"$0" -e0.25 --discard-untrimmed"}' > intermediate_script.sh

sh intermediate_script.sh

# 
