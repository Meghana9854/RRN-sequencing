# Using filtlong to QC the reads 

# counting the number of reads pre quality filtering #
for i in $(cat samplenames.txt)
do
echo "$i"
awk '{s++}END{print s/4}' ./1.data/"$i"_qls.fastq
done

#gzip ./1.data/*fastq

## filtering out reads that that are greater than Q12 and from 3500 to 5000bp using filtlong## since Q10 is 90% and Q20 is 99% so Q12 should be around 92 ##
## NOTE: since PacBio read quality was already hight than Q12 so we didnt loose any reads for the quality threshold, we only loose reads for those lower than##### 3500bp for the filtlong used for the RRN reads ###
#module load filtlong/0.2.0

#mkdir ./3.minimap/1.qual_filtered_reads

#for i in ./1.data/demux_QualLessThan_80/*_qls.fastq.gz
#do
#filtlong --min_length 3500 --min_mean_q 92 "$i" | gzip > ./3.minimap/1.qual_filtered_reads/"$(basename "$i" _qls.fastq.gz)".fastq.gz
#done 

#module unload filtlong/0.2.0

# counting the number of reads post quality filtering #

#for i in $(cat samplenames.txt)
#do
#echo "$i"
#awk '{s++}END{print s/4}' ./3.minimap/1.qual_filtered_reads/"$i".fastq
#done

## nanoplot those reads ##
#mkdir ./3.minimap/1.qual_filtered_reads/nanoplots

#cat ./3.minimap/1.qual_filtered_reads/*.fastq > ./3.minimap/1.qual_filtered_reads/onefile.fastq

module load nanoplot/1.28.2
NanoPlot --fastq ./3.minimap/1.qual_filtered_reads/onefile.fastq -o ./3.minimap/1.qual_filtered_reads/nanoplots
module unload nanoplot/1.28.2
