# Using filtlong to QC the reads # Performed on both ONT and PacBio data

## counting the number of reads pre quality filtering
for i in $(cat samplenames.txt)
do
echo "$i"
awk '{s++}END{print s/4}' ./1.data/"$i"_qls.fastq
done
gzip ./1.data/*_qls.fastq

## length filtering for size 3500-5000bp using cutadapt 
module load cutadapt/2.6
mkdir ./3.minimap/1.length_trimmed
for i in ./1.data/*_qls.fastq.gz
do
cutadapt -m 3500 -M 5000 -o ./3.minimap/1.length_trimmed/"$(basename "$i" _qls.fastq.gz)".fastq.gz "$i"
done
#module unload cutadapt/2.6

## quality filtering reads for >Q12 using filtlong
module load filtlong/0.2.0
mkdir ./3.minimap/2.qual_filtered_reads
for i in ./3.minimap/1.length_trimmed/*_qls.fastq.gz
do
filtlong --min_mean_q 92 "$i" | gzip > ./3.minimap/2.qual_filtered_reads/"$(basename "$i" _qls.fastq.gz)".fastq.gz
done 
module unload filtlong/0.2.0

## counting the number of reads post quality filtering
gunzip ./3.minimap/2.qual_filtered_reads/*.fastq.gz
for i in $(cat samplenames.txt)
do
echo "$i"
awk '{s++}END{print s/4}' ./3.minimap/2.qual_filtered_reads/"$i".fastq
done

## nanoplot to check if QC worked
mkdir ./3.minimap/2.qual_filtered_reads/nanoplots
cat ./3.minimap/2.qual_filtered_reads/*.fastq > ./3.minimap/2.qual_filtered_reads/onefile.fastq
module load nanoplot/1.28.2
NanoPlot --fastq ./3.minimap/2.qual_filtered_reads/onefile.fastq -o ./3.minimap/2.qual_filtered_reads/nanoplots
module unload nanoplot/1.28.2

#zip the files again
gzip ./3.minimap/2.qual_filtered_reads/*.fastq
