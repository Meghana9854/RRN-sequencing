# Demultiplexing perfomed only for ONT data using nanoplexer 

# First demultiplex based on the M13 barcoders (outer sequences)
mkdir demultiplexed_files
~/nanoplexer -b ./demultiplex_nanoplexer/barcode_M13.fa \           #provide list of M13 barcode sequences in .fa file 
-d ./demultiplex_nanoplexer/dual_barcode_pair_13.txt \              #provide dual barcode combinations in .txt file 
-p ./demultiplexed_files/ ./gridion_fastq_pass_unzipped/*.fastq     #provide output file nad input file destinations

# Then demultiplexing based on the RRN primers (inner sequences)
~/nanoplexer -b ./demultiplex_nanoplexer/barcode_primers.fa \
-d ./demultiplex_nanoplexer/dual_barcode_pair_primers.txt \
-p ./demultiplex_nanoplexer/demultiplexed_files/primers_demultiplexed/ ./demultiplex_nanoplexer/demultiplexed_files/*.fastq
