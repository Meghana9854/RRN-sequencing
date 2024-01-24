#### I. CHIMERA REMOVAL #####
## first aligning the sequencing data with our mock_DB containing genomes of all the species/representative species in the mock communities, using minimap2 ##

### Step 1 running minimap2 ###
module load minimap2/2.17-r974
mkdir ./3.minimap
mkdir ./3.minimap/2.chimera_removal
mkdir ./3.minimap/2.chimera_removal/1.minimap
for i in ./3.minimap/1.qual_filtered_reads/*.fastq
do
minimap2 -cx map-ont ./mock_DB/mock_DB.fasta "$i" \
         -z 70 > ./3.minimap/2.chimera_removal/1.minimap/"$(basename "$i" .fastq)".paf
done
module unload minimap2/2.17-r974

### Step 2: chimera removal with yacrd ###
# yacrd first detects the chimeras and then filters them out from the length trimmed fastq sequencing files # 
module load yacrd/1.0
mkdir ./3.minimap/chimera_removal/2.yacrd
for i in $(cat samplenames.txt)
do 
yacrd --input ./3.minimap/chimera_removal/1.minimap/"$i".paf \
      --output ./3.minimap/2.chimera_removal/2.yacrd/"$i"_report.yacrd \
      -c 4 -n 0.4 filter --input ./3.minimap/1.qual_filtered_reads/"$i".fastq \
      --output ./3.minimap/2.chimera_removal/2.yacrd/"$i"_chimeras_removed.fastq
done
module unload yacrd/1.0

### Step 3: counting the number of reads per file post chimera removal ###
gunzip ./minimap2/chimera_removal/2.yacrd/*_chimeras_removed.fastq.gz
for i in $(cat samplenames.txt)
do
echo "$i"
awk '{s++}END{print s/4}' ./3.minimap/2.chimera_removal/2.yacrd/"$i"_chimeras_removed.fastq
done
gzip ./minimap2/chimera_removal/2.yacrd/*_chimeras_removed.fastq

###### II. ASSIGNING TAXONOMY #####
### Step 4 taxonomy assignment using minimap2 against FANGORN ###
module load minimap2/2.17-r974
mkdir ./3.minimap/3.tax_assign/

## option 1.i: with FANGORN GTDB NR nrRep ##
mkdir ./3.minimap/3.tax_assign/GTDB
for i in $(cat samplenames.txt)
do
minimap2 -t 32 -cx map-ont ./Fangorn_db/GTDB_NR_207/nrRep.fasta \
         ./3.minimap/2.chimera_removal/2.yacrd/"$i"_chimeras_removed.fastq \
         -z 70 -f1000 > ./3.minimap/3.tax_assign_pb/GTDB/"$i"_GTDB_nrRep.paf
done

## option 1.ii: with FANGORN RefSeq NR nrRep ##
mkdir ./3.minimap/3.tax_assign/RefSeq
for i in $(cat samplenames.txt)
do
minimap2 -t 32 -cx map-ont ./Fangorn_db/RefSeq_NR_207/nrRep.fasta \
         ./3.minimap/2.chimera_removal/2.yacrd/"$i"_chimeras_removed.fastq \
         -z 70 -f1000 > ./3.minimap/3.tax_assign/RefSeq/"$i"_RefSeq_nrRep.paf
done

## option 2: with rrn_DBv2 ##
mkdir ./3.minimap/3.tax_assign/rrn_DBv2
for i in $(cat samplenames.txt)
do
minimap2 -t 32 -cx map-ont ./rrn_databases/rrn_DBv2.fasta \
         ./3.minimap/2.chimera_removal/2.yacrd/"$i"_chimeras_removed.fastq \
         -z 70 -f1000 > ./3.minimap/3.tax_assign/rrn_DBv2/"$i"_rrn_DBv2.paf
done

## option 3: with MIrRoR ##
mkdir ./3.minimap/3.tax_assign/mirror
for i in $(cat samplenames.txt)
do
minimap2 -t 32 -cx map-ont ./rrn_databases/MIrROR_DBDIR/MIrROR_DB_r01.mmi \
         ./3.minimap/2.chimera_removal/2.yacrd/"$i"_chimeras_removed.fastq \
         -z 70 -f1000 > ./3.minimap/3.tax_assign/mirror/"$i"_mirror.paf
done

module unload minimap2/2.17-r974
