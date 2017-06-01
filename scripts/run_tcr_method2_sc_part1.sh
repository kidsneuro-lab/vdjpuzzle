#!/bin/bash
#This script execut the RNA-Seq alingnamet pipeline sequentially, multithreading is implemented within the scripts (Trimmomatic, tophat and cufflinks)
#This script should be executed in parallel on different cells indicating the right directory ($1)

#PBS -W group_list=va1

#PBS -q workq

#PBS -l select=1:ncpus=1:mem=20G,walltime=72:00:00

#parameters
#parameter one needs to be the ABSOLUTE path where cell sequences are located WITHOUT /
if [[ -n "$P3" ]]; then
	param1=$P1 # fastq dir
	param2=$P2 # organism
	param3=$P3 # nth file
	param4=$P4 # pwd
	param5=$P5 # transcriptomic quantification?
	param6=$P6 # trim?
	param7=$P7 # (TCRA, TCRB) or (IGH, IGK, IGL)
	param8=$P8 # (a, b) or (h, k, l)
	param9=$P9 # threads
else
	param1=$1
	param2=$2
	param3=$3
	param4=$4
	param5=$5
	param6=$6
	param7=$7
	param8=$8
	param9=$9
fi

CELL_PATH=$param1

Q1="${CELL_PATH}/${param2}1.fastq.gz"
Q2="${CELL_PATH}/${param2}2.fastq.gz"
Q3=$param4/VDJ_p1_$param2
Q4=$param4
CHAIN_ARRAY=($param7)
CHAIN_PREFIX_ARRAY=($param8)
NTHREADS=($param9)

rm $Q3/overlapping_reads*
for prefix in "${CHAIN_PREFIX_ARRAY[@]}"
do
	rm $Q3/out1${prefix}.fastq
	rm $Q3/out2${prefix}.fastq
done
mkdir $Q3
mkdir $Q3/out

echo "P1: $param1
P2: $param2
P3: $param3
P4: $param4
P5: $param5
P6: $param6
P7: $param7
P8: $param8
P9: $param9"

echo "$ENSEMBL ${CHAIN_ARRAY[*]} ${CHAIN_PREFIX_ARRAY[*]}"

if [ "$param6" -ge 1 ]; then
	PAIR_1="${CELL_PATH}/PAIRED_${param2}1.fastq.gz"
	PAIR_2="${CELL_PATH}/PAIRED_${param2}2.fastq.gz"
	UNPAIR_1="${CELL_PATH}/UNPAIRED_${param2}1.fastq.gz"
	UNPAIR_2="${CELL_PATH}/UNPAIRED_${param2}2.fastq.gz"

	trimmomatic PE -phred33 $Q1 $Q2 $PAIR_1 $UNPAIR_1 $PAIR_2 $UNPAIR_2 ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:$LEADING TRAILING:$TRAILING SLIDINGWINDOW:$WINDOW_LEN:$WINDOW_QUAL MINLEN:$MINLEN > $CELL_PATH/log_trimmometric.txtfi
	tophat -o $Q3/out/tophat_run1 -p $NTHREADS $ENSEMBL $PAIR_1
	tophat -o $Q3/out/tophat_run2 -p $NTHREADS $ENSEMBL $PAIR_2
else
	tophat -o $Q3/out/tophat_run1 -p $NTHREADS $ENSEMBL $Q1
	tophat -o $Q3/out/tophat_run2 -p $NTHREADS $ENSEMBL $Q2
fi

if [ "$param5" -ge 1 ]; then
	cuffquant -o $CUFFOUTPUT/$param2 $ANNOTATION  $Q3/out/tophat_run1/accepted_hits.bam
	cuffquant -o $CUFFOUTPUT/$param2 $ANNOTATION  $Q3/out/tophat_run2/accepted_hits.bam
fi

index=0
for chain in "${CHAIN_ARRAY[@]}"
do
	if [[ "$param5" -ge 1 ]]; then
		bedtools -wa -abam $Q3/out/tophat_run1/accepted_hits.bam -b ${!chain} > $Q3/out/tophat_run1/overlapping_reads.bam
		bedtools -wa -abam $Q3/out/tophat_run2/accepted_hits.bam -b ${!chain} > $Q3/out/tophat_run2/overlapping_reads.bam
		samtools view -h -o $Q3/out/tophat_run1/overlapping_reads.sam $Q3/out/tophat_run1/overlapping_reads.bam
		samtools view -h -o $Q3/out/tophat_run2/overlapping_reads.sam $Q3/out/tophat_run2/overlapping_reads.bam
	else
		samtools view -h -o $Q3/out/tophat_run1/overlapping_reads.sam $Q3/out/tophat_run1/accepted_hits.bam
		samtools view -h -o $Q3/out/tophat_run2/overlapping_reads.sam $Q3/out/tophat_run2/accepted_hits.bam
	fi

	cat $Q3/out/tophat_run1/overlapping_reads.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > $Q3/overlapping_reads_1.fq
	cat $Q3/out/tophat_run2/overlapping_reads.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > $Q3/overlapping_reads_2.fq
	cat $Q3/overlapping_reads_1.fq | awk 'NR%4==1{printf ">%s\n", substr($0,2)}NR%4==2{print}' > $Q3/overlapping_reads_1.fa
	cat $Q3/overlapping_reads_2.fq | awk 'NR%4==1{printf ">%s\n", substr($0,2)}NR%4==2{print}' > $Q3/overlapping_reads_2.fa

	# find overlapping reads
	grep ">" $Q3/overlapping_reads_1.fa $Q3/overlapping_reads_2.fa | sed 's\>\\g' > $Q3/overlapping_readsIDs.txt
	zcat $Q1 | grep -f $Q3/overlapping_readsIDs.txt -A 3 -F > $Q3/out1${CHAIN_PREFIX_ARRAY[$index]}.fastq
	zcat $Q2 | grep -f $Q3/overlapping_readsIDs.txt -A 3 -F > $Q3/out2${CHAIN_PREFIX_ARRAY[$index]}.fastq

	Trinity --left $Q3/out1${CHAIN_PREFIX_ARRAY[$index]}.fastq --right $Q3/out2${CHAIN_PREFIX_ARRAY[$index]}.fastq --seqType fq --max_memory 10G --output $Q3/trinity_out_dir

	mv $Q3/trinity_out_dir/Trinity.fasta $Q3/${chain}.fa
	rm -r $Q3/trinity_out_dir
	rm $Q3/overlapping_reads*

	index=$((index+1))
done

mkdir $param4/summary

for chain in "${CHAIN_ARRAY[@]}"
do
	migmap -S $param3 -R ${chain//C} $Q3/${chain}.fa $param4/summary/${chain//C}_$param2
done

rm $CELL_PATH/merged*
gzip $CELL_PATH/*.fastq
