#!/bin/bash
#This script execut the RNA-Seq alingnamet pipeline sequentially, multithreading is implemented within the scripts (Trimmomatic, tophat and cufflinks)
#This script should be executed in parallel on different cells indicating the right directory ($1)

#PBS -P va1

#PBS -q normal

#PBS -l walltime=48:00:00
#PBS -l ncpus=8
#PBS -l mem=20G

#PBS -l wd

set -x -e -o pipefail # echo on, command fails causes script to exit, pipes fail

#parameters
#parameter one needs to be the ABSOLUTE path where cell sequences are located WITHOUT /
if [[ -n "$P3" ]]; then
	param1=$P1
	param2=$P2
	param3=$P3
	param4=$P4
	param5=$P5
	param6=$P6
	param7=$P7
	param8=$P8
	param9=$P9
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

# if PATH_PARAM has been passed to the script, then set PATH
# this is because sometimes PATH gets overwritten on slave nodes, even when using -V
if [ ! -z ${PATH_PARAM+x} ]; then
	export PATH=$PATH_PARAM
fi

CELL_PATH=$param1

FNAME1=`find -L ${CELL_PATH} -name "*1.fastq.gz" | grep -v "PAIRED" | xargs basename` # ${CELL_PATH}/${param2}1.fastq.gz"
FNAME2=`find -L ${CELL_PATH} -name "*2.fastq.gz" | grep -v "PAIRED" | xargs basename` #"${CELL_PATH}/${param2}2.fastq.gz"
Q1=${CELL_PATH}/$FNAME1
Q2=${CELL_PATH}/$FNAME2
Q3=$param4/VDJ_p1_$param2
Q4=$param4
CHAIN_ARRAY=($param7)
CHAIN_PREFIX_ARRAY=($param8)

# export PATH=/apps/bowtie2/2.1.0/:$PATH
# export PATH=/short/va1/fzl561/scRNAseq/Tools/bowtie/bowtie-1.1.2/:$PATH
# export PATH=/apps/java/jdk1.7.0_25/bin/:$PATH
# export PATH=/short/va1/fzl561/scRNAseq/Tools/igblastwrapper_linux64/bin/:$PATH
# module unload samtools/0.1.18
# module load tophat/2.0.7
# module unload samtools/1.3.1
# module load samtools/0.1.18
# module load bowtie2/2.1.0
# module load cufflinks/2.2.1
# module load bedtools/2.26.0

rm -f $Q3/overlapping_reads*
for prefix in "${CHAIN_PREFIX_ARRAY[@]}"
do
	rm -f $Q3/out1${prefix}.fastq
	rm -f $Q3/out2${prefix}.fastq
done
mkdir -p $Q3
mkdir -p $Q3/out

echo "P1: $param1
P2: $param2
P3: $param3
P4: $param4
P5: $param5
P6: $param6
P7: $param7
P8: $param8
P9: $param9"

echo "$TOPHAT $BEDTOOLS $SAMTOOLS $trinitypath $BOWTIE_INDEX ${CHAIN_ARRAY[*]} ${CHAIN_PREFIX_ARRAY[*]}"

if [ "$param6" -ge 1 ]; then
	PAIR_1="${CELL_PATH}/PAIRED_${FNAME1}"
	PAIR_2="${CELL_PATH}/PAIRED_${FNAME2}"
	UNPAIR_1="${CELL_PATH}/UNPAIRED_${FNAME1}"
	UNPAIR_2="${CELL_PATH}/UNPAIRED_${FNAME2}"

	$TRIMMOMATIC PE -phred33 $Q1 $Q2 $PAIR_1 $UNPAIR_1 $PAIR_2 $UNPAIR_2 ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:$LEADING TRAILING:$TRAILING SLIDINGWINDOW:$WINDOW_LEN:$WINDOW_QUAL MINLEN:$MINLEN > $CELL_PATH/log_trimmometric.txtfi
	$TOPHAT -o $Q3/out/tophat_run1 -p $param9 $BOWTIE_INDEX $PAIR_1
	$TOPHAT -o $Q3/out/tophat_run2 -p $param9 $BOWTIE_INDEX $PAIR_2
else
	$TOPHAT -o $Q3/out/tophat_run1 -p $param9 $BOWTIE_INDEX $Q1
	$TOPHAT -o $Q3/out/tophat_run2 -p $param9 $BOWTIE_INDEX $Q2
fi

if [ "$param5" -ge 1 ]; then
	cuffquant -o $CUFFOUTPUT/$param2 $ANNOTATION  $Q3/out/tophat_run1/accepted_hits.bam
	cuffquant -o $CUFFOUTPUT/$param2 $ANNOTATION  $Q3/out/tophat_run2/accepted_hits.bam
fi

# module unload samtools/0.1.18
# module load samtools/1.2

index=0
for chain in "${CHAIN_ARRAY[@]}"
do
	if [[ "param5" -ge 1 ]]; then
		$BEDTOOLS -wa -abam $Q3/out/tophat_run1/accepted_hits.bam -b ${!chain} > $Q3/out/tophat_run1/overlapping_reads.bam
		$BEDTOOLS -wa -abam $Q3/out/tophat_run2/accepted_hits.bam -b ${!chain} > $Q3/out/tophat_run2/overlapping_reads.bam
		$SAMTOOLS view -h -o $Q3/out/tophat_run1/overlapping_reads.sam $Q3/out/tophat_run1/overlapping_reads.bam
		$SAMTOOLS view -h -o $Q3/out/tophat_run2/overlapping_reads.sam $Q3/out/tophat_run2/overlapping_reads.bam
	else
		$SAMTOOLS view -h -o $Q3/out/tophat_run1/overlapping_reads.sam $Q3/out/tophat_run1/accepted_hits.bam
		$SAMTOOLS view -h -o $Q3/out/tophat_run2/overlapping_reads.sam $Q3/out/tophat_run2/accepted_hits.bam
	fi

	cat $Q3/out/tophat_run1/overlapping_reads.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > $Q3/overlapping_reads_1.fq
	cat $Q3/out/tophat_run2/overlapping_reads.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > $Q3/overlapping_reads_2.fq
	cat $Q3/overlapping_reads_1.fq | awk 'NR%4==1{printf ">%s\n", substr($0,2)}NR%4==2{print}' > $Q3/overlapping_reads_1.fa
	cat $Q3/overlapping_reads_2.fq | awk 'NR%4==1{printf ">%s\n", substr($0,2)}NR%4==2{print}' > $Q3/overlapping_reads_2.fa

	grep ">" $Q3/overlapping_reads_1.fa > $Q3/overlapping_readsID_1.fa
	grep ">" $Q3/overlapping_reads_2.fa > $Q3/overlapping_readsID_2.fa
	sed 's\>\\g' $Q3/overlapping_readsID_1.fa > $Q3/overlapping_readsID3.txt
	sed 's\>\\g' $Q3/overlapping_readsID_2.fa >> $Q3/overlapping_readsID3.txt

	# find fastq entries containing overlapping read IDs from either raw or trimmed fastq files
	# get rid of pesky -- lines which appear for some reason using "^\-\-$"
	if [ "$param6" -ge 1 ]; then # we are using trimmed reads
		zcat $PAIR_1 | grep -f $Q3/overlapping_readsID3.txt -A 3 -F | egrep -v "^\-\-$" > $Q3/out1${CHAIN_PREFIX_ARRAY[$index]}.fastq
		zcat $PAIR_2 | grep -f $Q3/overlapping_readsID3.txt -A 3 -F | egrep -v "^\-\-$" > $Q3/out2${CHAIN_PREFIX_ARRAY[$index]}.fastq
	else
		zcat $Q1 | grep -f $Q3/overlapping_readsID3.txt -A 3 -F > $Q3/out1${CHAIN_PREFIX_ARRAY[$index]}.fastq
		zcat $Q2 | grep -f $Q3/overlapping_readsID3.txt -A 3 -F > $Q3/out2${CHAIN_PREFIX_ARRAY[$index]}.fastq
	fi

	# rebuild the trinity index
	$trinitypath --left $Q3/out1${CHAIN_PREFIX_ARRAY[$index]}.fastq --right $Q3/out2${CHAIN_PREFIX_ARRAY[$index]}.fastq --seqType fq --max_memory 10G --output $Q3/trinity_out_dir

	mv $Q3/trinity_out_dir/Trinity.fasta $Q3/${chain}.fa
	rm -rf $Q3/trinity_out_dir
	rm -f $Q3/overlapping_reads*

	index=$((index+1))
done

# module unload java/jdk1.7.0_25
# module load java/jdk1.8.0_60
# module load blast/2.2.28+

mkdir -p $param4/summary

for chain in "${CHAIN_ARRAY[@]}"
do
	$MIGMAP -S $param3 -R ${chain//C} $Q3/${chain}.fa $param4/summary/${chain//C}_$param2
done

rm -f $CELL_PATH/merged*
gzip $CELL_PATH/*.fastq
