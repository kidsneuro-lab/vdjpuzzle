#!/bin/bash
#This script execut the RNA-Seq alingnamet pipeline sequentially, multithreading is implemented within the scripts (Trimmomatic, tophat and cufflinks)
#This script should be executed in parallel on different cells indicating the right directory ($1)

#PBS -W group_list=va1

#PBS -q workq

#PBS -l select=1:ncpus=1:mem=20G,walltime=96:00:00

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
else
	param1=$1
	param2=$2
	param3=$3
	param4=$4
	param5=$5
	param6=$6
	param7=$7
	param8=$8
fi

CELL_PATH=$param1

Q1="${CELL_PATH}/${param2}1.fastq.gz"
Q2="${CELL_PATH}/${param2}2.fastq.gz"
Q3=$param4/VDJ_p3_$param2
Q4=$param4
CHAIN_ARRAY=($param6)
CHAIN_PREFIX_ARRAY=($param7)

echo "P1: $param1
P2: $param2
P3: $param3
P4: $param4
P5: $param5
P6: $param6
P7: $param7
P8: $param8"

rm $Q3/overlapping_reads*
for prefix in "${CHAIN_PREFIX_ARRAY[@]}"
do
	rm $Q3/out1${prefix}.fastq
	rm $Q3/out2${prefix}.fastq
done
mkdir $Q3
mkdir $Q3/out

if [ "$param5" -ge 1 ]; then
	Q1="${CELL_PATH}/PAIRED_${param2}1.fastq.gz"
	Q2="${CELL_PATH}/PAIRED_${param2}2.fastq.gz"
fi


for chain in "${CHAIN_ARRAY[@]}"
do
	$BOWTIE --no-unal -p $param8 -k 1 --np 0 --rdg 1,1 --rfg 1,1 -x $Q4/assembledTCR_genome/${chain} -1 $Q1 -2 $Q2 --al-conc $Q3/reads_${chain}_%.fastq -S $Q4/${chain}.sam
done

	
for chain in "${CHAIN_ARRAY[@]}"
do
	$trinitypath --left $Q3/reads_TCRB_1.fastq --right $Q3/reads_TCRB_2.fastq --seqType fq --max_memory 10G --output $Q3/trinity_out_dir
	mv $Q3/trinity_out_dir/Trinity.fasta $Q3/${chain}.fa
	rm -r $Q3/trinity_out_dir
done

mkdir $param4/summary

echo "Running MIGMAP"
for chain in "${CHAIN_ARRAY[@]}"
do
	$JAVA18 -jar $MIGMAP -S $param3 -R ${chain//C} $Q3/${chain}.fa $param4/summary/${chain//C}_$param2
done
