#!/bin/bash
#This script execut the RNA-Seq alingnamet pipeline sequentially, multithreading is implemented within the scripts (tophat and cufflinks)
#This script should be executed in parallel on different cells indicating the right directory ($1)


#PBS -W group_list=va1

#PBS -q workq

#PBS -l select=1:ncpus=1:mem=20G,walltime=24:00:00

#parameters
#parameter one needs to be the ABSOLUTE path where cell sequences are located WITHOUT /
if [[ -n "$P3" ]]; then
	param1=$P1
	param2=$P2
	param3=$P3
	param4=$P4
	param5=$P5
	param6=$P6
else
	param1=$1
	param2=$2
	param3=$3
	param4=$4
	param5=$5
	param6=$6
fi

CELL_PATH=$param1
CHAIN_ARRAY=($param5)
CHAIN_PREFIX_ARRAY=($param6)

#unzip all files in a dir

echo "P1: $param1
P2: $param2
P3: $param3
P4: $param4
P5: $param5
P6: $param6"

index=0
for chain in "${CHAIN_ARRAY[@]}"
do
	$JAVA18 -jar $MIGMAP -S $param2 -R ${chain//C} --by-read $CELL_PATH/${chain}.fa $CELL_PATH/${chain}.out
	tail -n+2 $CELL_PATH/${chain}.out > $CELL_PATH/${chain}.tmp
	cut -f1 $CELL_PATH/${chain}.tmp -d " " > $CELL_PATH/reads_${CHAIN_PREFIX_ARRAY[$index]}.txt
	cut -c 2- $CELL_PATH/reads_${CHAIN_PREFIX_ARRAY[$index]}.txt | xargs -n 1 samtools faidx $CELL_PATH/${chain}.fa > $param3/matching_reads_${CHAIN_PREFIX_ARRAY[$index]}_$param4.fa

	index=$((index+1))
done
