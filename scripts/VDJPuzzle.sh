#!/bin/bash
#Source code implemented by Simone Rizzetto (UNSW, Australia), for enquiries and documentation please refer to https://github.com/Simo-88/VDJPuzzle
#LATEST VERSION

function printHelp {
	echo "VDJ Puzzle v1.0"
	echo "Usage: ./VDJPuzzle.sh directory_name [option]"
	echo "--qsub\texecutes on the cluster (overrides --CPU flag)"
	echo "--type=(t|b)\tspecifies tcell (default) or bcell analysis"
	echo "--CPU=n\truns with n processes"
	echo "--THR=n\truns bowtie and tophat with n threads (default 8)"
	echo "--species=(human|mouse)\tspecified human (default) or mouse organism"
	echo "--only-statistics\texecutes only summary statistics script"
	echo "--no-statistics\tdo not execute summary statistics script"
	echo "--transcriptomic\tBLAHHH"
	echo "--trim\tTrims reads using Trimmomatic"
	echo "--counts\tBLAHHH"
}

function testArg {
	if [[ "$1" =~ $2 ]]; then
		local result='true'
	else
		local result='false'
	fi
	echo "$result"
}

function waitForProcess {
	sleep_time=120
	status=`qstat -u $USER | grep "p${1}_${2}"`
	while [ -n "$status" ] # while $status is not empty
	do
		sleep $sleep_time
		status=`qstat -u $USER | grep "p${1}_${2}"`
	done
}

#CONFIGURE SCRIPT PATHS
export PATH=~/binf3111/VDJPuzzle/Tools/bowtie2/:$PATH
export PATH=~/binf3111/VDJPuzzle/Tools/bowtie/:$PATH
export PATH=~/binf3111/VDJPuzzle/Tools/samtools/bin/:$PATH
export PATH=~/binf3111/VDJPuzzle/Tools/java1.7/bin/:$PATH
export PATH=~/binf3111/VDJPuzzle/Tools/Blast/usr/bin/:$PATH
export PATH=~/binf3111/VDJPuzzle/Tools/igblastwrapper_linux64/bin/:$PATH
export PATH=~/binf3111/VDJPuzzle/Tools/cufflinks-2.2.1.Linux_x86_64/:$PATH
export MIGMAP=~/binf3111/VDJPuzzle/Tools/migmap-0.9.7/migmap-0.9.7.jar
export trinitypath=~/binf3111/VDJPuzzle/Tools/trinityrnaseq-2.0.6/Trinity
export TOPHAT=~/binf3111/VDJPuzzle/Tools/tophat/2.1.0/tophat2
export BOWTIE=~/binf3111/VDJPuzzle/Tools/bowtie2/bowtie2
export BOWTIEBUILD=~/binf3111/VDJPuzzle/Tools/bowtie2/bowtie2-build
export SAMTOOLS=~/binf3111/VDJPuzzle/Tools/samtools/bin/samtools
export BEDTOOLS=~/binf3111/VDJPuzzle/Tools/bedtools/bin/intersectBed
export JAVA18=~/binf3111/VDJPuzzle/Tools/java1.8/bin/java
export TRIMMOMATIC=~/binf3111/VDJPuzzle/Tools/Trimmomatic-0.36/trimmomatic-0.36.jar
export featureCounts=~/binf3111/VDJPuzzle/Tools/subread-1.5.1-Linux-x86_64/bin/featureCounts

### TRIMMING PARAMETERS ###
export ADAPTERS=$(echo $TRIMMOMATIC | egrep -o ".*\/")adapters/NexteraPE-PE.fa
export LEADING=3 
export TRAILING=3
export WINDOW_LEN=4
export WINDOW_QUAL=15
export MINLEN=50
export AVGQUAL=20

#SET DEFAULT PARAMETERS
PARALLEL=1
NTHREADS=1
COUNTS=0
RUNTCR=1
RUNR=1
ORGANISM=human
TRIM=0
QSUB=0
TRANSCRIPTOMIC=0
TYPE='t'

echo `date`

#PARSING PARAMETERS
for i in "$@"
do
	case $i in
		-h|--help)
		printHelp
		exit 0
		shift
		;;

		--qsub)
		QSUB=1
		shift
		;;

		--type=*)
		TYPE="${i#*=}"
		if [[ "$(testArg $TYPE '^[tb]$')" != "true" ]]; then
			printHelp
			exit
		fi
		shift
		;;

		--CPU=*)
		PARALLEL="${i#*=}"
		if [[ "$(testArg $PARALLEL '^[0-9]+$')" != "true" ]]; then
			printHelp
			exit
		fi
		shift
		;;

		--THR=*)
		NTHREADS="${i#*=}"
		if [[ "$(testArg $NTHREADS '^[0-9]+$')" != "true" ]]; then
			printHelp
			exit
		fi
		shift
		;;

		--only-statistics)
		RUNTCR=0
		shift
		;;

		--species=*)
		ORGANISM="${i#*=}"
		if [[ "$(testArg $ORGANISM '(mouse|human)')" != "true" ]]; then
			printHelp
			exit
		fi
		shift
		;;

		--transcriptomic)
		TRANSCRIPTOMIC=1
		shift
		;;

		--trim)
		TRIM=1
		shift
		;;

		--counts)
		COUNTS=1
		shift
		;;

		--no-statistics)
		RUNR=0
		shift
		;;

		*)
		if [[ ! -d "$i" ]]; then
			echo "Invalid SC_RNA_SEQ directory: $i"
			printHelp
			exit
		fi
		echo "SC_RNA_SEQ Directory $i"
		RNADIR=$i
		shift
		;;
	esac
done

if [[ -z "$RNADIR" ]]; then
	echo "No inputted cell directory!"
	printHelp
	exit
fi

#CONFIGURE REFERENCE PATHS

if [ "$TYPE" == "t" ]; then
	CHAIN_ARRAY=('TCRA' 'TCRB')
	CHAIN_PREFIX_ARRAY=('a' 'b')
	cell_type='tcell_'
else
	CHAIN_ARRAY=('IGH' 'IGK' 'IGL')
	CHAIN_PREFIX_ARRAY=('h' 'k' 'l')
	cell_type='bcell_'
fi

if [[ $TRANSCRIPTOMIC -ge 1 ]]; then
	cell_type=''
fi

export ENSEMBL=~/binf3111/VDJPuzzle/Tools/refGenome/$ORGANISM/Bowtie2Index/${cell_type}genome
export ANNOTATION=~/binf3111/VDJPuzzle/Tools/refGenome/$ORGANISM/Genes/genes.gtf

for chain in "${CHAIN_ARRAY[@]}"
do
	export $chain=~/binf3111/VDJPuzzle/Tools/refGenome/BED_files/${chain//C}_$ORGANISM.bed
done

export CUFFOUTPUT=$PWD/CuffQuant
export COUNTOUTPUT=$PWD/ReadCounts

echo "parallel = $PARALLEL
qsub = $QSUB
threads = $NTHREADS
counts = $COUNTS
organism = $ORGANISM
runtcr = $RUNTCR
runr = $RUNR
trim = $TRIM
type = $TYPE
transcriptomic = $TRANSCRIPTOMIC
chain_array = ${CHAIN_ARRAY[*]}
chain_prefix_array = ${CHAIN_PREFIX_ARRAY[*]}"

#MAKE SURE TO HAVE THE ABSOLUTE PATH WITHOUT / AT THE END
CURDIR="$PWD"
cd $RNADIR
RNADIR="$PWD"
cd $CURDIR

if [ "$TRANSCRIPTOMIC" -ge 1 ]; then
	mkdir  $CUFFOUTPUT
	echo "Transcriptomic quantification"
	echo -e "sample_name\tgroup" > $CUFFOUTPUT/sample_sheet.txt 
fi

#PART1 - GENERATE PRELIMINARY VDJ REPERTOIRE
if [[ $QSUB -ge 1 ]] ; then
	N=1
	for d in $RNADIR/*
	do
		echo $d
		di=${d%%/};
		echo $di
		nfiles=$(find $di/*fastq* -type f | wc -l)
		if [ "$nfiles" -ge 2 ]; then
			filename=$(basename $di)
			qsub -N "p1_${TYPE}_$N" -V -v P1=$di,P2=$filename,P3=$ORGANISM,P4=$PWD,P5=$TRANSCRIPTOMIC,P6=$TRIM,P7="${CHAIN_ARRAY[*]}",P8="${CHAIN_PREFIX_ARRAY[*]}",P9=$NTHREADS run_tcr_method2_sc_part1.sh
		fi
		N=$((N+1))
	done
	waitForProcess 1 $TYPE

else
	N=0
	for d in $RNADIR/*
	do
		echo $d
		di=${d%%/};
		echo $di
		nfiles=$(find $di/*fastq* -type f | wc -l)
		if [ "$nfiles" -ge 2 ]; then
			filename=$(basename $di)
			if (( $N % $PARALLEL == 0 )) ; then
				wait
			fi
			./run_tcr_method2_sc_part1.sh $di $filename $ORGANISM $PWD $TRANSCRIPTOMIC $TRIM "${CHAIN_ARRAY[*]}" "${CHAIN_PREFIX_ARRAY[*]}" $NTHREADS &
		fi
		N=$((N+1))
	done
	wait
fi

if [ "$TRANSCRIPTOMIC" -ge 1 ]; then
	find "$CUFFOUTPUT"/* -type d -exec echo -ne "{}/abundances.cxb\t" \; -exec sh -c "echo {} | sed 's/.*\///'" \; >> $CUFFOUTPUT/sample_sheet.txt
	cuffnorm --use-sample-sheet -o CuffNorm $ANNOTATION $CUFFOUTPUT/sample_sheet.txt
fi

rm -r preliminary
mkdir preliminary
mv VDJ_p1* preliminary/

if [ "$COUNTS" -ge 1 ]; then
	mkdir $COUNTOUTPUT
	find preliminary -name 'accepted_hits.bam' -type f > $COUNTOUTPUT/bam_files.txt
	mapfile -t <$COUNTOUTPUT/bam_files.txt
	$featureCounts --primary -a $ANNOTATION -o "$COUNTOUTPUT/featureCounts.csv" "${MAPFILE[@]}"
	sed 1d "$COUNTOUTPUT/featureCounts.csv" | cut -f1,7- | sed s/Geneid/id/ > "$COUNTOUTPUT/featureCounts_formatted.csv"
fi

echo "Preliminary TCR repertoire reconstructed."


#PART2 - BUILD A NEW REFERENCE GENOME BASED ON THE PRELIMINARY REPERTOIRE
for chain_prefix in "${CHAIN_PREFIX_ARRAY[@]}"
do
	rm all_sequences_${chain_prefix}.fa
done

if [[ $QSUB -ge 1 ]] ; then
	N=1
	for d in $PWD/preliminary/*
	do
		echo $d
		di=${d%%/};
		echo $di
		nfiles=$(find $di/*fastq* -type f | wc -l)
		if [ "$nfiles" -ge 2 ]; then
			filename=$(basename $di)
			qsub -N "p2_${TYPE}_$N" -V -v P1=$di,P2=$ORGANISM,P3=$PWD,P4=$N,P5="${CHAIN_ARRAY[*]}",P6="${CHAIN_PREFIX_ARRAY[*]}" run_tcr_method2_sc_part2.sh
		fi
		N=$((N+1))
	done
	waitForProcess 2 $TYPE

else
	N=0
	for d in $PWD/preliminary/*
	do
		echo $d
		di=${d%%/};
		echo $di
		nfiles=$(find $di/*fastq* -type f | wc -l)
		if [ "$nfiles" -ge 2 ]; then
			filename=$(basename $di)
			if (( $N % $PARALLEL == 0 )) ; then
				wait
			fi
			./run_tcr_method2_sc_part2.sh $di $ORGANISM $PWD $N "${CHAIN_ARRAY[*]}" "${CHAIN_PREFIX_ARRAY[*]}" &
			N=$((N+1))
		fi
	done
	wait
fi


for chain_prefix in "${CHAIN_PREFIX_ARRAY[@]}"
do
	cat matching_reads_${chain_prefix}_* > all_sequences_${chain_prefix}.fa
done

index=0
for chain in "${CHAIN_ARRAY[@]}"
do
	$BOWTIEBUILD -f all_sequences_${CHAIN_PREFIX_ARRAY[$index]}.fa $chain
	index=$((index+1))
done

mkdir assembledTCR_genome
mv *.bt2 assembledTCR_genome/

mkdir summary2
mv summary/* summary2/


#PART3 - 
if [[ $QSUB -ge 1 ]] ; then
	N=1
	for d in $RNADIR/*
	do
		echo $d
		di=${d%%/};
		echo $di
		nfiles=$(find $di/*fastq* -type f | wc -l)
		if [ "$nfiles" -ge 2 ]; then
			filename=$(basename $di)
			qsub -N "p3_${TYPE}_$N" -V -v P1=$di,P2=$filename,P3=$ORGANISM,P4=$PWD,P5=$TRIM,P6="${CHAIN_ARRAY[*]}",P7="${CHAIN_PREFIX_ARRAY[*]}",P8=$NTHREADS run_tcr_method2_sc_part3.sh
			N=$((N+1))
		fi
	done
	waitForProcess 3 $TYPE

else
	N=0
	for d in $RNADIR/*
	do
		echo $d
		di=${d%%/};
		echo $di
		nfiles=$(find $di/*fastq* -type f | wc -l)
		if [ "$nfiles" -ge 2 ]; then
			filename=$(basename $di)
			if (( $N % $PARALLEL == 0 )) ; then
				wait
			fi
			./run_tcr_method2_sc_part3.sh $di $filename $ORGANISM $PWD $TRIM "${CHAIN_ARRAY[*]}" "${CHAIN_PREFIX_ARRAY[*]}" $NTHREADS &
			N=$((N+1))
		fi
	done
	wait
fi

rm -r TCRsequences
mkdir TCRsequences
mv VDJ_p3* TCRsequences/

# UPDATE SUMMARY FOLDER WITH LATEST
echo '#!/bin/bash' > mvTCR.sh
for chain in "${CHAIN_ARRAY[@]}"
do
	find summary2 -name "${chain//C}_*" -exec awk 'END { if (NR > 1) print "mv " FILENAME " summary/" }' {} \; >> mvTCR.sh
done
chmod a+x mvTCR.sh
./mvTCR.sh

if [ "$RUNR" -ge 1 ]; then

	index=0
	for chain in "${CHAIN_ARRAY[@]}"
	do
		cat summary/${chain//C}* > summary/h${index}.txt
		head -1 summary/h${index}.txt > summary/header${index}.txt
		rm summary/${chain//C}tmp.txt
		rm summary/${chain//C}_cells.txt
		rm summary/${chain//C}.csv

		for tcr in summary/${chain//C}*;
		do
			if [ "$(cat $tcr | wc -l)" -ge 2 ]; then
				number="$(cat $tcr | wc -l)"
				number="$(expr $number - 1)"
				for i in `seq $number`; do
					echo ${tcr##*/} >> summary/${chain//C}_cells.txt
				done
				tail -n +2 $tcr >> summary/${chain//C}tmp.txt
			fi
		done

		cat summary/header${index}.txt summary/${chain//C}tmp.txt > summary/single_cells_${chain//C}.txt
		echo "CellID" > summary/single_cells_ID_${chain//C}.txt
		cat summary/${chain//C}_cells.txt >> summary/single_cells_ID_${chain//C}.txt
		paste -d"\t" summary/single_cells_ID_${chain//C}.txt summary/single_cells_${chain//C}.txt > summary/${chain//C}.csv
		rm summary/h${index}.txt
		rm summary/header${index}.txt
		rm summary/${chain//C}tmp.txt
		rm summary/${chain//C}_cells.txt
		rm summary/single_cells_${chain//C}.txt
		rm summary/single_cells_ID_${chain//C}.txt
		index=$((index+1))
	done
fi

# ERROR CORRECTION
# Requires fasta file of reference (reconstructed) sequences and paired end reads
# Will run migmap, take identified sequences, align reads and produce the corrected consensus sequence

export PATH=/apps/bwa/0.7.12/:$PATH
export PATH=/home/fzl561/scRNAseq/Tools/bcftools-1.3.1/bin/:$PATH

for chain in "${CHAIN_ARRAY[@]}"
do
	R1="reads_${chain}_1.fastq"
	R2="reads_${chain}_2.fastq"

	# Set up fasta filename containing all reconstructed sequences
	if [[ $TYPE == 't' ]] ; then
		# tcr_[ab].fa
		ORIG_REF="$(echo $chain | tr '[:upper:]' '[:lower:]' | sed -e 's/^.\{3\}/&_/').fa"
	else
		# ig_[hkl].fa
		ORIG_REF="$(echo $chain | tr '[:upper:]' '[:lower:]' | sed -e 's/^.\{2\}/&_/').fa"
	fi

	# create directory for storing results
	OUTPUT=$CURDIR/"$(echo ${chain//c} | tr '[:upper:]' '[:lower:]')"_err_corr

	echo "CURDIR:\t$CURDIR\nR1:\t$R1\nR2:\t$R2\nORIG_REF:\t$ORIG_REF\nOUTPUT:\t$OUTPUT\n"
	exit

	if [[ -e $OUTPUT ]]
	then
		echo "Output directory already exits, exiting ..."
		exit 1
	fi
	mkdir $OUTPUT

	for cell in $CURDIR/TCRsequences/*
	do 
		# read original data
		ALL_FASTA=$cell$ORIG_REF
		READS_1=$cell$R1
		READS_2=$cell$R2
		
		# skip the case/cell where a directory does not have required files
		if [ -e $ALL_FASTA -a -e $READS_1 -a -e $READS_2 ]
		then
			echo "[current] Reference file and reads found."
			sleep 2
		else
			echo "[current] Reference file and reads were not found ... moving to next cell."
			continue
		fi
		
		# create subdirectory for each cell analysed
		mkdir $OUTPUT/$(basename $cell)
		cd $OUTPUT/$(basename $cell)
		
		# construct path for final output file
		CORR_SEQ=$OUTPUT/$(basename $cell)/new_$ORIG_REF
		rm -vf accepted_seq_names.txt

		# orange.intersect
		if [[ $(hostname) == 'orange' ]] ; then
			$JAVA18 -jar $MIGMAP -S human --by-read -R $CHAIN_TYPE $ALL_FASTA - |
		else
			# sudo required for write permission on cluster
			sudo $JAVA18 -jar $MIGMAP -S human --by-read -R $CHAIN_TYPE $ALL_FASTA - |
		fi

		# first tab delimited entry from migmap is the sequence name
		# create a file with names of identified sequences 
		cut -f1 | while read x ; do
				if [[ "$x" == ">"* ]] ; then 
					echo "$x" >> accepted_seq_names.txt; 
				fi;
			done 

		rm -vf $CORR_SEQ
			
		# go through identified sequences and use as reference for alignment
		while read seq_name; do
			# create a fasta reference file with a single sequence, index it as required
			python $CURDIR/extract_fasta.py "$seq_name" $ALL_FASTA # creates reference.fasta
			samtools faidx reference.fasta
			bwa index reference.fasta
			
			# samtools v1.3.1
			bwa mem reference.fasta $READS_1 $READS_2 | samtools view -b | samtools sort -o alignment.bam
			
			samtools index alignment.bam
			echo "[current] samtools index finished"
			
			samtools mpileup -uf reference.fasta alignment.bam | 
				
			if [[ $(hostname) == 'orange' ]] ; then
				bcftools call -c | python $CURDIR/mpileup2cons.py reference.fasta
			else 
				# cluster 21
				bcftools view -cg - | python $CURDIR/mpileup2cons.py reference.fasta
			fi
			
			# http://samtools.sourceforge.net/mpileup.shtml
			# alternative but produces truncated consensus and ambiguous bases
			# use vcfutils, skip first line which is the sequence ID, convert to a single line, stop after the + character (now have the consensus sequence)
			# samtools mpileup -uf reference.fasta alignment.bam | bcftools view -cg - |
				# vcfutils.pl vcf2fq | tail -n +2 | tr -d '\n' | cut -f1 -d +
			
			echo "[current] writing sequence to file ..."
			cat new_reference.fasta >> $CORR_SEQ
			echo -e "[current] Output located at:\n"$CORR_SEQ
			
			# Even if we don't remove, they will be overwritten when the enxt sequence is processed.
			rm -vf new_reference.fasta
			rm -vf reference.fasta*
			rm -vf alignment.bam* 
			
			echo "[current] sleeping for 10 then moving to next sequence ..." && sleep 10
			
		done <accepted_seq_names.txt
		rm -vf accepted_seq_names.txt
		echo "[current] sleeping for 10 then moving to next cell ..." && sleep 10
	done
done


if [[ $TYPE == 't' ]]; then
	echo "TCR identification complete! Check your results in the summary directory."
else
	echo "Ig identification complete! Check your results in the summary directory."
fi
echo `date`
