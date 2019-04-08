#!/bin/bash
#This script execut the RNA-Seq alingnamet pipeline sequentially, multithreading is implemented within the scripts (Trimmomatic, tophat and cufflinks)
#This script should be executed in parallel on different cells indicating the right directory ($1)

#PBS -P va1

#PBS -q normal

#PBS -l walltime=24:00:00
#PBS -l ncpus=8
#PBS -l mem=10G

#PBS -l wd

set -x # echo on, command fails causes script to exit, pipes fail

#parameters
#parameter one needs to be the ABSOLUTE path where cell sequences are located WITHOUT /
if [[ -n "$P3" ]]; then
	CELL_PATH=$P1
	CELL_ID=$P2 #filename aka cellID
	ORGANISM=$P3
	WORKING_DIR=$P4
	TRANSCRIPTOMIC=$P5
	TRIM=$P6
	CHAIN_ARRAY=($P7)
	CHAIN_PREFIX_ARRAY=($P8)
	NTHREADS=$P9
	USE_MIGMAP_DETAILS=$P10
	TRIMG_FLAG=$P11
	SRA_FLAG=$P12
	INPUTBAM=$P13 # it is a [1,0] paramter indicating if bam files are used as input
	BAMFILES_PATH=$P14 # indicates bam files location
	CELL_NUMBER=$P15
	
else
	CELL_PATH=$1
	CELL_ID=$2
	ORGANISM=$3
	WORKING_DIR=$4
	TRANSCRIPTOMIC=$5
	TRIM=$6
	CHAIN_ARRAY=($7)
	CHAIN_PREFIX_ARRAY=($8)
	NTHREADS=$9
	USE_MIGMAP_DETAILS=${10}
	TRIMG_FLAG=${11}
	SRA_FLAG=${12}
	INPUTBAM=${13}
	BAMFILES_PATH=${14}
	CELL_NUMBER=${15}
		
	
fi

# if PATH_PARAM has been passed to the script, then set PATH
# this is because sometimes PATH gets overwritten on slave nodes, even when using -V
if [ ! -z ${PATH_PARAM+x} ]; then
	export PATH=$PATH_PARAM
fi

FNAME1=`find -L ${CELL_PATH} -name "*fastq.gz" | egrep ".+_(R1_001|R1|1)\.fastq\.gz" | grep -v "PAIRED" | xargs basename` #"${CELL_PATH}/${CELL_ID}1.fastq.gz"
FNAME2=`find -L ${CELL_PATH} -name "*fastq.gz" | egrep ".+_(R2_001|R2|2)\.fastq\.gz" | grep -v "PAIRED" | xargs basename` #"${CELL_PATH}/${CELL_ID}2.fastq.gz"
Q1=${CELL_PATH}/$FNAME1
Q2=${CELL_PATH}/$FNAME2
Q3=$WORKING_DIR/VDJ_p1_$CELL_ID

rm -f $Q3/overlapping_reads*
for prefix in "${CHAIN_PREFIX_ARRAY[@]}"
do
	rm -f $Q3/out1${prefix}.fastq
	rm -f $Q3/out2${prefix}.fastq
done
mkdir -p $Q3
mkdir -p $Q3/out
if [ "$ALIGNER" == "star" ]; then
	mkdir -p $Q3/out/star_both
fi

echo "P1: $CELL_PATH
P2: $CELL_ID
P3: $ORGANISM
P4: $WORKING_DIR
P5: $TRANSCRIPTOMIC
P6: $TRIM
P7: $CHAIN_ARRAY
P8: $CHAIN_PREFIX_ARRAY
P9: $NTHREADS
P10: $USE_MIGMAP_DETAILS
P11: $TRIMG_FLAG
P12: $SRA_FLAG
P13: $INPUTBAM
P14: $BAMFILES_PATH
P15: $CELL_NUMBER"

echo "$BOWTIE_INDEX $FASTA $ALIGNER $TOPHAT $BEDTOOLS $SAMTOOLS $trinitypath ${CHAIN_ARRAY[*]} ${CHAIN_PREFIX_ARRAY[*]}"

#
if [ "$INPUTBAM" -ge 1 ]; then
	#alignment is already performed before VDJPuzzle and input files are provided
	echo "fiding bam files"
	BAMFILE=$(find $BAMFILES_PATH -name '*.bam' -type f | grep $CELL_ID) #find all aligments and grep the one corresponding to the cell
	echo "$BAMFILE"

elif [ "$TRIM" -ge 1 ]; then
        PAIR_1="${CELL_PATH}/PAIRED_${FNAME1}"
        PAIR_2="${CELL_PATH}/PAIRED_${FNAME2}"
        UNPAIR_1="${CELL_PATH}/UNPAIRED_${FNAME1}"
        UNPAIR_2="${CELL_PATH}/UNPAIRED_${FNAME2}"

        trimmomatic PE -phred33 $Q1 $Q2 $PAIR_1 $UNPAIR_1 $PAIR_2 $UNPAIR_2 ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:$LEADING TRAILING:$TRAILING SLIDINGWINDOW:$WINDOW_LEN:$WINDOW_QUAL MINLEN:$MINLEN > $CELL_PATH/log_trimmometric.txtfi

	if [[ "$ALIGNER" == "star" ]]; then
		#/data/STAR-master/source/STAR --runThreadN $NTHREADS --runMode genomeGenerate --genomeDir $STAR_INDEX --genomeFastaFiles  $FASTA --sjdbGTFfile $ANNOTATION  --sjdbOverhang 99
		cd $Q3/out/star_both 
		STAR --genomeDir $STAR_INDEX --runThreadN $NTHREADS --readFilesCommand zcat --readFilesIn $PAIR_1 $PAIR_2 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterType BySJout  --outSAMmode Full --outSAMunmapped Within --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbGTFfile $ANNOTATION --sjdbGTFfeatureExon exon --quantMode TranscriptomeSAM
		cd $WORKING_DIR
		BAMFILE=$Q3/out/star_both/Aligned.sortedByCoord.out.bam
		
	else
		tophat -o $Q3/out/tophat_both -p $NTHREADS $BOWTIE_INDEX $PAIR_1 $PAIR_2
		BAMFILE=$Q3/out/tophat_both/accepted_hits.bam
		#rm $Q3/out/tophat_both/unmapped.bam
	fi

elif [ "$TRIMG_FLAG" -ge 1 ]; then
        echo "Trimming with trim-galore"
        filename1="${FNAME1%.*}"
        filename1="${filename1%.*}"
        filename2="${FNAME2%.*}"
        filename2="${filename2%.*}"
        echo $filename
        PAIR_1="${CELL_PATH}/${filename1}_val_1.fq.gz"
        PAIR_2="${CELL_PATH}/${filename2}_val_2.fq.gz"
        echo $PAIR_1
        trim_galore --paired -o "${CELL_PATH}" $Q1 $Q2 #this is hardcoded
	if [[ "$ALIGNER" == "star" ]]; then
		#/data/STAR-master/source/STAR --runThreadN $NTHREADS --runMode genomeGenerate --genomeDir  $WORKING_DIR/star_index --genomeFastaFiles  $FASTA --sjdbGTFfile $ANNOTATION  --sjdbOverhang 99
		cd $Q3/out/star_both
		STAR --genomeDir $STAR_INDEX --runThreadN $NTHREADS --readFilesCommand zcat --readFilesIn $PAIR_1 $PAIR_2 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterType BySJout  --outSAMmode Full --outSAMunmapped Within --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbGTFfile $ANNOTATION --sjdbGTFfeatureExon exon --quantMode TranscriptomeSAM 
		cd $WORKING_DIR
		BAMFILE=$Q3/out/star_both/Aligned.sortedByCoord.out.bam
	else
		tophat -o $Q3/out/tophat_both -p $NTHREADS $BOWTIE_INDEX $PAIR_1 $PAIR_2
		BAMFILE=$Q3/out/tophat_both/accepted_hits.bam
		#rm $Q3/out/tophat_both/unmapped.bam
	fi
else
        if [[ "$ALIGNER" == "star" ]]; then
		#/data/STAR-master/source/STAR --runThreadN $NTHREADS --runMode genomeGenerate --genomeDir  $WORKING_DIR/star_index --genomeFastaFiles  $FASTA --sjdbGTFfile $ANNOTATION  --sjdbOverhang 99	
		cd $Q3/out/star_both
		STAR --genomeDir $STAR_INDEX --runThreadN $NTHREADS --readFilesCommand zcat --readFilesIn $Q1 $Q2 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterType BySJout  --outSAMmode Full --outSAMunmapped Within --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbGTFfile $ANNOTATION --sjdbGTFfeatureExon exon --quantMode TranscriptomeSAM
		cd $WORKING_DIR
		BAMFILE=$Q3/out/star_both/Aligned.sortedByCoord.out.bam
	else
		tophat -o $Q3/out/tophat_both -p $NTHREADS $BOWTIE_INDEX $Q1 $Q2
		BAMFILE=$Q3/out/tophat_both/accepted_hits.bam
	#rm $Q3/out/tophat_both/unmapped.bam
	fi
fi

if [ "$TRANSCRIPTOMIC" -ge 1 ]; then
	if [ "$ALIGNER" == "star" ]; then
		cuffquant -o $CUFFOUTPUT/$CELL_ID --library-type fr-firststrand $ANNOTATION $BAMFILE
	else
		cuffquant -o $CUFFOUTPUT/$CELL_ID $ANNOTATION $BAMFILE
    fi
fi


index=0
for chain in "${CHAIN_ARRAY[@]}"
do
	echo $chain
	echo ${!chain}
	if [ "$ALIGNER" == "star" ]; then
		intersectBed -wa -abam $BAMFILE -b ${!chain} > $Q3/overlapping_reads.bam
		samtools view -h $Q3/overlapping_reads.bam | grep -av "^@" | awk '{print "@"$1"\n"$10"\n+\n"$11}' > $Q3/overlapping_reads.fq
		#cat $Q3/overlapping_reads.fq | awk 'NR%4==1{printf ">%s\n", substr($0,2)}NR%4==2{print}' > $Q3/overlapping_reads.fa
	else
		intersectBed -wa -abam $BAMFILE -b ${!chain} > $Q3/overlapping_reads.bam
        samtools view -h $Q3/overlapping_reads.bam | grep -av "^@" | awk '{print "@"$1"\n"$10"\n+\n"$11}' > $Q3/overlapping_reads.fq
	fi

        cat $Q3/overlapping_reads.fq | awk 'NR%4==1{printf ">%s\n", substr($0,2)}NR%4==2{print}' > $Q3/overlapping_reads.fa	
	
	if [[ $SRA_FLAG -ge 1 ]] ; then	
		# Because SRA has basic read ids, we need to use the whole read id to reduce the chance of that read id
		# occuring in the quality score of a read.
		grep -a ">" $Q3/overlapping_reads.fa | sed 's\>\@\g' | sed 's\$\/1\' > $Q3/overlapping_readsID.txt
		grep -a ">" $Q3/overlapping_reads.fa | sed 's\>\@\g' | sed 's\$\/2\' >> $Q3/overlapping_readsID.txt
		grep_x_param='-x'
	else
		grep -a ">" $Q3/overlapping_reads.fa | sed 's\>\\g' > $Q3/overlapping_readsID.txt
		grep_x_param=''
	fi

	# find fastq entries containing overlapping read IDs from either raw or trimmed fastq files
	# get rid of pesky -- lines which appear for some reason using "^\-\-$"
	if [ "$TRIM" -ge 1 ]; then # we are using trimmed reads
		zcat $PAIR_1 | grep -af $Q3/overlapping_readsID.txt -A 3 -F $grep_x_param | egrep -av "^\-\-$" > $Q3/out1${CHAIN_PREFIX_ARRAY[$index]}.fastq
		zcat $PAIR_2 | grep -af $Q3/overlapping_readsID.txt -A 3 -F $grep_x_param | egrep -av "^\-\-$" > $Q3/out2${CHAIN_PREFIX_ARRAY[$index]}.fastq
	else
		zcat $Q1 | grep -af $Q3/overlapping_readsID.txt -A 3 -F $grep_x_param | egrep -av "^\-\-$" > $Q3/out1${CHAIN_PREFIX_ARRAY[$index]}.fastq
		zcat $Q2 | grep -af $Q3/overlapping_readsID.txt -A 3 -F $grep_x_param | egrep -av "^\-\-$" > $Q3/out2${CHAIN_PREFIX_ARRAY[$index]}.fastq
	fi

	# rebuild the trinity index
	Trinity --left $Q3/out1${CHAIN_PREFIX_ARRAY[$index]}.fastq --right $Q3/out2${CHAIN_PREFIX_ARRAY[$index]}.fastq --seqType fq --max_memory 10G --output $Q3/trinity_out_dir

	mv $Q3/trinity_out_dir/Trinity.fasta $Q3/${chain}.fa
	rm -rf $Q3/trinity_out_dir

	index=$((index+1))
done

mkdir -p $WORKING_DIR/summary

for chain in "${CHAIN_ARRAY[@]}"
do
	if [[ $USE_MIGMAP_DETAILS -ge 1 ]] ; then
		migmap -S $ORGANISM -R ${chain//C} --data-dir=$CONDA_PREFIX/share/igblast --details fr1nt,cdr1nt,fr2nt,cdr2nt,fr3nt,fr4nt $Q3/${chain}.fa $WORKING_DIR/summary/${chain//C}_$CELL_ID
	else
		migmap -S $ORGANISM -R ${chain//C} --data-dir=$CONDA_PREFIX/share/igblast $Q3/${chain}.fa $WORKING_DIR/summary/${chain//C}_$CELL_ID
	fi
done

rm -f $CELL_PATH/merged*
gzip $CELL_PATH/*.fastq


# PART 2
index=0
for chain in "${CHAIN_ARRAY[@]}"
do
	migmap -S $ORGANISM -R ${chain//C} --by-read --data-dir=$CONDA_PREFIX/share/igblast $Q3/${chain}.fa $Q3/${chain}.out
	tail -n+2 $Q3/${chain}.out > $Q3/${chain}.tmp
	cut -f1 $Q3/${chain}.tmp -d " " > $Q3/reads_${CHAIN_PREFIX_ARRAY[$index]}.txt
	cut -c 2- $Q3/reads_${CHAIN_PREFIX_ARRAY[$index]}.txt | xargs -n 1 $SAMTOOLS faidx $Q3/${chain}.fa > $WORKING_DIR/matching_reads_${CHAIN_PREFIX_ARRAY[$index]}_$CELL_NUMBER.fa

	index=$((index+1))
done
