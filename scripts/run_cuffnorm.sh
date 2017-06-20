#!/bin/bash

#PBS -P va1

#PBS -q normal

#PBS -l walltime=48:00:00
#PBS -l ncpus=8
#PBS -l mem=128G

#PBS -l wd

#module load cufflinks/2.2.1

#ANNOTATION=/short/va1/fzl561/scRNAseq/refGenome/human/Genes/genes.gtf
#CUFFOUTPUT=`pwd`/CuffQuant

# if PATH_PARAM has been passed to the script, then set PATH
# this is because sometimes PATH gets overwritten on slave nodes, even when using -V
if [ ! -z ${PATH_PARAM+x} ]; then
	export PATH=$PATH_PARAM
fi

cuffnorm -p 8 --use-sample-sheet -o CuffNorm $ANNOTATION $CUFFOUTPUT/sample_sheet.txt
