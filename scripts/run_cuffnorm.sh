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

cuffnorm --use-sample-sheet -o CuffNorm $ANNOTATION $CUFFOUTPUT/sample_sheet.txt

