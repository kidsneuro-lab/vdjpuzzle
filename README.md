## VDJPuzzle 

## Addendum on 28 July 2023

![docker-image.yml](https://github.com/kidsneuro-lab/vdjpuzzle/actions/workflows/docker-image.yml/badge.svg?branch=main)

This is based on the original author's commit [c695ac](https://bitbucket.org/kirbyvisp/vdjpuzzle/commits/c695ac331151cc7a286bbf90c844c8c2948a84ec). It's unsure if the original maintainer is updating this package or not (last update was back in Mar 2020). This repo is purely intended to get `vdjpuzzle` to a working state with a Dockerhub image so that it's easier to run.

**Pull image from Dockerhub**
```
docker pull schnknc/vdjpuzzle:latest
```

**Run vdjpuzzle**
```
docker run --rm schnknc/vdjpuzzle
```

TCR and BCR reconstruction from scRNA-seq data


## Latest version update

* added option to start from bam files
* TCR and BCR expression quantification
* kallisto quantification
* Isotype Identification 
* membrane bound vs secreted isoforms detection

## Setup

Download this repository and move into the VDJPuzzle directory. 

To create the conda environment run `conda env create -f environment.yml`. If you don't have conda installed, you can find it [here](https://conda.io/docs/user-guide/install/index.html). It is advised to use [miniconda3](https://conda.io/miniconda.html)

Type `conda activate vdjpuzzle` to activate the environment (Note: if you have an earlier version of conda, you will need to type `source activate vdjpuzzle` instead).

Symlink or copy bin/vdjpuzzle into your path by copying the following command in the .bashrc file in your home directory and substituting path_to_vdjpuzzle_dir with the absolute VDJPuzzle directory

`export PATH=/path_to_vdjpuzzle_dir/bin:$PATH`

VDJPuzzle requires the [Ensembl reference genome](https://ccb.jhu.edu/software/tophat/igenomes.shtml). The latest human genome reference (v38), with its associated annotation file, has been indexed and can be downloaded [here](https://unsw-my.sharepoint.com/:f:/g/personal/z5168329_ad_unsw_edu_au/EvKq-aVsVDlFq3Imjjr80qsB-ukZHbDIP5F6Xc3YKe-3mg). 

run an example with `nohup vdjpuzzle Example --bowtie-index=path_to_bt2_index/genome --gtf=path_to_gene_annotations.gtf > LOG.txt &` from the VDJPuzzle directory, you can run it on a different directory but make sure that "Example" is pointing to the Example directory in this repository.

This command will take approximaly 30 minutes to complete. You will find the output in the summary_corrected directory.

## Run VDJPuzzle with a different reference genome
VDJPuzzle uses a BED file to locate the position of the VDJ genes in the genome. The BED files provided are built for the Ensembl reference genome.
If you would like to use a different reference genome, you can generate a new BED file using this [Python script](https://bitbucket.org/kirbyvisp/marmo/src/7cfeada825fb9a00d07ebe89a7e8599550b709f1/scripts/extract_receptors.py?at=master&fileviewer=file-view-default).

## Execution and parameters

Usage: `vdjpuzzle rna_seq_directory_name --bowtie-index=path_to_bt2_index/genome_prefix --gtf=path_to_gene_annotations.gtf [option]`

Note that --bowtie-index and --gtf parameters are mandatory. 

rna_seq_directory_name contains the fastq files organized by single cell (i.e. one sub-directory for each cell that include the fastq files from that cell, check the structure of the Example directory). All fastq files need to be zipped e.g. fastq.gz and paired data needs to be specified using \_1 and \_2 in the file name.

|parameter|description|
| ------------- |-------------|
|--help|show the help|
|--qsub|executes on the cluster (overrides --CPU flag)|
|--type=(t⎮b)|specifies tcell (default) or bcell analysis|
|--CPU=n|runs with n processes|
|--THR=n|runs bowtie and tophat with n threads (default 8)|
|--species=(human⎮mouse)|specified human (default) or mouse organism|
|--no-err-corr|Do not perform final error correction on consensus sequence|
|--only-statistics|Executes only summary statistics script|
|--no-statistics|Do not execute summary statistics script|
|--transcriptomic|Enable cuffquant/cuffnorm gene quantification measurements|
|--trim|Trim reads using Trimmomatic|
|--counts|Enable featureCounts gene quantification measurements|
|--bowtie-index=path\_to\_bt2\_index|Location of the bowtie index files including the prefix (e.g. /path/to/bt2/genome)|
|--star-index=path\_to\_star_index| star index directory 
|--align=star\_if star aligner is used instead of tophat|
|--gtf=path\_to\_gtf/|Location of the GTF annotation file|
|--bam=path\_to/bam_files/|Location of bam files. Bam Files should contain the CellIDs of the input RNA-seq data, and be organized in folders with the CellID as the name. E.g. for cells "CellA1" and "CellA2" the contents of the bam_files folder should be ./CellA1/CellA1.bam and ./CellA1/CellA1.bam.|

An additional script to plot gene expression as an heatmap annotated with mutation rates and other phenotype data is provided in scripts/mutation\_gene\_expression\_analysis.R

|parameter|description|
| ------------- |-------------|
|--help|show the help|
|-g file|Gene annotation in the CuffNorm directory|
|-f file|CuffNorm FPKM matrix|
|-a file|Annotation file for each cell. First column contains cell ID|

## Dockerfile file is also available for VDJPuzzle 

## Run VDJPuzzle on a cluster
VDJPuzzle support the execution on a system with PBS scheduler by adding the --qsub option. Every system has different parameters, thus make sure to change these parameters at the beginning of the .sh files in the script directory. 

## VDJPuzzle Output visualization on VDJView 
The output of VDJPuzzle in the final_receptor_results directory(TCR/BCR data) along with the gene expression data in the CuffNorm directory can be uploaded/visualized in VDJView. A complete guideline is provided on the [VDJView page](https://bitbucket.org/kirbyvisp/vdjview/src/master/). VDJView intergates multiple single cell visualization and analysis tools into a single R Shiny App.

## Citation

VDJPuzzle:
Simone Rizzetto, David NP Koppstein, Jerome Samir, Mandeep Singh, Joanne H Reed, Curtis H Cai, Andrew R Lloyd, Auda A Eltahla, Christopher C Goodnow, and Fabio Luciani. B-cell receptor reconstruction from single-cell RNA-seq with VDJPuzzle. Bioinformatics, Volume 34, Issue 16, 15 August 2018, Pages 2846–2847, https://doi.org/10.1093/bioinformatics/bty203

