# VDJPuzzle2 - README #

TCR and BCR reconstruction from scRNA-seq data

## Setup

Run `conda create --name vdjpuzzle --file environment-explicit.txt` to create the conda environment. 

Type `source activate vdjpuzzle` to activate the environment. 

Symlink or copy bin/vdjpuzzle into your path, and run `nohup ./vdjpuzzle Example > LOG.txt &`.

## Execution and parameters

Usage: `./vdjpuzzle directory_name [option]`

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
|--gtf=path\_to\_gtf\|Location of the GTF annotation file|

An additional script to plot gene expression as an heatmap annotated with mutation rates and other phenotype data is provided in scripts/mutation\_gene\_expression\_analysis.R

|parameter|description|
| ------------- |-------------|
|--help|show the help|
|-g file|Gene annotation used for CuffNorm|
|-f file|CuffNorm FPKM matrix|
|-a file|Annotation file for each cell. First column contains cell ID|


## Citation

Manuscript is under review, in the meantime you can cite our biorxiv:
Simone Rizzetto, David NP Koppstein, Jerome Samir, Mandeep Singh, Joanne H Reed, Curtis H Cai, Andrew R Lloyd, Auda A Eltahla, Christopher C Goodnow, and Fabio Luciani. B-cell receptor reconstruction from single-cell RNA-seq with VDJPuzzle. Biorxiv, 2017
