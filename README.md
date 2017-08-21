# VDJPuzzle2 - README #

## Setup

Run `conda create --name vdjpuzzle --file environment-explicit.txt` to create the conda environment. 

Type `source activate vdjpuzzle` to activate the environment. 

Symlink or copy bin/vdjpuzzle into your path, and run `nohup ./vdjpuzzle Example > LOG.txt &`.

## Execution and parameters

Usage: `./vdjpuzzle directory_name [option]`

|--help|show the help|
|--qsub|executes on the cluster (overrides --CPU flag)|
|--type=(t|b)|specifies tcell (default or bcell analysis|
|--CPU=n|runs with n processes|
|--THR=n|runs bowtie and tophat with n threads (default 8)|
|--species=(human|mouse)|specified human (default) or mouse organism|
|--no-err-corr|Do not perform final error correction on consensus sequence|
|--only-statistics|Executes only summary statistics script|
|--no-statistics|Do not execute summary statistics script|
|--transcriptomic|Enable cuffquant/cuffnorm gene quantification measurements|
|--trim|Trim reads using Trimmomatic|
|--counts|Enable featureCounts gene quantification measurements|

## Citation

Manuscript is under review, in the meantime you can cite the first version of VDJPuzzle:

Auda Eltahla*, Simone Rizzetto*, Mehdi Rasoli*, Brigid Betz-Stablein, Vanessa Venturi, Katherine Kedzierska, Andrew R Lloyd, Rowena A Bull and Fabio Luciani. Linking the T cell receptor to the single cell transcriptome in antigen-specific human T cells. Immunology and cell biology, 2016, doi: 10.1038/icb.2016.16G 