# README #

## Setup

Run `conda create --name vdjpuzzle --file environment-explicit.txt` to create the conda environment. 

Type `source activate vdjpuzzle` to activate the environment. 

Symlink or copy bin/vdjpuzzle into your path, and run `nohup ./vdjpuzzle Example > LOG.txt &`.

## Execution and parameters

Usage: `./vdjpuzzle directory_name [option]`

<table><tr><td>--help</td><td>show the help</td></tr>
<tr><td>--qsub</td><td>executes on the cluster (overrides --CPU flag)</td></tr>
<tr><td>--type=(t|b)</td><td>specifies tcell (default or bcell analysis</td></tr>
<tr><td>--CPU=n</td><td>runs with n processes</td></tr>
<tr><td>--THR=n</td><td>runs bowtie and tophat with n threads (default 8)</td></tr>
<tr><td>--species=(human|mouse)</td><td>specified human (default) or mouse organism</td></tr>
<tr><td>--no-err-corr</td><td>Do not perform final error correction on consensus sequence</td></tr>
<tr><td>--only-statistics</td><td>Executes only summary statistics script</td></tr>
<tr><td>--no-statistics</td><td>Do not execute summary statistics script</td></tr>
<tr><td>--transcriptomic</td><td>Enable cuffquant/cuffnorm gene quantification measurements</td></tr>
<tr><td>--trim</td><td>Trim reads using Trimmomatic</td></tr>
<tr><td>--counts</td><td>Enable featureCounts gene quantification measurements</td></tr>
</table>
## Citation

Manuscript is under review, in the meantime you can cite the first version of VDJPuzzle:<br>
Auda Eltahla*, Simone Rizzetto*, Mehdi Rasoli*, Brigid Betz-Stablein, Vanessa Venturi, Katherine Kedzierska, Andrew R Lloyd, Rowena A Bull and Fabio Luciani. Linking the T cell receptor to the single cell transcriptome in antigen-specific human T cells. Immunology and cell biology, 2016, doi: 10.1038/icb.2016.16G 