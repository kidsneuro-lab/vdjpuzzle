
mkdir $PWD/$kallisto
kallisto index -i chain.index  FASTA.fa
kallisto quant -i  chain.index -o TRA_out FASTQ1.fa.gz fastaQ2.fa.gz


