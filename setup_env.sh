conda create -n ctyper-env python=3.10
conda activate ctyper-env
conda install -y -c bioconda -c conda-forge htslib=1.21 eigen=3.4 snakemake=6.15.1 numpy pandas samtools blast bedtools minimap2 winnowmap gcc

make

