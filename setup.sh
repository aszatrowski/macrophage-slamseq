#!/bin/bash
export USE_CONDA_CACHE=0
echo "conda cache:" "$USE_CONDA_CACHE"
module load python/miniforge-25.3.0
echo "Loaded python/miniforge-25.3.0"

ENV_NAME="${1:-slamseq_v2}"

if [ -d "/project/lbarreiro/USERS/austin/envs/$ENV_NAME" ]; then
    source activate /project/lbarreiro/USERS/austin/envs/"$ENV_NAME"
    echo "Activated conda env $ENV_NAME."
    echo "Setup complete."
else
    echo "Error: Environment $ENV_NAME not found at /project/lbarreiro/USERS/austin/envs/$ENV_NAME"
    return 1 
fi

# snakemake --executor slurm --default-resources --jobs 1
# conda create -c bioconda --prefix=/project/lbarreiro/USERS/austin/envs/snakemake python=3.11 snakemake snakemake-executor-plugin-slurm
# conda create -c bioconda -c conda-forge --prefix=/project/lbarreiro/USERS/austin/envs/slamseq_v1 python=3.11 snakemake snakemake-executor-plugin-slurm fastqc trimmomatic star hisat-3n kallisto samtools bedtools openjdk
# conda create -c bioconda -c conda-forge --prefix=/project/lbarreiro/USERS/austin/envs/slamseq_fastp python=3.11 snakemake snakemake-executor-plugin-slurm fastqc fastp trimmomatic star hisat-3n kallisto samtools bedtools openjdk
# conda create -c conda-forge -c bioconda --prefix=/project/lbarreiro/USERS/austin/envs/slamseq_v2 snakemake snakemake-executor-plugin-slurm snakemake-storage-plugin-http fastp fastqc hisat-3n samtools bedtools kallisto=0.50.1 openjdk  

# source activate, not conda activate! per https://docs.rcc.uchicago.edu/software/apps-and-envs/python/
# source activate /project/lbarreiro/USERS/austin/envs/snakemake
