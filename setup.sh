export USE_CONDA_CACHE=0
echo "$USE_CONDA_CACHE"
module load python/miniforge-25.3.0
source activate /project/lbarreiro/USERS/austin/envs/slamseq_v1
source activate /project/lbarreiro/USERS/austin/envs/slamseq_fastp

snakemake --executor slurm --default-resources --jobs 1
# conda create -c bioconda --prefix=/project/lbarreiro/USERS/austin/envs/snakemake python=3.11 snakemake snakemake-executor-plugin-slurm
# conda create -c bioconda -c conda-forge --prefix=/project/lbarreiro/USERS/austin/envs/slamseq_v1 python=3.11 snakemake snakemake-executor-plugin-slurm fastqc trimmomatic star hisat-3n kallisto samtools bedtools openjdk
conda create -c bioconda -c conda-forge --prefix=/project/lbarreiro/USERS/austin/envs/slamseq_fastp python=3.11 snakemake snakemake-executor-plugin-slurm fastqc fastp trimmomatic star hisat-3n kallisto samtools bedtools openjdk

# source activate, not conda activate! per https://docs.rcc.uchicago.edu/software/apps-and-envs/python/
# source activate /project/lbarreiro/USERS/austin/envs/snakemake