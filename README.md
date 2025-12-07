# SLAM-seq rotation project
*Austin Szatrowski, based on a SLAMseq analysis pipeline by Jesse Lehman, Pai Lab @ UMass Med*

### Rulegraph:
![pipeline rulegraph](outputs/rulegraph.png)

## SETUP (START HERE)
To run this pipeline, you will need to (1) set up a conda environment, (2) build a apptainer/docker container for GRAND-SLAM, and (3) execute the pipeline.

1. Build the apptainer container: 
```bash
apptainer build /project/lbarreiro/USERS/austin/containers/gedi_1.0.6a.sif slamseq/gedi.def
```

2. Build the conda environment, and activate it:
```bash
conda create profile/environment.yaml > slamseq_v3 
source setup.sh
```

3. Execute the pipeline, using settings in `profile/config.yaml`
```bash
snakemake --workflow-profile ./profile
```
## Containers

**Make sure you do this before running the pipeline.**
* This is a container (like a conda environment), but it can specify any system configuration, including OS versions, memory, environment variables, and arbitrary software, rather than just packages on `conda-forge` or `bioconda`. Since GEDI/GRAND-SLAM is on neither, I've specified a container to pin its version and let it run in isolation, just as one would with conda.
    * The above command builds the container from the definition file (analogous to conda's `environment.yaml`) called `gedi.def`.
    * Build whereever is convenient (I use a personal directory, and make sure to set `container_path` in `config.yaml` to whereever you've built it, so snakemake can access it.)

## envs:
* `snakemake` is clean snakemake with the slurm execution plugin
* `slamseq_v2` is all the dependencies for Jesse's SLAM-seq pipeline with multiqc
    * major packages:
        * `snakemake`
            * `snakemake-executor-plugin-slurm`
        * `fastp`
        * `star`
        * `samtools`
        * `multiqc`

## Folders
* `.snakemake`: under-the-hood metadata for DAG creation, plus execution (`log/`) and rule (`slurm_logs/`) log files
* `.gedi` An annoying folder that GEDI (the software from GRAND-SLAM comes) creates. I'm afraid to touch it.
* `config/genomic` will contain the index files for GEDI. So far as I can tell, they have to be in a folder with this name, otherwise GEDI commands won't recognize them.
* `data`: symlinks to raw `.fastq` data in Hannah's personal folder, plus all intermediate outputs (before clearance by snakemake, see `temp()`) and space for temp files. Anything that is too large to be part of a GitHub repo, not human-readable, or isn't an interesting result goes in here
* `outputs`: results, metadata, and QC summary folders. Of note, a _summary_ of the fastp checks go here, but individual reports go in `data`.


## Global inputs
* `*.fastq.gz` compressed sequencing files
* `ref_genome_assembly` a `.fa` reference genome file against which to align. I have used the Barreiro shared hg38 one.
## Input file specification
* As it stands, each input `sample_id` (`/project/pi/path/to/data/{sample_id}.fastq.gz`) should be associated with a donor and a timepoint in `config.yaml` (not to be confused with `profile/config.yaml`, which is used for execution settings.)
    * `sample_id`


## Outputs
* GRAND-SLAM output fileset

## Configuration
In general, snakemake will submit each job (one rule, executing for one input) in the workflow as its own `sbatch` job. The SBATCH arguments are specified in `profile/config.yaml`:
```
executor: slurm # use slurm
cores: 32 # maximum cores during execution, shared across jobs
jobs: 8 # number of parallel jobs to run
use-singularity: true # activate singularity for container
singularity-args: "--bind /project/lbarreiro,/scratch/midway3/$USER" # give container access to all of lbarreiro & personal scratch
default-resources:
  slurm_account: pi-lbarreiro 
  slurm_partition: lbarreiro
```

### Rules
#### `cat_fastqs`
* Concatenates `.fastq.gz` files from the same sample (different lanes) together before preprocessing and alignment.
* Runs directly on login node

#### `process_fastp`
* Pulls in each sequencing file (R1 and R2 separately), and runs `fastp`, which runs a FastQC-style sequence QC, and adapter trimming
* Output: trimmed fastq sequence files, `html` report and `json` machine-readable metadata
    * saved in `data/trimmed/`

#### `star`
* explain copy to scratch
#### `multiqc`
* Reads `data/fastp_reports/` `logs/star/qc` for per-sample QC files, and generates a beautifully formatted HTML summary with interative plots
#### `index_bam`
#### `rename_with_donor_timepoint`
#### `mark_no4sU_samples`
#### `gedi_index_genome`
#### `bam_to_cit`
#### `grand_slam`
#### `plot_read_timeseries`