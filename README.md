# SLAM-seq rotation project
*January 2026*

*Austin Szatrowski [[email](mailto:aszatrowski@uchicago.edu)], based on a SLAMseq analysis pipeline by Jesse Lehman, Pai Lab @ UMass Med*

### Rulegraph:
![pipeline rulegraph](outputs/rulegraph.png)

## SETUP (START HERE)
To run this pipeline, you will need to (1) set up a conda environment, (2) build a apptainer/docker container for GEDI (a SLAM-seq quantification software package), and (3) execute the pipeline.

1. Build the apptainer container: 
```bash
apptainer build /project/your-pi/you/containers/gedi_R.sif slamseq/gedi.def
```
In `config.yaml`, set `container_path` to `/project/your-pi/you/containers/gedi_R.sif` so snakemake will know where to find it.

2. Build the conda environment, and activate it:
```bash
conda create profile/environment.yaml > /path/to/envs/slamseq_v3 
source setup.sh
```
* `source activate` may be different on your cluster.

3. Execute the pipeline, using settings in `profile/config.yaml`
```bash
snakemake --workflow-profile ./profile [-np]
```
* `-np` indicates a snakemake dry run, which prints out all the instances of all the rules to be executed and their corresponding shell commands. Gives an excellent overview of the work to be done, and is useful for debugging.
* `profile/config.yaml` is configured for slurm and UChicago's midway3 architecture; some revisions may be necessary.

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
* `data`: symlinks to raw `.fastq` data in Hannah's personal folder, plus all intermediate outputs (before clearance by snakemake, see `temp()`) and space for temp files. Anything that is too large to be part of a GitHub repo, not human-readable, or isn't an interesting result goes in here
* `outputs`: results, metadata, and QC summary folders. Of note, a _summary_ of the fastp checks go here, but individual reports go in `data`.


## Global inputs
* `*.fastq.gz` compressed sequencing files
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

## Approximate timing
* `cat_fastqs`: seconds
* `fastp`
    * 4 cores: ~45-60 min per sample
* `star`
    * 8 threads: ~1hr/10G of trimmed reads
* `bam_to_cit`: 18h for ~150G of aligned BAMs
* `grand_slam`: 
    * 16 threads: 5h for 76G CIT file