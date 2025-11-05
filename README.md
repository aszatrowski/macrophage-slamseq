# SLAM-seq rotation project
*Austin Szatrowski, based on a SLAMseq analysis pipeline by Jesse Lehman, Pai Lab @ UMass Med*

### Rulegraph:
![pipeline rulegraph](outputs/rulegraph_dag/rulegraph.png)

## envs:
* `snakemake` is clean snakemake with the slurm execution plugin
* `slamseq_v1` is all (I think) the dependencies for Jesse's SLAM-seq pipeline
* `slamseq_fastp` is all the dependencies plus fastp
* `slamseq_v2` is all the dependencies for Jesse's SLAM-seq pipeline with multiqc and featureCounts
    * major packages:
        * `snakemake`
            * `snakemake-executor-plugin-slurm`
            * `snakemake-storage-plugin-http`
        * `fastp`
        * `hisat-3n`
        * `samtools`
        * `multiqc`

## Folders
* `.snakemake`: under-the-hood metadata for DAG creation, plus execution (`log/`) and rule (`slurm_logs/`) log files
* `data`: symlinks to raw `.fastq` data in Hannah's personal folder, plus all intermediate outputs (before clearance by snakemake, see `temp()`) and space for temp files. Anything that is too large to be part of a GitHub repo, not human-readable, or isn't an interesting result goes in here
* `outputs`: results, metadata, and QC summary folders. Of note, a _summary_ of the fastp checks go here, but individual reports go in `data`.

## Pipeline
* `setup.sh` is a setup file (run before anything else) that loads the correct python version and conda environment for snakemake, its plugins (slurm executor) and all the packages.

### Global inputs
* `*.fastq.gz` compressed sequencing files
* `ref_genome_assembly` a `.fa` reference genome file against which to align. I have used the Barreiro shared hg38 one.

### Outputs
* Read counts (so far)

### Rules
#### `cat_fastqs`
* Concatenate `.fastq.gz` files from the same sample (different lanes) together before preprocessing and alignment.
* Runs directly on login node, no need for slurm

#### `process_fastp`
* Pulls in each sequencing file (R1 and R2 separately), and runs `fastp`, which **includes adapter trimming**
* Output: trimmed fastq sequence files, `html` report and `json` machine-readable metadata
    * saved in `data/trimmed/`

#### `build_hisat3n_index`
* Builds a `hisat-3n` index fileset from a reference `.fa` genome file
* I used the one in Barreiro lab's shared reference directory, seemed to work fine
* Inserts T>C and complementary G>A substitutions in the genome to match the SLAMseq substitutions
* Output: 16 total files
    * `data/hisat3n_indexes/hg38.3n.CT.[1-8].ht2`
    * `data/hisat3n_indexes/hg38.3n.GA.[1-8].ht2`

#### `align_hisat3n`
* Align `data/trimmed/*_R[n]_001.fastq` read files to hg38 using the hg38 hisat-3n reference fileset
    * First, copy the fastq files to a scratch folder `/scratch/midway3/aszatrowski/{sample_id}_{jobid}`, where the CPU-heavy alignment won't compete for IO with everything else. Empirically, this boosts CPU efficiency from ~0.23% to ~80%.
    * Run the alignment on the copied files, and write the output to log files
    * Pipe the results directly into `samtools` to convert the .sam output to .bam
    * Sort the alignments by genomic position, and write the result to `output.bam`
    * Copy output back to `/project/` with sample id
    * Wipe the `{sample_id}_{jobid}` directory
* Input: PAIRS of `data/trimmed/*_R[n]_001.fastq` files, since these represent the forward and reversed paired-end adapters used in the same sample
* Output: 
    * Per-sample aligned `.bam` files

#### `generate_tagvalues_file`
* `samtools view` counting requires string matching (per JL, it seems like there ought to be a better way), so this generates a `.txt` file containing the strings "2", "3", ..., "300" so it can call nascent transcripts.

#### `count_nascent_transcripts`
* Uses `samtools view` to read the `Yf:i:<str>` tag at the end of each transcript in `{sample_id}_aligned.bam`
* `Yf:i:<str>` contains the number of specified (T>C) substitutions present in the transcript, as determined by `hisat-3n`.
* If `<str>` appears in `tagvalues.txt` (from above), meaning the transcript contains ≥2 T>C substitutions, then that transcript is called as nascent
    * Per JL, there are slightly more careful ways to do this involving corrections for true SNPs and incomplete 4sU conversion, but this is fine for now.

#### `wget_hg38_gtf`
* Pulls in a gtf file for featureCounts. To be removed

#### `featureCounts_nascent`
* Run `featureCounts` over nascent transcript subset
* Going to be removed in favor of GRAND-SLAM for nascent transcript calling and quantification

#### `multiqc`
* Reads `data/fastp_reports/` `logs/hisat-3n/qc` `data/featurecounts/nascent/` for per-sample QC files, and generates a nicely formatted summary
