## CONFIG:
import os

# path to lab's shared hg38 .fa index for hisat3n indexing
assembly_path = '/project/lbarreiro/SHARED/REFERENCES/Homo_sapiens/GATK/GRCh38/GRCh38.primary_assembly.genome.fa'
# pi account for slurm
slurm_account = 'pi-lbarreiro'

# these are so lightweight that they can be run directly on the login node; no need for slurm
# will want to add figure generation to this
localrules: cat_fastqs, generate_tagvalues_file, multiqc


## CHOOSE FILES
# find all sample files in the folder, and remove _R1_001 & _R2_001 since paired-end reads will be processed together
# will need to update this to also get rid of lanes
# sample_ids = [f.removesuffix('_R1_001.fastq.gz').removesuffix('_R2_001.fastq.gz') 
#               for f in os.listdir('data/fastq_symlinks') 
#               if f.endswith('.fastq.gz')][0:2]

sample_ids = [
    # 'LB-HT-28s-HT-05_S5', # EVIL EVIL EVIL. WHY DOES IT RUN SO SLOWLY???? # HE WHO EATS OUR PRECIOUS SUs
    'LB-HT-28s-HT-02_S2',
    'LB-HT-28s-HT-03_S3', # topped out the memory at 80GB
    'LB-HT-28s-HT-06_S6',
    'LB-HT-28s-HT-07_S7',
    'LB-HT-28s-HT-08_S8',
    'LB-HT-28s-HT-09_S9',
    'LB-HT-28s-HT-10_S10',
    'LB-HT-28s-HT-12_S12',
    'LB-HT-28s-HT-16_S16',
    'LB-HT-28s-HT-17_S17', # smallest fileset (collectively R1 6.5KB + R2 6.4KB); use this as a test
    'LB-HT-28s-HT-18_S18', # second smallest fileset (collectively R1 15GB + R2 14GB)
    'LB-HT-28s-JL-01_S19', # also failed. EVIL.
    'LB-HT-28s-JL-02_S20',
    'LB-HT-28s-JL-04_S22',
    'LB-HT-28s-JL-05_S23',
    'LB-HT-28s-JL-06_S24',
    'LB-HT-28s-JL-07_S25',
    'LB-HT-28s-JL-08_S26'
]
LANES = [5, 6, 7, 8] # user-defined sequencing lanes

## OTHER USER-DEFINED SETTINGS
substitutions_min = 2 # minimum T>C substitutions for a transcript to be called 'nascent'

rule all:
    input: 
        expand(
            "data/nascent_counts/{sample_id}_nascent_counts.bam",
            sample_id = sample_ids
        ),
        'outputs/multiqc_report.html'

rule cat_fastqs:
    input: 
        lambda wildcards: expand(
            "data/fastq_symlinks/{{sample_id}}_L{seq_lane}_R{{end}}_001.fastq.gz",
            seq_lane=[f"{l:03d}" for l in LANES]
        )
    output: 
        temp('data/fastq_merged/{sample_id}_R{end}_001.fastq.gz')
    shell: 
        "cat {input} > {output}"

rule process_fastp:
    input: 
        r1 = 'data/fastq_merged/{sample_id}_R1_001.fastq.gz',
        r2 = 'data/fastq_merged/{sample_id}_R2_001.fastq.gz'
    output:
        r1 = temp('data/trimmed/{sample_id}_R1_001.fastq.gz'),
        r2 = temp('data/trimmed/{sample_id}_R2_001.fastq.gz'),
        html = "data/fastp_reports/{sample_id}.html",
        json = "data/fastp_reports/{sample_id}.json"
    threads: 4
    resources: 
        slurm_account = 'pi-lbarreiro',
        runtime = 75, # takes 45-50 min so far, adding a little buffer for bigger files
        mem = "3G"
    shell:
        (
            "fastp "
            "--in1 {input.r1} --in2 {input.r2} "
            "--out1 {output.r1} --out2 {output.r2} "
            "--html {output.html} "
            "--json {output.json} "
            "--thread {threads} "
            "--detect_adapter_for_pe " # dynamically detect adapter sequences
        )

rule build_hisat3n_index:
    input: 
        assembly = '/project/lbarreiro/SHARED/REFERENCES/Homo_sapiens/GATK/GRCh38/GRCh38.primary_assembly.genome.fa' # make variable
    output: 
        expand(
            'data/hisat3n_indexes/hg38.3n.{converted_bases}.{index_n}.ht2',
            converted_bases = ['CT', 'GA'], # using a C>T conversion + reverse complement
            index_n = range(1, 9)
        )
    resources: 
        slurm_account = 'pi-lbarreiro',
        runtime = 120,
        mem_mb = 24000
    shell:  
        (
            "hisat-3n-build "
            "--base-change T,C " # don't penalize T>C substitutions
            "{input.assembly} "
            "data/hisat3n_indexes/hg38 "
        )
    
rule align_hisat3n:
    input: 
        fastq_r1 = 'data/trimmed/{sample_id}_R1_001.fastq.gz',
        fastq_r2 = 'data/trimmed/{sample_id}_R2_001.fastq.gz',
        index = expand(
            'data/hisat3n_indexes/hg38.3n.{converted_bases}.{index_n}.ht2',
            converted_bases = ['CT', 'GA'], # T>C for SLAM-seq plus reverse complement
            index_n = range(1, 9) # 9 total index files
        )
    output:
        aligned_bam = 'data/aligned_bam/{sample_id}_aligned.bam'
    params:
        index_prefix = 'data/hisat3n_indexes/hg38',
        hisat_threads = 20,
        samtools_threads = 4,
        # turns out these operations are MASSIVELY IO limited, so we'll copy everything to scratch, where it won't compete for RW access
        scratch = lambda wildcards: f"/scratch/midway3/$USER/{wildcards.sample_id}_$SLURM_JOB_ID"
    log:
        qc = "logs/hisat-3n/qc/{sample_id}.log", # alignment quality
        log = "logs/hisat-3n/logs/{sample_id}.log" # timestamped log file
    threads: 24
    resources:
        job_name = lambda wildcards: f"{wildcards.sample_id}_align_hisat3n",
        slurm_account = 'pi-lbarreiro',
        mem_mb = "96G",
        runtime = 720    # 12 hours in minutes
    shell:
        """
        # Create local scratch directory
        mkdir -p {params.scratch}
        
        # Copy input files to local scratch
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying inputs to local scratch..." >> {log.log} 2>&1
        cp {input.fastq_r1} {params.scratch}/ 2>> {log.log}
        cp {input.fastq_r2} {params.scratch}/ 2>> {log.log}
        
        # Get basenames for local files
        R1_BASE=$(basename {input.fastq_r1})
        R2_BASE=$(basename {input.fastq_r2})
        
        # Run HISAT-3N and pipe directly to BAM
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting alignment..." >> {log.log} 2>&1
        hisat-3n \
            -x {params.index_prefix} \
            -1 {params.scratch}/$R1_BASE \
            -2 {params.scratch}/$R2_BASE \
            -p {params.hisat_threads} \
            --base-change T,C \
            --rna-strandness RF \
            -q \
            --new-summary \
            --summary-file {log.qc} \
            2>> {log.log} | \
        samtools view -@ {params.samtools_threads} -b | \
        samtools sort -@ {params.samtools_threads} -o {params.scratch}/output.bam \
        2>> {log.log}
        
        # Copy result back to /project/
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying output back..." >> {log.log} 2>&1
        cp {params.scratch}/output.bam {output.aligned_bam} 2>> {log.log}
        
        # Cleanup
        rm -rf {params.scratch}
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Alignment complete" >> {log.log} 2>&1
        """

rule generate_tagvalues_file:
    output: 
        tags = "data/tag_values.txt"
    run: 
        with open(output.tags, 'w') as f: # write sequence to .txt, one int per line
            for i in range(substitutions_min, 301):
                f.write(f'{i}\n')

rule count_nascent_transcripts:
    input: 
        bamfile = "data/aligned_bam/{sample_id}_aligned.bam",
        tags = "data/tag_values.txt" # file with integers 2-300 against which samtools will match Yf:
    output: 
        nascent_counts = "data/nascent_counts/{sample_id}_nascent_counts.bam"
    resources:
        slurm_account = slurm_account,
        runtime = 8,
        mem = "2G"
    shell: 
        (
            "samtools view "
            "--with-header " 
            "--bam "
            "-F 260 " # exclude 260 field—multiple alignments
            "-D Yf:{input.tags} " # only alignments with Yf: tag == STR where STR in input.tags
            "{input.bamfile} > {output.nascent_counts}"
        )

rule multiqc:
    input: 
        expand(
            "data/fastp_reports/{sample_id}.{filetype}",
            sample_id = sample_ids,
            filetype = ['json']
        ),
        expand(
            "logs/hisat-3n/qc/{sample_id}.log",
            sample_id = sample_ids
        ),
    output: 
        'outputs/multiqc_report.html'
    shell: 
        (
            'multiqc '
            'data/fastp_reports logs/hisat-3n/qc '
            '--force ' # overwrite existing report; otherwise it will attach a suffix that snakemake won't detect
            '--outdir outputs'
            # future: --ignore-samples for ones that failed to process
        )