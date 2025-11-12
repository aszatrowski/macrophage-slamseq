## CONFIG:
import os

# path to lab's shared hg38 .fa index for hisat3n and GEDI indexing
assembly_path = '/project/lbarreiro/SHARED/REFERENCES/Homo_sapiens/GATK/GRCh38/GRCh38.primary_assembly.genome.fa'
data_path = '/project/lbarreiro/DATA/SLAM-seq/pilot_2025.08.19/FastX'
container_path = '/project/lbarreiro/USERS/austin/containers/gedi_1.0.6a.sif'
# pi account for slurm

GEDI_INDEX_DIR = 'data/gedi_indexes'
# these are so lightweight that they can be run directly on the login node; no need for slurm
# will want to add figure generation to this
localrules: cat_fastqs, wget_hg38_gtf, multiqc, gedi_index_genome, make_bamlist


## CHOOSE FILES
# find all sample files in the folder, and remove _R1_001 & _R2_001 since paired-end reads will be processed together
# will need to update this to also get rid of lanes
# sample_ids = [f.removesuffix('_R1_001.fastq.gz').removesuffix('_R2_001.fastq.gz') 
#               for f in os.listdir('data/fastq_symlinks') 
#               if f.endswith('.fastq.gz')][0:2]

sample_ids = [
    'LB-HT-28s-HT-01_S1',
    'LB-HT-28s-HT-02_S2',
    'LB-HT-28s-HT-03_S3', # topped out the memory at 80GB. but runs now!
    'LB-HT-28s-HT-05_S5', # EVIL EVIL EVIL. WHY DOES IT RUN SO SLOWLY???? # HE WHO EATS OUR PRECIOUS SUs
    'LB-HT-28s-HT-06_S6',
    'LB-HT-28s-HT-07_S7',
    'LB-HT-28s-HT-08_S8',
    'LB-HT-28s-HT-09_S9',
    'LB-HT-28s-HT-10_S10',
    'LB-HT-28s-HT-12_S12',
    'LB-HT-28s-HT-16_S16',
    'LB-HT-28s-HT-17_S17',# smallest fileset (collectively R1 6.5KB + R2 6.4KB); use this as a test BUT has no nascent transcripts and fails featureCounts
    'LB-HT-28s-HT-18_S18', # second smallest fileset (collectively R1 15GB + R2 14GB)
    # 'LB-HT-28s-JL-01_S19', # also failed. EVIL. EVEN AFTER 20HOURS??
    'LB-HT-28s-JL-02_S20',
    'LB-HT-28s-JL-04_S22',
    'LB-HT-28s-JL-05_S23',
    'LB-HT-28s-JL-06_S24',
    'LB-HT-28s-JL-07_S25',
    'LB-HT-28s-JL-08_S26',
    'LB-HT-28s-JL-09_S27'
]

LANES = [5, 6, 7, 8]

SAMPLE_METADATA = {
    "LB-HT-28s-HT-01_S1": {"timepoint": 0, "donor": "1rep1"},
    "LB-HT-28s-HT-02_S2": {"timepoint": 15, "donor": "1rep1"},
    "LB-HT-28s-HT-03_S3": {"timepoint": 30, "donor": "1rep1"},
    "LB-HT-28s-HT-05_S5": {"timepoint": 60, "donor": "1rep1"},
    "LB-HT-28s-HT-06_S6": {"timepoint": 90, "donor": "1rep1"},
    "LB-HT-28s-HT-07_S7": {"timepoint": 105, "donor": "1rep1"},
    "LB-HT-28s-HT-08_S8": {"timepoint": 120, "donor": "1rep1"},
    "LB-HT-28s-HT-09_S9": {"timepoint": 120, "donor": "1rep1"},  # no4sU
    "LB-HT-28s-HT-10_S10": {"timepoint": 15, "donor": "2rep1"},
    "LB-HT-28s-HT-12_S12": {"timepoint": 45, "donor": "2rep1"},
    "LB-HT-28s-HT-16_S16": {"timepoint": 105, "donor": "2rep1"},
    "LB-HT-28s-HT-18_S18": {"timepoint": 120, "donor": "2rep1"},  # no4sU
    "LB-HT-28s-JL-01_S19": {"timepoint": 0, "donor": "2rep2"},
    "LB-HT-28s-JL-02_S20": {"timepoint": 15, "donor": "2rep2"},
    "LB-HT-28s-JL-03_S21": {"timepoint": 30, "donor": "2rep2"},
    "LB-HT-28s-JL-04_S22": {"timepoint": 45, "donor": "2rep2"},
    "LB-HT-28s-JL-05_S23": {"timepoint": 60, "donor": "2rep2"},
    "LB-HT-28s-JL-06_S24": {"timepoint": 75, "donor": "2rep2"},
    "LB-HT-28s-JL-07_S25": {"timepoint": 90, "donor": "2rep2"},
    "LB-HT-28s-JL-08_S26": {"timepoint": 105, "donor": "2rep2"},
    "LB-HT-28s-JL-09_S27": {"timepoint": 120, "donor": "2rep2"}
}

# DONORS = ["1rep1", "2rep1", "2rep2"]
DONORS = ["2rep2"]

def get_donor_samples(donor):
    return [sample_id for sample_id, meta in SAMPLE_METADATA.items() 
            if meta["donor"] == donor]

rule all:
    input: 
        expand(
            [
                f'{GEDI_INDEX_DIR}/{{input_basename}}.genes.tab',
                f'{GEDI_INDEX_DIR}/{{input_basename}}.index.cit',
                f'{GEDI_INDEX_DIR}/{{input_basename}}.transcripts.fasta',
                f'{GEDI_INDEX_DIR}/{{input_basename}}.transcripts.fi',
                f'{GEDI_INDEX_DIR}/{{input_basename}}.transcripts.tab',
                f'{GEDI_INDEX_DIR}/{{input_basename}}_metadata.oml'
            ],
            input_basename = ['ncbi_refseq_hg38.gtf']
        ),
        expand(
            "data/cit/{donor}.cit",
            donor = DONORS
        ),
        # 'outputs/multiqc_report.html'

rule cat_fastqs:
    input: 
        lambda wildcards: expand(
            f"{data_path}/{{sample_id}}_L{{seq_lane}}_R{{end}}_001.fastq.gz",  # ← double braces
            sample_id=wildcards.sample_id,
            end=wildcards.end,
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
            converted_bases = ['CT', 'GA'], # T>C for AM-seq plus reverse complement
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
        mem_mb = "96G",
        runtime = 1440    # 24 hours in minutes # please let this be enough
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

rule add_md_tags:
    input:
        bam = "data/aligned_bam/{sample_id}_aligned.bam",
        fasta = assembly_path
    output:
        bam = "data/md_tagged/{sample_id}.bam",
        bai = "data/md_tagged/{sample_id}.bam.bai"
    resources:
        runtime = 120,
        mem = "12G"
    threads: 8
    params:
        scratch = lambda wildcards: f"/scratch/midway3/$USER/{wildcards.sample_id}_$SLURM_JOB_ID",
    shell:
        """
        mkdir -p {params.scratch}
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying inputs to local scratch..."
        cp {input.bam} {params.scratch}/ 2>> {params.scratch}/samtools.log
        cp {input.fasta} {params.scratch}/ 2>> {params.scratch}/samtools.log

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting tagging..."
        # get basename 
        BASE=$(basename {input.bam})
        BASE_FA=$(basename {input.fasta})
        
        # run samtools on scratch-copied file, using scratch path and basename
        samtools calmd -b {params.scratch}/$BASE {params.scratch}/$BASE_FA 2> {params.scratch}/samtools.log | \
        samtools sort -@ {threads} -o {params.scratch}/output.bam 2>> {params.scratch}/samtools.log
        samtools index {params.scratch}/output.bam 2>> {params.scratch}/samtools.log

        # Copy result back to /project/
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying output back..."
        cp {params.scratch}/output.bam {output.bam} 2>> {params.scratch}/samtools.log
        cp {params.scratch}/output.bam.bai {output.bai} 2>> {params.scratch}/samtools.log
        

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Tagging complete."
        # clean up
        rm -rf {params.scratch}
        """

rule wget_hg38_gtf:
    output: 
        "data/ncbi_refseq_hg38/ncbi_refseq_hg38.gtf.gz"
    params:
        refseq_path = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz"
    shell: 
        (
            "wget {params.refseq_path} "
            "--output-document {output} "
            "--force-directories "
            "--show-progress "
        )

rule gedi_index_genome:
    input:
        fasta = assembly_path,
        gtf = "data/ncbi_refseq_hg38/{input_basename}.gz"
    output:
        genes_tab = f'{GEDI_INDEX_DIR}/{{input_basename}}.genes.tab',
        index_cit = f'{GEDI_INDEX_DIR}/{{input_basename}}.index.cit',
        transcripts_index = f'{GEDI_INDEX_DIR}/{{input_basename}}.transcripts.fasta',
        transcripts_fi = f'{GEDI_INDEX_DIR}/{{input_basename}}.transcripts.fi',
        transcripts_tab = f'{GEDI_INDEX_DIR}/{{input_basename}}.transcripts.tab',
        metadata = f'{GEDI_INDEX_DIR}/{{input_basename}}_metadata.oml'
    container:
        container_path
    params:
        output_dir = GEDI_INDEX_DIR
    log:
        "logs/gedi/index_{input_basename}.log"
    shell:
        f"""
            gedi \
            -e IndexGenome \
            -s {{input.fasta}} \
            -a {{input.gtf}} \
            -f {GEDI_INDEX_DIR} \
            -o {GEDI_INDEX_DIR}/{{wildcards.input_basename}}_metadata.oml \
            -n {{wildcards.input_basename}} \
            -nobowtie -nostar -nokallisto \
            2> {{log}}
        """

rule make_bamlist:
    input: 
        bams = lambda wildcards: expand(
            "data/md_tagged/{sample_id}.bam",
            sample_id = get_donor_samples(wildcards.donor)
        )
    output: 
        "data/bamlist_cit/{donor}.bamlist"
    shell: 
        "printf '%s\\n' {input} > {output}"

rule bam_to_cit:
    input:
        "data/bamlist_cit/{donor}.bamlist"
    output: 
        "data/cit/{donor}.cit"
    container:
        container_path
    resources:
       runtime = 30,
       mem = "2G",
    log:
        "logs/gedi/{donor}.log"
    shell: # converts aligned and tagged .bam files to GRAND-SLAM's custom CIT format
        """
        bamlist2cit {input} 2> {log}
        mv data/md_tagged/{wildcards.donor}.cit {output} 2> log
        """

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
            'data/fastp_reports logs/hisat-3n/qc ' # switch these to ALL for consistency in multiqc
            '--force ' # overwrite existing report; otherwise it will attach a suffix that snakemake won't detect
            '--outdir outputs'
            # future: --ignore-samples for ones that failed to process
        )