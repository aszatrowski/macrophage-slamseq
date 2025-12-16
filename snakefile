import os
## CONFIG:
configfile: "config.yaml"
# these are so lightweight that they can be run directly on the login node; no need for slurm
localrules: cat_fastqs, index_bam, rename_with_donor_timepoint, mark_no4sU_samples, multiqc

DONORS = ['donor1_rep2']
sample_ids = list(config['sample_ids'].keys())

rule all:
    input: 
        expand(
            "data/slam_quant/{donor}/grandslam.tsv.gz",
            donor = DONORS
        ),
        # 'outputs/multiqc_report.html'

rule cat_fastqs:
    """
    Joins fastq files together from multiple sequencing lanes. Lanes are specified in config['sequencing_lanes'].
    """
    input: 
        lambda wildcards: expand(
            f"{config['data_path']}/{{sample_id}}_L{{seq_lane}}_R{{end}}_001.fastq.gz",  # ← double braces
            sample_id=wildcards.sample_id,
            end=wildcards.end,
            seq_lane=[f"{l:03d}" for l in config['sequencing_lanes']]
        )
    output: 
        temp('data/fastq_merged/{sample_id}_R{end}_001.fastq.gz')
    shell: 
        "cat {input} > {output}"

rule fastp:
    """
    Runs fastp for quality control and adapter trimmming. fastp (https://github.com/OpenGene/fastp) is a faster alternative to fastqc + trimmomatic/cutadapt, and crucially only needs to be run once, rather than fastqc > trim fastqc. Though fastp can dynamically detect adapter sequences quite efficiently, in practice it misses bases on the ends, which introduces artificial mismatches which may reduce T>C substitution calling downstream, so manual specification is recommended.
    The parameters are configurd for paired-end sequencing, so you'll need to modify if you have single-end. See the fastp docs for details.
    Outputs an html webpage and a machine-readable json QC summary for each sample; multiqc will take the json as input for its own summary.
    """
    input: 
        r1 = 'data/fastq_merged/{sample_id}_R1_001.fastq.gz',
        r2 = 'data/fastq_merged/{sample_id}_R2_001.fastq.gz'
    output:
        r1 = temp('data/trimmed/{sample_id}_R1_001.fastq.gz'),
        r2 = temp('data/trimmed/{sample_id}_R2_001.fastq.gz'),
        html = "data/fastp_reports/{sample_id}.html",
        json = "data/fastp_reports/{sample_id}.json"
    params:
        adapter_sequence_r1 = config['adapter_sequence_r1_rc'],
        adapter_sequence_r2 = config['adapter_sequence_r2_rc']
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
            "--adapter_sequence={params.adapter_sequence_r1} "
            "--adapter_sequence_r2={params.adapter_sequence_r2}"
        )

rule star:
    """
    Uses STAR (https://github.com/alexdobin/STAR) to align the trimmed sequences.
    STAR is heavily I/O limited during its complex computations, so I have snakemake copy all the necessary files to the /scratch/ filesystem, where they won't compute with everyone else's I/O. According to my basic """
    input: 
        fastq_r1 = 'data/trimmed/{sample_id}_R1_001.fastq.gz',
        fastq_r2 = 'data/trimmed/{sample_id}_R2_001.fastq.gz',
        index = '/project/lbarreiro/SHARED/REFERENCES/Homo_sapiens/GATK/GRCh38/STAR'
    output:
        aligned_bam = 'data/aligned_bam/{sample_id}.bam',
        stats = 'logs/star/qc/{sample_id}.Log.final.out'
    params:
        # turns out these operations are MASSIVELY IO limited, so we'll copy everything to scratch, where it won't compete for RW access
        scratch = lambda wildcards: f"/scratch/midway3/$USER/{wildcards.sample_id}_$SLURM_JOB_ID",
        max_perread_sub_fraction = 0.25
    log:
        "logs/star/{sample_id}.log" # alignment quality
    benchmark:
        "benchmarks/{sample_id}.star_align.benchmark.txt"
    threads:
        8
    resources:
        job_name = lambda wildcards: f"{wildcards.sample_id}_align_star",
        mem = "32G",
        runtime = 240    # 4 hours in minutes
    shell:
        """
        # Create local scratch directory
        mkdir -p {params.scratch}
        
        # Copy input files to local scratch
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying inputs to local scratch..." >> {log} 2>&1
        cp {input.fastq_r1} {params.scratch}/ 2>> {log}
        cp {input.fastq_r2} {params.scratch}/ 2>> {log}
        cp -R {input.index} {params.scratch}/ 2>> {log}
        
        # Get basenames for local files
        R1_BASE=$(basename {input.fastq_r1})
        R2_BASE=$(basename {input.fastq_r2})
        INDEX_BASE=$(basename {input.index})
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting alignment..." >> {log} 2>&1
        STAR \
            --runMode alignReads \
            --genomeDir {params.scratch}/$INDEX_BASE \
            --readFilesIn {params.scratch}/$R1_BASE {params.scratch}/$R2_BASE \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --outFilterMismatchNmax 999 \
            --limitBAMsortRAM 16000000000 \
            --alignEndsType EndToEnd \
            --outFilterMismatchNoverReadLmax {params.max_perread_sub_fraction} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes MD NH AS \
            --outFileNamePrefix {params.scratch}/ \
            2> {log}
        
        # Copy results back to /project/
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying output back..." >> {log} 2>&1
        cp {params.scratch}/Aligned.sortedByCoord.out.bam {output.aligned_bam} 2>> {log}
        cp {params.scratch}/Log.final.out {output.stats} 2>> {log}
        cat {params.scratch}/Log.out >> {log}
        
        # Cleanup
        rm -rf {params.scratch}
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Alignment complete." >> {log} 2>&1
        """

rule gedi_index_genome:
    output:
        fasta = f"{config['gedi_index_dir']}/homo_sapiens.115.fasta",
        fi = f"{config['gedi_index_dir']}/homo_sapiens.115.fi",
        genes_tab = f"{config['gedi_index_dir']}/homo_sapiens.115.genes.tab",
        gtf = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf",
        index_cit = f"{config['gedi_index_dir']}/homo_sapiens.115.index.cit",
        transcripts_fasta = f"{config['gedi_index_dir']}/homo_sapiens.115.transcripts.fasta",
        transcripts_fi = f"{config['gedi_index_dir']}/homo_sapiens.115.transcripts.fi",
        transcripts_tab = f"{config['gedi_index_dir']}/homo_sapiens.115.transcripts.tab",
        oml = f"{config['gedi_index_dir']}/homo_sapiens.115.oml"
    container:
        config['container_path']
    resources:
        mem_mb= "8G",
        runtime = 5 
    params:
        output_dir = config['gedi_index_dir']
    shell:
        f"""
            mkdir -p config/genomic
            gedi -e IndexGenome -organism homo_sapiens -version 115 -f config/genomic -o {config['gedi_index_dir']}/homo_sapiens.115.oml -nomapping
        """

rule index_bam:
    input: 
        "data/aligned_bam/{sample_id}.bam"
    output: 
        "data/aligned_bam/{sample_id}.bam.bai"
    shell: 
        """
        samtools index --bai {input} -o {output}
        """

def get_sample_from_donor_timepoint(donor, timepoint):
    """
    Given a donor (as requested in rule all:), and a specific timepoint, retrieves the corresponding sample
    """
    for sample_id, info in config["sample_ids"].items():
        if info["donor"] == donor and info["timepoint"] == int(timepoint):
            return sample_id
    raise ValueError(
        f"No sample found for donor='{donor}' at timepoint='{timepoint}'. "
        f"Available timepoints for {donor}: {get_donor_timepoints(donor)}"
    )

rule rename_with_donor_timepoint:
    input:
        bam=lambda w: f"data/aligned_bam/{get_sample_from_donor_timepoint(w.donor, w.timepoint)}.bam",
        bai=lambda w: f"data/aligned_bam/{get_sample_from_donor_timepoint(w.donor, w.timepoint)}.bam.bai"
    output:
        bam="data/donor_timepoint_symlinks/{donor}/{timepoint}m.bam",
        bai="data/donor_timepoint_symlinks/{donor}/{timepoint}m.bam.bai"
    shell:
        """
        mkdir -p $(dirname {output.bam})
        ln -sf $(realpath {input.bam}) {output.bam}
        ln -sf $(realpath {input.bai}) {output.bai}
        """

rule mark_no4sU_samples:
    input:
        bam=lambda w: f"data/aligned_bam/{config['donor_controls'][w.donor]}.bam",
        bai=lambda w: f"data/aligned_bam/{config['donor_controls'][w.donor]}.bam.bai",
    output:
        bam="data/donor_timepoint_symlinks/{donor}/control_no4sU.bam",
        bai="data/donor_timepoint_symlinks/{donor}/control_no4sU.bam.bai"
    shell:
        """
        mkdir -p $(dirname {output.bam})
        ln -sf $(realpath {input.bam}) {output.bam}
        ln -sf $(realpath {input.bai}) {output.bai}
        """

def get_donor_timepoints(donor):
    """Get all timepoints that exist for a given donor"""
    timepoints = [info["timepoint"] for sample_id, info in config["sample_ids"].items() 
                  if info["donor"] == donor]
    return sorted(set(timepoints))

rule bam_to_cit:
    input:
        bams = lambda w: expand(
            "data/donor_timepoint_symlinks/{donor}/{timepoint}m.bam",
            donor = w.donor,
            timepoint = get_donor_timepoints(w.donor)
        ),
        bais = lambda w: expand(
            "data/donor_timepoint_symlinks/{donor}/{timepoint}m.bam.bai",
            donor = w.donor,
            timepoint = get_donor_timepoints(w.donor)
        ),
        # ADD IN CONTROL NO4SU FILE
        no4sU_bam="data/donor_timepoint_symlinks/{donor}/control_no4sU.bam",
        no4sU_bai="data/donor_timepoint_symlinks/{donor}/control_no4sU.bam.bai",
    output:
        cit_sample_set = "data/cit_sample_sets/{donor}.cit",
        cit_metadata = "data/cit_sample_sets/{donor}.cit.metadata.json",
    params:
        scratch = lambda wildcards: f"/scratch/midway3/$USER/{wildcards.donor}_$SLURM_JOB_ID",
        bam_basenames = lambda w, input: " ".join([os.path.basename(b) for b in input.bams]),
        java_xmx=lambda w, resources: int(resources.mem_mb * 0.875 / 1024),  # 87.5% in GB
        java_xms=lambda w, resources: max(4, int(resources.mem_mb * 0.125 / 1024))  # 12.5% in GB, min 4
    container:
        config["container_path"]
    resources:
        mem_mb = 32000,
        runtime = 1440 # 20 hours in minutes
    benchmark:
        "benchmarks/{donor}.bam_to_cit.benchmark.txt"
    shell:
        """
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying to scratch..."
        mkdir -p {params.scratch}
        cp {input.bams} {input.bais} {input.no4sU_bam} {input.no4sU_bai} {params.scratch}/
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copy complete."

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Generating output paths..."
        cit_output_path=$(realpath {output.cit_sample_set})
        metadata_output_path=$(realpath {output.cit_metadata})
        mkdir -p $(dirname $cit_output_path)
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Paths created."

        cd {params.scratch}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Setting Java max memory..."
        export _JAVA_OPTIONS="-Xmx{params.java_xmx}g -Xms{params.java_xms}g"

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Beginning CIT conversion..."
        gedi -e Bam2CIT -p output.cit \
            {params.bam_basenames} \
            $(basename {input.no4sU_bam})
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Conversion complete."

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying output back..."
        cp output.cit $cit_output_path
        cp output.cit.metadata.json $metadata_output_path
        rm -rf {params.scratch}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copy complete, scratch cleared."
        """

rule grand_slam:
    input: 
        cit_sample_set = "data/cit_sample_sets/{donor}.cit",
        cit_metadata = "data/cit_sample_sets/{donor}.cit.metadata.json",
        index_oml = rules.gedi_index_genome.output.oml
    output: 
        nascent_counts = "data/slam_quant/{donor}/grandslam.tsv.gz",
        sub_rates = "data/slam_quant/{donor}/grandslam.rates.tsv",
        mismatch_tsv = "data/slam_quant/{donor}/grandslam.mismatches.tsv",
        mismatch_plot = "data/slam_quant/{donor}/grandslam.mismatches.pdf",
        mismatch_positions = "data/slam_quant/{donor}/grandslam.mismatchpos.pdf",
        ntr_stats = "data/slam_quant/{donor}/grandslam.ntrstat.tsv",
    params:
        java_xmx=lambda w, resources: int(resources.mem_mb * 0.875 / 1024),  # 87.5% in GB
        java_xms=lambda w, resources: max(4, int(resources.mem_mb * 0.125 / 1024))  # 12.5% in GB, min 4
    container:
        config['container_path']
    resources:
       runtime = 480, # 8 hours in minutes
       mem_mb = 20000,
    threads:
        16,
    benchmark:
        "benchmarks/{donor}.grandslam.benchmark.txt"
    shell: 
        """
            export _JAVA_OPTIONS="-Xmx{params.java_xmx}g -Xms{params.java_xms}g"
            gedi -e Slam \
            -genomic {input.index_oml} \
            -reads {input.cit_sample_set} \
            -prefix data/slam_quant/{wildcards.donor}/grandslam \
            -introns \
            -no4sUpattern no4sU \
            -nthreads {threads} \
            -progress \
            -full \
            -plot \
            -progress
        """

rule multiqc:
    input: 
        expand(
            "data/fastp_reports/{sample_id}.{filetype}",
            sample_id = sample_ids,
            filetype = ['json']
        ),
        expand(
            'logs/star/qc/{sample_id}.Log.final.out',
            sample_id = sample_ids,
        )
    output: 
        'outputs/multiqc_report.html'
    shell: 
        (
            'multiqc '
            'data/fastp_reports logs/star ' # switch these to ALL for consistency in multiqc
            '--force ' # overwrite existing report; otherwise it will attach a suffix that snakemake won't detect
            '--outdir outputs'
        )