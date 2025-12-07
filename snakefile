## CONFIG:
configfile: "config.yaml"
# these are so lightweight that they can be run directly on the login node; no need for slurm
# will want to add figure generation to this
localrules: cat_fastqs, index_bam, rename_with_donor_timepoint, mark_no4sU_samples, multiqc

DONORS = ['donor1_rep2']
# sample_ids = [sample for donor in DONORS for sample in config['donor_sample_ids'][donor]]
sample_ids = list(config['sample_ids'].keys())
# print(sample_ids)


rule all:
    input: 
        # expand(
        #     "data/cit_sample_sets/{donor}.cit",
        #     donor = DONORS
        # ),
        expand(
            "data/slam_quant/{donor}/grandslam.tsv",
            donor = DONORS
        ),
        'outputs/multiqc_report.html'

rule cat_fastqs:
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

rule build_star_index:
    input:
        fasta = config["assembly_path"],
        gtf = "data/gencode_gtf/gencode.v49.primary_assembly.gtf"
    output:
        index=directory("data/star_genome_index")
    params:
        sjdbOverhang=config.get("read_length", 100) - 1,
        genomeSAindexNbases=config.get("genomeSAindexNbases", 14)
    threads: 8
    resources:
        mem_mb=35000,
        runtime = 180
    log:
        "logs/star/index_genome.log"
    shell:
        """
        mkdir -p {output.index}
        
        STAR --runMode genomeGenerate \
            --genomeDir {output.index} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.sjdbOverhang} \
            --genomeSAindexNbases {params.genomeSAindexNbases} \
            --runThreadN {threads} \
        >> {log} 2>&1
        """ 

rule star:
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
            --outFilterMismatchNoverReadLmax {params.max_perread_sub_fraction} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes nM NM MD AS \
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
        # GEDI INDEX FILES
        index = rules.gedi_index_genome.output
    output:
        cit_sample_set = "data/cit_sample_sets/{donor}.cit",
    params:
        scratch = lambda wildcards: f"/scratch/midway3/$USER/{wildcards.donor}_$SLURM_JOB_ID",
    container:
        config["container_path"]
    resources:
        mem = "20G",
        runtime = 540 # 9 hours in minutes
    benchmark:
        "benchmarks/{donor}.bam_to_cit.benchmark.txt"
    shell:
        # gedi -e Bam2CIT -p  data/cit/s17_s18_test.cit data/aligned_bam/LB-HT-28s-HT-17_S17.bam data/aligned_bam/LB-HT-28s-HT-18_S18.bam
        # does not run multithreaded; speeds are the same
        """
        echo "Copying to scratch..."
        cp {input.bams} {input.bais} {input.no4sU_bam} {input.no4sU_bai} {input.index} {params.scratch}/
        echo "Copy complete."
        cd {params.scratch}
        echo "Beginning CIT conversion..."
        gedi -e Bam2CIT -p output.cit $(basename {input.bams}) $(basename {input.no4sU_bam})
        echo "Conversion complete."
        echo "Copying output back..."
        cp output.cit {output.cit_sample_set}
        """

rule grand_slam:
    input: 
        cit_sample_set = "data/cit_sample_sets/{donor}.cit",
        index_oml = rules.gedi_index_genome.output.oml
    output: 
        nascent_counts = "data/slam_quant/{donor}/grandslam.tsv",
        mismatches_tsv = "data/slam_quant/{donor}/grandslam.mismatches.tsv",
        mismatches_pdf = "data/slam_quant/{donor}/grandslam.mismatches.pdf",
        mismatch_positions = "data/slam_quant/{donor}/grandslam.mismatchpos.pdf",
        ntr_stats = "data/slam_quant/{donor}/grandslam.ntrstat.tsv",
    container:
        config['container_path']
    resources:
       runtime = 360,
       mem = "16G",
    threads:
        16,
    benchmark:
        "benchmarks/{donor}.grandslam.benchmark.txt"
    shell: 
        """
            gedi -e Slam \
            -genomic {input.index_oml} \
            -reads {input.cit_sample_set} \
            -prefix data/slam_quant/{wildcards.donor}/grand_slam \
            -nthreads {threads} \
            -introns \
            -no4sUpattern no4sU \
            -progress \
            -full \
            -plot
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