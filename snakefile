## CONFIG:
configfile: "config.yaml"
# these are so lightweight that they can be run directly on the login node; no need for slurm
# will want to add figure generation to this
localrules: cat_fastqs, index_bam, mark_control_samples, multiqc

LANES = [5, 6, 7, 8]
DONORS = ["donor1_rep2"]
sample_ids = config['donor_sample_ids']['donor1_rep2']

rule all:
    input: 
        # lambda w: expand("data/aligned_bam/{sample_id}.bam",
        #                              sample_id=[s for s in config["donor_sample_ids"]["donor1_rep2"] 
        #                                     if s != config["control_sample_ids"]["donor1_rep2"]]),
        # lambda w: expand("data/aligned_bam/{sample_id}.bam.bai",
        #                              sample_id=[s for s in config["donor_sample_ids"]["donor1_rep2"] 
        #                                     if s != config["control_sample_ids"]["donor1_rep2"]]),
        # # Control sample_id (with .no4sU suffix)
        # lambda w: f"data/no4sU_tagged/{config['control_sample_ids']["donor1_rep2"]}.no4sU.bam",
        # lambda w: f"data/no4sU_tagged/{config['control_sample_ids']["donor1_rep2"]}.no4sU.bam.bai",
        expand(
            "data/cit_sample_sets/{donor}.cit",
            donor = DONORS
        )
        # 'outputs/multiqc_report.html'

rule cat_fastqs:
    input: 
        lambda wildcards: expand(
            f"{config['data_path']}/{{sample_id}}_L{{seq_lane}}_R{{end}}_001.fastq.gz",  # ← double braces
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

rule star_align:
    input: 
        fastq_r1 = 'data/trimmed/{sample_id}_R1_001.fastq.gz',
        fastq_r2 = 'data/trimmed/{sample_id}_R2_001.fastq.gz',
        index = '/project/lbarreiro/SHARED/REFERENCES/Homo_sapiens/GATK/GRCh38/STAR'
    output:
        aligned_bam = 'data/aligned_bam_{max_perread_sub_fraction}/{sample_id}.bam',
        stats = 'logs/star/qc/{sample_id}_{max_perread_sub_fraction}.Log.final.out'
    params:
        # turns out these operations are MASSIVELY IO limited, so we'll copy everything to scratch, where it won't compete for RW access
        scratch = lambda wildcards: f"/scratch/midway3/$USER/{wildcards.sample_id}_$SLURM_JOB_ID",
        # max_perread_sub_fraction = 0.25
    log:
        "logs/star/{sample_id}_{max_perread_sub_fraction}.log" # alignment quality
    benchmark:
        "benchmarks/{sample_id}.star_align_{max_perread_sub_fraction}.benchmark.txt"
    threads: 8
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
            --outFilterMismatchNoverReadLmax {wildcards.max_perread_sub_fraction} \
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

rule mark_control_samples:
    # while this rule appears general, it is only required for the no4sU samples as called in the second step
    input:
        bam="data/aligned_bam/{sample_id}.bam",
        bai="data/aligned_bam/{sample_id}.bam.bai"
    output:
        bam="data/no4sU_tagged/{sample_id}.no4sU.bam",
        bai="data/no4sU_tagged/{sample_id}.no4sU.bam.bai"
    shell:
        """
        mkdir -p data/no4sU_tagged
        ln -sf $(realpath {input.bam}) {output.bam}
        ln -sf $(realpath {input.bai}) {output.bai}
        """

rule bam_to_cit:
    input:
        bams_regular = lambda w: expand(
            "data/aligned_bam/{sample_id}.bam",
            sample_id=[s for s in config["donor_sample_ids"][w.donor] 
                if s != config["control_sample_ids"][w.donor]]
        ),
        bais_regular = lambda w: expand(
            "data/aligned_bam/{sample_id}.bam.bai",
            sample_id=[s for s in config["donor_sample_ids"][w.donor] 
                if s != config["control_sample_ids"][w.donor]]
        ),
        bam_control = lambda w: f"data/no4sU_tagged/{config['control_sample_ids'][w.donor]}.no4sU.bam",
        bai_control = lambda w: f"data/no4sU_tagged/{config['control_sample_ids'][w.donor]}.no4sU.bam.bai",
        # GEDI INDEX FILES
        index = rules.gedi_index_genome.output
    output:
        cit_sample_set = "data/cit_sample_sets/{donor}.cit"
    container:
        config["container_path"]
    resources:
        mem = "8G",
        runtime = 60
    shell:
        """
        Bam2CIT \
            {output.cit_sample_set} \
            {input.bams_regular} {input.bam_control} \
        """

rule grand_slam:
    input: 
        sample_cit = "data/cit/{sample_id}.cit",
        index = rules.gedi_index_genome.output
    output: 
        binom = "data/slam_quant/{sample_id}/grandslam.binom.tsv",
        doublehit = "data/slam_quant/{sample_id}/grandslam.doublehit.tsv",
        mismatches = "data/slam_quant/{sample_id}/grandslam.mismatches.tsv",
        # -full outputs
        # mismatches = "{sample_id}/slam_quant.mismatches.pdf",
        # double = "{sample_id}/slam_quant.double.pdf",
        # mismatchpos = "{sample_id}/slam_quant.mismatchpos.pdf",
        # mismatchposzoomed = "{sample_id}/slam_quant.mismatchposzoomed.pdf",
        # binomRates = "{sample_id}/slam_quant.binomRates.png",
        # mismatches = "{sample_id}/slam_quant.mismatches.pdf",
        # mismatches = "{sample_id}/slam_quant.mismatches.pdf",
        # mismatches = "{sample_id}/slam_quant.mismatches.pdf"
    container:
        config['container_path']
    resources:
       runtime = 30,
       mem = "8G",
    log:
        "logs/gedi/slam/{sample_id}.log"
    shell: 
        f"""
            gedi -e Slam \
            -genomic {{input.ref_oml}} \
            -reads data/aligned_bam/{{input.sample_cit}}
            -prefix data/slam_quant/LB-HT-28s-HT-18_S18/grandslam
            -nthreads {{threads}}
            -introns
            -progress
            # add -no4sU pattern
        """
rule multiqc:
    input: 
        expand(
            "data/fastp_reports/{sample_id}.{filetype}",
            sample_id = sample_ids,
            filetype = ['json']
        ),
        expand(
            'logs/star/qc/{sample_id}_{max_perread_sub_fraction}.Log.final.out',
            sample_id = sample_ids,
            max_perread_sub_fraction = [0.25, 0.1]
        )
    output: 
        'outputs/multiqc_report.html'
    shell: 
        (
            'multiqc '
            'data/fastp_reports logs/star ' # switch these to ALL for consistency in multiqc
            '--force ' # overwrite existing report; otherwise it will attach a suffix that snakemake won't detect
            '--outdir outputs'
            # future: --ignore-samples for ones that failed to process
        )