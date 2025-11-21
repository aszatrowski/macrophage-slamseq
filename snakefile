## CONFIG:
configfile: "config.yaml"
# these are so lightweight that they can be run directly on the login node; no need for slurm
# will want to add figure generation to this
localrules: cat_fastqs, index_bam, multiqc

## CHOOSE FILES
# find all sample files in the folder, and remove _R1_001 & _R2_001 since paired-end reads will be processed together
# will need to update this to also get rid of lanes
# sample_ids = [f.removesuffix('_R1_001.fastq.gz').removesuffix('_R2_001.fastq.gz') 
#               for f in os.listdir('data/fastq_symlinks') 
#               if f.endswith('.fastq.gz')][0:2]

sample_ids = [
    'LB-HT-28s-HT-01_S1',
    'LB-HT-28s-HT-02_S2',
    'LB-HT-28s-HT-03_S3',
    'LB-HT-28s-HT-05_S5',
    'LB-HT-28s-HT-06_S6',
    'LB-HT-28s-HT-07_S7',
    'LB-HT-28s-HT-08_S8',
    'LB-HT-28s-HT-09_S9',
    'LB-HT-28s-HT-10_S10',
    'LB-HT-28s-HT-12_S12',
    'LB-HT-28s-HT-16_S16',
    'LB-HT-28s-HT-17_S17',# smallest fileset (collectively R1 6.5KB + R2 6.4KB); use this as a test BUT has no nascent transcripts and fails featureCounts
    'LB-HT-28s-HT-18_S18',
    'LB-HT-28s-JL-01_S19',
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
        # expand(
        #     "data/cit/{sample_id}.cit",
        #     sample_id = sample_ids
        # ),
        expand(
            "{sample_id}/slam_quant.tsv", # {prefix}.tsv
            sample_id = 'LB-HT-28s-HT-17_S17'
        ),
        'outputs/multiqc_report.html'

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

rule wget_gencode_gtf:
    output: 
        temp("data/gencode_gtf/gencode.v49.primary_assembly.gtf")
    params:
        refseq_path = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz",
        folder = 'data/gencode_gtf'
    shell: 
        """
        wget {params.refseq_path} \
            --output-document {params.folder}/gencode.v49.primary_assembly.gtf.gz  \
            --force-directories \
            --show-progress
        gunzip {params.folder}/gencode.v49.primary_assembly.gtf.gz
        """

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
        aligned_bam = 'data/aligned_bam/{sample_id}.bam',
        stats = 'logs/star/qc/{sample_id}.Log.final.out'
    params:
        # turns out these operations are MASSIVELY IO limited, so we'll copy everything to scratch, where it won't compete for RW access
        scratch = lambda wildcards: f"/scratch/midway3/$USER/{wildcards.sample_id}_$SLURM_JOB_ID",
        max_perread_sub_fraction = 0.5
    log:
        "logs/star/{sample_id}.log" # alignment quality
    benchmark:
        "benchmarks/{sample_id}.star_align.benchmark.txt"
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
        gtf = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf",
        genes_tab = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf.genes.tab",
        index_cit = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf.index.cit",
        transcripts_tab = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf.transcripts.tab",
        transcripts_fi = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf.transcripts.fi",
        transcripts_fasta = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf.transcripts.fasta",
    container:
        config['container_path']
    resources:
        mem_mb= "8G",
        runtime = 5 
    params:
        output_dir = config['gedi_index_dir']
    log:
        "logs/gedi/index_homo_sapiens.log"
    shell:
        """
            mkdir -p config/genomic
            gedi \
                -e IndexGenome \
                -organism homo_sapiens \
                -version 115 \
                -p \
                -f config/genomic \
                -nobowtie -nostar -nokallisto 
            2&1> {log}
            # mv homo_sapiens.115.fasta homo_sapiens.115.fi homo_sapiens.115.gtf homo_sapiens.115.gtf.genes.tab homo_sapiens.115.gtf.index.cit homo_sapiens.115.gtf.transcripts.tab homo_sapiens.115.gtf.transcripts.fi homo_sapiens.115.gtf.transcripts.fasta config/genomic
        """

# rule make_bamlist:
#     input: 
#         bams = lambda wildcards: expand(
#             "data/aligned_bam/{sample_id}.bam",
#             sample_id = get_donor_samples(wildcards.donor)
#         )
#     output: 
#         "data/bamlist_cit/{donor}.bamlist"
#     shell: 
#         "printf '%s\\n' {input} > {output}"

rule index_bam:
    input: 
        "data/aligned_bam/{sample_id}.bam"
    output: 
        "data/aligned_bam/{sample_id}.bam.bai"
    shell: 
        """
        samtools index --bai {input} -o {output}
        """

rule bam_to_cit:
    input:
        bam = "data/aligned_bam/{sample_id}.bam",
        bai = "data/aligned_bam/{sample_id}.bam.bai",
        # GEDI index files
        fasta = f"{config['gedi_index_dir']}/homo_sapiens.115.fasta",
        fi = f"{config['gedi_index_dir']}/homo_sapiens.115.fi",
        gtf = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf",
        genes_tab = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf.genes.tab",
        index_cit = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf.index.cit",
        transcripts_tab = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf.transcripts.tab",
        transcripts_fi = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf.transcripts.fi",
        transcripts_fasta = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf.transcripts.fasta",
    output: 
        "data/cit/{sample_id}.cit"
    container:
        config['container_path']
    resources:
       runtime = 15,
       mem = "8G",
    shell: # converts aligned and tagged .bam files to GRAND-SLAM's custom CIT format
        """
        gedi -e Bam2CIT -p {output} {input}
        """

rule grand_slam:
    input: 
        sample_cit = "data/cit/{sample_id}.cit",
        fasta = f"{config['gedi_index_dir']}/homo_sapiens.115.fasta",
        fi = f"{config['gedi_index_dir']}/homo_sapiens.115.fi",
        gtf = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf",
        genes_tab = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf.genes.tab",
        index_cit = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf.index.cit",
        transcripts_tab = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf.transcripts.tab",
        transcripts_fi = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf.transcripts.fi",
        transcripts_fasta = f"{config['gedi_index_dir']}/homo_sapiens.115.gtf.transcripts.fasta",
    output: 
        slam_counts = "{sample_id}/slam_quant.tsv" # {prefix}.tsv
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
        """
        gedi -e Slam \
            -genomic {input.index_cit} \
            {input.sample_cit} \
            -progress \
            -prefix {wildcards.sample_id}/slam_quant
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
            sample_id = sample_ids
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