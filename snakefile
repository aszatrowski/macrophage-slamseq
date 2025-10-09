data_path = 'data'
# sample_ids = [f.removesuffix('.fastq.gz') for f in os.listdir(f'{data_path}/') if f.endswith('.fastq.gz')]
sample_ids = ['LB-HT-28s-HT-10_S10_L007', 'LB-HT-28s-HT-02_S2_L006', 'LB-HT-28s-JL-08_S26_L007']
n_tag = ['001']
rule all:
    input: 
        expand(
            "outputs/fastp_reports/{sample_id}_{n}.html",
            sample_id = sample_ids,
            n = n_tag
        ),
        expand(
            "outputs/fastp_reports/{sample_id}_{n}.json",
            sample_id = sample_ids,
            n = n_tag
        ),
        expand(
            'data/hisat3n_indexes/hg38.3n.{converted_bases}.{index_n}.ht2',
            converted_bases = ['CT', 'GA'], # using a C>T conversion + reverse complement
            index_n = range(1, 9)
        ),
        expand(
            'data/aligned_bam/{sample_id}_aligned.bam',
            sample_id = sample_ids
        )

rule process_fastp:
    input: # n is necessary so all parts have the same wildcards 
        r1 = f'{data_path}/fastq_symlinks/{{sample_id}}_R1_{{n}}.fastq.gz',
        r2 = f'{data_path}/fastq_symlinks/{{sample_id}}_R2_{{n}}.fastq.gz'
    output:
        r1 = f'{data_path}/trimmed/{{sample_id}}_R1_{{n}}.fastq.gz',
        r2 = f'{data_path}/trimmed/{{sample_id}}_R2_{{n}}.fastq.gz',
        html = "outputs/fastp_reports/{sample_id}_{n}.html",
        json = "outputs/fastp_reports/{sample_id}_{n}.json"
    log:
        "logs/fastp/{sample_id}_{n}.log"
    threads: 4
    params:
        extra = "--detect_adapter_for_pe"
    resources: 
        slurm_account = 'pi-lbarreiro',
        runtime = 30
    shell:
        """
        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --html {output.html} \
            --json {output.json} \
            --thread {threads} \
            {params.extra} \
            2> {log}
        """

rule build_hisat3n_index:
    input: 
        assembly = '/project/lbarreiro/SHARED/REFERENCES/Homo_sapiens/GATK/GRCh38/GRCh38.primary_assembly.genome.fa'
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
        f"""
        hisat-3n-build \
        --base-change T,C \
        {input.assembly} \
        {data_path}/hisat3n_indexes/hg38 \
        """
    
rule align_hisat3n:
    input: 
        fastq_r1 = 'data/trimmed/{sample_id}_R1_001.fastq.gz',
        fastq_r2 = 'data/trimmed/{sample_id}_R2_001.fastq.gz',
        index = expand(
            'data/hisat3n_indexes/hg38.3n.{converted_bases}.{index_n}.ht2',
            converted_bases = ['CT', 'GA'], # using a C>T conversion + reverse complement
            index_n = range(1, 9)
        )
    output:
        aligned_sam = 'data/aligned_sam/{sample_id}_aligned.sam'
    threads: 8
    resources: 
        slurm_account = 'pi-lbarreiro',
        runtime = 120,
        mem_mb = 24000
    shell: 
        (
            "hisat-3n "
            "-p {threads} " # num threads
            "-x data/hisat3n_indexes/hg38 "
            "--base-change T,C "
            "--rna-strandness RF " # not sure what this means, but it's in JL's pipeline
            "-q " # for .fastq, not .fasta
            "-1 {input.fastq_r1} -2 {input.fastq_r2} "
            "-S  {output.aligned_sam}"
        )

rule sam_to_bam:
    input: 
        samfile = "data/aligned_sam/{sample_id}_aligned.sam"
    output: 
        bamfile = "data/aligned_bam/{sample_id}_aligned.bam"
    resources: 
        slurm_account = 'pi-lbarreiro',
        runtime = 20,
        mem_mb = 16000
    shell: 
        """
        samtools view -h -bS -F 260 {input.samfile} | samtools sort > {output.bamfile}
        """