import os
storage:
    provider="http",
    retrieve=False

localrules: decompress_kallisto_index_tar, wget_kallisto_index_tar

slurm_account = 'pi-lbarreiro'
assembly_path = '/project/lbarreiro/SHARED/REFERENCES/Homo_sapiens/GATK/GRCh38/GRCh38.primary_assembly.genome.fa'
kallisto_idx_path = 'https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/human_index_nac.tar.xz'
data_path = 'data'
# sample_ids = [f.removesuffix('.fastq.gz') for f in os.listdir(f'{data_path}/') if f.endswith('.fastq.gz')]
# sample_ids = ['LB-HT-28s-HT-10_S10_L007', 'LB-HT-28s-HT-02_S2_L006', 'LB-HT-28s-JL-08_S26_L007']
sample_ids = ['LB-HT-28s-HT-10_S10_L007']
n_tag = ['001']

os.makedirs('data/samtools_temp', exist_ok=True) # for some reason samtools refuses to create its own dirs

rule all:
    input: 
        expand(
            "outputs/fastp_reports/{sample_id}_{n}.{filetype}",
            sample_id = sample_ids,
            n = n_tag,
            filetype = ['html', 'json']
        ),
        expand(
            'data/aligned_bam/{sample_id}_aligned.bam',
            sample_id = sample_ids
        ),
        expand(
            'outputs/kallisto/{sample_id}_counts.tsv',
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
        aligned_sam = temp('data/aligned_sam_temp/{sample_id}_aligned.sam')
    threads: 16
    resources: 
        slurm_account = 'pi-lbarreiro',
        runtime = 180,
        mem_mb = 24000
    shell: 
        (
            "hisat-3n "
            "-p {threads} " # num threads inherited from threads:
            "-x data/hisat3n_indexes/hg38 " # hisat-3n index to align to
            "--base-change T,C " # don't penalize T>C substitutions from SLAM-seq
            "--rna-strandness RF " # not sure what this means, but it's in JL's pipeline
            "-q " # for .fastq, not .fasta
            "-1 {input.fastq_r1} -2 {input.fastq_r2} " # r1 and r2 paired end files
            "-S  {output.aligned_sam}"
        )

rule sam_to_bam:
    input: 
        samfile = "data/aligned_sam_temp/{sample_id}_aligned.sam"
    output: 
        bamfile = "data/aligned_bam/{sample_id}_aligned.bam"
    threads: 4
    resources: 
        slurm_account = 'pi-lbarreiro',
        runtime = 30, # can bring down to ~30
        mem_mb = 16000
    shell: 
        (
            "samtools view "
            "--with-header " # with header
            "--bam " # output to .bam; -S is meaningless legacy
            "--exclude-flags 260 " # QC filter—excludes multimapped reads I believe
            "-@ {threads} " # num threads inherited from threads:
            "{input.samfile} "
            "| samtools sort "
            "-T data/samtools_temp/ "
            "> {output.bamfile}"
        )

rule wget_kallisto_index_tar:
    input:
        path = storage.http('https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/human_index_nac.tar.xz')
    output: 
        kallisto_index_tar = 'data/homo_sapiens_kallisto_index.tar.xz'
    shell:
        "wget -v -O {output.kallisto_index_tar} {input.path} "

rule decompress_kallisto_index_tar:
    input:
        kallisto_index_tar = 'data/homo_sapiens_kallisto_index.tar.xz'
    output: 
        kallisto_cdna = 'data/homo_sapiens_kallisto_index/cdna.txt',
        kallisto_index = 'data/homo_sapiens_kallisto_index/index.idx',
        kallisto_nascent = 'data/homo_sapiens_kallisto_index/nascent.txt',
        kallisto_t2g = 'data/homo_sapiens_kallisto_index/t2g.txt'
    shell:
        "tar -xvf {input.kallisto_index_tar} -C data/homo_sapiens_kallisto_index"

rule kallisto_quant:
    input:
        fastq_r1 = 'data/trimmed/{sample_id}_R1_001.fastq.gz',
        fastq_r2 = 'data/trimmed/{sample_id}_R2_001.fastq.gz',
        kallisto_index = 'data/homo_sapiens_kallisto_index/index.idx'
    output:
        kallisto_counts = 'outputs/kallisto/{sample_id}_counts.tsv'
    threads: 4
    resources: 
        slurm_account = slurm_account,
        runtime = 60,
        mem_mb = 16000
    shell: 
        (
            "kallisto quant "
            "-i {input.kallisto_index} "
            "-o {output.kallisto_counts} "
            "{input.fastq_r1} {input.fastq_r2} "
            "--plaintext " # output to human-readable .tsv not HDF5 binary
            "--threads {threads} "
            "--rf-stranded"
        )