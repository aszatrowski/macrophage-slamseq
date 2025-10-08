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
            f'{data_path}/hisat3n_indexes/{{sample_id}}_{{R}}_{{n}}.3n.{{converted_bases}}.{{index_n}}.ht2',
            R = ['R1', 'R2'],
            index_n = range(8),
            converted_bases = ['CT', 'GA'], # using a C>T conversion + reverse complement
            sample_id = sample_ids,
            n = n_tag
        )

rule process_fastp:
    input: # n is necessary so all parts have the same wildcards 
        r1 = f'{data_path}/{{sample_id}}_R1_{{n}}.fastq.gz',
        r2 = f'{data_path}/{{sample_id}}_R2_{{n}}.fastq.gz'
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
        slurm_account='pi-lbarreiro',
        runtime=60
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
        f'{data_path}/trimmed/{{sample_id}}_{{R}}_{{n}}.fastq.gz',
    output: 
        f'{data_path}/hisat3n_indexes/{{sample_id}}_{{R}}_{{n}}.3n.{{converted_bases}}.{{index_n}}.ht2',
    resources: 
        slurm_account='pi-lbarreiro',
        runtime=60,
        mem_mb=16000,
    shell:  
        f"""
        hisat-3n-build --base-change T,C {{input}} {data_path}/hisat3n_indexes/{{wildcards.sample_id}}_{{wildcards.R}}_{{wildcards.n}}
        """
    
rule align_hisat3n:
    input: 
    output: 
    resources: 
        slurm_account='pi-lbarreiro',
        runtime=60
    shell: 
        """
        hisat-3n \
        -x $indexfile \
        -p 5 \
        -q \ # for .fastq, not .fasta
        -1 $R1 -2 $R2 \
        -S  "$dir_out"/"$file_out_sam" --base-change T,C --rna-strandness RF
        """
