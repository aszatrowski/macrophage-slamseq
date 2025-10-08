data_path = 'data'
# sample_ids = [f.removesuffix('.fastq.gz') for f in os.listdir(f'{data_path}/') if f.endswith('.fastq.gz')]
sample_ids = ['LB-HT-28s-HT-10_S10_L007', 'LB-HT-28s-HT-02_S2_L006', 'LB-HT-28s-JL-08_S26_L007']
n_tag = ['001']
rule all:
    input: 
        expand(
            f'{data_path}/trimmed/{{sample_id}}_R1_{{n}}.fastq.gz',
            sample_id = sample_ids,
            n = n_tag
        ),
        expand(
            f'{data_path}/trimmed/{{sample_id}}_R2_{{n}}.fastq.gz',
            sample_id = sample_ids,
            n = n_tag
        ),
        expand(
            "outputs/fastp_reports/{sample_id}_{n}.html",
            sample_id = sample_ids,
            n = n_tag
        ),
        expand(
            "outputs/fastp_reports/{sample_id}_{n}.json",
            sample_id = sample_ids,
            n = n_tag
        )

rule process_fastp:
    input:  
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
