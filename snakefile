data_path = 'data'
# sample_ids = [f.removesuffix('.fastq.gz') for f in os.listdir(f'{data_path}/') if f.endswith('.fastq.gz')]
sample_ids = ['LB-HT-28s-HT-10_S10_L007_R1_001', 'LB-HT-28s-HT-02_S2_L006_R1_001', 'LB-HT-28s-JL-08_S26_L007_R1_001']
rule all:
    input: 
        expand(
            'outputs/fastqc_reports/{sample_id}_fastqc.zip',
            sample_id = sample_ids
        ),
        expand(
            'outputs/fastqc_reports/{sample_id}_fastqc.html',
            sample_id = sample_ids
        ),
        'outputs/fastqc_summary/fastqc_summary.csv',
        'outputs/fastqc_summary/fastqc_fails_warnings.csv'

rule fastqc:
    input:  f'{data_path}/{{sample_id}}.fastq.gz'
    output:
        'outputs/fastqc_reports/{sample_id}_fastqc.zip',
        'outputs/fastqc_reports/{sample_id}_fastqc.html'
    resources: 
        slurm_account='pi-lbarreiro',
        runtime=18
    shell: 
        """
        fastqc -o outputs/fastqc_reports/ -f fastq {input}
        """

rule summarize_fastqc:
    input: 
        expand(
            'outputs/fastqc_reports/{sample_id}_fastqc.html',
            sample_id = sample_ids
        ) 
    output: 
        summary = 'outputs/fastqc_summary/fastqc_summary.csv',
        fails_warnings = 'outputs/fastqc_summary/fastqc_fails_warnings.csv'
    resources: 
        slurm_account='pi-lbarreiro',
        runtime='20'
    script: "scripts/fastqc_summary.R"

rule test_cluster:
    output:
        "test_output.txt"
    resources: 
        slurm_account='pi-lbarreiro',
        runtime='1'
    shell:
        """
        echo "Cluster execution works!" > {output}
        echo "Hostname: $(hostname)" >> {output}
        echo "Date: $(date)" >> {output}
        echo "Job completed successfully" >> {output}
        """