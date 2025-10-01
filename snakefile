data_path = 'data'
run_ids = ['LB-HT-28s-HT-01_S1_L005_R1_001']
fastqc_exts = ['html', 'zip']

rule all:
    input: 
        expand(
            "outputs/fastqc_reports/fastqc_{run_id}.qc.{ext}",
            run_id = run_ids,
            ext = fastqc_exts
        ) 


rule fastqc:
    input:  f"{data_path}/{{run_id}}.fastq.gz"
    output: "outputs/fastqc_reports/fastqc_{run_id}.qc.{ext}"
    resources: 
        slurm_account='pi-lbarreiro',
        runtime='5'
    shell: 
        """
        fastqc -o outputs/fastqc_reports/ -f fastq {input}
        """

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