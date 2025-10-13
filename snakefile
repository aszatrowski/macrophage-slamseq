import os
storage: # configure HTTP call to pre-built hg38 kallisto index
    provider="http",
    retrieve=False,
    keep_local=False
kallisto_idx_path = 'https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/human_index_nac.tar.xz'

# path to lab's shared hg38 .fa index for hisat3n indexing
assembly_path = '/project/lbarreiro/SHARED/REFERENCES/Homo_sapiens/GATK/GRCh38/GRCh38.primary_assembly.genome.fa'
# pi account for slurm
slurm_account = 'pi-lbarreiro'

# these are so lightweight that they can be run directly on the login node; no need for slurm
# will want to add figure generation to this
localrules: decompress_kallisto_index_tar, wget_kallisto_index_tar, generate_tagvalues_file
# find all sample files in the folder, and remove _R1_001 & _R2_001 since paired-end reads will be processed together
sample_ids = [f.removesuffix('_R1_001.fastq.gz').removesuffix('_R2_001.fastq.gz') 
              for f in os.listdir('data/fastq_symlinks') 
              if f.endswith('.fastq.gz')][0:40]
print(len(sample_ids))

os.makedirs('data/samtools_temp', exist_ok=True) # for some reason samtools refuses to create its own dirs

rule all:
    input: 
        expand(
            "data/fastp_reports/{sample_id}_001.{filetype}",
            sample_id = sample_ids,
            filetype = ['html', 'json']
        ),
        expand(
            'data/kallisto/{sample_id}/abundance.{ext}',
            sample_id = sample_ids,
            ext = ['tsv', 'h5']
        ),
        expand(
            'data/kallisto/{sample_id}/run_info.json',
            sample_id = sample_ids
        ),
        expand(
            "data/nascent_counts/{sample_id}_nascent_counts.bam",
            sample_id = sample_ids
        ),
rule process_fastp:
    input: 
        r1 = 'data/fastq_symlinks/{sample_id}_R1_001.fastq.gz',
        r2 = 'data/fastq_symlinks/{sample_id}_R2_001.fastq.gz'
    output:
        r1 = 'data/trimmed/{sample_id}_R1_001.fastq.gz',
        r2 = 'data/trimmed/{sample_id}_R2_001.fastq.gz',
        html = "data/fastp_reports/{sample_id}_001.html",
        json = "data/fastp_reports/{sample_id}_001.json"
    log:
        "logs/fastp/{sample_id}_001.log"
    threads: 4
    resources: 
        slurm_account = 'pi-lbarreiro',
        runtime = 30
    shell:
        (
            "fastp "
            "--in1 {input.r1} --in2 {input.r2} "
            "--out1 {output.r1} --out2 {output.r2} "
            "--html {output.html} "
            "--json {output.json} "
            "--thread {threads} "
            "--detect_adapter_for_pe " # dynamically detect adapter sequences
            "2> {log}"
        )

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
        (
            "hisat-3n-build "
            "--base-change T,C "
            "{input.assembly} "
            "data/hisat3n_indexes/hg38 "
        )
    
rule align_hisat3n:
    input: 
        fastq_r1 = 'data/trimmed/{sample_id}_R1_001.fastq.gz',
        fastq_r2 = 'data/trimmed/{sample_id}_R2_001.fastq.gz',
        index = expand(
            'data/hisat3n_indexes/hg38.3n.{converted_bases}.{index_n}.ht2',
            converted_bases = ['CT', 'GA'], # using a T>C conversion + reverse complement
            index_n = range(1, 9)
        )
    output:
        aligned_sam = temp('data/aligned_sam_temp/{sample_id}_aligned.sam')
    threads: 12
    resources: 
        slurm_account = 'pi-lbarreiro',
        runtime = 210,
        mem = "32G"
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
        runtime = 45, # can bring down to ~30
        mem = "16G"
    shell: 
        (
            "samtools view "
            "--with-header " # with header
            "--bam " # output to .bam; -S is meaningless legacy
            "--exclude-flags 260 " # QC filter—excludes multimapped reads I believe
            "-@ {threads} " # num threads inherited from threads:
            "{input.samfile} "
            "| samtools sort " # pipe to samtools sort
            "-T data/samtools_temp/ " # use the created directory for temp files (many and large, ew)
            "> {output.bamfile}"
        )

rule wget_kallisto_index_tar:
    input:
        path = storage.http('https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/human_index_nac.tar.xz')
    output: 
        kallisto_index_tar = temp('data/homo_sapiens_kallisto_index.tar.xz')
    shell:
        "wget -v -O {output.kallisto_index_tar} {input.path}"

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
        kallisto_counts = 'data/kallisto/{sample_id}/abundance.tsv',
        kallisto_h5 = 'data/kallisto/{sample_id}/abundance.h5',
        kallisto_run_info = 'data/kallisto/{sample_id}/run_info.json'
    threads: 5
    resources: 
        slurm_account = slurm_account,
        runtime = 40,
        mem = "32G"
    shell: 
        (
            "kallisto quant "
            "-i {input.kallisto_index} "
            "-o data/kallisto/{wildcards.sample_id} "
            "{input.fastq_r1} {input.fastq_r2} "
            "--threads {threads} "
            "--rf-stranded"
        )

rule generate_tagvalues_file:
    output: 
        tags = "data/tag_values.txt"
    run: 
        with open(output.tags, 'w') as f: # write sequence to .txt, one int per line
            for i in range(2, 301):
                f.write(f'{i}\n')

rule count_nascent_transcripts:
    input: 
        bamfile = "data/aligned_bam/{sample_id}_aligned.bam",
        tags = "data/tag_values.txt" # file with integers 2-300 against which samtools will match Yf:
    output: 
        nascent_counts = "data/nascent_counts/{sample_id}_nascent_counts.bam"
    resources:
        slurm_account = slurm_account,
        runtime = 5,
        mem = "8G"
    shell: 
        (
            "samtools view "
            "--with-header " 
            "--bam "
            "-F 260 " # exclude 260 field—multiple alignments
            "-D Yf:{input.tags} " # only alignments with Yf: tag == STR where STR in input.tags
            "{input.bamfile} > {output.nascent_counts}"
        )