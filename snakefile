## CONFIG:
configfile: "config.yaml"
# these jobs are so lightweight that they can be run directly on the login node; no need for slurm or compute nodes
localrules: cat_fastqs, index_bam, rename_with_donor_timepoint, multiqc, calc_nascent_total_reads, merge_reads_across_donors, dge, map_ensg_genesymbol, volcano_plot, timecourse_dge_plot, effect_size_correlations_plot, venn_diagram_plot

DONORS = ['donor1_rep1', 'donor1_rep2']
sample_ids = list(config['sample_ids'].keys())

rule all:
    input: 
        expand(
            "outputs/venn_diagrams/venn_{comparison}.png",
            comparison = ['15_vs_0m','30_vs_0m','60_vs_0m','90_vs_0m','105_vs_0m','120_vs_0m']
        ),
        expand(
            "outputs/effect_size_correlations/corr_{comparison}.pdf",
            comparison = ['15_vs_0m','30_vs_0m','60_vs_0m','90_vs_0m','105_vs_0m','120_vs_0m']
        ),
        expand(
            "outputs/timecourse_dge_plot_{readtype}.pdf",
            readtype = ['total', 'nascent'],
        ),
        expand(
            "outputs/volcano_plots/volcanoplot_{readtype}-{comparison}.pdf",
            readtype = ['total', 'nascent'],
            comparison = ['15_vs_0m','30_vs_0m','60_vs_0m','90_vs_0m','105_vs_0m','120_vs_0m']
        ),
        expand(
            "outputs/dge_results/summary_stats_{readtype}.csv",
            readtype = ['total', 'nascent'],
        ),
        'outputs/multiqc_report.html'

rule cat_fastqs:
    """
    Joins fastq files together from multiple sequencing lanes. Lanes are specified in config['sequencing_lanes'].
    """
    input: 
        lambda wildcards: expand(
            f"{config['data_path']}/{{sample_id}}_L{{seq_lane}}_R{{end}}_001.fastq.gz",
            sample_id=wildcards.sample_id,
            end=wildcards.end,
            seq_lane=[f"{l:03d}" for l in config['sequencing_lanes']]
        )
    output: 
        temp('data/fastq_merged/{sample_id}_R{end}_001.fastq.gz')
    shell: 
        "cat {input} > {output}"

rule fastp:
    """
    Runs fastp for quality control and adapter trimming. fastp (https://github.com/OpenGene/fastp) is a faster alternative to fastqc + trimmomatic/cutadapt, and crucially only needs to be run once, rather than fastqc > trim > fastqc. Though fastp can dynamically detect adapter sequences quite efficiently, in practice it misses bases on the ends, which introduces artificial mismatches (adapter contamination), which may reduce T>C substitution calling downstream, so manual specification is recommended.
    The parameters are configurd for paired-end sequencing, so you'll need to modify if you have single-end. See the fastp docs for details.
    Outputs an html webpage and a machine-readable json QC summary for each sample; multiqc will take the json as input for its own summary.
    """
    input: 
        r1 = 'data/fastq_merged/{sample_id}_R1_001.fastq.gz',
        r2 = 'data/fastq_merged/{sample_id}_R2_001.fastq.gz'
    output:
        r1 = temp('data/trimmed/{sample_id}_R1_001.fastq.gz'),
        r2 = temp('data/trimmed/{sample_id}_R2_001.fastq.gz'),
        html = "data/fastp_reports/{sample_id}.html",
        json = "data/fastp_reports/{sample_id}.json"
    params:
        adapter_sequence_r1 = config['adapter_sequence_r1_rc'],
        adapter_sequence_r2 = config['adapter_sequence_r2_rc']
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
            "--adapter_sequence={params.adapter_sequence_r1} "
            "--adapter_sequence_r2={params.adapter_sequence_r2}"
        )

rule star:
    """
    Uses STAR (https://github.com/alexdobin/STAR) to align the trimmed sequences.
    STAR is heavily I/O limited during its complex alignment computations, so I have snakemake copy all the necessary files to the /scratch/ filesystem and run the alignment itself there, where they won't compute with everyone else's I/O. According to my basic benchmarking, this massively increases CPU efficiency (order ~30X), where efficiency = CPU active time / total clock time. When the alignment finishes, it copies the results back to their destination in /project/.../slamseq (the working directory) and deletes everything from /scratch/.
    
    Parameters used:
    --readFilesCommand zcat # unzips .fastq.gz files
    --runThreadN {threads} # number of threads. Flexible, 8 is fine
    --outFilterMismatchNmax 999 # maximal number of mismatches per transcript. Since this is SLAM-seq, some are expected, and we can't specify specific types (e.g. T>C), setting this to 999 effectively disables it.
    --limitBAMsortRAM 16000000000 # max RAM for BAM read sorting after alignment
    --alignEndsType EndToEnd # Don't soft-clip reads to better capture mismatches. Recommended by GRAND-SLAM docs.
    --outFilterMismatchNoverReadLmax {params.max_perread_sub_fraction} # see comment below
    --outSAMtype BAM SortedByCoordinate # Sort reads by genomic position when writing BAM output
    --outSAMattributes MD NH AS 
        # MD: string identifying mismatch/del positions in read
        # NH: total number of mismatches
        # AS: alignment score. Useful for alignment quality filtering, but not directly used here.
    """
    input: 
        fastq_r1 = 'data/trimmed/{sample_id}_R1_001.fastq.gz',
        fastq_r2 = 'data/trimmed/{sample_id}_R2_001.fastq.gz',
        index = config['star_index_path']
    output:
        aligned_bam = 'data/aligned_bam/{sample_id}.bam',
        stats = 'logs/star/qc/{sample_id}.Log.final.out'
    params:
        # define a folder in /scratch/ with the name of the sample and the slurm jobid, ensuring that it is unique.
        scratch = lambda wildcards: f"/scratch/midway3/$USER/{wildcards.sample_id}_$SLURM_JOB_ID",
        # maximum number of mismatches per transcript. A transcript with 25% T with maximal SLAM-seq substitution would be become 25% C; anything above this should be an error. In practice, substitution rates are never that high, and it doesn't make a difference
        max_perread_sub_fraction = 0.25
    log:
        "logs/star/{sample_id}.log"
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
            mkdir -p {params.scratch}

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying inputs to local scratch..." >> {log} 2>&1
            cp {input.fastq_r1} {params.scratch}/ >> {log} 2>&1
            cp {input.fastq_r2} {params.scratch}/ >> {log} 2>&1
            cp -R {input.index} {params.scratch}/ >> {log} 2>&1

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
                --alignEndsType EndToEnd \
                --outFilterMismatchNoverReadLmax {params.max_perread_sub_fraction} \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMattributes MD NH AS \
                --outFileNamePrefix {params.scratch}/ \
                >> {log} 2>&1

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying output back..." >> {log} 2>&1
            cp {params.scratch}/Aligned.sortedByCoord.out.bam {output.aligned_bam} >> {log} 2>&1
            cp {params.scratch}/Log.final.out {output.stats} >> {log} 2>&1
            cat {params.scratch}/Log.out >> {log}

            rm -rf {params.scratch}

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Alignment complete." >> {log} 2>&1
        """

rule index_bam:
    """
    Builds a .bai index for each BAM file. The Bam2CIT and GRAND-SLAM docs never mention requiring an index, but will fail if there is no file called $filename.bam.bai in the same folder as the provided input.
    """
    input: 
        "data/aligned_bam/{sample_id}.bam"
    output: 
        "data/aligned_bam/{sample_id}.bam.bai"
    shell: 
        """
        samtools index --bai {input} -o {output}
        """

def get_donor_timepoints(donor):
    """
    Helper function: retrieve all timepoints that exist for a given donor in the config file.
    """
    timepoints = [info["timepoint"] for sample_id, info in config["sample_ids"].items() 
                  if info["donor"] == donor]
    if donor == 'donor1_rep2':
        timepoints.append('no4sU') # donor1_rep2 is missing no4sU; this + get_donor_sample_id will redirect it to borrow S9, whihc is donor1_rep1's no4sU control
    return timepoints

def get_sample_from_donor_timepoint(donor, timepoint):
    """
    Given a donor (as requested in rule all:), and a specific timepoint, retrieves the corresponding sample id from config.yaml.
    - convert to str if not already str(15) -> '15'; str('no4sU') -> 'no4sU'
    - search through yaml dict for donor + timepoint_str
    - return if match
    - otherwise print error
    """
    for sample_id, info in config["sample_ids"].items():
        # donor1_rep1 and donor1_rep2 share the same control 4sU, so redirect the search.
        if donor == 'donor1_rep2' and timepoint == 'no4sU':
            donor = 'donor1_rep1'
            # print(f"donor id changed for {donor}.") 
        if info["donor"] == donor and str(info["timepoint"]) == str(timepoint):
            return sample_id
    print(
        f"No sample found for donor='{donor}' at timepoint='{timepoint}'. "
    )
    return None

rule rename_with_donor_timepoint:
    """
    When GRAND-SLAM creates its QC plots by sample, it uses their filenames, and complicated names are unwieldy and uninterpretable, so create a symlink to each sample with just its timepoint as a name.
    """
    input:
        bam=lambda w: f"data/aligned_bam/{get_sample_from_donor_timepoint(w.donor, w.timepoint_str)}.bam",
        bai=lambda w: f"data/aligned_bam/{get_sample_from_donor_timepoint(w.donor, w.timepoint_str)}.bam.bai"
    output:
        bam="data/donor_timepoint_symlinks/{donor}/{timepoint_str}.bam",
        bai="data/donor_timepoint_symlinks/{donor}/{timepoint_str}.bam.bai",
    shell:
        """
        mkdir -p $(dirname {output.bam})
        ln -sf $(realpath {input.bam}) {output.bam}
        ln -sf $(realpath {input.bai}) {output.bai}
        """

rule gedi_index_genome:
    """
    Build an index for GEDI's GRAND-SLAM (https://github.com/erhard-lab/gedi/wiki/GRAND-SLAM) to use. Produces a ton of difficult to interpret files, but the most important is the OML that tells GRAND-SLAM where everything else is. It being in a directory called ./config/genomic seems to be hard-coded into GEDI tools, so there it must go, unless there's some way around that.
    
    By default, it fetches the FASTA and GTF files from ensembl and builds the reference based on those, but if incompatiblity issues arise, you can manually provide your own versions of each, see docs at: https://github.com/erhard-lab/gedi/wiki/Preparing-genomes
    """
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
            mkdir -p {config['gedi_index_dir']}
            gedi -e IndexGenome -organism homo_sapiens -version 115 -f {config['gedi_index_dir']} -o {config['gedi_index_dir']}/homo_sapiens.115.oml -nomapping
        """

def parse_bam_to_cit_progress(progress_bool):
    """
    Map config['bam_to_cit_progress'] bool to '-p' argument or '' (empty)
    """
    if progress_bool == True:
        return "-p"
    elif progress_bool == False:
        return ''
    raise ValueError('Invalid bam_to_cit progress setting. Should be set to True or False.')
rule bam_to_cit:
    """
    GRAND-SLAM runs much faster on its custom Centered Interval Trees (CIT) format (https://github.com/erhard-lab/gedi/wiki/Mapped-reads). When all the samples in a timecourse are rolled together into a single CIT file, SNPs can be called jointly (they will be shared across samples, whereas 4sU-introduced T>C reference mismatches will not), and the whole process runs more efficiently. 
    This rule retrieves all the timepoints for a given donor using the helper function above, and generates a CIT file. Takes about 18 hours for 6 timepoints whose bams are ~20GB each.
    Makes use of the same copy to scratch trick as read alignment with STAR.
    Outputs one very large CIT file (though smaller than the sum of the BAM file sizes), and a metadata file that GRAND-SLAM will read to determine what the original files were called.
    """
    input:
        bams = lambda w: expand(
            "data/donor_timepoint_symlinks/{donor}/{timepoint_str}.bam",
            donor = w.donor,
            timepoint_str = get_donor_timepoints(w.donor),
        ),
        bais = lambda w: expand(
            "data/donor_timepoint_symlinks/{donor}/{timepoint_str}.bam.bai",
            donor = w.donor,
            timepoint_str = get_donor_timepoints(w.donor),
        ),
    output:
        cit_sample_set = "data/cit_sample_sets/{donor}.cit",
        cit_metadata = "data/cit_sample_sets/{donor}.cit.metadata.json",
    params:
        progress = parse_bam_to_cit_progress(config['bam_to_cit_progress']),
        scratch = lambda wildcards: f"/scratch/midway3/$USER/{wildcards.donor}_$SLURM_JOB_ID",
        bam_basenames = lambda w, input: " ".join([os.path.basename(b) for b in input.bams]),
        # GEDI is a Java program, and by default Java only allocates 16GB of RAM for any running programs. For large samples, this will not be enough (job crashed several times), so we set an environment variable for this allocation (java_xmx) to 87.5% of the total RAM allocation for the job. The remaining 12.5% (or 4GB, whichever is larger) is reserved for the system.
        java_xmx=lambda w, resources: int(resources.mem_mb * 0.875 / 1024),  # 87.5% in GB
        java_xms=lambda w, resources: max(4, int(resources.mem_mb * 0.125 / 1024)), # 12.5% in GB, min 4
    container:
        # Path to the apptainer container for GEDI
        config["container_path"]
    resources:
        mem_mb = 32000,
        runtime = 1440 # 20 hours in minutes
    benchmark:
        "benchmarks/{donor}.bam_to_cit.benchmark.txt"
    log:
        "logs/bam_to_cit/{donor}.log"
    shell:
        """
        log_path=$(realpath {log})
        echo "Logging to $log_path."
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying to scratch..." > $log_path 2>&1
        mkdir -p {params.scratch} >> $log_path 2>&1
        cp {input.bams} {input.bais} {params.scratch}/ >> $log_path 2>&1
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copy complete." >> $log_path 2>&1

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Generating output paths..." >> $log_path 2>&1
        cit_output_path=$(realpath {output.cit_sample_set}) >> $log_path 2>&1
        metadata_output_path=$(realpath {output.cit_metadata}) >> $log_path 2>&1
        mkdir -p $(dirname $cit_output_path) >> $log_path 2>&1
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Paths created." >> $log_path 2>&1

        cd {params.scratch} >> $log_path 2>&1
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Setting Java max memory..." >> $log_path 2>&1
        export _JAVA_OPTIONS="-Xmx{params.java_xmx}g -Xms{params.java_xms}g"

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Beginning CIT conversion..." >> $log_path 2>&1
        gedi -e Bam2CIT {params.progress} output.cit \
            {params.bam_basenames} \
             >> $log_path 2>&1
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Conversion complete." >> $log_path 2>&1

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying output back..." >> $log_path 2>&1
        cp output.cit $cit_output_path >> $log_path 2>&1
        cp output.cit.metadata.json $metadata_output_path >> $log_path 2>&1
        rm -rf {params.scratch} >> $log_path 2>&1
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copy complete, scratch cleared." >> $log_path 2>&1
        """

rule grand_slam:
    """
    Calls and quantifies nascent transcripts in each CIT file, using the GRAND-SLAM binomial mixture model (essentially, is it more likely that the T>C mismatches on this read arose from 4sU labeling or by sequencing error, after SNP correction?) Details in the GRAND-SLAM paper here: https://doi.org/10.1093/bioinformatics/bty256
    More details about individual steps are in README.md.
    Produces total read counts and nascent read fractions (grandslam.tsv.gz), in addition a dizzying array of QC data and plots, which I am only just beginning to understand. "Explanations" here: https://github.com/erhard-lab/gedi/wiki/GRAND-SLAM
    """
    input: 
        cit_sample_set = "data/cit_sample_sets/{donor}.cit",
        cit_metadata = "data/cit_sample_sets/{donor}.cit.metadata.json",
        index_oml = rules.gedi_index_genome.output.oml
    output: 
        # Total read counts, new-to-total read ratios (NTRs), and NTR distributions for all samples
        read_table = "data/slam_quant/{donor}/grandslam.tsv.gz",
        # Substitution rates per sample
        sub_rates = "data/slam_quant/{donor}/grandslam.rates.tsv",
        # Summary of mismatch types in each sample
        mismatch_summary_tsv = "data/slam_quant/{donor}/grandslam.mismatches.tsv",
        # TSV of every mismatch observed across samples
        every_mismatch_tsv = "data/slam_quant/{donor}/grandslam.mismatchdetails.tsv",
        # Plot of mismatch rates in each sample
        mismatch_plot = "data/slam_quant/{donor}/grandslam.mismatches.pdf",
        # Mismatch rates across all samples, with positions in the read. Often, there are artifactual mismatch spikes at the ends of reads that can be excluded with -trim5p and -trim3p
        mismatch_positions = "data/slam_quant/{donor}/grandslam.mismatchpos.pdf",
        # Same plot, but with narrower y-axis so spikes don't flatten everything else out
        mismatch_positions_zoomed = "data/slam_quant/{donor}/grandslam.mismatchposzoomed.pdf",
        # Same as mismatch_plot, but only for mismatches with read1 + read2 coverage 
        mismatches_double_plot = "data/slam_quant/{donor}/grandslam.double.pdf",
        # Same as mismatch_summary_tsv, but only for mismatches with read1 + read2 coverage 
        mismatches_double_tsv = "data/slam_quant/{donor}/grandslam.doublehit.tsv",
        # Number of T>C mismatches, and nascent/old RNA fractions in each sample
        ntr_stats = "data/slam_quant/{donor}/grandslam.ntrstat.tsv",
        # I don't know what these are. I'm so sorry.
        binom = "data/slam_quant/{donor}/grandslam.binom.tsv",
        binom_est = "data/slam_quant/{donor}/grandslam.binomEstimated.tsv",

        binom_overlap = "data/slam_quant/{donor}/grandslam.binomOverlap.tsv",
        binom_est_overlap = "data/slam_quant/{donor}/grandslam.binomOverlapEstimated.tsv",

        # ratios of exonic to intronic reads plot
        exon_intron_ratios = "data/slam_quant/{donor}/grandslam.exonintron.pdf",

        # per-sample Bayesian estimates of the NTR (as MAP), plus 0.95CIs and beta approximations to the posterior distribution
        beta_dist_approximations = "data/slam_quant/{donor}/grandslam.ext.tsv",
        # run parameters
        run_params = "data/slam_quant/{donor}/grandslam.param",
        # execution times for each step of GRAND-SLAM
        runtime = "data/slam_quant/{donor}/grandslam.runtime",
        # inferred strandness of sequencing experiment
        strandness = "data/slam_quant/{donor}/grandslam.strandness"
    params:
        # GEDI is a Java program, and by default Java only allocates 16GB of RAM for any running programs. For large samples, this will not be enough (job crashed several times), so we set an environment variable for this allocation (java_xmx) to 87.5% of the total RAM allocation for the job. The remaining 12.5% (or 4GB, whichever is larger) is reserved for the system.
        java_xmx=lambda w, resources: int(resources.mem_mb * 0.9 / 1024),  # 90% in GB
        java_xms=lambda w, resources: max(4, int(resources.mem_mb * 0.1 / 1024)),  # 10% in GB, min 4
        trim5p = 12,
        trim3p = 12,
    log:
        "logs/grandslam/{donor}.log"
    container:
        # Path to apptainer container for GEDI
        config['container_path']
    resources:
       runtime = 480, # 8 hours in minutes
       mem_mb = 42000,
    threads:
        14, # bump up to 24
    benchmark:
        "benchmarks/{donor}.grandslam.benchmark.txt"
    shell: 
        """
            echo "Setting Java memory..." > {log}
            export _JAVA_OPTIONS="-Xmx{params.java_xmx}g -Xms{params.java_xms}g" >> {log} 2>&1
            gedi -e Slam \
            -genomic {input.index_oml} \
            -reads {input.cit_sample_set} \
            -prefix data/slam_quant/{wildcards.donor}/grandslam \
            -trim5p {params.trim5p} \
            -trim3p {params.trim3p} \
            -no4sUpattern no4sU \
            -nthreads {threads} \
            -introns \
            -progress \
            >> {log} 2>&1
        """

rule multiqc:
    """
    Multiqc (https://docs.seqera.io/multiqc) produces a nice, interactive, all-in-one html dashboard of the QC data for numerous bioinformatics tools, using their standard QC output files. Here, we generate a report from fastp (adapter content % and post-filter per-base quality) and STAR (alignment rate, alignment quality, etc.) I might build a custom module for GRAND-SLAM.
    """
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
            'data/fastp_reports logs/star/qc '
            '--force ' # overwrite existing report; otherwise it will attach a suffix that snakemake won't detect
            '--outdir outputs'
        )

rule calc_nascent_total_reads:
    """
    Calculate nascent and total read counts for each donor, across timepoints. Reads in total read counts and maximum a posteriori (MAP) estimates of nascent-to-total ratio (NTR), and calculates nascent count as total read count * MAP. Writes gene x timepoint dataframes for total and nascent for each donor as output.
    """
    input: 
        read_table = "data/slam_quant/{donor}/grandslam.tsv.gz",
    output: 
        total_counts = "data/processed_reads/{donor}_reads_total.csv",
        nascent_counts = "data/processed_reads/{donor}_reads_nascent.csv"
    script: "scripts/calc_nascent_total_readcounts.R"

rule map_ensg_genesymbol:
    """
    Using GRAND-SLAM's output files that contain both gene alphanumeric names and ENSG IDs, create a file mapping one to other. Since gene names are non-unique, ENSGs will be used for merging and identification operations, but common names will be displayed on plots. This output file allows the former to be converted to the latter when necessary.
    """
    input: 
        expand(
            "data/slam_quant/{donor}/grandslam.tsv.gz",
            donor = DONORS[0] # only need to run once
        )
    output: 
        symbol_ensg_mapping = "data/processed_reads/ensg_genesymbol_mapping.csv"
    script: "scripts/map_ensg_symbol.R"

rule merge_reads_across_donors:
    """
    Merge read counts across donors into two CSV files, one for nascent counts and the other for total counts.
    """
    input: 
        read_counts = lambda wildcards: expand(
            "data/processed_reads/{donor}_reads_{readtype}.csv",
            donor = DONORS,
            readtype = wildcards.readtype
        )
    output: 
        merged_counts = "outputs/readcounts/merged_counts_{readtype}.csv",
    # the R script needs access to the donor list, but DONORS isn't passed as a wildcard (since it's passed in expand()), so we provide it explicitly here.
    params:
        donors = DONORS
    script: "scripts/merge_reads_across_donors.R"

rule dge:
    """
    Runs differential gene expression (DGE) analysis for the nascent and total (wildcard {readtype}) subsets of the data, and outputs summary statistics (logFC, pvalue, FDR) for each as a CSV, to be accessed for downstream analysis.
    """
    input: 
        merged_counts = "outputs/readcounts/merged_counts_{readtype}.csv",
        symbol_ensg_mapping = "data/processed_reads/ensg_genesymbol_mapping.csv",
    output: 
        dge_summary_stats = "outputs/dge_results/summary_stats_{readtype}.csv"
    script: "scripts/dge.R"

rule volcano_plot:
    """
    Generate a volcano plot given a {readtype} (total/nascent) and a time comparison (15m vs 0m, etc.,) using DGE summary statistics. Significance thresholds (drawn as lines on the plot) for FDR and logFC are passed as params.
    """
    input: 
        dge_summary_stats = "outputs/dge_results/summary_stats_{readtype}.csv"
    output: 
        volcano_plot = "outputs/volcano_plots/volcanoplot_{readtype}-{comparison}.pdf",
    params:
        fdr_threshold = 0.05,
        logFC_threshold = 1,
        palette = config['plot_color_palette'],
        plot_width = 8,
        plot_height = 12,
    script: "scripts/plot_volcano.R"

rule timecourse_dge_plot:
    """
    Generate a plot of gene fold change at each time point for genes (up to 15 by default, user-specified) that are significantly differentially expressed at at least n (2 by default, user-specified) time points.
    """
    input: 
        dge_summary_stats = "outputs/dge_results/summary_stats_{readtype}.csv"
    output: 
        timecourse_dge_plot = "outputs/timecourse_dge_plot_{readtype}.pdf",
    params:
        n_timepoints_dge_threshold = 2,
        fdr_threshold = 0.05,
        logFC_threshold = 1,
        max_genes_to_plot = 15,
        time_unit = "min",
    script: "scripts/plot_dge_timecourse.R"

rule effect_size_correlations_plot:
    """
    Plot the correlation of nascent and total fold change magnitudes, given a time point. Add a regression line and a y = x reference line.
    """
    input: 
        dge_summary_stats = expand(
            "outputs/dge_results/summary_stats_{readtype}.csv",
            readtype = ['nascent', 'total']
        )
    output: 
        corr_plot = "outputs/effect_size_correlations/corr_{comparison}.pdf"
    params:
        palette = config['plot_color_palette'],
    script: "scripts/plot_effect_size_correlations.R"

rule venn_diagram_plot:
    """
    Generate Venn diagram for genes that are DE in the nascent and the total subsets, given a tiempoint. As before, thresholds are specifiable.
    """
    input: 
        dge_summary_stats = expand(
            "outputs/dge_results/summary_stats_{readtype}.csv",
            readtype = ['nascent', 'total']
        )
    output: 
        venn_diagram = "outputs/venn_diagrams/venn_{comparison}.png"
    params:
        fdr_threshold = 0.05,
        logFC_threshold = 1,
        palette = config['plot_color_palette'],
    script: "scripts/plot_venn_diagram.R"