apptainer build containers/gedi_1.0.6a.sif slamseq/gedi.def

apptainer exec containers/gedi_1.0.6a.sif gedi --help

snakemake --use-singularity --cores 4