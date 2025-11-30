module load apptainer
apptainer build containers/gedi_1.0.6a.sif slamseq/gedi.def
apptainer shell --bind /project/lbarreiro/USERS/austin,/scratch/midway3/aszatrowski /project/lbarreiro/USERS/austin/containers/gedi_1.0.6a.sif