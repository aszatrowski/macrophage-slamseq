module load apptainer
apptainer build /project/lbarreiro/USERS/austin/containers/gedi_R.sif gedi.def
apptainer shell --bind /project/lbarreiro/USERS/austin,/scratch/midway3/aszatrowski /project/lbarreiro/USERS/austin/containers/gedi_R.sif