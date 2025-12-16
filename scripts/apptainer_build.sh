module load apptainer
apptainer build /project/lbarreiro/USERS/austin/containers/gedi_R.sif gedi.def
apptainer shell --bind /project/lbarreiro/USERS/austin,/scratch/midway3/aszatrowski /project/lbarreiro/USERS/austin/containers/gedi_R.sif

export _JAVA_OPTIONS="-Xmx17g -Xms4g"                                                                       
gedi -e Slam -genomic config/genomic/homo_sapiens.115.oml -reads data/cit_sample_sets/donor1_rep2.cit -prefix data/slam_quant/donor1_rep2_old/grand_slam -introns -no4sUpattern no4sU -nthreads 16 -progress -full -plot

gedi -e Slam -genomic /project/lbarreiro/USERS/austin/slamseq/config/genomic/homo_sapiens.115.oml -reads s27_no4sU_test.cit -prefix s27_test_output/grandslam -introns -progress