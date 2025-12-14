# before cd to scratch, get destination folder
cit_output_path=$(realpath data/cit_sample_sets/donor1_rep2.cit)
metadata_output_path=$(realpath data/cit_sample_sets/donor1_rep2.cit)
mkdir -p $(dirname output_path)
# gedi -e Slam ...
cd /scratch/midway3/aszatrowski/donor1_rep2_43174494
cp output.cit $output_path
cp output.cit.metadata.json $metadata_output_path