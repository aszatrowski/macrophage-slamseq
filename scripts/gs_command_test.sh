gedi -e IndexGenome -organism homo_sapiens -version 115 -n ensembl_hg38 -f config/genomic -o config/genomic/hg38.oml -D -p -nomapping

# WORKS FOR CIT BATCH CREATION. something about -id causes it to fail with "can only keep integer read names"
gedi -e Bam2CIT -p data/cit/s27_no4sU_test.cit data/aligned_bam/LB-HT-28s-JL-09_S27.bam data/no4sU_tagged/LB-HT-28s-HT-09_S9.no4sU.bam
gedi -e Bam2CIT -p cit/s27_no4sU_test.cit LB-HT-28s-JL-09_S27.bam LB-HT-28s-HT-09_S9.no4sU.bam

gedi -e Slam -genomic config/genomic/homo_sapiens.115.oml data/cit/LB-HT-28s-HT-17_S17.cit -progress -prefix data/slam_quant/LB-HT-28s-HT-17_S17/grandslam -nthreads 4 -introns



gedi -e Slam -genomic config/genomic/homo_sapiens.115.oml -reads data/aligned_bam/LB-HT-28s-HT-18_S18.bam -progress -prefix data/slam_quant/LB-HT-28s-HT-18_S18/grandslam -nthreads 4 -introns # WORKS up to modeling
gedi -e Slam -genomic config/genomic/homo_sapiens.115.oml -reads data/cit/LB-HT-28s-HT-17_S17.cit -progress -prefix data/slam_quant/LB-HT-28s-HT-17_S17/grandslam -nthreads 4 -introns

gedi -e Slam -genomic config/genomic/homo_sapiens.115.oml -reads data/cit/test.bamlist -progress -prefix data/slam_quant/donor1_rep2_test/grandslam -nthreads 4 -introns -no4sUpattern no4sU # works !!! might be slow though with only bams.
gedi -e Slam -genomic /project/lbarreiro/USERS/austin/slamseq/config/genomic/homo_sapiens.115.oml -reads cit/s27_no4sU_test.cit -prefix slam_quant/s27_no4sU_test/grandslam -nthreads 12 -introns -no4sUpattern no4sU -progress # WORKS, running now