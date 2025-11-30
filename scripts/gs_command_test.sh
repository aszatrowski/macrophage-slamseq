gedi -e IndexGenome -organism homo_sapiens -version 115 -n ensembl_hg38 -f config/genomic -o config/genomic/hg38.oml -D -p -nomapping
gedi -e Slam -genomic config/genomic/homo_sapiens.115.oml data/cit/LB-HT-28s-HT-17_S17.cit -progress -prefix data/slam_quant/LB-HT-28s-HT-17_S17/grandslam -nthreads 4 -introns



gedi -e Slam -genomic config/genomic/homo_sapiens.115.oml -reads data/aligned_bam/LB-HT-28s-HT-18_S18.bam -progress -prefix data/slam_quant/LB-HT-28s-HT-18_S18/grandslam -nthreads 4 -introns # WORKS up to modeling
gedi -e Slam -genomic config/genomic/homo_sapiens.115.oml -reads data/cit/LB-HT-28s-HT-17_S17.cit -progress -prefix data/slam_quant/LB-HT-28s-HT-17_S17/grandslam -nthreads 4 -introns