#!/bin/bash

for LANE in {2..5}
do
    echo $LANE

    Rscript /media/MPEdge16/MM137/sc/pre_scripts/retrieve_barcodes.R cd "/media/MPEdge16/MM137/sc/10x_cloud_dl/10872-MM-${LANE}_standard/"
    
    cd /media/MPEdge16/MM137/sc/

    nohup /media/MPEdge16/MM137/sc/singularity exec -B "/media/MPEdge16/MM137/sc/10x_cloud_dl/10872-MM-${LANE}_standard" /media/MPEdge16/MM137/sc/souporcell_latest.sif souporcell_pipeline.py -i "/media/MPEdge16/MM137/sc/10x_cloud_dl/10872-MM-${LANE}_standard/possorted_genome_bam.bam" -b "/media/MPEdge16/MM137/sc/10x_cloud_dl/10872-MM-${LANE}_standard/barcodes_R.tsv" -f /media/MPEdge16/MM137/sc/10x_cloud_dl/hg38.fa -t 40 -o "/media/MPEdge16/MM137/sc/souporcell_outs/10872-MM-${LANE}_soc" -k 3 &
done

exit 0
