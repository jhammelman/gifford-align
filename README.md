# Scripts for basic bioinformatics submissions to sge scheduler
-----------------------------------------------------------------------------------------------------------------------------
## Making TrackHubs
-----------------------------------------------------------------------------------------------------------------------------

python make-trackhub.py examples/trackhub-template.txt trackFolder
usage: make-trackhub.py [-h] [-norm] [-bs BINSIZE] experiment_template trackname
make-trackhub.py: error: the following arguments are required: experiment_template, trackname

trackFolder can be an existing or new folder in 00_UCSC_tracks

If trackFolder doesn’t exist, it will make one with all the appropriate stuff for a UCSC trackhub

If trackFolder does exist, it will append the lines for the new data to the trackDB.txt

This works for bed, bedpe, and bam. 

-----------------------------------------------------------------------------------------------------------------------------
## Alignment
-----------------------------------------------------------------------------------------------------------------------------

python align.py examples/rna-seq-template.txt START -rsem -rmdup -trim_adaptors
usage: align.py [-h] [-t NTHREADS] [-trim_adaptors] [-rmdup] [-quality QUALITY] [-rsem]
                [-cufflinks]
                experiment_template {STAR,bwa,bowtie2}
align.py: error: the following arguments are required: experiment_template, aligner

-----------------------------------------------------------------------------------------------------------------------------
## Peak calling
-----------------------------------------------------------------------------------------------------------------------------
python call_accessible_regions.py examples/peaks-template.txt
usage: call_accessible_regions.py [-h] [-t NTHREADS] [-p PVAL] [-mapq MAPQ] experiment_template
call_accessible_regions.py: error: the following arguments are required: experiment_template
