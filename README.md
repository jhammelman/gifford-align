-----------------------------------------------------------------------------------------------------------------------------
# Making TrackHubs
-----------------------------------------------------------------------------------------------------------------------------

python make-trackhub.py examples/trackhub-template.txt trackFolder

trackFolder can be an existing or new folder in 00_UCSC_tracks

If trackFolder doesn’t exist, it will make one with all the appropriate stuff for a UCSC trackhub

If trackFolder does exist, it will append the lines for the new data to the trackDB.txt

This works for bed,bedpe, and bam. 
