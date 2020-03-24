#!/bin/env python

import argparse
import subprocess
import os

'''
define constants
'''
UCSC_PATH='/archive/gl/shared/projects/wichterleMN/www/00_UCSC_tracks/'
BAM2BW='/archive/gl/shared/projects/wichterleMN/tools/bamCoverage'
BED2BB='/archive/gl/shared/projects/wichterleMN/tools/bedToBigBed'
interact='/archive/gl/shared/projects/wichterleMN/tools/interact.as'
BEDPEINTERACT='/archive/gl/shared/projects/wichterleMN/pipelines/bedpe2interact.py'
SAMTOOLS='/archive/gl/shared/projects/wichterleMN/tools/samtools'

chromSizes = {'mm10':'/archive/gl/shared/user/jhammelm/old_cluster/jhammelm/genomes/atac_mm10/mm10/mm10.chrom.sizes',
              'hg38':'/archive/gl/shared/user/jhammelm/old_cluster/jhammelm/genomes/atac_human/hg38/hg38.chrom.sizes'}

def ensure_dir(file_path):
    if not os.path.exists(file_path):
        os.makedirs(file_path)

def write_bigwig(filepath,filename,genome_build,opts):
    cmd = []
    cmd.append(SAMTOOLS+' index '+filepath)
    cmd.append(' '.join([BAM2BW +' --ignoreDuplicates --centerReads',
                     '--bam',filepath,
                     '--outFileName',UCSC_PATH+opts.trackname+'/'+genome_build+'/'+filename+'.bw',
                     '--binSize',str(opts.binSize),
                     '--outFileFormat bigwig']))
    if opts.normBam:
        cmd[-1] += '--normalizeUsing CPM'
        
    return cmd
    

def write_bigbed(filepath,filename,genome_build,opts):
    cmd = []
    cmd.append('sort -k1,1 -k2,2n '+filepath+' | cut -f1,2,3 > /tmp/'+filename+'.sorted')
    cmd.append(' '.join([BED2BB,'-type=bed3+ /tmp/'+filename+'.sorted',chromSizes[genome_build],
                         UCSC_PATH+opts.trackname+'/'+genome_build+'/'+filename+'.bb']))
    return cmd
    

def write_bigInteract(filepath,filename,genome_build,opts):
    cmd = []
    cmd.append(' '.join(['python',BEDPEINTERACT,filepath,'> /tmp/'+filename+'.bed']))
    cmd.append('sort -k1,1 -k2,2n /tmp/'+filename+'.bed > /tmp/'+filename+'.sorted')
    cmd.append(' '.join([BED2BB,'-as='+interact,'-type=bed5+13 /tmp/'+filename+'.sorted',chromSizes[genome_build],
                         UCSC_PATH+opts.trackname+'/'+genome_build+'/'+filename+'.bb']))
    return cmd
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('experiment_template')
    parser.add_argument('trackname')
    parser.add_argument('-norm','--normBam',action='store_true',default=False)
    parser.add_argument('-bs','--binSize',type=int,default=5)
    opts = parser.parse_args()

    genomes = []

    ensure_dir(UCSC_PATH+opts.trackname)


    with open(UCSC_PATH+opts.trackname+'/hub.txt','w') as f:
        f.write('\n'.join(['hub '+opts.trackname,
                           'shortLabel '+opts.trackname,
                           'longLabel '+opts.trackname,
                           'genomesFile genomes.txt',
                           'email jhammelm@mit.edu'])+'\n')
    hasheader=True
    for line in open(opts.experiment_template):
        if hasheader:
            hasheader=False
            header=line.strip().split(',')
            assert('samplepath' in header)
            assert('genome_build' in header)
            assert('sampletype' in header)
            continue
        
        data = line.strip().split(',')

        
        filepath = data[header.index('samplepath')]
        if 'samplename' not in header:
            filename = filepath.split('/')[-1]
            filename = '.'.join(filename.split('.')[:-1])
        else:
            filename = data[header.index('samplename')]
        
        filetype = data[header.index('sampletype')]
        genome_build = data[header.index('genome_build')]

        assert(filetype in ['bam','bed','bedpe'])
    
        ensure_dir(UCSC_PATH+opts.trackname+'/'+genome_build)
        if genome_build not in genomes:
            genomes.append(genome_build)

        if filetype == 'bam':
            cmds = write_bigwig(filepath,filename,genome_build,opts)
            with open(UCSC_PATH+opts.trackname+'/'+genome_build+'/trackDb.txt','a') as f:
                f.write('\n'.join(['track '+filename,
                                   'windowingFunction maximum',
                                   'bigDataUrl '+filename+'.bw',
                                   'shortLabel '+filename,
                                   'longLabel '+filename,
                                   'type bigWig'])+'\n\n')
        elif filetype == 'bed':
            cmds = write_bigbed(filepath,filename,genome_build,opts)
            
            with open(UCSC_PATH+opts.trackname+'/'+genome_build+'/trackDb.txt','a') as f:
                f.write('\n'.join(['track '+filename,
                                   'bigDataUrl '+filename+'.bb',
                                   'shortLabel '+filename,
                                   'longLabel '+filename,
                                   'type bigBed'])+'\n\n')            
        elif filetype == 'bedpe':
            cmds = write_bigInteract(filepath,filename,genome_build,opts)
            with open(UCSC_PATH+opts.trackname+'/'+genome_build+'/trackDb.txt','a') as f:
                f.write('\n'.join(['track '+filename,
                                   'bigDataUrl '+filename+'.bb',
                                   'shortLabel '+filename,
                                   'longLabel '+filename,
                                   'useScore on',
                                   'visibility full',
                                   'interactDirectional false',
                                   'type bigInteract'])+'\n\n')

        
        with open(UCSC_PATH+opts.trackname+'/commands.log','a') as f:
            f.write('\n'.join(cmds)+'\n')

        with open(UCSC_PATH+opts.trackname+'/'+filename+'.sh','w') as f:
            f.write('\n'.join(cmds))

        #subprocess.run(['qsub  -v PATH=$PATH -wd $PWD -N run_'+filename+' '+UCSC_PATH+opts.trackname+'/'+filename+'.sh'],shell=True)
        
with open(UCSC_PATH+opts.trackname+'/genomes.txt','w') as f:
    for genome in genomes:
        f.write('\n'.join(['genome '+genome,
                           'trackDb '+genome+'/trackDb.txt']
                          )+'\n')
                           
                
