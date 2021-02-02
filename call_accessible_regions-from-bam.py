#!/bin/env python

import argparse
import subprocess
import os
from collections import defaultdict


genome_size = {'mm10':1.87e9,
               'hg38':2.7e9}
genome_blacklist = {"mm10":"/archive/gl/shared/user/jhammelm/old_cluster/jhammelm/genomes/atac_mm10/mm10/mm10.blacklist.bed.gz",
                    "hg38":"/archive/gl/shared/user/jhammelm/old_cluster/jhammelm/genomes/atac_human/hg38/hg38.blacklist.bed.gz"}

def ensure_dir(file_path):
    if not os.path.exists(file_path):
        os.makedirs(file_path)

def process_bam(bam_info,args):
    cmds = []
    bam = bam_info['bamfile']+'.bam'
    bam_qname  = bam_info['experiment']+'/bams/'+bam_info['samplename']+'-qnamesort.bam'
    bam_filt = bam_info['experiment']+'/bams/'+bam_info['samplename']+'-filt1.bam'
    bam_fixmate = bam_info['experiment']+'/bams/'+bam_info['samplename']+'-fixmate.bam'
    final_bam=bam_info['experiment']+'/bams/'+bam_info['samplename']+'-filtered.bam'
    bam_plusstrand = bam_info['experiment']+'/bams/'+bam_info['samplename']+'-plus.bam'
    bam_minusstrand = bam_info['experiment']+'/bams/'+bam_info['samplename']+'-minus.bam'
    bamstats = bam_info['experiment']+'/stats/'+bam_info['samplename']+'.stats'
    #Filt Bam (mito, unmapped, mapq)
    if bam_info['readtype'] == 'se':
        cmds.append(["samtools sort -n","-@",str(args.nthreads),bam,"-o ",bam_qname])
        cmds.append(["samtools view  -F 1796 -u",bam_qname,"| samtools sort /dev/stdin -o",final_bam])
        cmds.append(["rm -f ",bam_qname])
        cmds.append(["samtools index",final_bam])
    elif bam_info['assay'] == 'dnase':
        cmds.append(["samtools sort -n",bam,"-o",bam_qname])
        cmds.append(["samtools fixmate -r",bam_qname,bam_fixmate])
        cmds.append(["rm -f",bam_qname])
        cmds.append(["samtools view -F 1804 -f 2 -u",bam_fixmate,
                     "| samtools sort /dev/stdin -o",final_bam])
        cmds.append(["rm -f",bam_fixmate])
        cmds.append(["samtools index",final_bam])
        cmds.append(["samtools flagstat",final_bam,">",bamstats])
    else:
        #technically, should shift bam reads but I'm having a mental breakdown because I had the commands flipped and
        #I don't want to re-run everything
        #cmds.append(["samtools view -f 13 -u",bam," | awk '{$4+=4;print $0}' | samtools sort /dev/stdin -n -o",bam_plusstrand])
        #cmds.append(["samtools view -f 33 -u",bam," | awk '{$4-=5;print $0}' | samtools sort /dev/stdin -n -o",bam_minusstrand])
        #cmds.append(["samtools merge ",bam_qname,bam_plusstrand,bam_minusstrand])
        #cmds.append(["rm -f",bam_plusstrand])
        #cmds.append(["rm -f",bam_minusstrand])
        cmds.append(["samtools sort -n ",bam,"-o ",bam_qname])
        cmds.append(["samtools fixmate -r",bam_qname,bam_fixmate])
        cmds.append(["rm -f",bam_qname])
        cmds.append(["samtools view -F 1804 -f 2 -u",bam_fixmate,
                     "| samtools sort /dev/stdin -o",final_bam])
        cmds.append(["rm -f",bam_fixmate])
        cmds.append(["samtools index",final_bam])
        cmds.append(["samtools flagstat",final_bam,">",bamstats])
    return cmds

def tag2region(bam_info,args):
    cmds=[]
    NPEAKS=300000

    bam=bam_info['experiment']+'/bams/'+bam_info['samplename']+'-filtered.bam'

    prefix = bam_info['experiment']+'/regions/'+bam_info['samplename']+'-threshold'+str(args.pval)
    #Macs2 peaks
    if bam_info['assay'] == 'atac':
        smooth_window = 73
    else:
        smooth_window = 150
    shiftsize = round(-smooth_window/2)
    if 'mm10' == bam_info['genome_build']:
        genome='mm'
    else:
        genome='hs'
    if bam_info['readtype'] == 'se':
        cmds.append(["macs2 callpeak -t",bam,"-f BAM","-n",prefix,"-g",genome,
                 "-p",str(args.pval),"--shift",str(shiftsize),'--extsize',str(smooth_window),'--nomodel',
                 '--keep-dup all','--call-summits'])
    elif bam_info['readtype'] == 'pe':
        cmds.append(["macs2 callpeak -t",bam,"-f BAMPE","-n",prefix,"-g",genome,
                 "-p",str(args.pval),"--shift",str(shiftsize),'--extsize',str(smooth_window),'--nomodel',
                 '--keep-dup all','--call-summits'])
    
    cmds.append(["sort -k 8gr,8gr",prefix+"_peaks.narrowPeak | awk 'BEGIN{OFS=\"\\t\"}{$4=\"Peak_\"NR ; print $0}' | gzip -nc >",prefix+".narrowPeak.gz"])
    cmds.append(["rm -f",prefix+"_peaks.narrowPeak"])
    cmds.append(["rm -f",prefix+"_peaks.xls"])
    cmds.append(["rm -f",prefix+"_summits.bed"])
    cmds.append(["bedtools intersect -v -a ",prefix+".narrowPeak.gz","-b",
                genome_blacklist[bam_info['genome_build']],
                "| awk 'BEGIN{OFS=\"\\t\"} {if ($5>1000) $5=1000; print $0}'",
                "| grep -P 'chr[\dXY]+[ \\t]'  | gzip -nc >",prefix+'.filt.narrowPeak.gz'])
    return cmds
                 
def process_csv(args):
    hasheader=True
    expts = defaultdict(list)
    for line in open(args.experiment_template):
        if hasheader:
            hasheader=False
            header=line.strip().split(',')
            assert('pe/se' in header)
            assert('assay' in header)
            assert('samplename' in header)
            assert('bamfile' in header)
            assert('genome_build' in header)
            continue
        
        data = line.strip().split(',')
        exptname = data[header.index('experiment')]        

        assay = data[header.index('assay')]
        assert(assay in ['atac','dnase'])
        
        readtype = data[header.index('pe/se')]
        assert(readtype in ['pe','se'])
               
        expts[exptname].append({'experiment':exptname,
                                'samplename':data[header.index('samplename')],
                                'bamfile':data[header.index('bamfile')],
                                'assay':data[header.index('assay')],
                                'genome_build':data[header.index('genome_build')],
                                'readtype':data[header.index('pe/se')]})
    return expts

def pool_bam(expt,bams_info):
    #returns list of samtools merge command
    #        pooled bam info dict
    
    cmd = ['samtools merge '+expt+'/bams/'+expt+'_pooled_reps-filtered.bam']
    for bam in bams_info:
        cmd.append(expt+'/bams/'+bam['samplename']+'-filtered.bam')
    return cmd,{'experiment':expt,
                'samplename':expt+'_pooled_reps',
                'assay':bam['assay'],
                'readtype':bam['readtype'],
                'genome_build':bam['genome_build']}

if __name__ == "__main__":
#Steps
    parser = argparse.ArgumentParser()
    parser.add_argument('experiment_template')
    parser.add_argument('-t','--nthreads',type=int,default=2)
    parser.add_argument('-p','--pval',default=0.01,type=float)
    #parser.add_argument('-mapq','--mapq',default=255,type=int)
    opts = parser.parse_args()
    
    expt_data = process_csv(opts)

    #pool bams, call on pooled
    for expt in expt_data.keys():
        cmds = []
        ensure_dir(expt)
        ensure_dir(expt+'/bams')
        ensure_dir(expt+'/stats')
        ensure_dir(expt+'/regions')
        for bam_info in expt_data[expt]:
            cmds.extend(process_bam(bam_info,opts))

        pool_cmd,pool_data = pool_bam(expt,expt_data[expt])
        cmds.append(pool_cmd)
        expt_data[expt].append(pool_data)
        for bam_info in expt_data[expt]:
            cmds.extend(tag2region(bam_info,opts))

            
        with open(expt+'/call_accessible.sh','w') as f:
            for cmd in cmds:
                f.write(' '.join(cmd)+'\n')
        subprocess.run(["qsub -m e -M jhammelm@mit.edu -pe slots.pe "+str(opts.nthreads)+" -v PATH=$PATH -wd $PWD -N call_accessible-"+expt+" ./"+expt+"/call_accessible.sh"],shell=True)
        
#run IDR on replicates? Not really necessary since we
#never use these for anything
#IDR
