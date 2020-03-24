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
    bam_qname = bam_info['experiment']+'/bams/'+bam_info['samplename']+'-qnamesort.bam'
    bam_filt = bam_info['experiment']+'/bams/'+bam_info['samplename']+'-filt1.bam'
    bam_fixmate = bam_info['experiment']+'/bams/'+bam_info['samplename']+'-fixmate.bam'
    final_bam=bam_info['experiment']+'/bams/'+bam_info['samplename']+'-filtered.bam'
    #Filt Bam (mito, unmapped, mapq)
    if bam_info['readtype'] == 'se':
        cmds.append(["samtools sort -n","-@",str(args.nthreads),bam,"-o ",bam_qname])
        cmds.append(["samtools view  -F 1804 -u ","-@",str(args.nthreads),"-q",str(args.mapq),bam_qname,"| samtools sort /dev/stdin -o",final_bam])
        cmds.append(["rm -f",bam_qname])
        cmds.append(["samtools index",final_bam])
    else:
        cmds.append(["samtools sort -n ",bam,"-o ",bam_qname])
        cmds.append(["samtools fixmate -r",bam_qname,bam_fixmate])
        cmds.append(["rm -f",bam_qname])
        cmds.append(["samtools view -f 2 -u",bam_fixmate,
                     "| samtools sort /dev/stdin -o",final_bam])
        cmds.append(["samtools view -F 1804 -f 2 -u",bam_fixmate,
                     "| samtools sort /dev/stdin -o",final_bam])
        cmds.append(["rm -f",bam_fixmate])
        cmds.append(["samtools index",final_bam])
    return cmds
    
def bam2tag(bam_info):
    cmds = []
    bam=bam_info['experiment']+'/bams/'+bam_info['samplename']+'-filtered.bam'
    tag = bam_info['experiment']+'/bams/'+bam_info['samplename']+'.tagAlign.gz'
    shifted_tag = bam_info['experiment']+'/bams/'+bam_info['samplename']+'.tn5.tagAlign.gz'
    outtmp = bam_info['experiment']+'/bams/'+bam_info['samplename']+'.srt.tmp.bam'

    cmds.append(["samtools sort -n ",bam," -o",outtmp])
    #Filt Bam -> Bedpe
    if bam_info['readtype'] == 'pe':
        #pair end
        cmds.append(["bedtools bamtobed -bedpe -mate1 -i", outtmp,"| awk 'BEGIN{OFS=\"\\t\"}{printf \"%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n\",$1,$2,$3,$9,$4,$5,$6,$10}' | gzip -nc >",tag])
    else: 
        #single end
        cmds.append(["bedtools bamtobed -i", bam,"| awk 'BEGIN{OFS=\"\\t\"}{$4=\"N\";$5=\"1000\";print $0}' | gzip -nc >",tag])
    #Tn5 shift if atac-seq
    if bam_info['assay'] == 'atac':
        cmds.append(["zcat ",tag," | awk 'BEGIN {OFS = FS = \"\t\"}{ if ($6 == \"+\") {$2 = $2 + 4} else if ($6 == \"-\") {$3 = $3 - 5} print $0}' | gzip -nc > ",shifted_tag])
    cmds.append(['rm ',outtmp])
    
    return cmds


def tag2region(bam_info,args):
    NPEAKS=300000
    if bam_info['assay'] == 'atac':
        tag = bam_info['experiment']+'/bams/'+bam_info['samplename']+'.tn5.tagAlign.gz'
    else:
        tag = bam_info['experiment']+'/bams/'+bam_info['samplename']+'.tagAlign.gz'
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
    cmds.append(["macs2 callpeak -t",tag,"-f BED","-n",prefix,"-g",genome,
                 "-p",str(args.pval),"--shift",str(shiftsize),'--extsize',str(smooth_window),'--nomodel',
                 '--SPMR','--keep-dup all','--call-summits'])
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
        cmd.append(expt+'/bams/'+bam['experiment']+'-filtered.bam')
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
    parser.add_argument('-mapq','--mapq',default=255,type=int)
    opts = parser.parse_args()
    
    expt_data = process_csv(opts)
    
    #pool bams, call on pooled
    for expt in expt_data.keys():
        cmds = []
        ensure_dir(expt)
        ensure_dir(expt+'/bams')
        ensure_dir(expt+'/regions')
        for bam_info in expt_data[expt]:
            cmds.extend(process_bam(bam_info,opts))
        pool_cmd,pool_data = pool_bam(expt,expt_data[expt])
        cmds.append(pool_cmd)
        expt_data[expt].append(pool_data)
        for bam_info in expt_data[expt]:
            cmds.extend(bam2tag(bam_info))
            cmds.extend(tag2region(bam_info,opts))
            
        with open(expt+'/call_accessible.sh','w') as f:
            for cmd in cmds:
                f.write(' '.join(cmd)+'\n')
                
        subprocess.run(["qsub -m e -M jhammelm@mit.edu -pe slots.pe "+str(opts.nthreads)+" -v PATH=$PATH -wd $PWD -N call_accessible-"+expt+" ./"+expt+"/call_accessible.sh"],shell=True)
        
#run IDR on replicates? Not really necessary since we
#never use these for anything
#IDR
