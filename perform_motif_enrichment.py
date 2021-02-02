#!/bin/env python


import argparse
import subprocess
import gzip
import os

distal_file = "/cluster/yuchun/projects/shared_data/mm10_hgTables.TSS2.5kb.coordSorted.bed"
mm10_fasta = "/archive/gl/shared/genomes/mm10/mm10.fa"
ame_path  = "/archive/gl/shared/projects/wichterleMN/tools/meme/bin/ame"
meme_path = "/archive/gl/shared/projects/wichterleMN/tools/meme/bin/dreme"
kmac_path= "/archive/gl/shared/software/gem/gem.jar"
homer_path = "/data/gl/g1/jhammelm/software/HOMER/bin/"
ame_database_all = "/archive/gl/shared/user/jhammelm/shared_data/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme"
ame_database_consensus = "/archive/gl/shared/user/jhammelm/shared_data/consensus_HOCOMOCOv11_core_MOUSE-affinityprop-handannotated.meme"

parser = argparse.ArgumentParser()
parser.add_argument('bedfile')
parser.add_argument('-o','--outfile',required=True)
parser.add_argument('-p','--percent',default=0,type=int)
parser.add_argument('-n','--topN',default=10000,type=int)
parser.add_argument('-db','--database',default='all',choices=['all','consensus'])
parser.add_argument('-bg','--background',default='shuffled')
parser.add_argument('-c','--sortcol',type=int,default=8)
parser.add_argument('-ame','--ame',default=False,action='store_true')
parser.add_argument('-meme','--meme',default=False,action='store_true')
parser.add_argument('-kmac','--kmac',default=False,action='store_true')
parser.add_argument('-homer','--homer',default=False,action='store_true')
opts = parser.parse_args()

if '.gz' in opts.bedfile:
    #print('is gzip')
    subprocess.call(['zcat '+opts.bedfile+' | sort -k1,1 -k2,2n | bedtools groupby -i '+opts.bedfile+' -grp 1,2,3 -c '+str(opts.sortcol)+' -o max > '+opts.outfile+'.tmp'],shell=True)
    #regions = [line.decode('utf8').strip().split() for line in gzip.open(opts.bedfile)]
else:
    subprocess.call(['cat '+opts.bedfile+' | sort -k1,1 -k2,2n | bedtools groupby -i '+opts.bedfile+' -grp 1,2,3 -c '+str(opts.sortcol)+' -o max > '+opts.outfile+'.tmp'],shell=True)

regions = [line.strip().split() for line in open(opts.outfile+'.tmp')]

if opts.percent != 0:
    N = int(len(regions)*(opts.percent/100.0))
else:
    N = opts.topN


#higher is more significant
regions_sorted = sorted(regions, key=lambda x: float(x[3]))[::-1]

top_n_regions = regions_sorted[:N]


with open(opts.outfile+'-select'+str(N)+'.bed','w') as f:
    f.write('\n'.join(['\t'.join(r) for r in top_n_regions]))

subprocess.call(['bedtools getfasta -fi '+mm10_fasta+' -bed '+opts.outfile+'-select'+str(N)+'.bed'+ ' -fo ' +opts.outfile+'-select'+str(N)+'.fa'],shell=True)

if opts.database == 'all':
    db = ame_database_all
else:
    db = ame_database_consensus
    
if opts.ame:
    if opts.background == 'shuffled':
        subprocess.call([ame_path + ' --control --shuffle-- -oc '+opts.outfile+'-ame '+opts.outfile+'-select'+str(N)+'.fa '+db],shell=True)
    elif opts.background == 'homer-matched':
        if not os.path.exists(opts.outfile+'-homer/background.fa'):
            subprocess.call([homer_path + 'findMotifsGenome.pl '+opts.outfile+'-select'+str(N)+'.bed mm10 ' + opts.outfile+'-homer -size given -dumpFasta'],shell=True)
        subprocess.call([ame_path + ' --control '+opts.outfile+'-homer/background.fa -oc '+opts.outfile+'-ame '+opts.outfile+'-select'+str(N)+'.fa '+db],shell=True)
    else:
        subprocess.call([ame_path + ' --control '+opts.background+' -oc '+opts.outfile+'-ame '+opts.outfile+'-select'+str(N)+'.fa '+db],shell=True)

if opts.meme:
    if opts.background == 'shuffled':
        subprocess.call(['python2 '+meme_path + ' -p ' +opts.outfile+'-select'+str(N)+'.fa '+ ' -oc '+opts.outfile+'-meme '],shell=True)
    elif opts.background == 'homer-matched':
        if not os.path.exists(opts.outfile+'-homer/background.fa'):
            subprocess.call([homer_path + 'findMotifsGenome.pl '+opts.outfile+'-select'+str(N)+'.bed mm10 ' + opts.outfile+'-homer -size given -dumpFasta'],shell=True)
        subprocess.call(['python2 '+meme_path + ' -p '+opts.outfile+'-select'+str(N)+'.fa -n '+opts.outfile+'-homer/background.fa -oc '+opts.outfile+'-meme '],shell=True)
    else:
        subprocess.call(['python2 '+meme_path + ' -p '+opts.outfile+'-select'+str(N)+'.fa -n '+opts.background+' -oc '+opts.outfile+'-meme '],shell=True)
        
if opts.kmac:
    cmd = 'java -Xmx8G -jar '+kmac_path+' KMAC --pos_seq '+opts.outfile+'-select'+str(N)+'.fa --k_win 100 --k_min 4 --k_max 13 --t 1 --k_seqs 10000 --k_top 10 --gap 4 --out_name '+ opts.outfile+'-kmac'
    if opts.background == 'homer-matched':
        if not os.path.exists(opts.outfile+'-homer/background.fa'):
            subprocess.call([homer_path + 'findMotifsGenome.pl '+opts.outfile+'-select'+str(N)+'.bed mm10 ' + opts.outfile+'-homer -size given -dumpFasta'],shell=True)
        cmd += ' --neg_seq '+opts.outfile+'-homer/background.fa'
    elif opts.background != 'shuffled':
        cmd += ' --neg_seq '+opts.background
        
    subprocess.call([cmd],shell=True)
    
if opts.homer:
    if opts.background != 'shuffled':
        cmd = homer_path + 'findMotifs.pl '+opts.outfile+'-select'+str(N)+'.fa fasta ' + opts.outfile+'-homer -fasta ' + opts.background
    else:
        cmd = homer_path + 'findMotifsGenome.pl '+opts.outfile+'-select'+str(N)+'.bed mm10 ' + opts.outfile+'-homer -size given'
    subprocess.call([cmd],shell=True)
