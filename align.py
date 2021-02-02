#!/bin/env python

import argparse
import subprocess
import os

def is_fq(filename):
    if '.fq' in filename:
        return True
    if '.fastq' in filename:
        return True
    
def strip_fq(name):
    name = name.split('.gz')[0]
    name = name.split('.fq')[0]
    name = name.split('.fastq')[0]
    return name

def get_ext(name):
    ext=''
    if '.fq' in name:
        ext += '.fq'
    else:
        ext += '.fastq'
    if '.gz' == name[-3:]:
        ext += '.gz'
    return ext

def ensure_dir(file_path):
    if not os.path.exists(file_path):
        os.makedirs(file_path)
        
parser=argparse.ArgumentParser()
parser.add_argument('experiment_template')
parser.add_argument('aligner',choices=['STAR','bwa','bowtie2'])
parser.add_argument('-t','--nthreads',default=4,type=int)
parser.add_argument('-trim_adaptors','--trim_adaptors',
                    action='store_true',default=False)
parser.add_argument('-multimap','--keepmultimappers',action='store_true',default=False)
parser.add_argument('-rmdup','--rmdup',
                    action='store_true',default=False)
parser.add_argument('-quality','--quality',type=int,
                    default=20)
parser.add_argument('-rsem','--rsem',
                    action='store_true',default=False)
parser.add_argument('-cufflinks','--cufflinks',action='store_true',default=False)
parser.add_argument('-local','--local',action='store_true',default=False)
opts=parser.parse_args()

hasheader=True
header=None

trimadaptorpath="/cluster/jhammelm/tools/TrimGalore/trim_galore"
STARpath="/cluster/software/STAR-2.5.2b/bin/Linux_x86_64/"
bwapath="/cluster/jhammelm/tools/bwa-0.7.17/"
bowtie2path="/cluster/software/bowtie2-2.2.9/"
rsempath="/cluster/software/RSEM/RSEM-1.3.0/rsem-calculate-expression"
cufflinkspath="/archive/gl/shared/software/cufflinks-2.2.1.Linux_x86_64/cufflinks"
fqcpath="/data/cgs/jhammelm/software/FastQC/fastqc"
genomepaths = {"mm10":{"STAR":"/data/cgs/jhammelm/genomes/mm10/ref/",
                       #"bowtie2":"/data/cgs/jhammelm/genomes/mm10/ref/mm10-bowtie2",
                       "bowtie2":"/archive/gl/shared/genomes/mm10/mm10",
                       "bwa":"/data/cgs/jhammelm/genomes/mm10/ref/mm10.fa"},
               "hg38":{"bwa":"/archive/gl/shared/user/jhammelm/old_cluster/jhammelm/genomes/atac_human/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
                       "bowtie2":"/archive/gl/shared/user/jhammelm/old_cluster/jhammelm/genomes/atac_human/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"},
            "hg19":{"bwa":"/cluster/genomes/GRCh37/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa",
                       "bowtie2":"/cluster/genomes/GRCh37/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"}}
rsemref = {"mm10":"/archive/gl/shared/projects/wichterleMN/shared_data/genomes/gencode-vm24-mm10/ref/mm10-STAR",
           "mm10-ERCC":"/archive/gl/shared/projects/wichterleMN/shared_data/genomes/gencode-vm24-mm10-ERCC/ref/mm10-STAR-ERCC",
           "mm10-TE":"/archive/gl/shared/projects/wichterleMN/transposable_elements/reference/mm10-STAR"}
cufflinksref = {"mm10":"/archive/gl/shared/projects/wichterleMN/shared_data/mm10_ncbi_refseq-2.gtf"}
for line in open(opts.experiment_template):
    if hasheader:
        hasheader=False
        header=line.strip().split(',')
        assert('samplename' in header)
        assert('genome_build' in header)
        assert('description' in header)
        continue
    data = line.strip().split(',')

    fastqpath = data[header.index('samplename')]
    description = data[header.index('description')]
    genome_build = data[header.index('genome_build')]
    if 'strand' in header:
        strand = data[header.index('strand')]
                        
    #assert(genome_build in genomepaths)
    
    if '/' in fastqpath:
        path = '/'.join(fastqpath.split('/')[:-1])
        name = fastqpath.split('/')[-1]
    else:
        path = '.'
        name = fastqpath
        
    fqfiles = sorted([f for f in os.listdir(path) if name in f and is_fq(f)])
    print(fqfiles)
    assert(len(fqfiles) <= 2)
    
    exptfolder=description+'_'+opts.aligner
    ensure_dir(exptfolder)
    
    ensure_dir(exptfolder+'/fastqc')
    ensure_dir(exptfolder+'/trimgalore')
    with open(exptfolder+'/run_align.sh','w') as f:
        f.write('#!/bin/bash\n')
        #f.write('conda activate align-pipeline\n')
        f.write(fqcpath+' '+' '.join(fqfiles)+' -o '+exptfolder+'/fastqc\n')
        if opts.trim_adaptors:
            if len(fqfiles) == 2:
                f.write(trimadaptorpath+' --quality '+str(opts.quality)+' --paired '+' '.join(fqfiles) + '\n')
            else:
                f.write(trimadaptorpath+' --quality '+str(opts.quality)+' '+fqfiles[0] + '\n')
            if len(fqfiles)==2:
                if '.gz' in fqfiles[0]:
                    fqfiles = ['$PWD/'+strip_fq(fqfiles[0])+'_val_1.fq.gz','$PWD/'+strip_fq(fqfiles[1])+'_val_2.fq.gz']
                else:
                      fqfiles = ['$PWD/'+strip_fq(fqfiles[0])+'_val_1.fq','$PWD/'+strip_fq(fqfiles[1])+'_val_2.fq']                  
            else:
                if '.gz' in fqfiles[0]:
                    fqfiles = ['$PWD/'+strip_fq(f)+'_trimmed.fq.gz' for f in fqfiles]
                else:     
                      fqfiles = ['$PWD/'+strip_fq(fqfiles[0])+'_trimmed.fq']                 
        #write alignment pipeline
        
        assert(opts.aligner in ['STAR','bwa','bowtie2'])
        if not opts.rsem:
            alignment = []
            if opts.aligner == 'STAR':
                alignment.append(STARpath+'STAR --runMode alignReads')
                alignment.extend(['--readFilesIn']+fqfiles)
                if '.gz' in fqfiles[0]:
                    alignment.append('--readFilesCommand zcat')
                alignment.append('--outStd SAM')
                alignment.append('--genomeDir '+genomepaths[genome_build][opts.aligner])
                alignment.append('--runThreadN '+str(opts.nthreads))
                alignment.append('--sjdbGTFfile '+str(cufflinksref[genome_build]))
            elif opts.aligner == 'bwa':
                alignment.append(bwapath + 'bwa mem')
                alignment.append('-t '+str(opts.nthreads))
                if opts.keepmultimappers:
                    alignment.append('-a')
                alignment.append(genomepaths[genome_build][opts.aligner])
                alignment.extend(fqfiles)
            elif opts.aligner == 'bowtie2':
                alignment.append(bowtie2path+'bowtie2')
                if opts.keepmultimappers:
                    alignment.append('-k')
                alignment.append('-p '+str(opts.nthreads))
                alignment.append('-x '+genomepaths[genome_build][opts.aligner])
                if len(fqfiles) == 2:
                    alignment.append(' -1 '+fqfiles[0])
                    alignment.append(' -2 '+fqfiles[1])
                else:
                    alignment.append(' -U '+fqfiles[0])
            alignment.append('| samtools view -bS - > '+exptfolder+'/'+description+'.bam')
            f.write(' '.join(alignment)+'\n')
            f.write('samtools index '+exptfolder+'/'+description+'.bam\n')
            f.write('samtools flagstat '+exptfolder+'/'+description+'.bam > '+exptfolder+'/'+description+'.bamstats\n')
            f.write('samtools sort -n  '+exptfolder+'/'+description+'.bam  > '+exptfolder+'/'+description+'-namesort.bam\n')
            f.write('samtools fixmate -m -O bam '+exptfolder+'/'+description+'-namesort.bam '+exptfolder+'/'+description+'-fixmate.bam\n')
            f.write('samtools sort  '+exptfolder+'/'+description+'-fixmate.bam  > '+exptfolder+'/'+description+'-sorted.bam\n')
            f.write('samtools markdup -r '+exptfolder+'/'+description+'-sorted.bam '+exptfolder+'/'+description+'-rmdup.bam\n')
            f.write('samtools sort '+exptfolder+'/'+description+'-rmdup.bam > '+exptfolder+'/'+description+'-final.bam\n')
                
            f.write('rm '+exptfolder+'/'+description+'-rmdup.bam\n')
            f.write('rm '+exptfolder+'/'+description+'-namesort.bam\n')
            f.write('rm '+exptfolder+'/'+description+'-sorted.bam\n')
            f.write('rm '+exptfolder+'/'+description+'-fixmate.bam\n')
        else:
            rsem = [rsempath]
            if len(fqfiles) > 1:
                rsem.append('--paired-end')
            rsem.append('--output-genome-bam')
            rsem.extend(['-p',str(opts.nthreads)])
            rsem.extend(fqfiles)
            if opts.aligner == "STAR":
                rsem.extend(['--star','--star-path',STARpath])
                if '.gz' in fqfiles[0]:
                    rsem.append('--star-gzipped-read-file')
            elif opts.aligner == "bowtie2":
                rsem.extend(['--bowtie2','--bowtie2-path',bowtie2path])
            rsem.append(rsemref[genome_build])
            rsem.append(exptfolder+'/'+description+'-rsem')
            rsem.append('--append-names')
            f.write(' '.join(rsem)+'\n')
        if opts.cufflinks:
            cl = [cufflinkspath]
            cl.extend(['-g',cufflinksref[genome_build],'--GTF',
                       exptfolder+'/'+description+'-final.bam'])
            
            counts =  ['htseq-count']
            if strand == 'forward':
                cl.extend(['--library-type','fr-firststrand'])
                counts.extend([exptfolder+'/'+description+'-final.bam '+cufflinksref[genome_build],'-s yes'])
            elif strand == 'reverse':
                cl.extend(['--library-type','fr-secondstrand'])
                counts.extend([exptfolder+'/'+description+'-final.bam '+cufflinksref[genome_build],'-s reverse'])
            else:
                cl.extend(['--library-type','unstranded'])
                counts.extend([exptfolder+'/'+description+'-final.bam '+cufflinksref[genome_build],'-s no'])
                
            counts.append(' > '+exptfolder+'/'+description+'.counts')
            cl.append(exptfolder+'/'+description+'.cufflinks')
            f.write(' '.join(cl)+'\n')
            f.write(' '.join(counts)+'\n')
        f.write('mv '+fastqpath+'*trimming_report.txt '+exptfolder+'/trimgalore\n')
        f.write('rm '+fastqpath+'*_val*\n')

    if opts.local:
        subprocess.run(['./'+exptfolder+'/run_align.sh'],shell=True)
    else:
        subprocess.run(["qsub -m e -M jhammelm@mit.edu -pe slots.pe "+str(opts.nthreads)+" -v PATH=$PATH -wd $PWD -N align_"+description+'_'+opts.aligner+' ./'+exptfolder+'/run_align.sh'],shell=True)
    
