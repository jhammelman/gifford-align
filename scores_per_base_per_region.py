#!/bin/env python
import pyBigWig
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('bed')
parser.add_argument('bigwig')
parser.add_argument('outprefix')
parser.add_argument('-sort','--sort',action='store_true',default=False)
parser.add_argument('-bin','--binsize',default=1,type=int)
parser.add_argument('-trim','--trim_to',default=0,type=int)
opts = parser.parse_args()

bw = pyBigWig.open(opts.bigwig)
all_regions = []
all_beds = []
for line in open(opts.bed):
    cols = line.strip().split('\t')
    if opts.trim_to != 0:
        mid = int((int(cols[2]) - int(cols[1]))/2)
        start = int(mid - opts.trim_to/2) + int(cols[1])
        end = int(mid + opts.trim_to/2) + int(cols[1])
    else:
        start = int(cols[1])
        end = int(cols[2])
    vals = np.array(bw.values(cols[0],start,end))
    if opts.binsize != 1:
        all_regions.append(np.mean(vals.reshape((-1,opts.binsize)),axis=1))
    else:
        all_regions.append(vals)
    all_beds.append('\t'.join([cols[0],str(start),str(end)]))
all_regions = np.array(all_regions)
if opts.sort:
    all_regions_sorted = all_regions[np.argsort(np.nansum(all_regions,axis=1))[::-1],:]
    all_beds = [all_beds[i] for i in np.argsort(np.nansum(all_regions,axis=1))[::-1]]
    with open(opts.outprefix+'.values_sorted.bed', 'w') as f:
        f.write('\n'.join(all_beds))
    np.savetxt(opts.outprefix+'.values_sorted.txt',all_regions_sorted)
    np.savetxt(opts.outprefix+'.values_sum.txt',np.nansum(all_regions_sorted,axis=0)/all_regions_sorted.shape[0])
else:
    np.savetxt(opts.outprefix+'.values.txt',all_regions)
    np.savetxt(opts.outprefix+'.values_sum.txt',np.nansum(all_regions,axis=0)/all_regions.shape[0])
    
