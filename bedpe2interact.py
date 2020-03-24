#!/bin/bash

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('bedfile')
opts = parser.parse_args()

for line in open(opts.bedfile):
    data = line.strip().split()
    minstart = str(min(int(data[1]),int(data[4])))
    maxend = str(max(int(data[2]),int(data[5])))
    if data[0] == data[3]:
        print('\t'.join([data[0],
                        minstart,
                        maxend,
                        '.',
			'1000',
			'1000',
                        '.',
			'0',
                        data[0],
                        data[1],
                        data[2],
                        '.',
			'.',
                        data[3],
                        data[4],
                        data[5],
                        '.',
			'.']))
