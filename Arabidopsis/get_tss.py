### script to create tss file

import sys

data = open(sys.argv[1])

with open('tss_arabidopsis.txt', 'w') as tss:
    tss.write("{}\t{}\t{}\t{}\t{}\t{}\n".format('CHROM','FRAG_START','FRAG_STOP', 'STRAND', 'DESCRIPTION', 'MODEL'))
    for line in data.readlines():
        chrom, _, kind, start, stop, _, strain, _,  name = line.strip().split('\t')
        if kind == "protein":
            model = name.split(';')[1].split("=")[1]
            tss.write("{}\t{}\t{}\t{}\t\"{}\"\t{}\n".format(chrom, start, stop, strain,name,model))

