### script to create 4tss.txt file  for Drosophila melanogaster (analogous to 4tts also)

import sys

data = open(sys.argv[1]) #annotation in gff format

with open('../data/outputs/4tss_dm.txt', 'w') as tss:
    tss.write("{}\t{}\t{}\t{}\t{}\t{}\n".format('CHROM','FRAG_START','FRAG_STOP', 'STRAND', 'DESCRIPTION', 'MODEL'))
    for line in data.readlines():
        if line[0] == '>':
            break
        if line[0] != '#':
            chrom, _, kind, start, stop, _, strain, _,  name = line.strip().split('\t')
            if kind == "mRNA":
                model = name.split(';')[0].split("=")[1]
                #if model[:2] == 'CG': # this indicates coding genes
                tss.write("{}\t{}\t{}\t{}\t\"{}\"\t{}\n".format(chrom, start, stop, strain, name, model))

