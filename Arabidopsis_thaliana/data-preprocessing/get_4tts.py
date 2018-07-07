## script to create 4tts file
import sys

data = open(sys.argv[1])

with open('./data/4tts_arab.txt', 'w') as tts:
    tts.write("{}\t{}\t{}\t{}\t{}\t{}\n".format('CHROM', 'FRAG_START', 'FRAG_STOP', 'STRAND', 'DESCRIPTION', 'MODEL'))
    for line in data.readlines():
        chrom, _, kind, _, start, stop, _, strain, _, name = line.strip().split('\t')
        if kind == 'mRNA':
            model = name.split(';')[1].split('=')[1]
            tts.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, stop, start, strain, name, model))
