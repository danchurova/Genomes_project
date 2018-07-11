### script to create 4tss file (analogous to 4tts also) for homo sapiens release 28

import sys

data = open(sys.argv[1])

with open('./data/4tss_transcripts_hs.txt', 'w') as tss:
    tss.write("{}\t{}\t{}\t{}\t{}\t{}\n".format('CHROM','FRAG_START','FRAG_STOP', 'STRAND', 'DESCRIPTION', 'MODEL'))
    for line in data.readlines():
        if line[0] == '#':
            continue
        else:
            chrom, _, kind, start, stop, _, strain, _,  name = line.strip().split('\t')
        if kind == "transcript":
            if name.split(';')[4].split('=')[1] == 'protein_coding':
                model = name.split(';')[1].split("=")[1]
                tss.write("{}\t{}\t{}\t{}\t\"{}\"\t{}\n".format(chrom, start, stop, strain, name, model))

