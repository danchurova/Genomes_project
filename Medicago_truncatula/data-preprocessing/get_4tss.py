### script to create 4tss file (analogous to 4tts also)
## medicago truncatula
import sys

data = open(sys.argv[1])

with open('../data/outputs/4tss_mrnas_mt.txt', 'w') as tss:
    tss.write("{}\t{}\t{}\t{}\t{}\t{}\n".format('CHROM','FRAG_START','FRAG_STOP', 'STRAND', 'DESCRIPTION', 'MODEL'))
    for line in data.readlines():
        if line[0] != '#':
            chrom, _, kind, start, stop, _, strain, _,  name = line.strip().split('\t')
            if kind == "mRNA":
                model = name.split(';')[1].split("=")[1]
                tss.write("{}\t{}\t{}\t{}\t\"{}\"\t{}\n".format(chrom, start, stop, strain, name, model))

