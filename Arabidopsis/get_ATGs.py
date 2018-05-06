### script to create ATGs file

import sys

data = open(sys.argv[1])

with open('ATGs_arabidopsis.txt', 'w') as atgs:
    atgs.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t{}\t{}\n".format('CHROM','FRAG_START','FRAG_STOP', 'STRAND', 'DESCRIPTION', 'MODEL', 'ATG', 'INCLUDE?','HAS TSS', 'IN'))
    for line in data.readlines():
        chrom, _, kind, start, stop, _, strain, _,  name = line.strip().split('\t')
        if kind == "protein":
            model = name.split(';')[1].split("=")[1]
            prot_name = name
        if kind == "five_prime_UTR":
            if strain == "+":
                atg = start
            else:
                atg = stop
            atgs.write("{}\t{}\t{}\t{}\t\"{}\"\t{}\t{}\t{}\t\t{}\t{}\n".format(chrom, start, stop, strain,prot_name,model,atg,"1", name, "1"))

