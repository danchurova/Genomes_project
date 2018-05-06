#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 15:26:19 2018

@author: mac
"""
import sys


with open (sys.argv[1]) as vcf:
    #with open('SNPS.txt', 'w') as SNPS:
    counter = 0
    print('"global"\t"coords"\t\t"Chr_snp"\t"from"\t"to"\t"maf"\t"nchrobs"\t"cov"\t"sum"')
    for line in vcf:
        if line[0:1] == '#':
            continue
        else:
            line_list = line.split('\t')
            if line_list[4] != "." and len(line_list[4]) == 1:
                counter += 1
                print('{}\t{}\t{}\t{}\t"{}"\t"{}"'.format(counter, line_list[1], line_list[1], line_list[0],line_list[3], line_list[4]) + '\t1'*4)

