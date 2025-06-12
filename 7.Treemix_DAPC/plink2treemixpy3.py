#!/usr/bin/env python3

import sys
import os
import gzip

if len(sys.argv) < 3:
    print("plink2treemix.py [gzipped input file] [gzipped output file]")
    print("ERROR: improper command line")
    exit(1)

infile = gzip.open(sys.argv[1], mode='rt', encoding='utf-8', errors='replace')
outfile = gzip.open(sys.argv[2], mode='wt', encoding='utf-8')

pop2rs = dict()
rss = []
rss2 = set()

# Skip the first header line
line = infile.readline()
line = infile.readline()

while line:
    line = line.strip().split()
    rs = line[1]
    pop = line[2]
    mc = line[6]
    total = line[7]
    if rs not in rss2:
        rss.append(rs)
    rss2.add(rs)
    if pop not in pop2rs:
        pop2rs[pop] = dict()
    if rs not in pop2rs[pop]:
        pop2rs[pop][rs] = f"{mc} {total}"
    line = infile.readline()

pops = list(pop2rs.keys())
outfile.write(" ".join(pops) + "\n")

for rs in rss:
    for pop in pops:
        tmp = pop2rs[pop][rs].split()
        c1 = int(tmp[0])
        c2 = int(tmp[1])
        c3 = c2 - c1
        outfile.write(f"{c1},{c3} ")
    outfile.write("\n")

infile.close()
outfile.close()
