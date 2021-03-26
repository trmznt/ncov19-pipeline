#!/usr/bin/env python3

import sys, gzip

depths = []
with gzip.open(sys.argv[1], 'r') as fin:
    for line in fin:
        d = int(line.strip().split()[-1])
        if d > 5:
            depths.append(d)

total = sum(depths)
if total == 0: average = 0
else: average = total / len(depths)
print('%s %d %d' % (sys.argv[2], average, len(depths) ))

