#!/usr/bin/python
description='''Also, merge any adjacent blocks that are separated by no more than the distance gap_threshold.
Usage: python remove-gaps-fibd.py gap_threshold file.fibd.gz (pipe output to selected file)
Output is one line per tract per pair of individuals, with positionsl of strt and end of tract.
'''
# highly modified by plr from process-ibd.py by Browning

import sys, gzip
from math import isnan
import re
import coalpedigree as coal
from optparse import OptionParser

parser = OptionParser(description=description)
parser.add_option("-b","--ibdfile",dest="ibdfile",help="name of file to read ibd from (or '-' for stdin)",default="-")
parser.add_option("-o","--outfile",dest="outfile",help="name of file to output merged ibd from (or '-' for stdout)",default="-")
parser.add_option("-g","--gaplen",dest="gaplen",help="merge blocks separated by a block no longer than this long", default=0.0)
parser.add_option("-c","--chromfile",dest="chromfile",help="file to parse chromosome lengths from (e.g. logfile for that run)", default=None)
(options,args) =  parser.parse_args()

# size of maximum gap to merge
gapthresh = float(options.gaplen)

# find chromosome lengths
try:
    if options.chromfile is None:
        chromfile = coal.fileopt( re.sub("\.fibd\.gz$",".log",options.ibdfile), "r")
    else:
        chromfile = coal.fileopt(options.chromfile,"r")
except:
    raise ValueError("Invalid chromfile specification.")

for line in chromfile:
    if line.startswith("chromosome ending positions:"):
        chrends = eval( line[28:].strip() )
        break
try:
    chrstartends = zip( [0.0]+chrends[:-1], chrends )
except NameError:
    raise ValueError("Can't find 'chromosome ending positions:' in chromfile.")
    

results = {}
infile = coal.fileopt(options.ibdfile,"r")
outfile = coal.fileopt(options.outfile,"w")
for line in infile:
    # note: what is read is NOT the "length" but it is different (by one position) than the output.  
    # Beagle reports the starting (inclusive) and ending (exclusive) positions.
    (id1, id2, start, end) = line.split()
    if id1 == "id1":
        # oops this is the header
        continue
    start = float(start)
    end = float(end)
    if (id1,id2) in results:
        currentlist = results[(id1,id2)]
        markremove = [False for x in currentlist]
        for j in range(len(currentlist)):
             x = currentlist[j]
             ovlap = x[0]<=end and start<=x[1]
             # is a gap?
             gap = ( x[0] <= end+gapthresh ) and ( start <= x[1] + gapthresh )
             # is shorter than an adjacent segment?
             # gaplen is max( start1-end2, start2-end1 ) since the min is negative
             gap = gap and ( max( x[0]-end, start-x[1] ) < max( end-start, x[1]-x[0] ) )
             if ovlap or gap:
                 # overlap
                 markremove[j] = True
                 start = min(x[0],start)
                 end = max(x[1],end)
        if sum(markremove):
             results[(id1,id2)] = [x for x,r in zip(currentlist,markremove) if not r]
        results[(id1,id2)].append([start,end])
    else:
        results[(id1,id2)] = [[start,end]]

# Note: changed 'minscore' to 'score'.
outfile.write("id1 id2 start end\n")
for x in results:
    # id1 = x[0]; id2 = x[1]
    # on the same chromosome?
    for y in results[x]:
        # start,end = y
        breakpoints = [ y[0] ] + [ z for z in chrends if y[0] < z and z < y[1] ] + [ y[1] ]
        for k in xrange(len(breakpoints)-1):
            outfile.write( " ".join(map(str,x) + map(str,y[k:(k+2)])) + "\n" )

outfile.close()
