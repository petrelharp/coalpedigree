#!/usr/bin/python
description='''Take a file of IBD blocks and run it through something like the detection procudure:
    merge together nearby blocks and retain only those above a certain length.

    Example usage:
        python winnow.py -b this-ibdfile.fibd.gz -o this-ibdfile-winnowed.fibd.gz -l this-ibdfile-winnowed.log -g .05 -m .02
'''

from optparse import OptionParser
import coalpedigree as cp

parser = OptionParser(description=description)
parser.add_option("-b","--ibdfile",dest="ibdfile",help="name of file to read ibd from (or '-' for stdout)",default="-")
parser.add_option("-o","--outfile",dest="outfile",help="name of output ibd file (or '-' for stdin)")
parser.add_option("-l","--logfile",dest="logfile",help="name of log file (or '-' for stdout)",default="-")
parser.add_option("-n","--minminlen",dest="minminlen",help="totally ignore any blocks shorter than this length", default=0.0)
parser.add_option("-g","--gaplen",dest="gaplen",help="merge blocks separated by a block no longer than this long", default=0.0)
parser.add_option("-m","--minlen",dest="minlen",help="only keep any blocks (including merged ones) at least this long", default=0.0)
(options,args) =  parser.parse_args()

minminlen = float(options.minminlen)
gaplen = float(options.gaplen)
minlen = float(options.minlen)

ibdfile = cp.fileopt(options.ibdfile, "r")
outfile = cp.fileopt(options.outfile, "w")
logfile = cp.fileopt(options.logfile, "w")

results = {}
# results is a dict indexed by pairs of ids
# with entries a list of [start,end,nsegs]

logfile.write("winnow.py: "+options.ibdfile+"\n")
logfile.write("writing to " + options.outfile + "\n")
logfile.write("totally ignoring blocks below: " + str(minminlen) + "\n")
logfile.write("merging closer than gaplen: " + str(gaplen) + "\n")
logfile.write("outputting merged blocks longer than: " + str(minlen) + "\n")

nin = 0
nskip = 0
nout = 0

header = ibdfile.readline().split()

for line in ibdfile:
    (id1, id2, start, end) = line.split()
    start = float(start)
    end = float(end)
    if end-start < minminlen:
        nskip = nskip + 1
        continue
    nin = nin + 1
    if (id1,id2) in results:
        currentlist = results[(id1,id2)]
        markremove = [False for x in currentlist]
        for j in range(len(currentlist)):
            x = currentlist[j]
            ovlap = x[0]<=end and start<=x[1]
            # is a gap?
            gap = ( ( x[0]<=end+gaplen ) and ( start<=x[1]+gaplen ) )
            # is shorter than an adjacent segment?
            # gaplen is max( start1-end2, start2-end1 ) since the min is negative
            gap = ( gap and ( max( x[0]-end, start-x[1] ) < max( end-start, x[1]-x[0] ) ) )
            if ovlap or gap:
                # overlap
                markremove[j] = True
                start = min(x[0],start)
                end = max(x[1],end)
        nsegs = 1 + sum( [ x[2] for x,r in zip(currentlist,markremove) if r] )
        if sum(markremove):
            results[(id1,id2)] = [x for x,r in zip(currentlist,markremove) if not r]
        results[(id1,id2)].append([start,end,nsegs])
    else:
        results[(id1,id2)] = [[start,end,1]]

outfile.write("id1 id2 start end nsegs\n")
for x in results:
    id1 = x[0]; id2 = x[1]
    for y in results[x]:
        nout = nout+1
        outfile.write( " ".join(map(str,x) + map(str,y))+"\n" )

logfile.write("Done merging " + str(nin) + " blocks into " + str(nout) + " blocks; totally omitted " + str(nskip) + " blocks.\n" )
