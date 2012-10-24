#!/usr/bin/python
description = '''Simulate IBD segments in a diploid population.
    infile should contain code defining
       sampsizes -- a dict whose values are the sample sizes
       ancnefn(pop,t) -- a function returning a dict whose keys are population names and whose values are effective population sizes in generation t
       migprobs(pop,t) -- a function returning a dict whose keys are pairs (tuples) of population names (x,y) 
           and whose values are backwards migration rates from x to y in generation t
    Here pop is a dict whose keys give the names of the populations (the same as the keys of sampsizes)
       and if ancnefn or migprobs can also be dicts rather than functions (i.e. constant).
    Some of these can be passed in on the command line.

    Example:
          python sim-ibd-pedigree.py -i sim-demographics-2.py -t 10 -b test.ibd.gz -l test.log
    '''

from optparse import OptionParser
import coalpedigree as coal
import re
import time
import subprocess, os, sys
import signal
# import pdb


parser = OptionParser(description=description)
# parser.add_option("-c","--coalfile",dest="coalfile",help="name of file to write final coalescent info to (or '-' for stdout)",default="-")
parser.add_option("-b","--ibdfile",dest="ibdfile",help="name of file to write final ibd blocks to (or '-' for stdout)",default="-")
parser.add_option("-l","--logfile",dest="logfile",help="name of log file (or '-' for stdout)",default="-")
parser.add_option("-i","--infile",dest="infile",help="name of input file to get parameters from (or '-' for stdin)")
parser.add_option("-t","--ngens",dest="ngens",help="total number of generations to simulate",default="10")
parser.add_option("-n","--nesize",dest="nesize",help="effective population size")
parser.add_option("-m","--migprob",dest="migprob",help="migration probability")
parser.add_option("-s","--samplesizes",dest="sampsizes",help="sample sizes")
parser.add_option("-e","--minlen",dest="minlen",help="minimum length of IBD block to record IN MORGANS (default value 0.005M = 0.5cM)",default=None)
parser.add_option("-g","--gaplen",dest="gaplen",help="gap length, IN MORGANS: blocks closer together than this will be recorded even if shorter than minlen (default value 0.5M = 50cM)",default=None)
(options,args) =  parser.parse_args()


# read in parameters etc
if options.infile is not None:
    infile = coal.fileopt(options.infile,"r")
    inparams = infile.read()
    exec(inparams)

# command line options supercede statements in infile
ibdfile = coal.fileopt(options.ibdfile, "w")
# coalfile = coal.fileopt(options.coalfile, "w")
logfile = coal.fileopt(options.logfile, "w")
ngens = int(options.ngens)

if options.nesize is not None:
    ancnefn = lambda pop,t: {}.fromkeys(pop.keys(),int(options.nesize))
else:
    try:
        if type(ancnefn) == type({}):
            ancnefn = lambda pop,t: ancnefn
    except NameError:
        print("Effective pop sizes (ancne) needs to be specified on the command line (-n) or in the input file (-i).")
        raise

if options.migprob is not None:
    migprobs = lambda pop,t: {}.fromkeys([(x,y) for x in pop.keys() for y in pop.keys()],float(options.migprob))
else:
    try:
        if type(migprobs) == type({}):
            migprobs = lambda pop,t: migprobs
    except NameError:
        print("Migration rates (migprobs) need to be specified on the command line (-m) or in the input file (-i).")
        raise


if options.sampsizes is not None:
    sampsizes = map(int,re.split("[, ]*",options.sampsizes))
    sampsizes = dict( zip(range(len(sampsizes)),sampsizes) )
else:
    try:
        if not type(sampsizes) ==type({}):
            raise TypeError("sampsizes should be a dict!")
    except NameError:
        raise TypeError("Sampsizes must be set in infile or on command line.")

if options.minlen is not None:
    minlen = float( options.minlen )
try:
    if not type(minlen) == type(0.0):
        raise TypeError("minlen is not defined: " + str(minlen))
except NameError:
    # default value
    minlen = 0.005

if options.gaplen is not None:
    gaplen = float( options.gaplen )
try:
    if not type(gaplen) == type(0.0):
        raise TypeError("gaplen is not defined: " + str(gaplen))
except NameError:
    # default value
    gaplen = 0.5

# catch ctrl-c gracefully
_exitnow = []

def catch_int(signal,frame):
    _exitnow.append(True)
    if len(_exitnow)>1:
        logfile.write("Caught SIGINT twice, terminating immediately.\n")
        logfile.flush()
        raise SystemExit
    else:
        logfile.write("Caught SIGINT, exiting after this generation.  SIGINT again to terminate.\n")
        logfile.flush()
        pass

signal.signal( signal.SIGINT, catch_int )

# initialize
pop = coal.initpop(sampsizes)

# sanity checks
mignames = reduce( lambda x,y: x+y, [ [u,v] for (u,v) in migprobs(pop,t=1).keys() ] )
ancnenames = ancnefn(pop,t=1).keys()
subpopnames = pop.keys()
if (not set(mignames) == set(ancnenames)) or (not set(subpopnames) == set(ancnenames) ):
    raise TypeError("Inconsistent population names -- using command-line samplesizes?")

# record "version number"
githash, giterr = subprocess.Popen(["git",'--git-dir='+os.path.abspath(os.path.dirname(sys.argv[0])+'/.git'),'rev-parse','HEAD'], stdout=subprocess.PIPE).communicate()
if giterr:
    githash = "(none available)"
logfile.write("sim-ibd-pedigree.py -- githash " + githash + time.strftime("%d %h %Y %H:%M:%S", time.localtime()) + "\n")
logfile.write("\n")
logfile.write("options "+str(options)+"\n")
logfile.write("\n")
# logfile.write("coal output: " + str(options.coalfile)+"\n")
logfile.write("ibd output: " + str(options.ibdfile)+"\n")
logfile.write("input: " + str(options.infile)+"\n")
logfile.write("--------------------------------\n")
logfile.write(inparams)
logfile.write("--------------------------------\n")
logfile.write("ngens: " + str(ngens)+"\n")
logfile.write("sampsizes: " + str(sampsizes)+"\n")
logfile.write("minlen: " + str(minlen)+"\n")
logfile.write("gaplen: " + str(gaplen)+"\n")
logfile.write("Ne at 1: " + str(ancnefn(pop,t=1))+"\n")
logfile.write("migprob at 1: " + str(migprobs(pop,t=1))+"\n")
logfile.write("chromosome ending positions: " + str(list(coal.chrpos)+[coal.chrlen]) + "\n")
logfile.write("\n")
logfile.write("Beginning ------------\n")

for t in xrange(ngens):
    logfile.write("  t="+str(t)+"\n")
    if t%10==0:
        logfile.write("    census (num indivs, num segments): " + str(coal.census(pop))+ "\n")
        logfile.flush()
    coal.parents(pop,ancne=ancnefn(pop,t),migprobs=migprobs(pop,t),t=t)
    if _exitnow:
        # there's been a ctrl-c; stop now.
        break

logfile.write("    census (num indivs, num segments): " + str(coal.census(pop))+ "\n")
logfile.write("Done with simulation at " + time.strftime("%d %h %Y %H:%M:%S", time.localtime()) + "; now writing out IBD info.\n" )

coal.writeibd(pop,minlen=minlen,gaplen=gaplen,outfile=ibdfile)
# writecoal(ibdict,outfile=coalfile)
# pdb.set_trace()

logfile.write("\n")
logfile.write("Closing ibd file...\n")

ibdfile.close()

logfile.write("All done at " + time.strftime("%d %h %Y %H:%M:%S", time.localtime()) + "\n" )
# coalfile.close()
logfile.close()
