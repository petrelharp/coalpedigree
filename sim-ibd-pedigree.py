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
    '''

from optparse import OptionParser
import coalpedigree as coal
import re
import time
import subprocess


parser = OptionParser(description=description)
parser.add_option("-o","--outfile",dest="outfile",help="name of output file (or '-' for stdout)",default="-")
parser.add_option("-p","--popfile",dest="popfile",help="name of file to write final population state to (or '-' for stdout)",default="-")
parser.add_option("-l","--logfile",dest="logfile",help="name of log file (or '-' for stdout)",default="-")
parser.add_option("-i","--infile",dest="infile",help="name of input file to get parameters from (or '-' for stdin)")
parser.add_option("-t","--ngens",dest="ngens",help="total number of generations to simulate",default="10")
parser.add_option("-n","--nesize",dest="nesize",help="default effective population size",default="10000")
parser.add_option("-m","--migprob",dest="migprob",help="default migration probability",default=".001")
parser.add_option("-s","--samplesizes",dest="sampsizes",help="sample sizes")
(options,args) =  parser.parse_args()


# read in parameters etc
if options.infile is not None:
    infile = coal.fileopt(options.infile,"r")
    inparams = infile.read()
    exec(inparams)

# command line options supercede statements in infile
outfile = coal.fileopt(options.outfile, "w")
popfile = coal.fileopt(options.popfile, "w")
logfile = coal.fileopt(options.logfile, "w")
ngens = int(options.ngens)

try:
    if type(ancnefn) == type({}):
        ancnefn = lambda pop,t: ancnefn
except NameError:
    ancnefn = lambda pop,t: {}.fromkeys(pop.keys(),int(options.nesize))

try:
    if type(migprobs) == type({}):
        migprobs = lambda pop,t: migprobs
except NameError:
    migprobs = lambda pop,t: {}.fromkeys([(x,y) for x in pop.keys() for y in pop.keys()],float(options.migprob))


if options.sampsizes is not None:
    sampsizes = map(int,re.split("[, ]*",options.sampsizes))
    sampsizes = dict( zip(range(len(sampsizes)),sampsizes) )
else:
    try:
        if not type(sampsizes) ==type({}):
            raise TypeError("sampsizes should be a dict!")
    except NameError:
        raise TypeError("Sampsizes must be set in infile or on command line.")

# initialize
pop = coal.initpop(sampsizes)

# sanity checks
mignames = reduce( lambda x,y: x+y, [ [u,v] for (u,v) in migprobs(pop,t=1).keys() ] )
ancnenames = ancnefn(pop,t=1).keys()
subpopnames = pop.keys()
if (not set(mignames) == set(ancnenames)) or (not set(subpopnames) == set(ancnenames) ):
    raise TypeError("Inconsistent population names -- using command-line samplesizes?")

# record "version number"
githash, giterr = subprocess.Popen(["git",'rev-parse','HEAD'], stdout=subprocess.PIPE).communicate()
if giterr:
    githash = ""
logfile.write("sim-ibd-pedigree.py -- " + githash + time.strftime("%d %h %Y %H:%M:%S", time.localtime()) + "\n")
logfile.write("\n")
logfile.write("output: " + str(options.outfile)+"\n")
logfile.write("final population state: " + str(options.popfile)+"\n")
logfile.write("input: " + str(options.infile)+"\n")
logfile.write("ngens: " + str(ngens)+"\n")
logfile.write("sampsizes: " + str(sampsizes)+"\n")
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
    pop = coal.parents(pop,t=t,ancne=ancnefn(pop,t),writeto=outfile)
    pop = coal.migrate(pop,migprobs=migprobs(pop,t))

popfile.write( str(pop) )

logfile.write("    census (num indivs, num segments): " + str(coal.census(pop))+ "\n")
logfile.write("\n")
logfile.write("All done at " + time.strftime("%d %h %Y %H:%M:%S", time.localtime()) + "\n" )
outfile.close()
popfile.close()
logfile.close()
