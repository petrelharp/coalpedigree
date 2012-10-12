#!/usr/bin/python
# coalpedigree, copyright Peter L. Ralph, 2012

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# simulate backwards through a pedigree
# individuals are encoded as
#   location (x), generation (t)
#   samples inheriting on which contiguous piece of chromosome (n,a,b)
# each generation, coalesce, recombine, migrate.

import random
import gzip
import sys
import bisect
from collections import defaultdict
import heapq
import itertools as it

# chrlen = 1.0 # chromosome length
# Default: human chromosome lengths (almost all of them)
chrlens = ( 2.6200830, 2.4403647, 2.0997332, 1.9748005, 1.8958412, 1.7470526, 1.7235610, 1.5652810, 1.4628861, 1.6040336, 1.4388007, 1.5640816, 1.1934065, 1.0483476, 1.1093523, 1.1868929, 1.1961966, 1.0574726, 0.9257134, 0.8362057, 0.5479314, 0.5559191 )  # default chromosome lengths (in Morgans)
chrpos = tuple( [ sum( chrlens[0:k] ) for k in range(1,len(chrlens)) ] )   # cumulative sum: position if lined up end-to-end
chrlen = sum(chrlens)  # the last one (total length)

# this is actually the number of parents... so it will work (?) for ploidy=4 but will not be biological.
ploidy = 2

def initpop(sampsizes):
    '''A pop is dict of: 
       location: [ dict of [ maternal chrom, paternal chrom ]s ]
        keys in each subpop dict are integers
        chromsomes are two lists, [ pos_i ] (in order) and [ anc_i ], with anc_i the ancestral chromosome of [pos_i,pos_i+1).
       Sampled diploid individual n has chromosomes numbers ploidy*n ... (ploidy+1)*n-1.
    '''
    pop = {}.fromkeys( sampsizes.keys() )
    k = 0
    for x in pop:
        pop[x] = {}.fromkeys( xrange(sampsizes[x]) )
        for ind in pop[x]:
            pop[x][ind] =  [ [ [ 0.0 ],  [ k+j ] ] for j in range(ploidy) ]
            k += ploidy
    return pop


def census(pop):
    '''Return number of individuals and number of breakpoints
    in each subpopulation.
    '''
    nindivs = [ len(pop[x]) for x in pop ]
    nchroms = [ sum( map( lambda u: len(u[0][0])+len(u[1][0]), pop[x].values() ) ) for x in pop ]
    return zip( nindivs, nchroms )


def sanity(pop,print_details=False):
    '''Do some sanity checks.
    '''
    errors = []
    for x in pop:
        for ind in pop[x]:
            try:
                # diploid?
                if len(pop[x][ind]) != ploidy:
                    errors.append( (x,ind) )
                    if print_details:
                        print "Oops!  Are we not diploid?"
                        print x, ind, pop[x][ind]
                # all chromosomes of the proper form?
                for chrom in pop[x][ind]:
                    if len(chrom)!=2 or len(chrom[0])!=len(chrom[1]) or any( [ (chrom[0][k+1]<=chrom[0][k]) for k in xrange(len(chrom[0])-1) ] ):
                        errors.append((x,ind))
                        if print_details:
                            print "Malformed chromosomes?"
                            print x, ind, pop[x][ind]


            except:
                print "Something bad just happened with"
                print pop[x][ind]
                raise
    if not errors:
        print "All good!", census(pop)
    else:
        print len(errors), "errors, oh dear!"
        return errors


def combine(chroms,ibd):
    '''Given a list of diploid individuals, combine them into a single one.
       For each of the maternal/paternal chromosomes,
       move through all contributing chromosomes in parallel,
       updating the current state for each, and combining them into the output state.

       When at least two nonempty states are combined,
       this is a coalescence;
       record the collection of states and the position of the segment in the list 'ibd'.
    '''
    if len(chroms)==1:
        return chroms[0]
    matpat = [[],[]]
    for ii in range(ploidy):
        whereat = [ 0 for chrom in chroms ]  # current locations on each chromosome
        nextones = [ ( chroms[k][ii][0] if len(chroms[k][ii])>0  else None ) for k in range(len(chroms)) ]
        thisstate = [ set() for u in nextones ]
        while any( nextones ):
            nextpos = min( [ u[0] for u in nextones if u ] )
            for j in range(len(whereat)):
                if nextones[j]:
                    if (nextones[j][0] == nextpos):
                        thisstate[j] = nextones[j][1]
                        whereat[j] += 1
                        try:
                            nextones[j] = chroms[j][ii][whereat[j]] 
                        except IndexError:
                            # i.e. we've got to the end of chroms[j][ii]
                            nextones[j] = None
            combining = [ x for x in thisstate if x ]
            combstate = set().union(*combining)
            if ibd and len( combining ) > 1:
                nextnextpos = min( [ u[0] for u in nextones if u ] )
                ibd.append( (nextpos,nextnextpos,combining) )
            matpat[ii].append( (nextpos, combstate) )
    return matpat

def getrecombs():
    # assign crossover locations -- begin with random subset of breaks between chromosomes
    # crossovers will be a list of crossover locations in increasing order
    crossovers = [ x for x in chrpos ]
    z = random.expovariate(1)
    while z < chrlen :
        crossovers.append(z)
        z += random.expovariate(1)
    crossovers.sort()
    return crossovers

def parents(pop,ancne,migprobs,t=0,ibdict=None,writeto=None):
    '''Reallocate each chromosome to the ancestors,
    and then recombine within each individual to resolve the mat/pat chromosomes.
    Optionally, write out to writeto rather than recording in ibdict.
      Note that combine() actually does the recording of events of interest in ibdict.
    '''
    cumne = [ ancne[x] for x in pop ]
    cumne = dict( zip( pop.keys(), [ sum(cumne[:k]) for k in range(len(cumne)) ] ) )

    parentdict = defaultdict()  # this is of the form chrom: parent, with parent diploid
    recombdict = defaultdict(getrecombs)  # this is of the form parent: [recomb locs]
    for x in pop:
        if len(pop)>1:
            cumprobs = [ migprobs.get((x,y),0) for y in pop.keys() ] 
            cumprobs = [ sum( cumprobs[:(k+1)] ) for k in range(len(cumprobs)) ]
            loclabs = pop.keys() + [x]

            def pickparent():
                y = loclabs[ bisect.bisect_left( cumprobs, random.random() )]
                return cumne[y] + random.sample(xrange(ancne[y]),1)[0]
        else:
            def pickparent():
                return random.sample(xrange(ancne.values()[0]),1)[0]

        parentdict.default_factory = pickparent

        for ind in pop[x]:
            for ii in xrange(ploidy):
                pos,anc = pop[x][ind][ii]
                newpos = []
                newanc = []
                for k in xrange(len(pos)):
                    mapa = parentdict[ anc[k] ]
                    recombs = recombdict[mapa]
                    # recombs[0] < ... < recombs[whichseg-1] < pos[k] <= recombs[whichseg]
                    whichseg = bisect.bisect_left( recombs, pos[k] )
                    if whichseg == len(recombs) or pos[k] != recombs[whichseg]:
                        newpos.append( pos[k] )
                        newanc.append( mapa*ploidy + (whichseg%ploidy) )
                    while whichseg < len(recombs) and ( k == (len(pos)-1) or recombs[whichseg] < pos[k+1] ):
                        newpos.append( recombs[whichseg] )
                        whichseg += 1
                        newanc.append( mapa*ploidy + (whichseg%ploidy) )
                pop[x][ind][ii] = newpos,newanc

def writecoal(ibdict,filename="coalpedigree.coal.gz",outfile=None,closeafter=True,writeheader=True):
    '''Write out all info in ibdict into a file.
    The file is of the form:
        time, unique identifying number, chromsome position start, end, population label of ancestral individual, and coalescing ids,
    where ids is the string representation of the list of sets of coalescing ids
    e.g. '[[2,3],[9],[1,10]]' means that 1,2,3,9,10 share this segment in a common ancestor for the first time,
    but that 2 and 3 had "already" inherited the segment in a more recent ancestor, as did 1 and 10.
    '''
    if outfile is None:
        outfile = fileopt(filename,"w")
    blockid = 0
    if writeheader:
        header = ["t", "blockid", "start", "end", "location", "ids"]
        outfile.write(" ".join(header)+"\n") 
    for t in ibdict:
        for x in ibdict[t]:
            for start,end,ids in ibdict[t][x]:
                outfile.write(" ".join(map(str,[t,blockid,start,end,x])+[str(map(list,ids)).translate(None," ")])+"\n")
                blockid += 1
    if closeafter:
        outfile.close()
    print blockid, "blocks written"


def writeibd(pop,minlen=0.0,gaplen=0.0,filename="coalpedigree.ibd.gz",simplify=True,outfile=None):
    '''Write out pairwise IBD info from the output of collectibd,
    restricting to segments at least minlen long OR at least as close as gaplen to another segment,
    and optionally simplifying by pasting together adjacent blocks shared with the same individual.
    Note that id1,id2 are in sorted order to avoid redundancy.

    Do this by stepping along the genome in parallel along all chromosomes.
    '''
    if outfile is None:
        outfile = fileopt(filename,"w")
        newfile = True
    else:
        newfile = False
    header = ["id1", "id2", "start", "end"]
    outfile.write(" ".join(header)+"\n") 
    chromlist = []
    for x in pop:
        for ind in pop[x]:
            chromlist.extend( [ it.izip(pos,anc) for pos,anc in pop[x][ind] ] )
    shared = {}.fromkeys(xrange(len(chromlist)))  # ind : { other:pos } where began sharing with other at pos
    shortones = {}.fromkeys(xrange(len(chromlist)))  # pending short ones that might be appended if there's another soon
    for ii in shared:
        shared[ii] = {}
        shortones[ii] = defaultdict(lambda : [])
    thestack = []
    for j in xrange(len(chromlist)):
        heapq.heappush( thestack, (chromlist[j].next(),j) )
    # initialize current ancestral individual states of each chromosome
    curstate = [ None ]*len(chromlist)
    while thestack:
        (pos,anc),curi = heapq.heappop( thestack )
        curstate[curi] = anc
    # sanity
    if filter( lambda x: x is None, curstate ):
        raise ValueError("writeibd: Wrong number of 0.0s found?!?")
    # reinit
    for j in xrange(len(chromlist)):
        try:
            heapq.heappush( thestack, (chromlist[j].next(),j) )
        except StopIteration:
            # oops, nothing more in chromlist[j]
            pass
    # initialize sharing segments
    for i,j in it.ifilter( lambda x: curstate[x[0]]==curstate[x[1]], it.product(xrange(len(curstate)),xrange(len(curstate))) ):
        if i != j:
            shared[i][j] = shared[j][i] = 0.0
    while thestack:
        # if anc has changed, must end segments j participates in at pos and begin new ones,
        #   updating curstate and shared
        (pos,anc),curi = heapq.heappop(thestack)
        if curstate[ curi ] != anc:
            for j,start in shared[curi].items():
                if shortones[curi][j]:
                    dogaps = [ shend > start-gaplen for shstart,shend in shortones[curi][j] ]
                    for shstart,shend in it.compress(shortones[curi][j],dogaps):
                        outfile.write( " ".join(map(str,[curi,j,shstart,shend])) + "\n" )
                    del shortones[j][curi]
                    del shortones[curi][j]
                else:
                    dogaps = [False]
                if any(dogaps) or pos-start > minlen:
                    outfile.write( " ".join(map(str,[curi,j,start,pos])) + "\n" )
                else:
                    shortones[curi][j].append( (start,pos) )
                    if not curi in shortones[j]:
                        shortones[j][curi] = shortones[curi][j]
                del shared[curi][j]
                del shared[j][curi]
            for k in [ k for k,x in it.izip(xrange(len(curstate)),curstate) if x==anc ]:
                shared[curi][k] = shared[k][curi] = pos
            curstate[curi] = anc
        try:
            heapq.heappush( thestack, ( chromlist[curi].next(), curi ) )
        except StopIteration:
            pass
    # finish off dangling blocks
    for i in shared.keys():
        for j,start in shared[i].items():
            if i<j:
                if shortones[i][j]:
                    dogaps = [ shend > start-gaplen for shstart,shend in shortones[i][j] ]
                    for shstart,shend in it.compress(shortones[i][j],dogaps):
                        outfile.write( " ".join(map(str,[i,j,shstart,shend])) + "\n" )
                else:
                    dogaps = [False]
                if any(dogaps) or pos-start > minlen:
                    outfile.write( " ".join(map(str,[i,j,start,chrlen])) + "\n" )
    # all done!
    if newfile:
        outfile.close()


def fileopt(fname,opts):
    '''Return the file referred to by fname, open with options opts;
    if fname is "-" return stdin/stdout; if fname ends with .gz run it through gzip.
    '''
    if fname == "-":
        if opts == "r":
            fobj = sys.stdin
        elif opts == "w":
            fobj = sys.stdout
        else:
            print "Something not right here."
    elif fname[len(fname)-3:len(fname)]==".gz":
        fobj = gzip.open(fname,opts)
    else:
        fobj = open(fname,opts)
    return fobj
