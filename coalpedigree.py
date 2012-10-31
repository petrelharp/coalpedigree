#!/usr/bin/python
# coalpedigree, copyright Peter L. Ralph, 2012

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but without any warranty; without even the implied warranty of
#    merchantability or fitness for a particular purpose.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


import random
import gzip
import sys
import bisect
from collections import defaultdict, Counter
import heapq
import itertools as it

# chrlen = 1.0 # chromosome length
# Default: human chromosome lengths (almost all of them)
chrlens = ( 2.6200830, 2.4403647, 2.0997332, 1.9748005, 1.8958412, 1.7470526, 1.7235610, 1.5652810, 1.4628861, 1.6040336, 1.4388007, 1.5640816, 1.1934065, 1.0483476, 1.1093523, 1.1868929, 1.1961966, 1.0574726, 0.9257134, 0.8362057, 0.5479314, 0.5559191 )  # default chromosome lengths (in Morgans)
chrpos = tuple( [ sum( chrlens[0:k] ) for k in range(1,len(chrlens)) ] )   # cumulative sum: position if lined up end-to-end
chrlen = sum(chrlens)  # the last one (total length)

# this is actually the number of parents... so it will work (?) for ploidy=4 but will not be biological.
ploidy = 2

# this number should be larger than Ne will ever be
maxne = 10**8  # since maximum integer size on some machines is 2**31, this restricts us to at most 21 populations?

def initpop(sampsizes):
    '''A pop is list of [ maternal chrom, paternal chrom ]s
          chromsomes are two lists, [ pos_i ] (in order) and [ anc_i ], with anc_i the ancestral chromosome of [pos_i,pos_i+1).
       Sampled diploid individual n in population k has chromosomes numbers k*maxne + ploidy*n ... k*maxne + (ploidy+1)*n-1.
    '''
    ordlabs = sorted(sampsizes.keys())  # ensure consistent order
    if type(sampsizes)==type({}):
        sampsizes = [ sampsizes[x] for x in ordlabs ]
    diploids = [ k*maxne+j for k in xrange(len(sampsizes)) for j in xrange(sampsizes[k]) ]
    pop = [ [ [ [ 0.0 ],  [ maxne*(k//maxne)+ploidy*(k%maxne)+j ] ] for j in xrange(ploidy) ] for k in diploids ]
    return pop


def census(pop,sampsizes=None):
    '''Return number of individuals and number of breakpoints
    in each subpopulation.
    '''
    nchroms = sum( map( lambda u: len(u[0][0])+len(u[1][0]), pop ) )
    subpops = []
    if sampsizes:
        ordlabs = sorted(sampsizes.keys())  # ensure consistent order
        k = 0
        for nsamps in [ sampsizes[x] for x in ordlabs ]:
            thispop = Counter()
            for j in xrange(k,k+nsamps):
                for ii in xrange(ploidy):
                    thispop.update( getsubpop(pop[j][ii][1],ordlabs) )
            k = k+nsamps
            subpops.append( thispop )
    return nchroms, subpops


def getsubpop(inds,ordlabs):
    '''Which subpopulation is a chromosome in?'''
    return [ ordlabs[ ind // maxne ] for ind in inds ]


def sanity(pop,sampsizes=None,ancne=None,print_details=False):
    '''Do some sanity checks.
    '''
    errors = []
    for ind in xrange(len(pop)):
        try:
            # diploid?
            if len(pop[ind]) != ploidy:
                errors.append( ind )
                if print_details:
                    print "Oops!  Are we not diploid?"
                    print ind, pop[ind]
            # all chromosomes of the proper form?
            for chrom in pop[ind]:
                if len(chrom)!=2 or len(chrom[0])!=len(chrom[1]) or any( [ (chrom[0][k+1]<=chrom[0][k]) for k in xrange(len(chrom[0])-1) ] ):
                    errors.append(ind)
                    if print_details:
                        print "Malformed chromosomes?"
                        print ind, pop[ind]
        except:
            print "Something bad just happened with"
            print pop[ind]
            raise
    if not errors:
        print "All good!", census(pop,sampsizes)
    else:
        print len(errors), "errors, oh dear!"
        return errors


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


def parentfactory(ancne,migprobs):
    '''Return a function that will pick a random (diploid) parent for a given (diploid) individual.
    Note that migration probabilities are *reverse-time* migration probabilities,
    i.e. migprobs[(x,y)] is the probability that a parent of someone in x is from y.
    '''
    if len(ancne)>1:
        ordlabs = sorted(ancne.keys())  # ensure consistent order
        cumprobs = [ [ migprobs.get((x,y),0) for y in ordlabs ] for x in ordlabs ]
        # assign missing probabilities to diagonal
        totprobs = map( sum, cumprobs )
        if any([p>1 for p in totprobs]):
            raise ValueError("Migration probabilities sum to >1, oops!")
        for k in xrange(len(cumprobs)):
            cumprobs[k][k] += 1-totprobs[k]
            cumprobs[k] = [ sum(cumprobs[k][:j+1]) for j in xrange(len(cumprobs[k])) ]
        def pickparent(ind):
            y = bisect.bisect_left( cumprobs[ind//maxne], random.random() )
            return y*maxne + random.sample(xrange(ancne[ordlabs[y]]),1)[0]
    else:
        def pickparent(ind):
            return random.sample(xrange(ancne.values()[0]),1)[0]
    return pickparent


def parents(pop,ancne,migprobs,t=0,ibdict=None,writeto=None):
    '''Reallocate each chromosome to the ancestors,
    and then recombine within each individual to resolve the mat/pat chromosomes.
    Optionally, write out to writeto rather than recording in ibdict.
    '''
    parentdict = {}  # this is of the form chrom: diploid parent
    pickparent = parentfactory(ancne,migprobs)  # function to choose a parent
    recombdict = defaultdict(getrecombs)  # this is of the form diploid indiv: [recomb locs]
    for ind in xrange(len(pop)):
        for ii in xrange(ploidy):
            pos,anc = pop[ind][ii]
            newpos = []
            newanc = []
            for k in xrange(len(pos)):
                try:
                    # anc[k] is the current ancestral chromosome; mapa is the diploid parent that chromosome came from
                    mapa = parentdict[ anc[k] ]
                except KeyError:
                    mapa = parentdict[ anc[k] ] = pickparent(maxne*(anc[k]//maxne)+(anc[k]%maxne)//ploidy)  # pickparent wants a diploid index
                # each haploid (e.g. anc[k]) is the product of a unique meiosis (between chromosomes of mapa)
                recombs = recombdict[ anc[k] ]
                # recombs[0] < ... < recombs[whichseg-1] < pos[k] <= recombs[whichseg]
                whichseg = bisect.bisect_left( recombs, pos[k] )
                if whichseg == len(recombs) or pos[k] != recombs[whichseg]:
                    newpos.append( pos[k] )
                    newanc.append( maxne*(mapa//maxne) + (mapa%maxne)*ploidy + (whichseg%ploidy) )
                while whichseg < len(recombs) and ( k == (len(pos)-1) or recombs[whichseg] < pos[k+1] ):
                    newpos.append( recombs[whichseg] )
                    whichseg += 1
                    newanc.append( maxne*(mapa//maxne) + (mapa%maxne)*ploidy + (whichseg%ploidy) )
            pop[ind][ii] = newpos,newanc
    # all done!
    return None



def writeibd(pop,minlen=0.0,gaplen=0.0,filename="coalpedigree.ibd.gz",simplify=True,outfile=None):
    '''Write out pairwise IBD info from the output of collectibd,
    restricting to segments at least minlen long OR at least as close as gaplen to another segment,
    and optionally simplifying by pasting together adjacent blocks shared with the same individual.

    Do this by stepping along the genome in parallel along all chromosomes.
    '''
    if outfile is None:
        outfile = fileopt(filename,"w")
        newfile = True
    else:
        newfile = False
    header = ["id1", "id2", "start", "end"]
    outfile.write(" ".join(header)+"\n") 
    chromlist = []  # list of iterators along each chromosome
    for ind in xrange(len(pop)):
        chromlist.extend( [ it.izip(pos,anc) for pos,anc in pop[ind] ] )
    shared = {}.fromkeys(xrange(len(chromlist)))  # ind : { other:pos } where began sharing with other at pos
    shortones = {}.fromkeys(xrange(len(chromlist)))  # pending short ones that might be appended if there's another soon
    for ii in shared:
        shared[ii] = {}
        shortones[ii] = defaultdict(lambda : [])
    thestack = [] # next positions for each, sorted
    pending = {}  # next events occurring at the same junction position
    for j in xrange(len(chromlist)):
        heapq.heappush( thestack, (chromlist[j].next(),j) )
    # initialize current ancestral individual states of each chromosome
    curstate = [ None ]*len(chromlist)
    while thestack:
        (pos,anc),curi = heapq.heappop( thestack )
        curstate[curi] = anc
    # sanity
    if any( [ x is None for x in curstate ]):
        raise ValueError("writeibd: Wrong number of 0.0s found?!?")
    # refill thestack
    for j in xrange(len(chromlist)):
        try:
            heapq.heappush( thestack, (chromlist[j].next(),j) )
        except StopIteration:
            # ah, nothing more in chromlist[j]
            pass
    # initialize sharing segments
    for i,j in it.ifilter( lambda x: curstate[x[0]]==curstate[x[1]], it.product(xrange(len(curstate)),xrange(len(curstate))) ):
        if i != j:
            shared[i][j] = shared[j][i] = 0.0
    pos = 0.0 # where we're at now
    while pending or thestack:
        # load all events occurring at the same junction position into pending;
        #   then process each, not dealing with relationships to others still pending.
        # for each, if anc has changed, must end segments j participates in at pos and begin new ones,
        #   updating curstate and shared
        if (not pending) or (thestack and (thestack[0][0][0] == pos)):
            (pos,anc),curi = heapq.heappop(thestack)
            pending[curi] = anc
            try:
                heapq.heappush( thestack, ( chromlist[curi].next(), curi ) )
            except StopIteration:
                pass
        else:
            while pending:
                curi,anc = pending.popitem()
                if (pos in chrpos) or (curstate[ curi ] != anc):
                    for j,start in shared[curi].items():
                        if j not in pending and ((pos in chrpos) or (curstate[j]!=anc)):
                            # finish off blocks that are not pending
                            #   and no longer have the same state
                            #   ... but finish blocks regardless if we're at a chromosome boundary
                            if j in shortones[curi]:
                                dogaps = [ shend > start-gaplen for shstart,shend in shortones[curi][j] ]
                                for shstart,shend in it.compress(shortones[curi][j],dogaps):
                                    outfile.write( " ".join(map(str,[curi,j,shstart,shend])) + "\n" )
                                del shortones[j][curi]
                                del shortones[curi][j]
                            else:
                                dogaps = [False]
                            if any(dogaps) or pos-start > minlen:
                                # if any( [ (start-z)*(z-pos)>0 for z in chrpos ] ):
                                #     import pdb; pdb.set_trace()
                                outfile.write( " ".join(map(str,[curi,j,start,pos])) + "\n" )
                            else:
                                shortones[curi][j].append( (start,pos) )
                                if not curi in shortones[j]:
                                    shortones[j][curi] = shortones[curi][j]
                            del shared[curi][j]
                            del shared[j][curi]
                    for k in [ k for k in xrange(len(curstate)) if (k != curi) and (k not in pending) and (k not in shared[curi]) and (curstate[k]==anc) ]:
                        # begin new blocks with other individuals, only if they 
                        #  don't have a pending state change in the current set of operations,
                        #  aren't already sharing a block, and have the same state,
                        shared[curi][k] = shared[k][curi] = pos
                    # ok, now reset the current state
                    curstate[curi] = anc
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
