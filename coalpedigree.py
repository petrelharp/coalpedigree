#!/usr/bin/python

# simulate backwards through a pedigree
# individuals are encoded as
#   location (x), generation (t)
#   samples inheriting on which contiguous piece of chromosome (n,a,b)
# each generation, coalesce, recombine, migrate.

import random
import gzip
import sys
import bisect

# chrlen = 1.0 # chromosome length
chrlens = ( 2.6200830, 2.4403647, 2.0997332, 1.9748005, 1.8958412, 1.7470526, 1.7235610, 1.5652810, 1.4628861, 1.6040336, 1.4388007, 1.5640816, 1.1934065, 1.0483476, 1.1093523, 1.1868929, 1.1961966, 1.0574726, 0.9257134, 0.8362057, 0.5479314, 0.5559191 )  # chromosome lengths (in Morgans)
chrpos = tuple( [ sum( chrlens[0:k] ) for k in range(1,len(chrlens)) ] )   # cumulative sum: position if lined up end-to-end
chrlen = sum(chrlens)  # the last one (total length)

def initpop(sampsizes):
    '''A pop is dict of: 
       location: [ dict of [ maternal chrom, paternal chrom ]s ]
        keys in each subpop dict are integers
        chromsomes are of the form [ (pos_i,S_i) , ... ],
          sorted by position, with S_i the set of sample chromosomes inheriting from [pos_i,pos_i+1).
          The final position always has S_n = {}.
       Sampled diploid individual n has chromosomes numbers 2n and 2n+1.
    '''
    pop = {}.fromkeys( sampsizes.keys() )
    k = 0
    for x in pop:
        pop[x] = {}.fromkeys( [ ii for ii in xrange(sampsizes[x]) ] )
        for ind in pop[x]:
            pop[x][ind] =  [ [(0.0,set([k])),(chrlen,set([]))], [(0.0,set([k+1])),(chrlen,set([]))] ]
            k += 2
    return pop

def effpopsize(pop,t):
    return {}.fromkeys(pop.keys(),300)

def census(pop):
    '''Return number of individuals and number of breakpoints
    in each subpopulation.
    '''
    nindivs = [ len(pop[x]) for x in pop ]
    nchroms = [ len(reduce(lambda u,v:u+(v[0]+v[1]),pop[x].values(),[])) for x in pop ]
    return zip( nindivs, nchroms )

def nsegs(pop):
    '''Return how many segments are in each individual.
    '''
    nsegs = {}
    for x in pop:
        nsegs[x] = {}
        for ind in pop[x]:
            nsegs[x][ind] = sum( map( len, pop[x][ind] ) )
    return nsegs

def lengthspectrum(pop):
    '''Return the set of segment lengths in each subpopulation.
    '''
    seglens = {}.fromkeys(pop.keys())
    for x in pop:
        seglens[x] = reduce( lambda u,v: u+v, map( lambda chroms: [ b-a for a,b,k in chroms[0]+chroms[1] ], pop[x].values() ) )
        seglens[x].sort()
    return seglens

def sanity(pop,print_details=False):
    '''Do some sanity checks.
    '''
    ploidy = []
    errors = []
    # accumulate the bits of chromosome from each sample here, to see if any have gone missing
    pieces = {}
    for x in pop:
        for ind in pop[x]:
            try:
                # diploid?
                if len(pop[x][ind]) != 2:
                    errors.append( (x,ind) )
                    if print_details:
                        print "Oops!  Are we not diploid?"
                        print x, ind, pop[x][ind]
                # all chromosomes of the proper form?
                if any( [any(map(lambda w: len(w)!= 2 or type(w[0])!=type(0.0) or type(w[1])!=type(set([])), z)) for z in pop[x][ind]] ):
                    errors.append( (x,ind) )
                    if print_details:
                        print "Oops: malformed individual!"
                        print x, ind, pop[x][ind]
                for chrom in range(len(pop[x][ind])):
                    # do we still have whole chromsomes?
                    for j in xrange(len(pop[x][ind][chrom])-1):
                        a,k = pop[x][ind][chrom][j]
                        b = pop[x][ind][chrom][j+1][0]
                        for sampind in k:
                            pieces[sampind] = pieces.get(sampind,[]) + [(a,b)]
            except:
                print "Something bad just happened with"
                print pop[x][ind]
                raise
    for ind in pieces:
        starts = [ a for a,b in pieces[ind] ]
        starts.sort()
        ends = [ b for a,b in pieces[ind] ]
        ends.sort()
        if not ( starts[0]==0 and ends[len(ends)-1]==chrlen and all( [ starts[ii+1]==ends[ii] for ii in xrange(len(starts)-1) ] ) ):
            if print_details:
                print "Ooops!  where did", ind, "go?"
                print pieces[ind]
            errors.append( ind )
    if not errors:
        print "All good!", census(pop)
    else:
        print len(errors), "errors, oh dear!"
        return errors

def migrate(pop,migprobs):
    '''Do migration:
        migprobs[(x,y)] gives the probability that an individual now in x has a parent in y.
          with missing entries implicitly zero.
    '''
    newpop = {}.fromkeys( pop.keys() )
    for x in newpop:
        newpop[x] = {}
    jj = 0  # make sure we have a unique integer key for each ancestor
    for x in pop:
        cumprobs = [ migprobs.get((x,y),0) for y in pop.keys() ] 
        for ii in range(1,len(cumprobs)):
            cumprobs[ii] = cumprobs[ii]+cumprobs[ii-1]
        whichlocs = [ bisect.bisect_left( cumprobs, random.random() ) for k in xrange(len(pop[x])) ]
        newlocs = [ (pop.keys()+[x])[k] for k in whichlocs ]
        for y in newpop:
            newpop[y].update( [ [kk,pop[x][ind]] for kk,ind,z in zip(xrange(jj,jj+len(newlocs)),pop[x],newlocs) if z == y ] )
            jj += len(newlocs)
    return newpop

def crossover(chrom):
    '''Given a list of segments (i.e. a single chromosome)
    return an assignment of subsegments to the parental chromosomes.
      and chrom is a (set of) chromosome(s), i.e. a list of pairs (location,samples)
         ordered by locations, with samples the sample ids inheriting from the next segment to the right
    '''
    if not chrom:
        # nothing to do
        return None
    minloc = chrom[0][0]
    maxloc = chrom[-1][0]  # recall this is where state returns to [] finally
    # assign crossover locations -- begin with random subset of breaks between chromosomes
    # crossovers will be a list of crossover locations in increasing order (begins with zero, or smallest location!)
    crossovers = [ x for x in chrpos if (x > minloc) and (x < maxloc) and (random.random() < .5) ]
    z = minloc + random.expovariate(1)
    while z < maxloc :
        crossovers.append(z)
        z += random.expovariate(1)
    crossovers.sort()
    crossovers.append( [maxloc + 1] )  # this has to be something larger than all breaks in the chromosomes
    # all done generating crossovers
    copyto = random.sample([0,1],1)[0]
    matpat = [[],[]]
    laststate = set([])
    j = 0  # where we are in the chromosome
    for k in xrange(len(crossovers)):
        if len(laststate)>0:
            # note that this won't happen for k==0
            matpat[copyto].append( (crossovers[k-1],laststate) )
        while (j < len(chrom)) and (chrom[j][0] < crossovers[k]):
            matpat[copyto].append( chrom[j] )
            laststate = chrom[j][1]
            j += 1
        if j < len(chrom):
            if len(laststate)>0:
                matpat[copyto].append( (crossovers[k],set([])) )
            # switch to other chromosome
            copyto = ((copyto+1) % 2)
    return matpat

def combine(chroms,ibd):
    '''Given a list of diploid individuals, combine them into a single one.
    '''
    if len(chroms)==1:
        return chroms[0]
    matpat = [[],[]]
    for ii in xrange(2):
        whereat = [ 0 for chrom in chroms ]  # current locations
        nextones = [ ( chroms[j][ii][0] if len(chroms[j][ii])>0  else None ) for j in xrange(len(chroms)) ]
        thisstate = [ set([]) for u in nextones ]
        while any( nextones ):
            nextpos = min( [ u[0] for u in nextones if u ] )
            for j in xrange(len(whereat)):
                if nextones[j]:
                    if (nextones[j][0] == nextpos):
                        thisstate[j] = nextones[j][1]
                        whereat[j] += 1
                        if whereat[j] < len(chroms[j][ii]):
                            nextones[j] = chroms[j][ii][whereat[j]] 
                        else:
                            nextones[j] = None
            combstate = reduce( lambda u,v: u|v, thisstate )
            if len( filter( lambda u: len(u)>0, thisstate ) ) > 1:
                nextnextpos = min( [ u[0] for u in nextones if u ] )
                ibd.append( (nextpos,nextnextpos,combstate) )
            matpat[ii].append( (nextpos, combstate ) )
    return matpat


def rearrange(chroms,ibd):
    '''
    OBSOLETE:
    
    Take a list of overlapping segments and replace it by the list of nonoverlapping segments
    with labels equal to the concatenation of the labels in the original, e.g.
      (a,a+2,[1,2]) , (a+1,a+3,[1,3]) -> (a,a+1,[1,2]) , (a+1,a+2,[1,1,2,3]) , (a+2,a+3,[1,3])

    Also, record such newly overlapping (IBD) segments here.
    '''
    for ii in range(len(chroms)):
        if len(chroms[ii])<=1:
            continue
        breakpoints = [ (a,+1,k) for (a,b,k) in chroms[ii] ] + [ (b,-1,k) for (a,b,k) in chroms[ii] ]
        breakpoints.sort()
        newchrom = []
        start = breakpoints[0][0]
        # this is greater than 1 iff we're in a new segment.  it should never be negative.
        cumsum = breakpoints[0][1]  
        labels = breakpoints[0][2].copy() 
        for x,bit,k in breakpoints[1:len(breakpoints)]:
            # process all changes at the same location first
            if not x==start:
                if labels:
                    newchrom.append( (start,x,labels.copy()) )
                    if cumsum>1:
                        ibd.append( (start,x,labels.copy()) )
                start = x
            if bit == +1:
                labels.update( k )
            else:
                labels.difference_update( k )
            cumsum += bit
        # # sanity check
        # if labels:
        #     print "rearrange error! left with some?", chroms, labels, breakpoints
        chroms[ii] = newchrom

def parents(pop,t=0,ibdict=None,ancne=None,writeto=None):
    '''Reallocate each chromosome to the ancestors,
    and then recombine within each individual to resolve the mat/pat chromosomes.
    Optionally, write out to writeto rather than recording in ibdict.
    '''
    if ancne is None:
        ancne = effpopsize(pop,t)
    newpop = {}.fromkeys( pop.keys() )
    if (ibdict is None) or writeto:
        ibdict = {}
    ibdict.setdefault(t,{})
    for x in pop:
        newpop[x] = {}
        ibd = ibdict[t].setdefault(x,[])
        for ind in pop[x]:
            # only redistribute nonempty chromosomes
            nonempties = filter( lambda x: len(x)>0, pop[x][ind] )
            parents = [ int( ancne[x]*random.random() ) for k in xrange(len(nonempties)) ]
            # ensure parents are unique
            while len(set(parents)) < len(nonempties):
                parents = [ int( ancne[x]*random.random() ) for k in xrange(len(nonempties)) ]
            for ii in range(len(nonempties)):
                # each ancestral individual gets a list of chromosomes they passed on to an offspring
                newpop[x].setdefault(parents[ii],[]).append( crossover( nonempties[ii] ) )
        # now reassort those chromosomes into the parental chromosomes
        for ind in newpop[x]:
            # split out overlapping bits
            newpop[x][ind] = combine( newpop[x][ind], ibd )
    if writeto:
        writeibd(ibdict,outfile=writeto,closeafter=False,writeheader=(t==0))
    return newpop

def writeibd(ibdict,filename="coalpedigree.ibd.gz",writeinfo=True,outfile=None,closeafter=True,writeheader=True):
    '''Output the ibd info in ibdict, with or without timing and location information.
    Records all overlaps greater than minlen, either pairwise or as a list.
    Since all n choose 2 combinations can get HUGE, only records a sufficient list of these.
    Identification of multiple IBD is helped along by the "blockid", unique to each block.
    '''
    if outfile is None:
        outfile = fileopt(filename,"w")
    blockid = 0
    if writeheader:
        header = ["id1", "id2", "start", "end"]
        if writeinfo:
            header += ["blockid", "t", "location"]
        outfile.write(" ".join(header)+"\n") 
    for t in ibdict:
        for x in ibdict[t]:
            if ibdict[t][x]:
                for start,end,ids in ibdict[t][x]:
                    blockid += 1
                    tmpids = ids.copy()
                    ida = tmpids.pop()
                    while( tmpids ):
                        idb = tmpids.pop()
                        output = [ida,idb,start,end]
                        if writeinfo:
                            output += [ blockid, t, x ]
                        outfile.write(" ".join(map(str,output))+"\n")
    if closeafter:
        outfile.close()
    print blockid, "blocks written"


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
