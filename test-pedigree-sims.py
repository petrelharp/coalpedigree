import sys
import time
import math, random
import pdb, copy
# import gc
# import objgraph, inspect

sys.path.append('/home/peter/projects/coalpedigree')
import coalpedigree as coal

reload(coal)

## Do one population, constant size

random.seed(1234)

sampsizes = dict( a=2 )
migprobs = dict( [ (('a','a'),0.0) ] )
ancne = dict(a=10)
coal.chrlens = ( 2.0, 1.0 )
coal.chrpos = tuple( [ sum( coal.chrlens[0:k] ) for k in range(1,len(coal.chrlens)) ] )   # cumulative sum: position if lined up end-to-end
coal.chrlen = sum(coal.chrlens)  # the last one (total length)

start = time.time()

pop = coal.initpop(sampsizes)
ibdict = {}
poplist = []
for t in xrange(20):
    coal.parents(pop,t=t,ibdict=ibdict,migprobs=migprobs,ancne=ancne)
    coal.sanity(pop,print_details=True)
    coal.writeibd(pop,minlen=0.0,gaplen=0.0,filename="test-fibd-"+("%(t)02d" % {'t':t})+".gz")
    poplist.append( copy.deepcopy(pop) )

print time.time()-start

# run this to look at the above.
rscript <- '''
source("~/projects/genome/ibd-blocks-fns.R")
fnames <- list.files(".","test-fibd-[0-9]*.gz")
.chrlens <- c(2,1)
.chrstarts <- c(0,2,3)

ibd <- lapply( fnames, read.table, header=TRUE )
for (k in seq_along(ibd)) {
    ibd[[k]]$chrom <- findInterval( ibd[[k]]$start, .chrstarts, rightmost.closed=TRUE )
    tmp <- length(.chrstarts) - findInterval( -ibd[[k]]$end, rev(-.chrstarts), rightmost.closed=TRUE )
    if ( any( tmp!= ibd[[k]]$chrom ) ) { stop("oops! blocks spanning chromosome gap in generation ",k) }
    ibd[[k]]$mapstart <- ibd[[k]]$start - .chrstarts[ibd[[k]]$chrom]
    ibd[[k]]$mapend <- ibd[[k]]$end - .chrstarts[ibd[[k]]$chrom]
}

pdf(file="blocks.pdf",width=7,height=4)
invisible( lapply( ibd, plotindivs, allids=0:3, chrspace=.1 ) )
dev.off()
'''


# coal.writecoal(ibdict,filename="test.coal.gz")
# objgraph.show_chain(objgraph.find_backref_chain(random.choice(objgraph.by_type('set')),inspect.ismodule),filename='chain.png')

## Do migration and expanding populations
## Three populations:
##        c splits from a 50 generations ago and grows over 20 generations
##        a and b are similar sized, constant

sampsizes = dict( zip(['a','b','c'],[1500,1000,500]) )

splittime = 50
growthtime = 30
splitsize = 10000
nesize = 300000

def ancnefn(pop,t):
    ancne = dict( a=nesize, b=nesize )
    if t>splittime:
        ancne['c'] = splitsize
    else:
        ancne['c'] = int( nesize - (nesize-splitsize)*math.exp( (splittime-t) * math.log(.1)/growthtime ) )
    return ancne

def migprobs(pop,t):
    migprobs = {}.fromkeys( [('a','b'), ('a','c'), ('b','a'), ('c','a')], .005  )
    if t>splittime:
        migprobs[('b','a')] = 1.0
        del migprobs[('a','b')]
    return migprobs

start = time.time()

# objgraph.show_growth(limit=3) 
# objgraph.show_most_common_types()

pop = coal.initpop(sampsizes)
for t in xrange(150):
    coal.parents(pop,t=t,migprobs=migprobs(pop,t),ancne=ancnefn(pop,t))
    print t

print time.time()-start

# objgraph.show_most_common_types()
 
nzeros = [ sum( map( lambda x: sum( map( lambda y: sum( [ ( 1 if len(z[1])==0 else 0 ) for z in y ] ), x ) ), d.values() ) ) for d in pop.values() ]  # number of zero-sets
nsets = [ sum( map( lambda x: sum( map( lambda y: sum( [ len(z[1]) for z in y ] ), x ) ), d.values() ) ) for d in pop.values() ]  # number of set elements
nbreaks = [ sum( map( lambda x: len(x[0])+len(x[1]), d.values() ) ) for d in pop.values() ]  # number of breakpoints
ninds = map( len, pop.values() )
[ (1.0*c/d,1.0*a/b,1.0*c/(b-a)) for a,b,c,d in zip( nzeros, nbreaks, nsets, ninds ) ]  # mean number of sets, mean proportion of sets that are empty and mean size of nonzero sets
# this is how many numbers?
sum(nbreaks) , sum(nsets)
newpop = {}.fromkeys(range(3))
ninds = [40896, 15736, 29079]
for x in range(3):
    newpop[x] = {}.fromkeys(range(ninds[x]))
    for y in newpop[x]:
        newpop[x][y] = [ [ ( 0.0, set(range(1)) ) for j in range(20) ], [ ( 0.0, set() ) for j in range(20) ] ]
    

del pop
gc.collect()
objgraph.show_most_common_types()
objgraph.show_growth()
showthis = random.choice(objgraph.by_type('set')); showthis; objgraph.show_chain(objgraph.find_backref_chain(showthis,inspect.ismodule),filename='chain.png')



####### testing
import coalpedigree as coal
import time

sampsizes = dict( zip([0,1,2],[500,1000,1500]) )
def ancne (t):
    if t<25:
        return dict( zip([0,1,2],[10000,10000,10000]) )
    else:
        return dict( zip([0,1,2],[10000,10000,100]) )

ancne = dict( zip([0,1,2],[10000,10000,10000]) )
pop = coal.initpop(sampsizes)
ibdict = {}
migprobs = {}.fromkeys( [ (x,y) for x in pop.keys() for y in pop.keys() ], .01 )
# outfile = file("/tmp/testing.gz","w")

starttime = time.time()

for t in xrange(30):
    pop = coal.parents(pop,t=t,ibdict=ibdict,ancne=ancne)
    # pop = coal.parents(pop,t=t,writeto=outfile)
    #coal.sanity(pop)
    pop = coal.migrate(pop,migprobs)
    #coal.sanity(pop)

time.time() - starttime

outfile.close()

######## time it better

import timeit

setup = '''
import coalpedigree
sampsizes = dict( zip([0,1,2],[500,1000,1500]) )
pop = coalpedigree.initpop(sampsizes)
ibdict = {}
'''

stmt = '''
for n in xrange(300):
    pop = coalpedigree.parents(pop,t=n,ibdict=ibdict)
    pop = coalpedigree.migrate(pop)
'''

t = timeit.Timer(stmt=stmt,setup=setup)
t.timeit(1)

t = timeit.Timer(stmt="x = random.sample(xrange(10000),2)",setup="import random")
t.timeit()
# 4.282135963439941
t = timeit.Timer(stmt='''
                 x = [ int(10000*random.random()) for k in xrange(2) ]
                 while len(set(x))<2:
                     x = [ int(10000*random.random()) for k in xrange(2) ]
''',setup="import random")
t.timeit()
# 1.1984889507293701


####### profile

import cProfile, pstats, trace

stmt = '''
sampsizes = dict( a=30 )
pop = coal.initpop(sampsizes)
outfile = file("test.fibd.gz","w")
for t in xrange(40):
    pop = coal.parents(pop,t=t,writeto=outfile,ancne=dict(a=10000))
'''

stmt = '''
sampsizes = dict( a=30 )
pop = coal.initpop(sampsizes)
ibdict = {}
outfile = file("test.fibd.gz","w")
for t in xrange(40):
    pop = coal.parents(pop,t=t,ibdict=ibdict,ancne=dict(a=10000))
'''

cProfile.run(stmt,"/tmp/simprof")
p = pstats.Stats("/tmp/simprof")
p.strip_dirs()
p.print_stats()

tracer = trace.Trace(trace=1,count=1)
tracer.run(stmt)
r = tracer.results()
r.write_results(show_missing=True, coverdir="/tmp")
