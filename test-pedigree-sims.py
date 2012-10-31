import sys
import time
import math, random
import pdb, copy
# import gc
# import objgraph, inspect

sys.path.append('/home/peter/projects/coalpedigree')
import coalpedigree as coal

reload(coal)

random.seed(1234)

sampsizes = dict( a=2, b=2 )
migprobs = dict( [ (('a','b'),0.01),(('b','a'),0.01) ] )
# ancne = dict(a=10,b=10)
def ancnefn(t):
    return dict(a=max(10,100-5*t),b=10)

coal.chrlens = ( 2.0, 1.0 )
coal.chrpos = tuple( [ sum( coal.chrlens[0:k] ) for k in range(1,len(coal.chrlens)) ] )   # cumulative sum: position if lined up end-to-end
coal.chrlen = sum(coal.chrlens)  # the last one (total length)

start = time.time()

pop = coal.initpop(sampsizes)
ibdict = {}
poplist = []
for t in xrange(20):
    print " ----", t
    coal.parents(pop,t=t,ibdict=ibdict,migprobs=migprobs,ancne=ancnefn(t))
    coal.sanity(pop,sampsizes=sampsizes,ancne=ancnefn(t),print_details=True)
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
    if ( any( tmp!= ibd[[k]]$chrom ) ) { warning("oops! blocks spanning chromosome gap in generation ",k) }
    ibd[[k]]$mapstart <- ibd[[k]]$start - .chrstarts[ibd[[k]]$chrom]
    ibd[[k]]$mapend <- ibd[[k]]$end - .chrstarts[ibd[[k]]$chrom]
}

pdf(file="blocks.pdf",width=7,height=4)
invisible( lapply( ibd, plotindivs, allids=0:7, chrspace=.1 ) )
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
    coalpedigree.parents(pop,t=n,ibdict=ibdict)
    coalpedigree.migrate(pop)
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
sampsizes = dict( a=10, b=10 )
ancne = dict(a=10000,b=10000)
pop = coal.initpop(sampsizes)
outfile = file("test.fibd.gz","w")
migprobs = dict( [ (('a','b'),0.01),(('b','a'),0.01) ] )
for t in xrange(40):
    coal.parents(pop,t=t,migprobs=migprobs,ancne=ancne)
'''

cProfile.run(stmt,"/tmp/simprof")
p = pstats.Stats("/tmp/simprof")
p.strip_dirs()
p.print_stats()

tracer = trace.Trace(trace=0,count=1)
tracer.run(stmt)
r = tracer.results()
r.write_results(show_missing=True, coverdir="/tmp")
