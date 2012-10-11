#!/usr/bin/python
'''One population, constant size 10,000.

   Run with e.g.
    BFIX=$(printf %05d $RANDOM); time nice -19 python ~/projects/coalpedigree/sim-ibd-pedigree.py -b ${BFIX}-growing-migration-1.fibd.gz -l ${BFIX}-growing-migration-1.log -t 150 -i ~/projects/coalpedigree/sim-demographics-1.py &

'''
import math

sampsizes = dict( a=200 )

# smaller Ne and smaller chromosomes
nesize = 10000

coal.chrlens = ( 2.0, 2.0, 1.0, 1.0 )
coal.chrpos = tuple( [ sum( coal.chrlens[0:k] ) for k in range(1,len(coal.chrlens)) ] )   # cumulative sum: position if lined up end-to-end
coal.chrlen = sum(coal.chrlens)  # the last one (total length)

def ancnefn(pop,t):
    ancne = dict( a=nesize )
    return ancne

