#!/usr/bin/python
'''Parameters for migration between expanding populations
   Three populations:
          c splits from a 50 generations ago and grows over 20 generations
          a and b are similar sized, constant

   Run with e.g.
    time nice -19 python ~/projects/genome/sim-ibd-pedigree.py -o growing-migration-2.fibd.gz -l growing-migration-2.log -t 150 -i ~/projects/genome/sim-demographics-2.py 

'''
import math

sampsizes = dict( zip(['a','b','c'],[1500,1000,500]) )

# smaller Ne and smaller chromosomes
splittime = 50
growthtime = 30
splitsize = 1000
nesize = 30000

coal.chrlens = ( 2.0, 1.0 )
coal.chrpos = tuple( [ sum( coal.chrlens[0:k] ) for k in range(1,len(coal.chrlens)) ] )   # cumulative sum: position if lined up end-to-end
coal.chrlen = sum(coal.chrlens)  # the last one (total length)

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

