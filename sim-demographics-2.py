#!/usr/bin/python
'''Parameters for migration between expanding populations
   Three populations:
          c splits from a 50 generations ago and grows over 20 generations
          a and b are similar sized, constant

   Run with e.g.
    BFIX=$(printf %05d $RANDOM); time nice -19 python ~/projects/coalpedigree/sim-ibd-pedigree.py -b ${BFIX}-growing-migration-2.fibd.gz -l ${BFIX}-growing-migration-2.log -t 150 -i ~/projects/coalpedigree/sim-demographics-2.py &

'''
import math

sampsizes = dict( zip(['a','b','c'],[150,100,50]) )

# smaller Ne and smaller chromosomes
splittime = 50
growthtime = 30
splitsize = 1000
nesize = 30000

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

