#!/usr/bin/python
'''Parameters for migration between expanding populations
   Two populations:
          No contact since splittime.
          a and b are similar sized, constant.

   Run with e.g.
    ( BFIX=$(printf %05d $RANDOM); time nice -19 python ~/projects/coalpedigree/sim-ibd-pedigree.py -b ${BFIX}-demographics-expansion.fibd.gz -l ${BFIX}-demographics-expansion.log -t 300 -i ~/projects/coalpedigree/sim-demographics-expansion.py -e .001 )&

'''
import math

sampsizes = dict( a=100 )

splittime = 50
growthtime = 30
splitsize = 30000
nesize = 3000000


def ancnefn(t):
    if t>splittime:
        ancne = dict( a=splitsize )
    else:
        ancne = dict( a=int( nesize - (nesize-splitsize)*math.exp( (splittime-t) * math.log(.1)/growthtime ) ) )
    return ancne

def migprobs(t):
    migprobs = {}.fromkeys([('a','a')],1.0)
    return migprobs
