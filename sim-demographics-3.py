#!/usr/bin/python
'''Parameters for migration between expanding populations
   Two populations:
          No contact since splittime.
          a and b are similar sized, constant.

   Run with e.g.
    ( BFIX=$(printf %05d $RANDOM); time nice -19 python ~/projects/coalpedigree/sim-ibd-pedigree.py -b ${BFIX}-growing-migration-3.fibd.gz -l ${BFIX}-growing-migration-3.log -t 300 -i ~/projects/coalpedigree/sim-demographics-3.py )&

'''
import math

sampsizes = dict( zip(['a','b'],[500,500]) )

# smaller Ne and smaller chromosomes
splittime = 60
nesize = 3000000

def ancnefn(pop,t):
    ancne = dict( a=nesize, b=nesize )
    return ancne

def migprobs(pop,t):
    if t>splittime:
        migprobs = {}.fromkeys( [('a','b'), ('b','a')], .005  )
    else:
        migprobs = {}.fromkeys( [('a','b'), ('b','a')], 0.0  )
    return migprobs

