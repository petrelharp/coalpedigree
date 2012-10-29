#!/usr/bin/python
'''Parameters for migration between expanding populations
   Two populations:
          No contact since splittime.
          a and b are similar sized, constant.

   Run with e.g.
    ( BFIX=$(printf %05d $RANDOM); time nice -19 python ~/projects/coalpedigree/sim-ibd-pedigree.py -b ${BFIX}-migration-3.fibd.gz -l ${BFIX}-migration-3.log -t 200 -i ~/projects/coalpedigree/sim-demographics-3.py -e .001 )&

'''
import math

sampsizes = dict( zip(['a','b'],[500,500]) )

splittime = 60
nesize = 3000000

def ancnefn(t):
    ancne = dict( a=nesize, b=nesize )
    return ancne

def migprobs(t):
    if t>splittime:
        migprobs = {}.fromkeys( [('a','b'), ('b','a')], .005  )
    else:
        migprobs = {}.fromkeys( [('a','b'), ('b','a')], 0.0  )
    return migprobs

