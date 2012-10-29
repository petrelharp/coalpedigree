#!/usr/bin/python
'''Simple simulation for testing purposes.
   Two populations:
          No contact since splittime.
          a and b are similar sized, constant.

   Run with e.g.
    ( BFIX=$(printf %05d $RANDOM); time nice -19 python ~/projects/coalpedigree/sim-ibd-pedigree.py -b ${BFIX}-migration-0.fibd.gz -l ${BFIX}-migration-0.log -t 200 -i ~/projects/coalpedigree/sim-demographics-0.py -e .0001 )&

'''
import math

sampsizes = dict( zip(['a','b'],[20,20]) )

splittime = 10
nesize = 1000

def ancnefn(t):
    ancne = dict( a=nesize, b=nesize )
    return ancne

def migprobs(t):
    if t>splittime:
        migprobs = {}.fromkeys( [('a','b'), ('b','a')], .05  )
    else:
        migprobs = {}.fromkeys( [('a','b'), ('b','a')], 0.0  )
    return migprobs

