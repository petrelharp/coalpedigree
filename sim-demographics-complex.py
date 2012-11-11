#!/usr/bin/python
'''Parameters for migration between expanding populations
   Two populations:
          No contact since splittime.
          a and b are similar sized, constant.

   Run with e.g.
    ( BFIX=$(printf %05d $RANDOM); time nice -19 python ~/projects/coalpedigree/sim-ibd-pedigree.py -b ${BFIX}-demographics-complex.fibd.gz -l ${BFIX}-demographics-complex.log -t 300 -i ~/projects/coalpedigree/sim-demographics-complex.py -e .001 )&

'''
import math

sampsizes = dict( a=100 )

recenttime = 5
flatsize = 30000
flattime = 15
decaytime = 30
decayrate = math.log(9)/11
bigsize = 800000
growthtime = 60
growthrate = math.log(29)/30
oldsize = 10000


def ancnefn(t):
    if t<recenttime:
        ancne = flatsize*math.exp(5-t)
    elif t<flattime:
        ancne = flatsize
    elif t<decaytime:
        ancne = (bigsize-flatsize) * math.exp( -(30-t)*decayrate )
    elif t<growthtime:
        ancne = bigsize - (bigsize-oldsize) * math.exp( (t-60)*growthrate )
    else:
        ancne = oldsize
    return dict( a=int(ancne) )

def migprobs(t):
    migprobs = {}.fromkeys([('a','a')],1.0)
    return migprobs
