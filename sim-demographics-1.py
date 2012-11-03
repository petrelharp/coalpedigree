#!/usr/bin/python
'''One population, constant size 10,000.

   Run with e.g.
    BFIX=$(printf %05d $RANDOM); time nice -19 python ~/projects/coalpedigree/sim-ibd-pedigree.py -b ${BFIX}-fixed-size-1.fibd.gz -l ${BFIX}-fixed-size-1.log -t 300 -i ~/projects/coalpedigree/sim-demographics-1.py -e .001 &

'''

sampsizes = dict( a=100 )

nesize = 100000

def ancnefn(t):
    ancne = dict( a=nesize )
    return ancne

def migprobs(t):
    migprobs = {}.fromkeys([('a','a')],1.0)
    return migprobs
