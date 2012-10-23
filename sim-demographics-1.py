#!/usr/bin/python
'''One population, constant size 10,000.

   Run with e.g.
    BFIX=$(printf %05d $RANDOM); time nice -19 python ~/projects/coalpedigree/sim-ibd-pedigree.py -b ${BFIX}-growing-migration-1.fibd.gz -l ${BFIX}-growing-migration-1.log -t 250 -i ~/projects/coalpedigree/sim-demographics-1.py -e .001 &

'''

sampsizes = dict( a=1000 )

nesize = 10000

def ancnefn(pop,t):
    ancne = dict( a=nesize )
    return ancne

def migprobs(pop,t):
    migprobs = {}.fromkeys([('a','a')],1.0)
    return migprobs
