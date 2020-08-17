# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 09:57:47 2020

@author: sjsty
"""

import numpy as np
import numpy.random as rnd
import matplotlib as mpl
import matplotlib.pyplot as plt

LRn222 = 3.825*24*60*60 # the half-life for Rn-222 (in seconds)
LPo218 = 3.05*60        # the half-life for Po-218
DC1HL = np.array([LRn222,LPo218])
DC1Lambda = np.log(2)/DC1HL # the decay constants for this first decay chain in units of 1/min
DC1Lambda

LRn220 = 54.5 # the half-life for Rn-220 (in seconds)
LPo216 = 0.158        # the half-life for Po-216
DC2HL = np.array([LRn220,LPo216])
DC2Lambda = np.log(2)/DC2HL # the decay constants for this first decay chain in units of 1/min
DC2Lambda

def gen_inputs(n, *args):
    inputs = [[0]*len(args)]*12
    for i in range(12):
        inputs[i][0] = np.exp(-args[0]*15*(i+1))
        for j in range(1,len(args)):
            inputs[i][j] = 0
            for r in range(j+1):
                tmp = 1
                for q in range(j+1):
                    if q == r:
                        continue
                    tmp *= args[q]/(args[q]-args[r])
                tmp *= args[r]
                tmp *= np.exp(-args[r]*15*(i+1))
                inputs[i][j] += tmp
            inputs[i][j] /= args[j]
    return inputs

def expcount(n, *args):
    countlist = np.array([0]*12)
    
    for _ in range(n):
        templist=[int(x)//12 for x in np.cumsum([rnd.exponential(a) for a in args])]
        
        for j in templist:
            if j <12:
                countlist[j] = countlist[j]+1
    return countlist


n=10000
    
counts = expcount(n,1/DC1Lambda) + expcount(n,1/DC2Lambda)
gen_inputs(n,1,2)

def get_inputs(*args):
    inputs = [[0]]*12
    times = [0,15,30,45,60,75,90,105,120,135,150,165,180]
    for i in range(12):
        inputs[i] = inputs[i]+(np.exp(-args[0]*times[i+1])-np.exp(-args[0]*times[i]))
        progeny_sum = 0
        for r in range(len(args)):
            tmp=1
            for q in range(len(args)):
                if q==r:
                    continue
                tmp *= args[q]/(args[q]-args[r])
            progeny_sum = progeny_sum + tmp*(np.exp(-args[r]*times[i+1])-np.exp(-args[r]*times[i]))
        inputs[i] = inputs[i]+progeny_sum
    return inputs
        
get_inputs(4,5)
