# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 09:57:47 2020

@author: sjsty
"""

import numpy as np
import numpy.random as rnd
import matplotlib as mpl
import matplotlib.pyplot as plt

LRn222 = 3.825 * 24 * 60 * 60  # the half-life for Rn-222 (in seconds)
LPo218 = 3.05 * 60  # the half-life for Po-218
DC1HL = np.array([LRn222, LPo218])
DC1Lambda = np.log(2) / DC1HL  # the decay constants for this first decay chain in units of 1/min
DC1Lambda

LRn220 = 54.5  # the half-life for Rn-220 (in seconds)
LPo216 = 0.158  # the half-life for Po-216
DC2HL = np.array([LRn220, LPo216])
DC2Lambda = np.log(2) / DC2HL  # the decay constants for this first decay chain in units of 1/min
DC2Lambda


def gen_inputs(*rates):
    proportions = [[0] * 3 for _ in range(12)]
    for i in range(12):
        for j in range(len(rates)):
            for r in range(j + 1):
                tmp = 1
                for q in range(j + 1):
                    if q != r:
                        tmp *= rates[q] / (rates[q] - rates[r])
                tmp *= np.exp(-rates[r] * 15 * i) - np.exp(-rates[r] * 15 * (i + 1))
                proportions[i][j] += tmp
    return [sum(p) for p in proportions]


def expcount(n, *args):
    countlist = np.array([0] * 12)

    for _ in range(n):
        templist = [int(x) // 15 for x in np.cumsum([rnd.exponential(a) for a in args])]

        for j in templist:
            if j < 12:
                countlist[j] = countlist[j] + 1
    return countlist


n = 10000
k = 2
counts = np.array(expcount(n, *(1 / DC2Lambda)))
for _ in range(k - 1):
    counts += np.array(expcount(n, *(1 / DC2Lambda)))  # + expcount(n,*(1/DC1Lambda))

counts = counts / k
print([t for t in zip(counts, gen_inputs(*DC2Lambda))])
