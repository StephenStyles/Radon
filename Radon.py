# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 09:57:47 2020

@author: sjsty
"""

import numpy as np
import numpy.random as rnd
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

LRn222 = 3.825 * 24 * 60 * 60  # the half-life for Rn-222 (in seconds)
LPo218 = 3.05 * 60  # the half-life for Po-218
LPb210 = (26.8 * 60) + (19.7 * 60) + (1.5e-4)  # the half-lives for Pb-214, Bi-214, and Po-214
DC1HL = np.array([LRn222, LPo218])
DC1Lambda = np.log(2) / DC1HL  # the decay constants for this first decay chain in units of 1/min

LRn220 = 54.5  # the half-life for Rn-220 (in seconds)
LPo216 = 0.158  # the half-life for Po-216
DC2HL = np.array([LRn220, LPo216])
DC2Lambda = np.log(2) / DC2HL  # the decay constants for this first decay chain in units of 1/min

time_to_run = 900;
sample_time = 15;
n_samples = time_to_run // sample_time;


def gen_inputs(*rates):
    proportions = [[0] * 3 for _ in range(n_samples)]
    for i in range(n_samples):
        for j in range(len(rates)):
            for r in range(j + 1):
                tmp = 1
                for q in range(j + 1):
                    if q != r:
                        tmp *= rates[q] / (rates[q] - rates[r])
                tmp *= np.exp(-rates[r] * i) - np.exp(-rates[r] * (i + 1))
                proportions[i][j] += tmp
    inputs = [sum(p) for p in proportions]
    return inputs  # /inputs[0]


def expcount(n, *args):
    countlist = np.array([0] * n_samples)

    for _ in range(n):
        templist = [int(x) for x in np.cumsum([rnd.exponential(a) for a in args])]

        for j in templist:
            if j < len(countlist):
                countlist[j] = countlist[j] + 1
    return countlist


time_to_run = 60;
sample_time = 15;
n_samples = time_to_run // sample_time;

k = 5
n1_mean = [[0] * 15 for _ in range(60)]
n1_stdv = [[0] * 15 for _ in range(60)]
n2_mean = [[0] * 15 for _ in range(60)]
n2_stdv = [[0] * 15 for _ in range(60)]

for i in range(60):
    sample_time = i + 1
    for j in range(15):
        time_to_run = 60 * (j + 1)
        n_samples = time_to_run // sample_time;
        n1_est = []
        n2_est = []
        for _ in range(k):
            n1 = 100000
            n2 = 50000
            counts = np.array(expcount(n1, *(1 / DC1Lambda))) + np.array(expcount(n2, *(1 / DC2Lambda)))
            input1 = gen_inputs(*DC1Lambda)
            input2 = gen_inputs(*DC2Lambda)

            inputs = np.transpose(np.vstack((input1, input2)))
            estimated_amounts = LinearRegression(fit_intercept=False).fit(inputs, counts)
            n1_est.append(estimated_amounts.coef_[0])
            n2_est.append(estimated_amounts.coef_[1])
        n1_mean[i][j] = np.mean(n1_est)
        n1_stdv[i][j] = np.std(n1_est)
        n2_mean[i][j] = np.mean(n2_est)
        n2_mean[i][j] = np.std(n2_est)
        print("With sample time {}s and run time {}s:".format(sample_time, time_to_run))
        print("\tThe mean for Rn222 is {0} with standard deviation {1}".format(np.mean(n1_est), np.std(n1_est)))
        print("\tThe mean for Rn220 is {0} with standard deviation {1}".format(np.mean(n2_est), np.std(n2_est)))
