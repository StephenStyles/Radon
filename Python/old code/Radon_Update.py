# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 11:31:38 2020

@author: sjsty
"""
import numpy as np
import numpy.random as rnd
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import threading
import concurrent.futures
import logging
import time
import os
import csv


def gen_inputs(sample_time, n_samples,offset, rates, counts=None):
    if counts is None:
        counts = [1] * len(rates)
    exp_rates = [[0] * len(rates) for _ in range(n_samples+offset)]
    for t_i in range(n_samples+offset):
        for i in range(len(rates)):
            for j in range(i,len(rates)):
                for r in range(i,j+1):
                    tmp = 1
                    for q in range(i,j+1):
                        if q != r:
                            tmp *= rates[q] / (rates[q] - rates[r])
                    tmp *= np.exp(-rates[r] * sample_time * t_i) - np.exp(-rates[r] * sample_time * (t_i + 1))
                    exp_rates[t_i][i] += tmp*counts[j]
    return exp_rates[offset:]


def expcount(n, sample_time, n_samples, *args,counts):
    countlist = np.array([0] * (n_samples+30))
    #countlist = np.array([0] * P)
    gen=rnd.default_rng()

    decays = np.cumsum(gen.exponential(args,(n,len(args))),1)
    indices = [i for i, x in enumerate(counts) if x == 1]
    decays = decays[:,[indices]]
    decays = decays//sample_time

    for i in decays.flatten():
        if int(i) < len(countlist):
            countlist[int(i)] += 1
    return countlist[30:]
    #return countlist


LRn222 = 3.8235 * 24 * 60 * 60  # the half-life for Rn-222 (in seconds)
LPo218 = 3.098 * 60  # the half-life for Po-218
LPb214 = 26.8 * 60 # the half-life for Pb-214
LBi214 = 19.9 * 60 # the half-life for Bi-214
LPo214 = 164.3e-6 # the half-life for Po-214
DC1HL = np.array([LRn222, LPo218,LPb214,LBi214,LPo214])
DC1AD = np.array([     1,      1,     0,     0,     1])
DC1Lambda = np.log(2) / DC1HL  # the decay constants for this first decay chain in units of 1/min


LRn220 = 55.6  # the half-life for Rn-220 (in seconds)
LPo216 = 0.145  # the half-life for Po-216
LPb212 = 10.64 * 60 * 60 # the half-life for Pb-212
LBi212 = 60.55 * 60 # the half-life for Bi-212
# Technically, Bi-212 can both alpha and beta decay, but the beta decay mode then alpha decays almost immediately
DC2HL = np.array([LRn220, LPo216,LPb212,LBi212])
DC2AD = np.array([     1,      1,     0,     1])
DC2Lambda = np.log(2) / DC2HL  # the decay constants for this first decay chain in units of 1/min

files = ['Radon_5min_Trial1.xlsx', 'Radon_5min_Trial2.xlsx','Radon_5min_Trial3.xlsx','Thoron_5min_Trial1.xlsx']

for file in files:
    RadonDF = pd.read_excel(file)

    st = 3
    tt = len(RadonDF.index)*3
    offset=30
    ns = tt//st
    in_rn222=np.array(gen_inputs(st,ns,offset,DC1Lambda,counts=DC1AD))
    in_rn222=in_rn222[:,0]
    #in_po218=gen_inputs(st,ns,*(DC1Lambda[1:]))
    #in_po214=gen_inputs(st,ns,DC1Lambda[4])
    
    in_rn220=np.array(gen_inputs(st,ns,offset,DC2Lambda,counts=DC2AD))
    in_rn220=in_rn220[:,0]
    # in_po216=gen_inputs(st,ns,*(DC2Lambda[1:]))
    # in_bi212=gen_inputs(st,ns,DC2Lambda[3])
    
    # lr=LinearRegression(fit_intercept=False).fit(np.transpose(np.vstack((in_rn222,in_po218,in_po214, in_rn220,in_po216,in_bi212))), RadonDF['Counts'])
    # rn222_est,po218_est,po214_est,rn220_est, po216_est,bi212_est = lr.coef_
    # lr=LinearRegression(fit_intercept=False).fit(np.array(in_rn222).reshape(-1,1), RadonDF['Counts'])
    # rn222_est = lr.coef_
    
    lr=LinearRegression(fit_intercept=False).fit(np.transpose(np.vstack((in_rn222,in_rn220))), RadonDF['Counts'])
    rn222_est,rn220_est = lr.coef_
    
    act_rn222 = [rn222_est*DC1Lambda[0]/3.7e-2/0.3, rn222_est*DC1Lambda[0]/0.3e-3]
    act_rn220 = [rn220_est*DC2Lambda[0]/3.7e-2/0.3, rn220_est*DC2Lambda[0]/0.3e-3]
    print("st: {}s, tt: {}s, Rn222 => Estimate: {:1f}, Activity: {:1f} pCi/L, {:1f} Bq/m^3".format(st,tt,rn222_est,act_rn222[0],act_rn222[1]))
    print("st: {}s, tt: {}s, Rn220 => Estimate: {:1f}, Activity: {:1f} pCi/L, {:1f} Bq/m^3".format(st,tt,rn220_est,act_rn220[0],act_rn220[1]))
    
    
    gendata = pd.DataFrame({'Counts':expcount(int(np.ceil(rn222_est)),st,ns,*(1/DC1Lambda),counts=DC1AD)+expcount(int(np.ceil(rn220_est)),st,ns,*(1/DC2Lambda),counts=DC2AD)})
    expected = pd.DataFrame({'Expected': int(np.ceil(rn222_est))*in_rn222+int(np.ceil(rn220_est))*in_rn220})
    
    plt.plot( RadonDF.index, RadonDF["Counts"], marker='o', markerfacecolor='blue', markersize=4, color='blue', linewidth=1, label ="Radon" )
    plt.plot( gendata.index, gendata, marker='o', markerfacecolor='red', markersize=4, color='red', linewidth=1,label = "Prediction")
    plt.plot( expected.index, expected, marker='', color='black', linewidth=3,label = "Expected Value")
    plt.legend()
    plt.show()

# RadonDF = pd.read_excel("Radon_Continuous_Trial1.xlsx")

# st = 3
# tt = len(RadonDF.index)*3
# offset=30
# ns = tt//st
# in_rn222=np.array(gen_inputs(st,ns,offset,DC1Lambda,counts=DC1AD))
# in_rn222=in_rn222[:,0]
# #in_po218=gen_inputs(st,ns,*(DC1Lambda[1:]))
# #in_po214=gen_inputs(st,ns,DC1Lambda[4])

# in_rn220=np.array(gen_inputs(st,ns,offset,DC2Lambda,counts=DC2AD))
# in_rn220=in_rn220[:,0]
# # in_po216=gen_inputs(st,ns,*(DC2Lambda[1:]))
# # in_bi212=gen_inputs(st,ns,DC2Lambda[3])

# # lr=LinearRegression(fit_intercept=False).fit(np.transpose(np.vstack((in_rn222,in_po218,in_po214, in_rn220,in_po216,in_bi212))), RadonDF['Counts'])
# # rn222_est,po218_est,po214_est,rn220_est, po216_est,bi212_est = lr.coef_
# # lr=LinearRegression(fit_intercept=False).fit(np.array(in_rn222).reshape(-1,1), RadonDF['Counts'])
# # rn222_est = lr.coef_

# lr=LinearRegression(fit_intercept=False).fit(np.transpose(np.vstack((in_rn222,in_rn220))), RadonDF['Counts'])
# rn222_est,rn220_est = lr.coef_

# act_rn222 = [rn222_est*DC1Lambda[0]/3.7e-2/0.3, rn222_est*DC1Lambda[0]/0.3e-3]
# act_rn220 = [rn220_est*DC2Lambda[0]/3.7e-2/0.3, rn220_est*DC2Lambda[0]/0.3e-3]
# print("st: {}s, tt: {}s, Rn222 => Estimate: {:1f}, Activity: {:1f} pCi/L, {:1f} Bq/m^3".format(st,tt,rn222_est,act_rn222[0],act_rn222[1]))
# #print("st: {}s, tt: {}s, Rn220 => Estimate: {:1f}, Activity: {:1f} pCi/L, {:1f} Bq/m^3".format(st,tt,rn220_est,act_rn220[0],act_rn220[1]))


# gendata = pd.DataFrame({'Counts':expcount(int(np.ceil(rn222_est)),st,ns,*(1/DC1Lambda),counts=DC1AD)+expcount(int(np.ceil(rn220_est)),st,ns,*(1/DC2Lambda),counts=DC2AD)})
# expected = pd.DataFrame({'Expected': int(np.ceil(rn222_est))*in_rn222+int(np.ceil(rn220_est))*in_rn220})

# #plt.plot(RadonDF.index, RadonDF["Counts"], '-o')
# plt.plot( RadonDF.index, RadonDF["Counts"], marker='o', markerfacecolor='blue', markersize=4, color='blue', linewidth=1, label ="Radon" )
# plt.plot( gendata.index, gendata, marker='o', markerfacecolor='red', markersize=4, color='red', linewidth=1,label = "Prediction")
# plt.plot( expected.index, expected, marker='', color='black', linewidth=3,label = "Expected Value")
# plt.legend()
# plt.show()

