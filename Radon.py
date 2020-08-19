#!/usr/bin/python3
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
import threading
import concurrent.futures
import logging
import time
import os
import csv

LRn222 = 3.825 * 24 * 60 * 60  # the half-life for Rn-222 (in seconds)
LPo218 = 3.05 * 60  # the half-life for Po-218
LPb210 = (26.8 * 60) + (19.7 * 60) + (1.5e-4)  # the half-lives for Pb-214, Bi-214, and Po-214
DC1HL = np.array([LRn222, LPo218])
DC1Lambda = np.log(2) / DC1HL  # the decay constants for this first decay chain in units of 1/min

LRn220 = 54.5  # the half-life for Rn-220 (in seconds)
LPo216 = 0.158  # the half-life for Po-216
DC2HL = np.array([LRn220, LPo216])
DC2Lambda = np.log(2) / DC2HL  # the decay constants for this first decay chain in units of 1/min

max_threads = len(os.sched_getaffinity(0))
thread_pool = threading.BoundedSemaphore(max_threads)

def gen_inputs(sample_time, n_samples, *rates):
    proportions = [[0] * 3 for _ in range(n_samples)]
    for i in range(n_samples):
        for j in range(len(rates)):
            for r in range(j + 1):
                tmp = 1
                for q in range(j + 1):
                    if q != r:
                        tmp *= rates[q] / (rates[q] - rates[r])
                tmp *= np.exp(-rates[r] * i * sample_time) - np.exp(-rates[r] * (i + 1) * sample_time)
                proportions[i][j] += tmp
    return [sum(p) for p in proportions]


def expcount(n, sample_time, n_samples, *args):
    countlist = np.array([0] * n_samples)

    for _ in range(n):
        templist = [int(x)//sample_time for x in np.cumsum([rnd.exponential(a) for a in args])]

        for j in templist:
            if j < len(countlist):
                countlist[j] = countlist[j] + 1
    return countlist

def runtrial_thread(args):
    with thread_pool:
        logging.info("Thread %s: starting", args[0])
        runtrial(*(args[1:]))

def runtrial(st,tt,i,j):
    ns = tt//st
    in_rn222=gen_inputs(st,ns,*DC1Lambda)
    in_rn220=gen_inputs(st,ns,*DC2Lambda)
    rn222_est=[0]*10
    rn220_est=[0]*10
    for k in range(10):
        out=np.array(expcount(100000,st,ns,*(1/DC1Lambda)))+np.array(expcount(50000,st,ns,*(1/DC2Lambda)))
        lr=LinearRegression(fit_intercept=False).fit(np.transpose(np.vstack((in_rn222,in_rn220))),out)
        rn222_est[k],rn220_est[k] = lr.coef_
    print("st: {}s, tt: {}s, Rn222 => mean: {:1f}, std: {:1f}".format(st,tt,np.mean(rn222_est),np.std(rn222_est)))
    print("st: {}s, tt: {}s, Rn220 => mean: {:1f}, std: {:1f}".format(st, tt, np.mean(rn220_est), np.std(rn220_est)))
    rn222_mean[i][j] = np.mean(rn222_est)
    rn222_stdv[i][j] = np.std(rn222_est)
    rn220_mean[i][j] = np.mean(rn220_est)
    rn220_stdv[i][j] = np.std(rn220_est)

if __name__ == "__main__":
    format = "%(asctime)s: %(message)s"
    logging.basicConfig(format=format, level=logging.INFO, datefmt="%H:%M:%S")

    rn222_mean = [[0] * 15 for _ in range(60)]
    rn222_stdv = [[0] * 15 for _ in range(60)]
    rn220_mean = [[0] * 15 for _ in range(60)]
    rn220_stdv = [[0] * 15 for _ in range(60)]

    jobs = [(i,i%60+1,60*(i//60+1),i%60,i//60) for i in range(60*15)]
    threads = list()
    for job in jobs:
        logging.info("Main\t: create and start thread %d", job[0])
        # if(len(threads) >= 11):
        #     threads[0].join()
        #     threads.pop(0)
        x = threading.Thread(target=runtrial_thread,args=(job,))
        threads.append(x)
        x.start()
    for t in threads:
        t.join()
    with open("rn222_mean.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rn222_mean)
    with open("rn222_std.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rn222_stdv)
    with open("rn220_mean.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rn220_mean)
    with open("rn220_std.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rn220_stdv)
