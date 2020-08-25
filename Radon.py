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

LRn222 = 3.8235 * 24 * 60 * 60  # the half-life for Rn-222 (in seconds)
LPo218 = 3.098 * 60  # the half-life for Po-218
LPb214 = 26.8 * 60  # the half-life for Pb-214
LBi214 = 19.9 * 60  # the half-life for Bi-214
LPo214 = 164.3e-6  # the half-life for Po-214
DC1HL = np.array([LRn222, LPo218, LPb214, LBi214, LPo214])
DC1AD = np.array([1, 1, 0, 0, 1])
DC1Lambda = np.log(2) / DC1HL  # the decay constants for this first decay chain in units of Hz

LRn220 = 55.6  # the half-life for Rn-220 (in seconds)
LPo216 = 0.145  # the half-life for Po-216
LPb212 = 10.64 * 60 * 60  # the half-life for Pb-212
LBi212 = 60.55 * 60  # the half-life for Bi-212
# Technically, Bi-212 can both alpha and beta decay, but the beta decay mode then alpha decays almost immediately
DC2HL = np.array([LRn220, LPo216, LPb212, LBi212])
DC2AD = np.array([1, 1, 0, 1])
DC2Lambda = np.log(2) / DC2HL  # the decay constants for this first decay chain in units of Hz

max_threads = 12 # len(os.sched_getaffinity(0))
thread_pool = threading.BoundedSemaphore(max_threads)


# This function generates the expected number of counts in each sample period, as a proportion of the initial
# amount of the radioactive substance, starting from each point in the decay chain
# Parameters
## sample_time: Duration of each sample period
## n_samples: Overall number of sample periods
## rates: The decay rates (λ) for each isotope in the decay chain
## α_decay: The number of alpha particles emitted when this isotope decays
def gen_inputs(sample_time, n_samples, rates, α_decay=None, offset=0):
    if α_decay is None:
        # If alpha decay information is omitted, assume all decays are alpha decays
        α_decay = [1] * len(rates)
    exp_rates = [[0] * len(rates) for _ in range(n_samples + offset)]
    for t_i in range(n_samples + offset):
        # t_i -> index of the sample period
        for i in range(len(rates)):
            # i -> starting isotope
            # The following all comes from Edward's DEQ solution
            for j in range(i, len(rates)):
                for r in range(i, j + 1):
                    tmp = 1
                    for q in range(i, j + 1):
                        if q != r:
                            tmp *= rates[q] / (rates[q] - rates[r])
                    tmp *= np.exp(-rates[r] * sample_time * t_i) - np.exp(-rates[r] * sample_time * (t_i + 1))
                    exp_rates[t_i][i] += tmp * α_decay[j]
    return exp_rates[offset:]


# This function generates sample data for a decay chain using the exponential distribution of each isotope in the chain
# Parameters
## n: number of particles of the first isotope in the chain. All other isotopes are assumed to not be present at time
##    t = 0
## sample_time: Duration of each sample period
## n_samples: Overall number of sample periods
## rates: The decay rates (λ) for each isotope in the decay chain
## α_decay: The number of alpha particles emitted when this isotope decays
def expcount(n, sample_time, n_samples, rates, α_decay=None):
    if α_decay is None:
        # If alpha decay information is omitted, assume all decays are alpha decays
        α_decay = [1] * len(rates)
    countlist = np.array([0] * n_samples)  # Stores the counted decays in each sample period
    gen = rnd.default_rng()  # newest method for using numpy.random
    indices = [i for i, c in enumerate(α_decay) if c == 1]  # Which columns correspond to alpha decays?
    decays = np.cumsum(gen.exponential(1 / rates, (n, len(rates))),
                       1)  # Determine what time each decay in the chain occurs
    # gen.exponential(1/rates) produces the decay intervals for a single particle. The second argument to
    # gen.exponential is a size tuple, which specifies that from each exponential distribution, n samples should be
    # drawn. This corresponds to having n particles of the starting isotope, then determining for each particle
    # the decay intervals along the decay chain. Doing a cumulative sum of these intervals gives the time at which each
    # decay occurs.
    decays = decays[:, indices] // sample_time
    # Discard the non-alpha decays, and determine which sample period each decay would occur in
    for i in decays.flatten():  # All decays remaining are alpha, and have no particular order
        if int(i) < len(countlist):
            # If the particle decays in one of our samples, then we count it. If it decays after we finish sampling, we
            # don't count it.
            countlist[int(i)] += 1
    return countlist


# This function generates sample data for a decay chain using the exponential distribution of each isotope in the chain
# This differs from the function above because it keeps track of which particles have decayed at each time step
# Parameters
## n: number of particles of the first isotope in the chain. All other isotopes are assumed to not be present at time
##    t = 0
## sample_time: Duration of each sample period
## n_samples: Overall number of sample periods
## rates: The decay rates (λ) for each isotope in the decay chain
## α_decay: The number of alpha particles emitted when this isotope decays
def exp_count(n, sample_time, n_samples, rates, α_decay=None):
    if α_decay is None:
        # If alpha decay information is omitted, assume all decays are alpha decays
        α_decay = [1] * len(rates)
    countlist = np.array([0] * n_samples)  # Stores the counted decays in each sample period
    state = [0] * len(rates)  # Stores the current number of each type of particle
    state[0] = n  # Starting with n particles of first isotope, and no progeny
    for i in range(n_samples):
        # For each sample period, determine the state at the end of the sample period and the counts seen during this
        # sample period.
        state, countlist[i] = exp_state(state, sample_time, rates, α_decay)
    return countlist


# This function generates a sample data for a single interval for a decay chain using the exponential
# distribution of each isotope in the chain.
# Parameters
## init_state: The initial state of the system (how many particles of each isotope
## interval: The length of the interval of interest
## rates: The decay rates (λ) for each isotope in the decay chain
## α_decay: The number of alpha particles emitted when this isotope decays
## Note init_state, rates, and α_decay must all be the same length and dimension
def exp_state(init_state, interval, rates, α_decay=None):
    if α_decay is None:
        # If alpha decay information is omitted, assume all decays are alpha decays
        α_decay = np.array([1] * len(init_state))
    if len(init_state) != len(rates) or len(init_state) != len(α_decay):
        # Raise an exception if the lengths are inconsistent
        raise Exception("Length of state, rates, and counts must match")
    count = np.array([0] * len(init_state))  # Track how many particles of each type have decayed
    gen = rnd.default_rng()  # Set up np.random
    for type in range(len(init_state)):
        # For each isotope in the chain, we regard it as the first isotope in its own chain, where there are no
        # predecessors, and there is no progeny. We then determine the decay intervals along the chain for each particle
        # we started with, and sum these intervals to determine the times of decay
        # The shortened chain may be regarded as starting with `init_state[type]` number of particles, and having the
        # decay rates `rates[type:]`, which includes this particle and everything after it
        decays = gen.exponential(1 / (rates[type:]), (init_state[type], len(rates) - type))
        count[type:] += np.sum(np.cumsum(decays, 1) < interval, 0)
        # We then add to our count of how many of each particle decayed, offset so that our determined decays match
        # with the count variable. The expression `np.cumsum(decays, 1) < interval` returns a boolean matrix, where each
        # row is a particle and each column is a possible decay. If the entry is True, that particle has undergone that
        # decay, because the time of decay is less than our interval. If it is false, that particle decays after our
        # interval has ended. Summing this boolean array over the rows gives the number of particles of each type that
        # decayed in the interval.
    state = np.array(init_state) - count  # Update the state with decay information
    state[1:] += count[0:-1]  # Each decay produces one of the next element in the chain
    return state, np.dot(count, α_decay)  # Return the state and the count of alpha decays


def runtrial_thread(args):
    with thread_pool:
        logging.info("Thread %s: starting", args[0])
        runtrial(*(args[1:]))


def runtrial(st, tt, i, j):
    ns = tt // st
    offset_n = int(np.ceil(90/st))
    in_rn222 = np.array(gen_inputs(st, ns, DC1Lambda, DC1AD, offset_n))
    in_rn222 = in_rn222[:,0]
    in_rn220 = np.array(gen_inputs(st, ns, DC2Lambda, DC2AD, offset_n))
    in_rn220 = in_rn220[:,0]
    rn222_est = [0.]*10
    rn220_est = [0.]*10
    for k in range(10):
        out = np.array(exp_count(1000000, st, ns + 30, DC1Lambda)) + np.array(
            exp_count(1000 // 6, st, ns + 30, DC2Lambda))
        out = out[30:]
        lr = LinearRegression(fit_intercept=False).fit(np.transpose(np.vstack((in_rn222,in_rn220))), out)
        rn222_est[k], rn220_est[k] = lr.coef_
    print("st: {}s, tt: {}s, Rn222 => mean: {:1.1f}, std: {:1.1f}".format(st, tt, np.mean(rn222_est),np.std(rn222_est)))
    print("st: {}s, tt: {}s, Rn220 => mean: {:1.1f}, std: {:1.1f}".format(st, tt, np.mean(rn220_est), np.std(rn220_est)))
    ratio_est = np.array(rn222_est)/np.array(rn220_est)
    print("st: {}s, tt: {}s, ratio => mean: {:1.1f}, std: {:1.1f}".format(st, tt, np.mean(ratio_est), np.std(ratio_est)))
    rn222_mean[i][j] = np.mean(rn222_est)
    rn222_stdv[i][j] = np.std(rn222_est)
    rn220_mean[i][j] = np.mean(rn220_est)
    rn220_stdv[i][j] = np.std(rn220_est)


if __name__ == "__main__":
    format = "%(asctime)s: %(message)s"
    logging.basicConfig(format=format, level=logging.INFO, datefmt="%H:%M:%S")

    n_period_grid = 5
    n_period_start = 3
    n_period_step = 3
    n_time_grid = 25
    n_time_start = 1 * 60
    n_time_step = 10

    rn222_mean = [[0] * n_time_grid for _ in range(n_period_grid)]
    rn222_stdv = [[0] * n_time_grid for _ in range(n_period_grid)]
    rn220_mean = [[0] * n_time_grid for _ in range(n_period_grid)]
    rn220_stdv = [[0] * n_time_grid for _ in range(n_period_grid)]

    jobs = [(i, (i % n_period_grid) * n_period_step + n_period_step, n_time_step * (i // n_period_grid) + n_time_start,
             i % n_period_grid, i // n_period_grid) for i in range(n_time_grid * n_period_grid)]
    threads = list()
    for job in jobs:
        logging.info("Main\t: create and start thread %d", job[0])
        x = threading.Thread(target=runtrial_thread, args=(job,))
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

    sdprod = np.multiply(rn222_stdv, rn220_stdv)
    i, j = np.unravel_index(np.argmin(sdprod), sdprod.shape)
    print("Minimum std product found with period {}s and total sample duration {}s", i * n_period_step + n_period_step, j * n_time_step + n_time_start)
    print("Radon-220 estimate at minimum: {:g} ± {:g}", rn220_mean[i][j], rn220_stdv[i][j])
    print("Radon-222 estimate at minimum: {:g} ± {:g}", rn222_mean[i][j], rn222_stdv[i][j])