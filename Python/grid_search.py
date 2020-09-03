import logging
import os
import threading
import csv

import numpy as np

from decay_chain_simulation import DecayChainSimulation
from sklearn.linear_model import LinearRegression
from radon_setup import DC1Lambda, DC1Mode, DC2Lambda, DC2Mode

max_threads = os.cpu_count()
thread_pool = threading.BoundedSemaphore(max_threads)

grid_period_num = 3
grid_period_first = 3
grid_period_step = 1
grid_duration_num = 3
grid_duration_first = 4 * 60
grid_duration_step = 30

runs_per_point = 25

Rn220_amount = 150
Rn222_amount = 750_000


def grid_point_thread(args):
    with thread_pool:
        logging.info("Thread %s: starting", args[0])
        grid_point(*(args[1:]))


def grid_point(ts: float, p: int, tf: float, i: int, j: int):
    global rn222_mean, rn222_stdv, rn220_mean, rn220_stdv
    Rn222_DCS = DecayChainSimulation(DC1Lambda, DC1Mode, n=Rn222_amount)
    Rn220_DCS = DecayChainSimulation(DC2Lambda, DC2Mode, n=Rn220_amount)
    x_222 = Rn222_DCS.expected_counts(ts, p, tf, counted=["a"])
    x_220 = Rn220_DCS.expected_counts(ts, p, tf, counted=["a"])
    rn222_est = [0.] * runs_per_point
    rn220_est = [0.] * runs_per_point
    for k in range(runs_per_point):
        Rn220_DCS.reset_simulation()
        Rn222_DCS.reset_simulation()
        out = (np.array(Rn222_DCS.simulate_counts(ts, p, tf, ["a"])) +
               np.array(Rn220_DCS.simulate_counts(ts, p, tf, ["a"])))

        lr = LinearRegression(fit_intercept=False).fit(np.transpose(np.vstack((x_222, x_220))), out)
        rn222_est[k], rn220_est[k] = lr.coef_
    fstring = "Sample period: {} s, number of samples: {}, Rn{} => mean: {:1.1f}, std: {:1.1f}, activity: {} pCi/L"
    print(fstring.format(ts, p, 222, np.mean(rn222_est), np.std(rn222_est),
                         np.mean(rn222_est) * Rn222_DCS.λ[0] / (0.037 * 0.3)))
    print(fstring.format(ts, p, 220, np.mean(rn220_est), np.std(rn220_est),
                         np.mean(rn220_est) * Rn220_DCS.λ[0] / (0.037 * 0.3)))
    rn222_mean[i][j] = np.mean(rn222_est)
    rn222_stdv[i][j] = np.std(rn222_est)
    rn220_mean[i][j] = np.mean(rn220_est)
    rn220_stdv[i][j] = np.std(rn220_est)


if __name__ == "__main__":
    fstr = "%(asctime)s: %(message)s"
    logging.basicConfig(format=fstr, level=logging.INFO, datefmt="%H:%M:%S")

    rn222_mean = [[0] * grid_duration_num for _ in range(grid_period_num)]
    rn222_stdv = [[0] * grid_duration_num for _ in range(grid_period_num)]
    rn220_mean = [[0] * grid_duration_num for _ in range(grid_period_num)]
    rn220_stdv = [[0] * grid_duration_num for _ in range(grid_period_num)]

    jobs = [
        (
            i,  # Thread ID
            (i % grid_period_num) * grid_period_step + grid_period_step,  # Sampling period
            grid_duration_step * (i // grid_period_num) + grid_duration_first,  # Sampling duration
            90,  # Filling time offset
            i % grid_period_num,  # Output storage row index
            i // grid_period_num,  # Output storage column index
        ) for i in range(grid_duration_num * grid_period_num)
    ]
    threads = list()
    for job in jobs:
        logging.info("Main\t: create and start thread %d", job[0])
        x = threading.Thread(target=grid_point_thread, args=(job,))
        threads.append(x)
        x.start()
    for t in threads:
        t.join()
    with open("../Data/rn222_mean.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rn222_mean)
    with open("../Data/rn222_std.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rn222_stdv)
    with open("../Data/rn220_mean.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rn220_mean)
    with open("../Data/rn220_std.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(rn220_stdv)
