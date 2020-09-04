import logging

from decay_chain_simulation import DecayChainSimulation
from radon_setup import DC1Lambda, DC1Mode, DC2Lambda, DC2Mode
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

files = ['../Data/Radon_5min_Trial1.xlsx',
         '../Data/Radon_5min_Trial2.xlsx',
         '../Data/Radon_5min_Trial3.xlsx',
         '../Data/Thoron_5min_Trial1.xlsx',
         ]
sample_times = [3, 3, 3, 3]
filling_times = [90, 90, 90, 90]

if __name__ == "__main__":

    Rn222_DCS = DecayChainSimulation(DC1Lambda, DC1Mode)
    Rn220_DCS = DecayChainSimulation(DC2Lambda, DC2Mode)
    assert (len(files) > 0) and (len(files) == len(sample_times)) and (len(files) == len(filling_times))
    sq = np.sqrt(len(files))
    nc = int(sq)
    if sq-int(sq) >= 0.5:
        nc += 1
    nr = int(np.ceil(len(files)/nc))
    fig, axes = plt.subplots(nr, nc)
    axes = [ax for row in axes for ax in row]
    fig.show()
    for file, sample_time, filling_time, ax in zip(files, sample_times, filling_times, axes[:len(files)]):
        RadonDF = pd.read_excel(file)
        number_samples = len(RadonDF.index)

        x_222 = Rn222_DCS.expected_counts(sample_time, number_samples, filling_time, counted=["a"])
        x_220 = Rn220_DCS.expected_counts(sample_time, number_samples, filling_time, counted=["a"])

        lr = LinearRegression(fit_intercept=False).fit(np.transpose(np.vstack((x_222, x_220))), RadonDF['Counts'])
        rn222_est, rn220_est = lr.coef_

        act_rn222 = [rn222_est * DC1Lambda[0] / 3.7e-2 / 0.3, rn222_est * DC1Lambda[0] / 0.3e-3]
        act_rn220 = [rn220_est * DC2Lambda[0] / 3.7e-2 / 0.3, rn220_est * DC2Lambda[0] / 0.3e-3]
        fstr = "st: {}s, tt: {}s, Rn{} => Estimate: {:1f}, Activity: {:1f} pCi/L, {:1f} Bq/m^3"
        print(fstr.format(sample_time, number_samples * sample_time, 222, rn222_est, act_rn222[0], act_rn222[1]))
        print(fstr.format(sample_time, number_samples * sample_time, 220, rn220_est, act_rn220[0], act_rn220[1]))

        Rn222_DCS.setup_simulation(int(np.ceil(rn222_est)))
        Rn220_DCS.setup_simulation(int(np.ceil(rn220_est)))
        Rn222_sim = Rn222_DCS.simulate_counts(sample_time, number_samples, filling_time, ["a"])
        Rn220_sim = Rn220_DCS.simulate_counts(sample_time, number_samples, filling_time, ["a"])
        gendata = pd.DataFrame({'Counts':
                                    (np.array(Rn222_sim) +
                                     np.array(Rn220_sim))
                                })

        expected = pd.DataFrame({'Expected':
                                     (int(np.ceil(rn222_est)) * np.array(x_222) +
                                      int(np.ceil(rn220_est)) * np.array(x_220))
                                 })

        ax.plot(RadonDF.index, RadonDF["Counts"], marker='o', markerfacecolor='blue', markersize=4,
                             color='blue',
                             linewidth=1, label="Radon")
        ax.plot(gendata.index, gendata, marker='o', markerfacecolor='red', markersize=4, color='red',
                             linewidth=1,
                             label="Prediction")
        ax.plot(expected.index, expected, marker='', color='black', linewidth=3, label="Expected Value")
        ax.legend()
        fig.canvas.draw()
    plt.show()
