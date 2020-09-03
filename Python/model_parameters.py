from decay_chain import DecayChain
from radon_setup import DC1Lambda, DC1Mode, DC2Lambda, DC2Mode
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    Rn222_DC = DecayChain(DC1Lambda, DC1Mode)
    Rn220_DC = DecayChain(DC2Lambda, DC2Mode)

    sample_time = 3  # in seconds
    sample_duration = 60 * 5  # in seconds
    filling_time = 90  # in seconds

    x_220 = Rn220_DC.expected_counts(sample_time, sample_duration // sample_time, filling_time, counted=["a"])
    x_222 = Rn222_DC.expected_counts(sample_time, sample_duration // sample_time, filling_time, counted=["a"])

    print("x_220: ")
    print(x_220)
    plt.plot(x_220)
    plt.show()
    print("x_222: ")
    print(x_222)
    plt.plot(x_222)
    plt.show()

    X = np.transpose(np.vstack((x_222, x_220)))

    print("X: ")
    print(X)

    X_pi = np.matmul(np.linalg.inv(np.matmul(np.transpose(X), X)), np.transpose(X))
    print("(X' X)\\X': ")
    print(X_pi)
