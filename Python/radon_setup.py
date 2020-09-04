import numpy as np


LRn222 = 3.8235 * 24 * 60 * 60  # the half-life for Rn-222 (in seconds)
LPo218 = 3.098 * 60  # the half-life for Po-218
LPb214 = 26.8 * 60  # the half-life for Pb-214
LBi214 = 19.9 * 60  # the half-life for Bi-214
LPo214 = 164.3e-6  # the half-life for Po-214
DC1HL = np.array([LRn222, LPo218, LPb214, LBi214, LPo214])
DC1Mode = ["a", "a", "b", "b", "a"]
DC1Lambda = np.log(2) / DC1HL  # the decay constants for this first decay chain in units of Hz

LRn220 = 55.6  # the half-life for Rn-220 (in seconds)
LPo216 = 0.145  # the half-life for Po-216
LPb212 = 10.64 * 60 * 60  # the half-life for Pb-212
LBi212 = 60.55 * 60  # the half-life for Bi-212
# Technically, Bi-212 can both alpha and beta decay, but the beta decay mode then alpha decays almost immediately
DC2HL = np.array([LRn220, LPo216, LPb212, LBi212])
DC2Mode = ["a", "a", "b", "a"]
DC2Lambda = np.log(2) / DC2HL  # the decay constants for this first decay chain in units of Hz
