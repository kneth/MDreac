# Calculate diffusion coefficient
# Preprocessing of raw simulation data:
# cat mdreac.diff  | awk '{print $1 "," $2}' > mdreac.diff.csv

import numpy as np
from numpy import genfromtxt
from scipy import stats

data = genfromtxt('mdreac.diff.csv', delimiter=',')
slope, intercept, r_value, p_value, std_err = stats.linregress(data[:, 0], data[:, 1])
print(0.25 * slope)
