# Calculate diffusion coefficient
# Preprocessing of raw simulation data:
# cat mdreac.diff  | awk '{print $1 "," $2}' > mdreac.diff.csv

import numpy as np
from numpy import genfromtxt

from sklearn.linear_model import LinearRegression

data = genfromtxt('mdreac.diff.csv', delimiter=',')
t = data[:, 0].reshape(-1, 1)
d = data[:, 1]

regressor = LinearRegression()
model = regressor.fit(t, d)
print(0.25 * model.coef_[0])

