# Calculate diffusion coefficient
# Preprocessing of raw simulation data:
# cat mdreac.diff  | awk '{print $1 "," $2}' > mdreac.diff.csv

import sys, getopt

import numpy as np
from numpy import genfromtxt

from sklearn.linear_model import LinearRegression

def usage(prog):
    print('Usage: ', prog, ' -N number')
    sys.exit(0)

def main(argv):
    N = 1024
    try:
        opts, args = getopt.getopt(argv, 'hd:N:')
    except getopt.GetoptError:
        usage(argv[0])
    for opt, arg in opts:
        if opt == '-h':
            usage(argv[0])
        elif opt == '-N':
            N = float(arg)


    data = genfromtxt('mdreac.diff.csv', delimiter=',')
    t = data[:, 0].reshape(-1, 1)
    d = data[:, 1]

    regressor = LinearRegression()
    model = regressor.fit(t, d)
    print(0.25 * N * model.coef_[0])

if __name__ == "__main__":
    main(sys.argv[1:])