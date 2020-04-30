# Calculate rate constants from event data
# Processing of raw simulation data:
# cat mdreac.reac | tr -s " " | cut -c2- | tr " " ","  | cut -f2- -d"," > mdreac.reac.csv

import sys, getopt

import numpy as np
from numpy import genfromtxt

def usage(prog):
    print('Usage: ', prog, ' -d number -N number')
    sys.exit(2)

def main(argv):
    rho = 0.8
    N = 1024

    try:
        opts, args = getopt.getopt(argv, 'hd:N:')
    except getopt.GetoptError:
        usage(argv[0])
    for opt, arg in opts:
        if opt == '-h':
            usage(argv[0])
        elif opt == '-d':
            rho = float(arg)
        elif opt == '-N':
            N = float(arg)

    data = genfromtxt('mdreac.reac.csv', delimiter=',')
    m = np.mean(data, 0)[0]
    print(2*m/(rho*N))

if __name__ == "__main__":
   main(sys.argv[1:])