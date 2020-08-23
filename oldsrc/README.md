The old MDreac
===================

In this directory you find the source code of the Molecular Dynamics
(MD) program I developed as a ph.d. student. It is used for studying a
system undergoing a phase transition while the three species are
participating in a set of chemical reactions.

The source code (`mdreac.F`) is a unified source for a number of
programs:

* `mdreac` is the primary MD program
* `geninit` can be used to generate initial configurations for
  simulations
* `parmdreac` and `pargeninit` are parallel versions of the above
  programs

The parallel versions (using PVM) haven't been built for many years,
and it is unlikely that they will work.

You can find the LaTeX source code of my ph.d. thesis at
https://github.com/kneth/phd-thesis

Building
---------

To build the MD program (`mdreac`) and the aux. program (`geninit`)
you are required to have GNU Fortran installed. Moreover, you will
need to install GNU Make. Building is as simple as:

    make

Usage
------

Both `mdreac` and `geninit` are using a simple text file named
`mdreac.param`. An example is:


```
nA              512
nX              256
nY              256
nZ              0
rho             0.8
rcut1           2.5
rcut2           1.122462
temperature     1.5
nsteps          250000
probability 1   0.001
probability 2   0.0011
probability 3   0.001
Rreaction       0.96116
```

The parameters are:

* `nA`, `nX`, `nY`, `nZ`: The number of particles of each species in
  the initial configuration.
* `rho`: the (number) density.
* `temperature`: the temperature.
* `rcut1`, `rcut2`: Cut-off distance for even and odd pair of particles.
* `nsteps`: Number of time steps to simulate.
* `probability 1`, `probability 2`, `probability 3`: The reaction probabilities.
* `Rreaction`: The distance where reactions can occur.

To perform a simulation, you must first generate the initial
configuration (positions and velocities of the particles). You use
`geninit` for that. No reaction will occur, and you should set
`nsteps` large enough to ensure that the configuration is relaxed.

Using initial configurating (`conf.init`), you can perform the actual
simulation using `mdreac`. The program will write a number of files
with data from the simulation. These output files are plain text
files, and for long simulation it is not uncommon tat `mdreac`
produces vast amount of data.

To analyze the data from the simulations, the following two Python utilities can be used:

* `diffcoeff.py`: Calculation of diffusion coefficient. Requires
  [NumPy](https://numpy.org/).
* `rateconst.py`: Calculation of rate constant. Requires
  [NumPy](https://numpy.org/) and [SciPy](https://scipy.org/).

Protocols
----------

In the directory `protocols` three small shell scripts can be
found. Each scripts run a series of simulations. All parameters except
temperature are fixed, and the series is a sweep over
temperature. Moreover, the scripts are calculating diffusion
coefficients and rate constants. The output can easily be plotted
using software like [GNUplot](http://www.gnuplot.info/).

The three scripts are:

* `large.sh`: Simulates a system of 65536 particles in a series of 16
  temperatures. Both oscillating reactions and phase transition will
  be observed. The total disk usage is 1.5 GB.
* `medium-density.sh`: Simulates a system of 1024 particles in a
  series of 42 temperatures. No phase transition will be
  observed. The total disk usage is 8.4 GB.
* `phd.sh`: Simulates a system of 1024 particles in a series of 42
  temperatures. Both phase transition and oscillating reactions will
  be observed. This is the protocol which is closed to the result
  presented in my Ph.D. thesis. The total disk usage is 8.4 GB.

All three scripts are expected to run for very long. `phd.sh` has
been observed to have the following durations:

* Intel Atom D2700 @ 2.13 GHz: 50.5 hours
* Intel Xeon E3-1241 @ 3.50 GHz: 5 hours (thanks to @fhoefling for
  contributing the number)
