#!/bin/bash

rho=0.8

mkdir density-${rho}
cd density-${rho}

for T in 0.4 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.2 3.4 3.6 3.8 4.0 4.33 4.66 5.0 5.33 5.66 6.0 6.5 7.0 7.5 8.0 9.0 10.0
do
    mkdir $T
    cd $T
    cat <<EOF >>mdreac.param
nA              512
nX              256
nY              256
nZ              0
rho             ${rho}
rcut1           2.5
rcut2           2.5
temperature     ${T}
nsteps          250000
probability 1   0.001
probability 2   0.0011
probability 3   0.001
Rreaction       0.96116
EOF
    cp ../../conf.init.1024 conf.init
    ../../../mdreac
    cat mdreac.diff  | awk '{print $1 "," $2}' > mdreac.diff.csv
    D=$(python3 ../../../diffcoeff.py)
    cat mdreac.reac | tr -s " " | cut -c2- | tr " " ","  | cut -f2- -d"," > mdreac.reac.csv
    k=$(python3 ../../../rateconst.py -N 1024 -d ${rho})
    echo "${T},${D},${k}"
    cd ..
done
