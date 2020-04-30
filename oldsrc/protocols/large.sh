#!/bin/bash

rho=0.8

mkdir large-${rho}
cd large-${rho}

for T in 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.6 3.0 3.5 4.0 5.0 6.0 8.0 10.0
do
    mkdir $T
    cd $T
    cat <<EOF >>mdreac.param
nA              32768
nX              16384
nY              16384 
nZ              0
rho             ${rho}
rcut1           2.5
rcut2           1.122462
temperature     ${T}
nsteps          25000
probability 1   0.001
probability 2   0.0011
probability 3   0.001
Rreaction       0.96116
EOF
    cp ../../conf.init.65536 conf.init
    ../../../mdreac
    cat mdreac.diff  | awk '{print $1 "," $2}' > mdreac.diff.csv
    D=$(python3 ../../../diffcoeff.py)
    cat mdreac.reac | tr -s " " | cut -c2- | tr " " ","  | cut -f2- -d"," > mdreac.reac.csv
    k=$(python3 ../../../rateconst.py -N 65536 -d ${rho})
    echo "${T},${D},${k}"
    cd ..
done
