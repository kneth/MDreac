#!/bin/bash

rho=0.8

mkdir density-${rho}
cd density-${rho}

for T in 0.4 0.6 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.33 3.66 4.0 5.0 6.0 7.0 8.0 9.0 10.0
do
    echo "#######################"
    echo "T = ${T}, rho = ${rho}"
    echo "#######################"

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
nsteps          500000
probability 1   0.001
probability 2   0.0011
probability 3   0.001
Rreaction       0.96116
EOF
    cp ../../conf.init.1024 conf.init
    ../../../mdreac
    cat mdreac.diff  | awk '{print $1 "," $2}' > mdreac.diff.csv
    D=$(python3 ../../../diffcoeff)
    cat mdreac.reac | tr -s " " | cut -c2- | tr " " ","  | cut -f2- -d"," > mdreac.reac.csv
    k=$(python3 ../../../rateconst.py -N 1024 -d ${rho})
    echo "${T},${D},${k}"
    cd ..
done