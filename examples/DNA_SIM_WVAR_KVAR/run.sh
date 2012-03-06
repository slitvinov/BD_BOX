#! /bin/bash

set -u

# Wi=240
vellist="0 2.88e-7 4.33e-6 1.15e-5 2.88e-4"
nbeads=60
klist="\
 0.0 \
 1e-5 \
 4e-5 \
 1e-4 \
 2e-4 \
 3e-4 \
 4e-4 \
 1e-3"

for vel in ${vellist}; do
    for k in ${klist}; do
	for ip in $(seq 1 1); do
            echo ${ip}
            echo ${vel}
	    echo ${k}
	done
    done 
done| parallel -j 16 -N3 ./runone.sh ${nbeads} {1} {2} {3}
