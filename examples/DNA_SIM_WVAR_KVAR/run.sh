#! /bin/bash

set -u

# Wi=240
#vellist=1.92e-05
vellist=0.0
klist=1e-5

nbeads=60

# list of all $k$
# klist="\
# 1e-5 \
# 4e-5 \
# 1e-4 \
# 2e-4 \
# 3e-4 \
# 4e-4 \
# 1e-3"

for vel in ${vellist}; do
    for k in ${klist}; do
	for ip in $(seq 1 1); do
            echo ${ip}
            echo ${vel}
	    echo ${k}
	done
    done 
done| ~/bin/parallel -j 100 -N3 --verbose ./runone.sh ${nbeads} {1} {2} {3}
