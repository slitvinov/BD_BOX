#! /bin/bash

set -e
set -u

np=1
N=60
seq 1 1 | parallel -j ${np} ./runone.sh {} ${N}


