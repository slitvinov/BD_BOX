#! /bin/bash

set -e
set -u

np=1
N=60
seq 1 1 | ${HOME}/bin/parallel -j ${np} --verbose ./runone.sh {}  ${N}


