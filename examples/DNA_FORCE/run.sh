#! /bin/bash

set -e
set -u

np=8
vel_grad=1.16e-5

seq 1 32 | ${HOME}/bin/parallel -j ${np} --verbose ./runone.sh {} ${vel_grad} 1e-3


