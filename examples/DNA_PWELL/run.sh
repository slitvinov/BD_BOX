#! /bin/bash

np=8
seq 1 80 | ~/bin/parallel  -j ${np} --verbose ./runone.sh {} 0.0
