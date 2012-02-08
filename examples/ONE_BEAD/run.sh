#! /bin/bash

np=8
seq 1 1 | parallel --verbose -j ${np} ./runone.sh {} 1

