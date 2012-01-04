#! /bin/bash

# 4e-6
seq 1 80 | parallel  --verbose ./runone.sh {} 0.0
seq 1 80 | parallel  --verbose ./runone.sh {} 1.0e-6
seq 1 80 | parallel  --verbose ./runone.sh {} 1.6e-5
