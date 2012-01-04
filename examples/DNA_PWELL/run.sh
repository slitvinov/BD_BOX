#! /bin/bash

# 4e-6
seq 1 80 | parallel  --verbose ./runone.sh {} 2e-6
seq 1 80 | parallel  --verbose ./runone.sh {} 4e-6
seq 1 80 | parallel  --verbose ./runone.sh {} 8e-6
