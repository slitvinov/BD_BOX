#! /bin/bash

# 4e-6
seq 1 800 | parallel  --verbose ./runone.sh {} 0.0
