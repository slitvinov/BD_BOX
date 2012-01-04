#! /bin/bash

seq 1 1000 | parallel  --verbose ./runone.sh
