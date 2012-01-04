#! /bin/bash

seq 1 800 | parallel  --verbose ./runone.sh
