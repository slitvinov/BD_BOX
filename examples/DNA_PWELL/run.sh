#! /bin/bash

seq 17 1000 | parallel  --verbose ./runone.sh
