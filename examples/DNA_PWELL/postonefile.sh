#! /bin/bash

set -e
set -u
xyz=$1
basename=${xyz%.*}
nskip=1000
printf "start processing: %s\n" ${basename}
awk -v cm=1 -f ../../scripts/xyz2punto.awk \
    ${basename}.xyz > ${basename}.punto

awk -v fr=1 -f ../../scripts/xyz2punto.awk ${basename}.xyz | \
    awk -v ns=${nskip} 'NR>ns' > ${basename}.fr

awk -f ../../scripts/autocorr.awk ${basename}.fr > \
    ${basename}.fr-corr-X

awk -v idx=2 -f ../../scripts/autocorr.awk ${basename}.fr > \
    ${basename}.fr-corr-Y

printf "finish processing: %s\n" ${basename}
