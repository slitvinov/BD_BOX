#! /bin/bash

xyz=$1
basename=${xyz%.*}
nksip=1000
printf "start processing: %s\n" ${basename}
awk -v cm=1 -f ../../scripts/xyz2punto.awk \
    ${basename}.xyz > ${basename}.punto

awk -v fr=1 -f ../../scripts/xyz2punto.awk ${basename}.xyz | \
    awk -v ns=${nskip} 'NR>ns' > ${basename}.fr

awk -f ../../scripts/autocorr.awk ${basename}.fr > \
    ${basename}.fr-corr-X

awk -v idx=2 -f ../../scripts/autocorr.awk ${basename}.fr > \
    ${basename}.fr-corr-Y
    #cp ${basename}.fr-corr-X ~/Dropbox/
    #cp ${basename}.fr-corr-Y ~/Dropbox/

printf "finish processing: %s\n" ${basename}
