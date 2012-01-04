#! /bin/bash

set -e
set -u

function getparlist() {
    local pat=$1
    ls *.xyz -1 | awk -v RS="." /${pat}/ | sed "s/${pat}//" | sort | uniq
}

ls *.xyz | parallel  --verbose ./postonefile.sh {}

for pat in $(getparlist vel_grad); do
    awk -f ../../scripts/av.awk dna_pwell.proc*vel_grad${pat}*.fr-corr-X  > dna_pwell.vel_grad${pat}fr-corr-X
    awk -f ../../scripts/av.awk dna_pwell.proc*vel_grad${pat}*.fr-corr-Y  > dna_pwell.vel_grad${pat}fr-corr-Y
    scp dna_pwell.vel_grad${pat}fr-corr-[XY] litvinov@kana:~/Dropbox/
done

