#! /bin/bash

set -e
set -u
bddir=../../

function getparlist() {
    local pat=$1
    ls *.xyz -1 | awk -v RS="." /${pat}/ | sed "s/${pat}//" | sort | uniq
}

ls *.xyz | parallel  --verbose ./postonefile.sh {}

for pat in $(getparlist vel_grad); do
    awk -f ${bddir}/scripts/av.awk dna_pwell.proc*vel_grad${pat}*.fr-corr-X  > dna_pwell.vel_grad${pat}fr-corr-X
    awk -f ${bddir}/scripts/av.awk dna_pwell.proc*vel_grad${pat}*.fr-corr-Y  > dna_pwell.vel_grad${pat}fr-corr-Y
    awk -f ${bddir}/scripts/av.awk dna_pwell.proc*vel_grad${pat}*.msd  > dna_pwell.vel_grad${pat}msd
    scp dna_pwell.vel_grad${pat}fr-corr-[XY] litvinov@kana:~/Dropbox/
    scp dna_pwell.vel_grad${pat}msd litvinov@kana:~/Dropbox/
done

