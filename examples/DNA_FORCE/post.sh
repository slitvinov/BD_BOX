#! /bin/bash

set -e
set -u
source vars.sh

function getparlist() {
    local pat=$1
    ls *.xyz -1 | xargs -I{} -n1 basename {} .xyz  | awk -v RS="^" 1 | awk "/${pat}/" | sed "s/${pat}//" | sort | uniq
}

pat=pwell_A
ls *4e-4*.xyz | parallel  --verbose ./postonefile.sh {}

for val in 4e-4; do
    awk -f ${bddir}/scripts/av.awk dna_pwell^proc*${pat}${val}*.fr-corr-X  > dna_pwell.${pat}${val}fr-corr-X
    awk -f ${bddir}/scripts/av.awk dna_pwell^proc*${pat}${val}*.fr-corr-Y  > dna_pwell.${pat}${val}fr-corr-Y
    awk -f ${bddir}/scripts/av.awk dna_pwell^proc*${pat}${val}*.msd  > dna_pwell.${pat}${val}msd
    awk -f ${bddir}/scripts/av.awk dna_pwell^proc*${pat}${val}*.msdcm  > dna_pwell.${pat}${pat}msdcm
    scp dna_pwell.${pat}${val}fr-corr-[XY] litvinov@kana:~/Dropbox/
    scp dna_pwell.${pat}${val}msd litvinov@kana:~/Dropbox/
done

