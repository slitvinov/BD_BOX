#! /bin/bash

set -e
set -u
source vars.sh

function getparlist() {
    local pat=$1
    ls *.xyz -1 | xargs -I{} -n1 basename {} .xyz  | awk -v RS="^" '/${pat}/' | sed 's/${pat}//g' | sort | uniq
}

ls *.xyz | parallel  --verbose ./postonefile.sh {}

pat=pwell_A
for val in $(getparlist ${pat}); do
    awk -f ${bddir}/scripts/av.awk dna_pwell.proc*${pat}${val}*.fr-corr-X  > dna_pwell.${pat}${val}fr-corr-X
    awk -f ${bddir}/scripts/av.awk dna_pwell.proc*${pat}${val}*.fr-corr-Y  > dna_pwell.${pat}${val}fr-corr-Y
    awk -f ${bddir}/scripts/av.awk dna_pwell.proc*${pat}${val}*.msd  > dna_pwell.${pat}${val}msd
    awk -f ${bddir}/scripts/av.awk dna_pwell.proc*${pat}${val}*.msdcm  > dna_pwell.${pat}${pat}msdcm
    scp dna_pwell.${pat}${val}fr-corr-[XY] litvinov@kana:~/Dropbox/
    scp dna_pwell.${pat}${val}msd litvinov@kana:~/Dropbox/
done

