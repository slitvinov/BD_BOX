#! /bin/bash

set -e
set -u
source vars.sh

function getparlist() {
    local pat=$1
    ls *.xyz -1 | xargs -I{} -n1 basename {} .xyz  | awk -v RS="^" 1 | awk "/${pat}/" | sed "s/${pat}//" | sort | uniq
}

pat=pwell_A
ns=1000

for val in $(getparlist ${pat}); do
    res=$(awk -v ns=${ns} -f ${bddir}/scripts/ext.awk dna_pwell^proc*${pat}${val}*.punto)
    echo ${val} ${res}
done

