#! /bin/bash
set -e
set -u
N=60
iproc=$1
id=dna_pwell.proc${iproc}
seed=$RANDOM
printf "seed: %s\n" ${seed} > "/dev/stderr"
printf "run: %s\n" ${id}
awk -v N=${N} -f ../../scripts/polygen.awk > ${id}.str
../../src/bd_box  dna_pwell.prm \
    --xyz_filename=${id}.xyz \
    --str_filename=${id}.str \
    --out_filename=${id}.out \
    --rand_seed=${seed} \
    --pwell_A=0.0
printf "stop: %s\n" ${id} > "/dev/stderr"
