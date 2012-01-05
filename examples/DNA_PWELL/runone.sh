#! /bin/bash
set -e
set -u
bddir=../../
N=60
iproc=$1
vel_grad=$2
id=dna_pwell.proc${iproc}.vel_grad${vel_grad}
seed=$RANDOM
printf "seed: %s\n" ${seed} > "/dev/stderr"
printf "run: %s\n" ${id}
awk -v N=${N} -f ${bddir}/scripts/polygen.awk > ${id}.str
${bddir}/src/bd_box dna_pwell.prm \
    --xyz_filename=${id}.xyz \
    --str_filename=${id}.str \
    --out_filename=${id}.out \
    --rand_seed=${seed} \
    --vel_grad_tensor="0.0 ${vel_grad} 0.0 0.0 0.0 0.0 0.0 0.0 0.0"

printf "stop: %s\n" ${id} > "/dev/stderr"
