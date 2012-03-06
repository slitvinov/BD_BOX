#! /bin/bash
set -e
set -u
set -x

source vars.sh
N=$1
iproc=$2
vel_grad=$3
pwell_A=$4

id=dna_free~proc${iproc}~nb${N}~vel${vel_grad}~pwell_A${pwell_A}
seed=$RANDOM
printf "seed: %s\n" ${seed} > "/dev/stderr"
printf "id: %s\n" ${id} > "/dev/stderr"
awk -v N=${N} -f ${bddir}/scripts/polygen.awk > ${id}.str

nice -n 19 ${bddir}/src/bd_box dna_free.prm \
    --str_filename=${id}.str \
    --out_filename=${id}.out \
    --dcd_filename=${id}.dcd \
    --rand_seed=${seed} \
    --vel_grad_tensor="0.0 ${vel_grad} 0.0 0.0 0.0 0.0 0.0 0.0 0.0" \
    --pwell_A=${pwell_A}

printf "stop: %s\n" ${id} > "/dev/stderr"
