#! /bin/bash
set -e
set -u
set -x

source vars.sh

# tau_60: 1.2479655e+7 [Ps]
# 
iproc=$1
N=$2
vel_grad=1.84e-5

id=dna_free~proc${iproc}~nb${N}~vel${vel_grad}
seed=$RANDOM
printf "seed: %s\n" ${seed} > "/dev/stderr"
printf "id: %s\n" ${id} > "/dev/stderr"
awk -v N=${N} -f ${bddir}/scripts/polygen.awk > ${id}.str

${bddir}/src/bd_box dna_free.prm \
    --str_filename=${id}.str \
    --out_filename=${id}.out \
    --dcd_filename=${id}.dcd \
    --rand_seed=${seed} \
    --vel_grad_tensor="0.0 ${vel_grad} 0.0 0.0 0.0 0.0 0.0 0.0 0.0"

printf "stop: %s\n" ${id} > "/dev/stderr"
