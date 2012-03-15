#! /bin/bash
set -e
set -u
set -x

source vars.sh

iproc=1
ktheta=$1
N=60

id=dna_free~proc${iproc}~nb${N}~ktheta${ktheta}
seed=$RANDOM
printf "seed: %s\n" ${seed} > "/dev/stderr"
printf "id: %s\n" ${id} > "/dev/stderr"
awk -v N=${N} \
    -v useangle=1 \
    -v theta0=180.0 \
    -v ktheta=${ktheta} \
    -f ${bddir}/scripts/polygen.awk > ${id}.str

${bddir}/src/bd_box dna_free.prm \
    --str_filename=${id}.str \
    --out_filename=${id}.out \
    --dcd_filename=${id}.dcd \
    --rst_filename=${id}.rst \
    --rand_seed=${seed} \

printf "stop: %s\n" ${id} > "/dev/stderr"