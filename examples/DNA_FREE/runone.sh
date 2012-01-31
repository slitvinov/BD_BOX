#! /bin/bash
set -e
set -u
set -x

source vars.sh

iproc=$1
N=$2

id=dna_free~proc.${iproc}
seed=$RANDOM
printf "seed: %s\n" ${seed} > "/dev/stderr"
printf "id: %s\n" ${id} > "/dev/stderr"
awk -v N=${N} -f ${bddir}/scripts/polygen.awk > ${id}.str

${bddir}/src/bd_box dna_free.prm \
    --str_filename=${id}.str \
    --out_filename=${id}.out \
    --dcd_filename=${id}.dcd \
    --rand_seed=${seed} \

printf "stop: %s\n" ${id} > "/dev/stderr"
