#! /bin/bash
set -e
set -u
set -x

source vars.sh
N=60
iproc=$1

id=dna_free
seed=$RANDOM
printf "seed: %s\n" ${seed} > "/dev/stderr"
echo -n "${id}" > "/dev/stderr"
awk -v N=${N} -f ${bddir}/scripts/polygen.awk > ${id}.str
${bddir}/src/bd_box dna_free.prm \
    --xyz_filename=${id}.xyz \
    --str_filename=${id}.str \
    --out_filename=${id}.out \
    --dcd_filename=${id}.dcd \
    --rand_seed=${seed} \

printf "stop: %s\n" ${id} > "/dev/stderr"
