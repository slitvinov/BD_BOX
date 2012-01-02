#! /bin/bash

N=20
#pwell_A=1e-5

function run() {
    local iproc=$1
    local id=dna_pwell.proc${iproc}
    local seed=$RANDOM
    echo "seed:" ${seed}
    awk -v N=${N} -f ../../scripts/polygen.awk > ${id}.str
    ../../src/bd_box  dna_pwell.prm \
	--xyz_filename=${id}.xyz \
	--dcd_filename=${id}.dcd \
	--enr_filename=${id}.enr \
	--str_filename=${id}.str \
	--out_filename=${id}.out \
	--pqr_filename=${id}.pqr \
	--rand_seed=${seed} \
	--pwell_A=0.0
}


nproc=8
for ip in $(seq 1 ${nproc}); do
    run ${ip} &
done
