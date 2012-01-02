#! /bin/bash

N=60
pwell_A=1e-5
awk -v N=${N} -f ../../scripts/polygen.awk > dna_pwell.str
../../src/bd_box  dna_pwell.prm \
--xyz_filename=dna_pwell.xyz \
--dcd_filename=dna_pwell.dcd \
--enr_filename=dna_pwell.enr \
--str_filename=dna_pwell.str \
--out_filename=dna_pwell.out \
--pqr_filename=dna_pwell.pqr \
--pwell_A=0.0
