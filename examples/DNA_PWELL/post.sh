#! /bin/bash

awk -v cm=1 -f ../../scripts/xyz2punto.awk \
    dna_pwell.xyz > dna_pwell.punto
awk -v cm=1 -v fr=1 -f ../../scripts/xyz2punto.awk dna_pwell.xyz > \
    dna_pwell.fr
awk -f ../../scripts/autocorr.awk dna_pwell.fr > \
       dna_pwell.fr-corr

