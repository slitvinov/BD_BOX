#! /bin/bash


ls *.xyz | parallel  --verbose ./postonefile.sh {}

awk -f ../../scripts/av.awk dna_pwell.proc*.fr-corr-X  > dna_pwell.fr-corr-X
awk -f ../../scripts/av.awk dna_pwell.proc*.fr-corr-Y  > dna_pwell.fr-corr-Y
scp dna_pwell.fr-corr-[XY] litvinov@kana:~/Dropbox/fujie
