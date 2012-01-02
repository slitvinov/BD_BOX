#! /bin/bash


nskip=10
for xyz in *.xyz; do
    basename=${xyz%.*}
    printf "processing: %s\n" ${basename}
    awk -v cm=1 -f ../../scripts/xyz2punto.awk \
	${basename}.xyz > ${basename}.punto

    awk -v fr=1 -f ../../scripts/xyz2punto.awk ${basename}.xyz | \
	awk -v ns=${nskip} 'NR>ns' > ${basename}.fr

    awk -f ../../scripts/autocorr.awk ${basename}.fr > \
	${basename}.fr-corr-X

    awk -v idx=2 -f ../../scripts/autocorr.awk ${basename}.fr > \
	${basename}.fr-corr-Y
    cp ${basename}.fr-corr-X ~/Dropbox/
    cp ${basename}.fr-corr-Y ~/Dropbox/
done

awk -f ../../scripts/av.awk dna_pwell.proc*.fr-corr-X  > dna_pwell.fr-corr-X
awk -f ../../scripts/av.awk dna_pwell.proc*.fr-corr-Y  > dna_pwell.fr-corr-Y
cp dna_pwell.fr-corr-[XY] ~/Dropbox/

