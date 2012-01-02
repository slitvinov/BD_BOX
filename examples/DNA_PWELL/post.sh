#! /bin/bash



for xyz in *.xyz; do
    basename=${xyz%.*}
    printf "processing: %s\n" ${basename}
    awk -v cm=1 -f ../../scripts/xyz2punto.awk \
	${basename}.xyz > ${basename}.punto

    awk -v fr=1 -f ../../scripts/xyz2punto.awk ${basename}.xyz > \
	${basename}.fr

    awk -f ../../scripts/autocorr.awk ${basename}.fr > \
	${basename}.fr-corr-X

    awk -v idx=2 -f ../../scripts/autocorr.awk ${basename}.fr > \
	${basename}.fr-corr-Y
done
