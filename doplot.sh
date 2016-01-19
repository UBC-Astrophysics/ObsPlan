#!/bin/bash
out=${1/.fits.gz/}
for i in  ; do
    python ObsPlan.py --gal-map 2MPZ.gz_0.03_0.04.fits.gz --cumprob 0.99 --no-savefigures ${1} $i > ${out}_0.03_0.04_${i}.dat
    python ObsPlan.py --gal-map 2MPZ.gz_0.03_0.04_smoothed.fits.gz --cumprob 0.99 --no-savefigures ${1} $i > ${out}_0.03_0.04_s_${i}.dat
    python ObsPlan.py  --cumprob 0.99 --no-savefigures ${1} $i > ${out}_raw_${i}.dat
done
ctioga2 --name $out --ylog \
	-c blue  ${out}_0.03_0.04_128.dat@5:6  \
	-c green ${out}_0.03_0.04_s_128.dat@5:6  \
	-c red   ${out}_raw_128.dat@5:6  \
    --line-style dashes \
    -c blue  ${out}_0.03_0.04_64.dat@5:6  \
    -c green ${out}_0.03_0.04_s_64.dat@5:6  \
    -c red   ${out}_raw_64.dat@5:6  \
    --line-style dots \
  -c blue   ${out}_0.03_0.04_32.dat@5:6  \
  -c green  ${out}_0.03_0.04_s_32.dat@5:6  \
  -c red    ${out}_raw_32.dat@5:6  \
    --line-style solid \
	-c blue  ${out}_0.03_0.04_16.dat@5:6  \
	-c green ${out}_0.03_0.04_s_16.dat@5:6  \
	-c red   ${out}_raw_16.dat@5:6  \
    -x 'Cumulative Probability' -y 'Number of Fields'
