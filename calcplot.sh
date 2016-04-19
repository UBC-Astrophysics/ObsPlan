#!/bin/bash
if [ $# -lt 2 ]; then
    echo Format:
    echo
    echo $0 LIGOSkyMap GalMAP
    exit
fi
out=${1/.fits.gz/}
galmap=${2/.fits.gz/}
if false
then
    for i in 16 64 256
    do
	python ObsPlan.py --gal-map $2 --cumprob 0.99 --no-savefigures ${1} $i > ${out}_${galmap}_${i}.dat
	python ObsPlan.py --gal-map ${galmap}_smoothed.fits.gz --cumprob 0.99 --no-savefigures ${1} $i > ${out}_${galmap}_s_${i}.dat
	python ObsPlan.py  --cumprob 0.99 --no-savefigures ${1} $i > ${out}_raw_${i}.dat
    done
fi
if false
then
    ctioga2 --name ${out}_${galmap} --ylog \
	    -c blue  ${out}_${galmap}_256.dat@5:6  \
	    -c green ${out}_${galmap}_s_256.dat@5:6  \
	    -c red   ${out}_raw_256.dat@5:6  \
	    --line-style dashes \
	    -c blue  ${out}_${galmap}_128.dat@5:6  \
	    -c green ${out}_${galmap}_s_128.dat@5:6  \
	    -c red   ${out}_raw_128.dat@5:6  \
	    --line-style dots \
	    -c blue  ${out}_${galmap}_64.dat@5:6  \
	    -c green ${out}_${galmap}_s_64.dat@5:6  \
	    -c red   ${out}_raw_64.dat@5:6  \
	    --line-style solid \
	    -c blue   ${out}_${galmap}_32.dat@5:6  \
	    -c green  ${out}_${galmap}_s_32.dat@5:6  \
	    -c red    ${out}_raw_32.dat@5:6  \
	    --line-style dashes \
	    -c blue  ${out}_${galmap}_16.dat@5:6  \
	    -c green ${out}_${galmap}_s_16.dat@5:6  \
	    -c red   ${out}_raw_16.dat@5:6  \
	    -x 'Cumulative Probability' -y 'Number of Fields'
fi
if true
then
    ctioga2 --name ${out}_${galmap} --ylog \
	    -c red --line-style solid ${out}_raw_256.dat@5:6  \
	    --line-style dashes ${out}_${galmap}_s_256.dat@5:6  \
	    --line-style dots ${out}_${galmap}_256.dat@5:6  \
	    -c green --line-style solid ${out}_raw_64.dat@5:6  \
	    --line-style dashes ${out}_${galmap}_s_64.dat@5:6  \
	    --line-style dots ${out}_${galmap}_64.dat@5:6  \
	    -c blue --line-style solid ${out}_raw_16.dat@5:6  \
	    --line-style dashes ${out}_${galmap}_s_16.dat@5:6  \
	    --line-style dots ${out}_${galmap}_16.dat@5:6  \
	    -x 'Cumulative Probability' -y 'Number of Fields'
fi    
if false
then
    ctioga2 --name ${out}_${galmap} --ylog \
	    -c red --line-style solid ${out}_raw_256.dat@5:6  \
	    --line-style dashes ${out}_${galmap}_s_256.dat@5:6  \
	    --line-style dots ${out}_${galmap}_256.dat@5:6  \
	    T125738_bayestar_prod_mask_256.dat@5:6 \
	    -c green --line-style solid ${out}_raw_64.dat@5:6  \
	    --line-style dashes ${out}_${galmap}_s_64.dat@5:6  \
	    --line-style dots ${out}_${galmap}_64.dat@5:6  \
	    T125738_bayestar_prod_mask_64.dat@5:6 \
	    -c blue --line-style solid ${out}_raw_16.dat@5:6  \
	    --line-style dashes ${out}_${galmap}_s_16.dat@5:6  \
	    --line-style dots ${out}_${galmap}_16.dat@5:6  \
	    T125738_bayestar_prod_mask_16.dat@5:6 \
	    -x 'Cumulative Probability' -y 'Number of Fields'
fi    
