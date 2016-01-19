#!/bin/bash
if [ $# -lt 2 ]; then
    echo Format:
    echo
    echo $0 Directory GalRange
    exit
fi
out=SkyMaps/${1}/SkyMap_OutFile${1}
ctioga2 --name ${1}_${2} --ylog \
	-c blue  ${out}_128_${2}.txt@5:6  \
	-c green ${out}_128_${2}smoothed.txt@5:6  \
	-c red   ${out}_128_NoGalMap.txt@5:6  \
    --line-style dashes \
	-c blue  ${out}_64_${2}.txt@5:6  \
	-c green ${out}_64_${2}smoothed.txt@5:6  \
	-c red   ${out}_64_NoGalMap.txt@5:6  \
    --line-style dots \
	-c blue  ${out}_32_${2}.txt@5:6  \
	-c green ${out}_32_${2}smoothed.txt@5:6  \
	-c red   ${out}_32_NoGalMap.txt@5:6  \
    -x 'Cumulative Probability' -y 'Number of Fields'
