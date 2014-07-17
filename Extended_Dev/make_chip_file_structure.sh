#!/bin/bash
BINDIR=/project/projectdirs/cmb/modules/carver/hpcports_gnu/wcstools-3.8.7_55c82b60-3.0/bin
file=$1

ccd=`${BINDIR}/gethead ccdid ${file}`

if [${#ccd} = 1]; then 
chip=C0${ccd}
mkdir -p $chip
fi
 
if [${#ccd} = 2]; then
chip=C${ccd}
mkdir -p $chip
fi

mv $1 $chip