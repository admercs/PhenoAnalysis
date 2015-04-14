#!/bin/sh

# This script performs geographical subsetting on all netcdf files
# in a directory and writes them into another directory

# input directory
#indir="/Users/stockli/ecmwf/era-interim"
indir="/store/msclim/stockli/ecmwf/era-interim"

# output directory
#outdir="/Users/stockli/ecmwf/era-interim-alps"
outdir="/store/msclim/stockli/ecmwf/era-interim-harvard"

# geographical subdomain of output files
lonmin=5.0
lonmax=10.0
latmin=45.0
latmax=50.0

lonmin=-75.0
lonmax=-70.0
latmin=40.0
latmax=45.0

# input file prefix
prefix="*"

# input file suffix
suffix=".nc"

# generate output directory
if [ ! -d ${outdir} ]
then
    mkdir ${outdir}
fi

# find list of files matching pattern
infiles=`find ${indir} -name "${prefix}*${suffix}"`

# extract geographical subset and write to output directory
for infile in ${infiles}
do
    f=`basename ${infile}`
    cdo -sellonlatbox,${lonmin},${lonmax},${latmin},${latmax} ${infile} ${outdir}/${f}

done



