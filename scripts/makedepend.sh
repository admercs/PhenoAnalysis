#!/bin/sh

# wrapper for makedepend.pl
# generates Fortran 90 dependencies including module depencencies for
# source directory specified below and adds the dependencies to the 
# Makefile.am file in the source directory

curdir=`pwd`

srcdir="../src"
file_src="sourcefiles.dat"
file_dir="sourcedirs.dat"
file_dep="dependencies.dat"

cd ${srcdir}

echo ${srcdir} > ${file_dir}
ls -1 *.F90 > ${file_src}

# remove old dependencies from Makefile.am
sed '1,/# DO NOT EDIT BELOW THIS LINE/p;/$/ d' Makefile.am > Makefile.tmp
mv Makefile.tmp Makefile.am

# generate new dependencies
${curdir}/makedepend.pl ${file_dir} ${file_src} > ${file_dep}

cat ${file_dep} >> Makefile.am

# remove temporary files
rm ${file_src}
rm ${file_dir}
rm ${file_dep}

