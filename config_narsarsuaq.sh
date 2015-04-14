#!/bin/sh

# This is a sample shell script that runs the configuration on OSX (fink+locally compiled netcdf / hdf libs)

#./configure --enable-debug --with-netcdf=/usr/local/netcdf3 --with-hdf4=/usr/local/hdf4 --with-jpeg=/sw
#./configure --enable-debug --with-netcdf=/usr/local/netcdf4 --with-hdf4=/usr/local/hdf4 --with-jpeg=/sw
./configure --enable-debug --enable-mpi --with-netcdf=/usr/local/netcdf4-parallel --with-hdf4=/usr/local/hdf4 --with-jpeg=/sw --with-mpi=/usr/local/openmpi
