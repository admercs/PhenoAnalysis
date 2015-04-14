#!/bin/sh

# This is a sample shell script that runs the configuration on the NCAR/UCAR Yellowstone Supercomputer

#Unload the NetCDF module before running this script
# module rm netcdf

HDF4_DIR=/glade/u/home/stockli
export CFLAGS=${INC_NCAR}
export FCFLAGS=${INC_NCAR}
export LDFLAGS="-L${HDF4_DIR}/lib ${LIB_NCAR}"
./configure --enable-debug --enable-mpi
