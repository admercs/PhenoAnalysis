#!/bin/sh
export CC=cc
export FC=ftn
export HDF4_DIR=/apps/albis/hdf4/4.2.6/gnu_434_static
./configure --enable-debug --enable-mpi --with-hdf4=${HDF4_DIR}
#./configure --enable-mpi --with-hdf4=${HDF4_DIR}
