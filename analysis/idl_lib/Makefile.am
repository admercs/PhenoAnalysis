# This script compiles the shared libraries like regular programs
# since LIBTOOL hasn't been stable on all systems

# modify install directories
bindir = @prefix@/lib

# programs (in this case shared libraries) to be built
bin_PROGRAMS = enkf_wrapper_r4.so enkf_wrapper_r8.so

# source files
enkf_wrapper_r4_so_SOURCES = enkf_wrapper_r4.c ../../src/enkf_module_r4.F90 enkf_driver_r4.F90
enkf_wrapper_r8_so_SOURCES = enkf_wrapper_r8.c ../../src/enkf_module_r8.F90 enkf_driver_r8.F90

# headers that are distributed but not installed
# noinst_HEADERS = mpfit.h

# compiler and linker flags follow
AM_CPPFLAGS = @AM_CPPFLAGS@ -fPIC
AM_LDFLAGS  = @SHLIBFLAGS@
LIBS = @MATH_LIBS@
#AM_CFLAGS =
#AM_FCFLAGS =
CC=@CC_NOMPI@
FC=@FC_NOMPI@

