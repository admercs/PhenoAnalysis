#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.50])
AC_INIT([phenoanalysis], [2.0.0], [reto.stockli@gmail.com])
AC_CONFIG_SRCDIR([src/phenoanalysis.F90])

# Checks for Operating System
AC_CANONICAL_HOST
#AC_CANONICAL_TARGET

# produce makefiles (turn on warnings on automake)
AM_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE([-Wall -Werror tar-ustar])

# Checks for programs and compilers
if test "x$CC" == x; then
    CCLIST="icc gcc xlc_r xlc pgcc cc"
else
    CCLIST=$CC
fi
OLD_CFLAGS=$CFLAGS
AC_PROG_CC($CCLIST)
CFLAGS=$OLD_CFLAGS

if test "x$FC" == x; then
    FCLIST="ifort gfortran xlf90_r xlf95_r xlf90 xlf95 pgf90 pgf95 f90 f95 ftn"
else
    FCLIST=$FC
fi
OLD_FCFLAGS=$FCFLAGS
AC_PROG_FC($FCLIST)
FCFLAGS=$OLD_FCFLAGS

# AM_PROG_CC_C_O

AC_CHECK_PROG(HAVE_sed,sed,yes,no)
if test "x$HAVE_sed" = xno; then
    AC_MSG_ERROR([please make sure that the <sed> program can be found])
fi

AC_CHECK_PROG(HAVE_perl,perl,yes,no)
if test "x$HAVE_perl" = xno; then
    AC_MSG_ERROR([please make sure that the <perl> program can be found])
fi

AC_CHECK_PROG(HAVE_POE,bsub,yes,no)
AC_CHECK_PROG(HAVE_PBS,qsub,yes,no)
AC_CHECK_PROG(HAVE_SLURM,sbatch,yes,no)

if test "x$HAVE_POE" = xyes; then
    AC_DEFINE(HAVE_POE,1,Define 1 to run with the POE job submission system)
    AC_MSG_NOTICE(The POE job submission system will be used)
else
    AC_DEFINE(HAVE_POE,0,Define 1 to run with the POE job submission system)
    AC_MSG_NOTICE(The code will be executed without the POE job submission system)
fi

if test "x$HAVE_PBS" = xyes && test "x$HAVE_POE" = xno; then
    AC_DEFINE(HAVE_PBS,1,Define 1 to run with the PBS job submission system)
    AC_MSG_NOTICE(The PBS job submission system will be used)
else
    AC_DEFINE(HAVE_PBS,0,Define 1 to run with the PBS job submission system)
    AC_MSG_NOTICE(The code will be executed without the PBS job submission system)
fi

if test "x$HAVE_SLURM" = xyes && test "x$HAVE_POE" = xno && test "x$HAVE_PBS" = xno; then
    AC_DEFINE(HAVE_SLURM,1,Define 1 to run with the SLURM job submission system)
    AC_MSG_NOTICE(The SLURM job submission system will be used)
else
    AC_DEFINE(HAVE_SLURM,0,Define 1 to run with the SLURM job submission system)
    AC_MSG_NOTICE(The code will be executed without the SLURM job submission system)
fi

# Check for Fortran free-form code and F90 extension
AC_LANG_PUSH(Fortran)
AC_FC_SRCEXT(F90)
AC_FC_FREEFORM
AC_LANG_POP(Fortran)

# Enable Debugging?
AC_MSG_CHECKING([whether to build with debug information])
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
        [enable compilation with debugging information (default=no)])],
    [ENABLE_DEBUG=yes],
    [ENABLE_DEBUG=no])
AC_MSG_RESULT([$ENABLE_DEBUG])

if test "x$ENABLE_DEBUG" = xyes; then
    AC_DEFINE(ENABLE_DEBUG,1,[Debugging Mode])
else
    AC_DEFINE(ENABLE_DEBUG,0,[Optimization Mode])
fi

# Enable MPI-type parallelization?
AC_MSG_CHECKING([whether to build with MPI support])
AC_ARG_ENABLE([mpi],
    [AS_HELP_STRING([--enable-mpi],
        [enable compilation with MPI support (default=no)])],
    [ENABLE_MPI=yes],
    [ENABLE_MPI=no])
AC_MSG_RESULT([$ENABLE_MPI])

# Checks for libraries.
# AC_CHECK_LIB([m],[ceil],,AC_MSG_NOTICE(math library libm is not present, but may be needed for HDF5))

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([assert.h glob.h stdio.h getopt.h math.h],,AC_MSG_ERROR([required system header file not found]))

# JPEG (needed for hdf4)
HAVE_JPEG=yes
AC_MSG_CHECKING([whether a location for the JPEG library was specified])
AC_ARG_WITH([jpeg],[AS_HELP_STRING([--with-jpeg=<directory>],
                    [Specify location the JPEG library. 
                     Configure will expect to find subdirs include and lib.])],
		   [AC_MSG_RESULT([$withval])
		    LDFLAGS="$LDFLAGS -L$withval/lib"
    		    CPPFLAGS="$CPPFLAGS -I$withval/include"],
    		   [AC_MSG_RESULT([no])])

AC_CHECK_HEADERS([jpeglib.h],,HAVE_JPEG=no)
AC_CHECK_LIB([jpeg],[jpeg_start_decompress],,HAVE_JPEG=no)

# ZLIB (needed for hdf4)
HAVE_ZLIB=yes
AC_MSG_CHECKING([whether a location for the zlib library was specified])
AC_ARG_WITH([zlib],[AS_HELP_STRING([--with-zlib=<directory>],
                    [Specify location the zlib library. 
                     Configure will expect to find subdirs include and lib.])],
		   [AC_MSG_RESULT([$withval])
		    LDFLAGS="$LDFLAGS -L$withval/lib"
    		    CPPFLAGS="$CPPFLAGS -I$withval/include"],
    		   [AC_MSG_RESULT([no])])

AC_CHECK_HEADERS([zlib.h],,HAVE_ZLIB=no)
AC_CHECK_LIB([z],[compress2],,HAVE_ZLIB=no)

# SZIP (SZLIB, alternative for ZLIB)
HAVE_SZLIB=yes
AC_MSG_CHECKING([whether a location for the szlib library was specified])
AC_ARG_WITH([szlib],[AS_HELP_STRING([--with-szlib=<directory>],
                    [Specify location the szlib library. 
                     Configure will expect to find subdirs include and lib.])],
		   [AC_MSG_RESULT([$withval])
		    LDFLAGS="$LDFLAGS -L$withval/lib"
    		    CPPFLAGS="$CPPFLAGS -I$withval/include"],
    		   [AC_MSG_RESULT([no])])

AC_CHECK_HEADERS([szlib.h],,HAVE_SZLIB=no)
AC_CHECK_LIB([sz],[SZ_BufftoBuffCompress],,HAVE_SZLIB=no)

# HDF4
HAVE_HDF4=yes
AC_MSG_CHECKING([whether a location for the HDF4 library was specified])
AC_ARG_WITH([hdf4],[AS_HELP_STRING([--with-hdf4=<directory>],
                    [Specify location the HDF4 library. 
                     Configure will expect to find subdirs include and lib.])],
		   [AC_MSG_RESULT([$withval])
		    LDFLAGS="$LDFLAGS -L$withval/lib"
    		    CPPFLAGS="$CPPFLAGS -I$withval/include"],
    		   [AC_MSG_RESULT([no])])

# No NetCDF symbols in HDF4 library allowed!
AC_CHECK_HEADERS([mfhdf.h hdf.h],,HAVE_HDF4=no)
AC_CHECK_LIB([df],[Hopen],,HAVE_HDF4=no)
AC_CHECK_LIB([mfhdf],[hdf_read_attrs],,HAVE_HDF4=no)
if test "x$HAVE_HDF4" = xno; then
   AC_MSG_ERROR([Can't find HDF4 library, or can't link a program to it.])
fi

# One of the math libraries below is needed for successful compilation
# On Darwin VECLIB will be automatically chosen.
# On the IBM for the XLF compiler the ESSL library is the preferred math library
# On AMD Processors the ACML math library should be used
# On Intel Processors the MKL math library should be used
# On any other host/cpu, please either link the open source Lapack and BLAS, or the ATLAS libraries

MATH_LIBS=""
HAVE_MATH=no
HAVE_VECLIB=no
HAVE_BLAS=no
HAVE_LAPACK=no
HAVE_ATLAS=no
HAVE_ACML=no
HAVE_MKL=no
HAVE_ESSL=no
HAVE_MASS=no

# predefine math libraries depending on host type
case "$host_os" in
     *darwin*)
     # VECLIB
     MATH_LIBS="-framework veclib"
     LIBS="$LIBS $MATH_LIBS"
     HAVE_VECLIB=yes
     HAVE_MATH=yes
     ;;
     *)
     ;;
esac

# predefine math libraries depending on compiler
case "$FC" in
     *ifort*)
     MATH_LIBS="-lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_rt"
     LIBS="$LIBS $MATH_LIBS"
     HAVE_MKL=yes
     ;;
     *xlf*)
     LDFLAGS="-q64 $LDFLAGS"
     MATH_LIBS="-lmass -lmassv -lessl -llapack"
     LIBS="$LIBS $MATH_LIBS"
     HAVE_MASS=yes
     HAVE_ESSL=yes
     HAVE_LAPACK=yes
     AC_MSG_CHECKING([whether a location for the LAPACK library was specified])
     AC_ARG_WITH([lapack],[AS_HELP_STRING([--with-lapack=<directory>],
                    [Specify location the LAPACK library. 
                     Configure will expect to find subdirs include and lib.])],
		   [AC_MSG_RESULT([$withval])
		    LDFLAGS="$LDFLAGS -L$withval/lib"
    		    CPPFLAGS="$CPPFLAGS -I$withval/include"],
    		   [AC_MSG_RESULT([no])])

      ;;
      *)
      ;;
esac

# MKL : only checked if the Intel Compiler is used
if test "x$HAVE_MKL" = xyes; then
   AC_MSG_CHECKING([whether a location for the MKL library was specified])
   AC_ARG_WITH([mkl],[AS_HELP_STRING([--with-mkl=<directory>],
                    [Specify location the MKL library. 
                     Configure will expect to find subdirs include and lib.])],
		   [AC_MSG_RESULT([$withval])
		    LDFLAGS="$LDFLAGS -L$withval/lib/intel64"
    		CPPFLAGS="$CPPFLAGS -I$withval/include/intel64/lp64"],
    		[AC_MSG_RESULT([no])])

# This check does not work right now
   AC_CHECK_LIB([mkl_lapack95_lp64],[dgesv],,HAVE_MKL=no)
fi

if test "x$HAVE_MKL" = xyes; then
   HAVE_MATH=yes
fi

# ESSL
if test "x$HAVE_ESSL" = xyes; then
   AC_MSG_CHECKING([whether a location for the ESSL library was specified])
   AC_ARG_WITH([essl],[AS_HELP_STRING([--with-essl=<directory>],
                    [Specify location the ESSL library. 
                     Configure will expect to find subdirs include and lib.])],
		   [AC_MSG_RESULT([$withval])
		    LDFLAGS="$LDFLAGS -L$withval/lib"
    		    CPPFLAGS="$CPPFLAGS -I$withval/include"],
    		   [AC_MSG_RESULT([no])])

# This check does not work right now
   AC_CHECK_LIB([essl],[dgesv],,HAVE_ESSL=no)
fi

# MASS (needed in combination with ESSL on IBM AIX)
if test "x$HAVE_MASS" = xyes; then
   AC_MSG_CHECKING([whether a location for the MASS library was specified])
   AC_ARG_WITH([mass],[AS_HELP_STRING([--with-mass=<directory>],
                    [Specify location the MASS library. 
                     Configure will expect to find subdirs include and lib.])],
		   [AC_MSG_RESULT([$withval])
		    LDFLAGS="$LDFLAGS -L$withval/lib"
    		    CPPFLAGS="$CPPFLAGS -I$withval/include"],
    		   [AC_MSG_RESULT([no])])

# This check does not work right now
   AC_CHECK_LIB([mass],[tsinh],,HAVE_MASS=no)
fi

# ESSL and MASS found: still need LAPACK, but not BLAS
if test "x$HAVE_ESSL" = xyes && test "x$HAVE_MASS" = xyes && test "x$HAVE_LAPACK" = xyes; then
   HAVE_MATH=yes
fi

# All other systems: check for ACML, BLAS and LAPACK or ATLAS

# ACML
if test "x$HAVE_MATH" = xno; then
   HAVE_ACML=yes
   AC_MSG_CHECKING([whether a location for the ACML library was specified])
   AC_ARG_WITH([acml],[AS_HELP_STRING([--with-acml=<directory>],
                    [Specify location the ACML library. 
                     Configure will expect to find subdirs include and lib.])],
		   [AC_MSG_RESULT([$withval])
		    LDFLAGS="$LDFLAGS -L$withval/lib"
    		    CPPFLAGS="$CPPFLAGS -I$withval/include"],
    		   [AC_MSG_RESULT([no])])

   # check ACML linking in Fortran
   AC_LANG_PUSH(Fortran)
   AC_CHECK_LIB([acml],[dgesv],,HAVE_ACML=no)
   AC_LANG_POP(Fortran)
fi

if test "x$HAVE_ACML" = xyes; then
   HAVE_MATH=yes
   MATH_LIBS="-lacml"
fi

# BLAS
if test "x$HAVE_MATH" = xno; then
   HAVE_BLAS=yes
   AC_MSG_CHECKING([whether a location for the BLAS library was specified])
   AC_ARG_WITH([blas],[AS_HELP_STRING([--with-blas=<directory>],
                    [Specify location the BLAS library. 
                     Configure will expect to find subdirs include and lib.])],
		   [AC_MSG_RESULT([$withval])
		    LDFLAGS="$LDFLAGS -L$withval/lib"
    		    CPPFLAGS="$CPPFLAGS -I$withval/include"],
    		   [AC_MSG_RESULT([no])])

   # check BLAS linking in Fortran
   AC_LANG_PUSH(Fortran)
   AC_CHECK_LIB([blas],[dgemv],,HAVE_BLAS=no)
   AC_LANG_POP(Fortran)
fi

# LAPACK (depends on the presence of BLAS)
if test "x$HAVE_MATH" = xno; then
   HAVE_LAPACK=yes
   AC_MSG_CHECKING([whether a location for the LAPACK library was specified])
   AC_ARG_WITH([lapack],[AS_HELP_STRING([--with-lapack=<directory>],
                    [Specify location the LAPACK library. 
                     Configure will expect to find subdirs include and lib.])],
		   [AC_MSG_RESULT([$withval])
		    LDFLAGS="$LDFLAGS -L$withval/lib"
    		    CPPFLAGS="$CPPFLAGS -I$withval/include"],
    		   [AC_MSG_RESULT([no])])

   # check LAPACK linking in Fortran
   AC_LANG_PUSH(Fortran)
   AC_CHECK_LIB([lapack],[dlamch],,HAVE_LAPACK=no)
   AC_LANG_POP(Fortran)
fi

# last resort: just use BLAS and LAPACK
if test "x$HAVE_LAPACK" = xyes && test "x$HAVE_BLAS" = xyes; then
   HAVE_MATH=yes
   MATH_LIBS="-llapack -lblas"
fi

# ATLAS
if test "x$HAVE_MATH" = xno; then
   HAVE_ATLAS=yes
   AC_MSG_CHECKING([whether a location for the ATLAS library was specified])
   AC_ARG_WITH([atlas],[AS_HELP_STRING([--with-atlas=<directory>],
                    [Specify location the ATLAS library. 
                     Configure will expect to find subdirs include and lib.])],
		   [AC_MSG_RESULT([$withval])
		    LDFLAGS="$LDFLAGS -L$withval/lib"
    		    CPPFLAGS="$CPPFLAGS -I$withval/include"],
    		   [AC_MSG_RESULT([no])])
   # check ATLAS linking in Fortran
#   AC_LANG_PUSH(Fortran)
   AC_CHECK_LIB([atlas],[ATL_dscal],,HAVE_ATLAS=no)
   AC_CHECK_LIB([lapack],[dgemm],,HAVE_ATLAS=no)
#   AC_LANG_POP(Fortran)
fi

if test "x$HAVE_ATLAS" = xyes; then
   HAVE_MATH=yes
   MATH_LIBS="-latlas"
fi

# Define C compiler flags by compiler
case "$CC" in
     # GNU
     *gcc*)
     if test "x$ENABLE_DEBUG" = xyes; then
	CFLAGS="$CFLAGS -g"
     else
	CFLAGS="$CFLAGS -O3 -msse3 -funroll-loops -ftree-vectorize"
     fi
#     SHLIBFLAGS="-shared -undefined dynamic_lookup"
     SHLIBFLAGS="-shared"
     ;;
     # Portland Group
     *pgc*)
     if test "x$ENABLE_DEBUG" = xyes; then
	CFLAGS="$CFLAGS -g"
     else
	CFLAGS="$CFLAGS -O3 -Mcache_align -fastsse -Munroll=n:4 -Mipa=fast,inline"
     fi
     SHLIBFLAGS="-shared"
     ;;
     # PathScale
     *path*)
     if test "x$ENABLE_DEBUG" = xyes; then
	CFLAGS="$CFLAGS -O0 -C -g3 -trapuv -zerouv -fPIC -fno-fast-stdlib"
     else
	CFLAGS="$CFLAGS -Ofast"
     fi
     SHLIBFLAGS="-shared"
     ;;
     # INTEL
     *icc*)
     if test "x$ENABLE_DEBUG" = xyes; then
	CFLAGS="$CFLAGS -g"
     else
	CFLAGS="$CFLAGS -O3"
     fi
     SHLIBFLAGS="-shared"
     ;;
     # IBM
     *xlc*)
     if test "x$ENABLE_DEBUG" = xyes; then
	CFLAGS="$CFLAGS -g"
     else
	CFLAGS="$CFLAGS -O3"
     fi
     SHLIBFLAGS="-shared"
     ;;
     # Cray: how to test for the underlying compiler in the cc wrapper?
     *cc*)
     if test "x$ENABLE_DEBUG" = xyes; then
	CFLAGS="$CFLAGS -g"
     else
	CFLAGS="$CFLAGS -O3"
     fi
     SHLIBFLAGS="-shared"
     ;;
     *) AC_MSG_ERROR([C Compiler $CC currently not supported])
     ;;
esac

# Define Fortran compiler flags by compiler
case "$FC" in
     # GNU
     *gfortran*)
     if test "x$ENABLE_DEBUG" = xyes; then
     	FCFLAGS="$FCFLAGS -g -W -Wall -Wextra -pedantic-errors -fimplicit-none -fbacktrace -fbounds-check -ffpe-trap=invalid,zero,overflow"
     else
	FCFLAGS="$FCFLAGS -O3 -fimplicit-none -msse3 -funroll-loops -ftree-vectorize"
     fi
     ;;
     # Portland Group
     *pgf*)
     ;;
     # PathScale
     *path*)
     if test "x$ENABLE_DEBUG" = xyes; then
     	FCFLAGS="$FCFLAGS -O0 -C -g3 -trapuv -zerouv -fPIC -fno-fast-stdlib"
     else
	FCFLAGS="$FCFLAGS -Ofast"
     fi
     ;;
     # INTEL
     *ifort*)
     if test "x$ENABLE_DEBUG" = xyes; then
     	FCFLAGS="$FCFLAGS -g -check all -warn all,nodec,interfaces -gen_interfaces -traceback -fpe0"
     else
	FCFLAGS="$FCFLAGS -O3 -fno-alias -ip"
     fi
     ;;
     # IBM
     *xlf*)
     if test "x$ENABLE_DEBUG" = xyes; then
     	FCFLAGS="$FCFLAGS -g -qessl -qstrict -qcheck -C -qsigtrap -qflttrap=ov:zero:en -qsuffix=cpp=F90"
     else
	FCFLAGS="$FCFLAGS -O3 -qessl -qmaxmem=-1 -qarch=auto -qtune=auto -qipa -qstrict -qsuffix=cpp=F90"
     fi
     ;;
     # Cray: how to test for the underlying compiler in the ftn wrapper?
     *ftn*)
     # assume: gnu
     if test "x$ENABLE_DEBUG" = xyes; then
     	FCFLAGS="$FCFLAGS -g -W -Wall -Wextra -pedantic-errors -fimplicit-none -fbacktrace -fbounds-check -ffpe-trap=invalid,zero,overflow"
     else
	FCFLAGS="$FCFLAGS -O3 -fimplicit-none -msse3 -funroll-loops -ftree-vectorize"
     fi
     ;;
     *) AC_MSG_ERROR([Fortran Compiler $FC currently not supported])
     ;;
esac

# save Non-MPI compiler selection for compilation of serial code
CC_NOMPI=$CC
FC_NOMPI=$FC

# MPI
if test	"x$ENABLE_MPI" = xyes; then
   	HAVE_MPI=yes

   # find MPI C and F90 compilers
   if test "x$FC" = xftn; then
	# Cray
	CC_MPI="cc"
	FC_MPI="ftn"
   else
	# all others
	CCLIST="mpcc_r mpcc mpiicc mpicc"
	FCLIST="mpxlf90_r mpxlf90 mpxlf95_r mpxlf95 mpxlf_r mpxlf mpiifort mpif90"

	AC_CHECK_PROGS(CC_MPI,$CCLIST,none)
	AC_CHECK_PROGS(FC_MPI,$FCLIST,none)
   fi

   # assign the MPI wrapper to the C or F90 compiler
   if test "x$MPICC" = xnone; then
      HAVE_MPI=no
   else
      CC=$CC_MPI
   fi
   if test "x$MPIFC" = xnone; then
      HAVE_MPI=no
   else
      FC=$FC_MPI
   fi

   AC_CHECK_LIB([pthread],[pthread_create])

   AC_MSG_CHECKING([whether a location for the MPI library was specified])
   AC_ARG_WITH([mpi],[AS_HELP_STRING([--with-mpi=<directory>],
               [Specify location the mpi library. 
                Configure will expect to find subdirs bin, include and lib.])],
		[AC_MSG_RESULT([$withval])
		 MPIPATH="$withval/bin"
		 LDFLAGS="$LDFLAGS -L$withval/lib"
#    		 CPPFLAGS="$CPPFLAGS -I$withval/include"
#		 FCFLAGS="$FCFLAGS -I$withval/include"
		],
    		[AC_MSG_RESULT([no])])

   AC_CHECK_LIB([mpi_r],[MPI_Init],,AC_CHECK_LIB([mpi],[MPI_Init],,AC_CHECK_LIB([mpich],[MPI_Init],,HAVE_MPI=no)))
   # Change to the Fortran 90 language
#   AC_LANG_PUSH(Fortran)

# Try link a simple MPI program.  If fail, try again with -lmpi.
#   AC_LINK_IFELSE([
#program main
#include 'mpif.h'
#integer :: ierr
#call MPI_INIT( ierr )
#end],,HAVE_MPI=no)

# Change to the C language
#   AC_LANG_POP(Fortran)

else
   HAVE_MPI=no
fi

# disable MPI statements in code during compile time
# if no MPI compiler environment was found
if test "x$HAVE_MPI" = xno; then
   CC=$CC_NOMPI
   FC=$FC_NOMPI
   CPPFLAGS="$CPPFLAGS -DHIDE_MPI"
fi

# define Automake Conditionals
#AM_CONDITIONAL([DEBUG], [test "x$ENABLE_DEBUG" = xyes])

# Check for HDF5 library
# We have to do that at the very end since the HDF5 library might have
# dependencies on the MPI compiler
#HAVE_HDF5=no
#HAVE_HDF5=yes
#AC_MSG_CHECKING([whether a location for the HDF5 library was specified])
#AC_ARG_WITH([hdf5],[AS_HELP_STRING([--with-hdf5=<directory>],
#                    [Specify location the HDF5 library. 
#                     Configure will expect to find subdirs include and lib.])],
#		   [AC_MSG_RESULT([$withval])
#		    LDFLAGS="$LDFLAGS -L$withval/lib"
#    		    CPPFLAGS="$CPPFLAGS -I$withval/include"],
#    		   [AC_MSG_RESULT([no])])

#AC_CHECK_HEADERS([hdf5.h],,HAVE_HDF5=no)
#AC_CHECK_LIB([hdf5],[H5Fopen],,HAVE_HDF5=no)
#AC_CHECK_LIB([hdf5_hl],[H5DSis_scale],,HAVE_HDF5=no)

# NetCDF
# We have to do that at the very end since the NetCDF library might have
# dependencies on the MPI compiler and on the HDF5 library
AC_MSG_CHECKING([whether a location for the NetCDF library was specified])
AC_ARG_WITH([netcdf],[AS_HELP_STRING([--with-netcdf=<directory>],
                    [Specify location the NetCDF library. 
                     Configure will expect to find subdirs include and lib.])],
		   [AC_MSG_RESULT([$withval])
		    LDFLAGS="$LDFLAGS -L$withval/lib"
    		    CPPFLAGS="$CPPFLAGS -I$withval/include"
    		    NC_BINPATH="$withval/bin"],
    		   [AC_MSG_RESULT([no])])

AC_ARG_WITH([netcdf-include],[AS_HELP_STRING([--with-netcdf-include=<directory>],
                    [Specify location the NetCDF include directory. 
                     Configure will expect to find netcdf.h and netcdf.mod there.])],
		   [AC_MSG_RESULT([$withval])
    		    CPPFLAGS="$CPPFLAGS -I$withval"],
    		   [AC_MSG_RESULT([no])])

AC_ARG_WITH([netcdf-lib],[AS_HELP_STRING([--with-netcdf-lib=<directory>],
                    [Specify location the NetCDF library itself. 
                     Configure will expect to find libnetcdf.a there.])],
		   [AC_MSG_RESULT([$withval])
		    LDFLAGS="$LDFLAGS -L$withval"]
    		   [AC_MSG_RESULT([no])])

#AC_ARG_WITH([netcdf-bin],[AS_HELP_STRING([--with-netcdf-bin=<directory>],
#                    [Specify location the NetCDF binaries itself. 
#                     Configure will expect to find ncdump or nc-config there.])],
#		   [AC_MSG_RESULT([$withval])
#		    NC_BINPATH="$withval"]
#    		   [AC_MSG_RESULT([no])])

#AC_CHECK_PROG( HAVE_NC_CONFIG, nc-config, [yes], [no], $NC_BINPATH )
#if test x$HAVE_NC_CONFIG = xyes; then
#    HAVE_NETCDF4=`$NC_BINPATH/nc-config --has-nc4`
#fi

# NetCDF: C library
HAVE_NETCDF=yes
AC_CHECK_HEADER(netcdf.h,, HAVE_NETCDF=no)
if test "x$HAVE_NETCDF" = xno; then
   AC_MSG_ERROR([Can't find NetCDF header file.])
fi

AC_CHECK_LIB([netcdf],[nc_open],,HAVE_NETCDF=no)
if test "x$HAVE_NETCDF" = xno; then
   AC_MSG_ERROR([Can't find NetCDF library, or can't link a program to it.])
fi

# check for separate fortran NetCDF library.
# We have to do that at the very end since any other
# library check above involving the C compiler could
# have problem with the fortran version of this library here
HAVE_NETCDF_F90=yes
AC_LANG_PUSH(Fortran)
AC_CHECK_LIB([netcdf],[nf_close],,HAVE_NETCDF_F90=no)
if test "x$HAVE_NETCDF_F90" = xno; then
   HAVE_NETCDF_F90=yes
   AC_CHECK_LIB([netcdff],[nf_close],,HAVE_NETCDF_F90=no)
   if test "x$HAVE_NETCDF_F90" = xno; then
      AC_MSG_ERROR([Can't find NetCDF Fortran 90 library, or can't link a program to it.])
   fi    
fi
AC_LANG_POP(Fortran)

# NetCDF-4
AC_CHECK_LIB([netcdf],[NC4_open],HAVE_NETCDF4=yes,HAVE_NETCDF4=no)

# NetCDF-4 with MPI
AC_CHECK_LIB([netcdf],[nc_use_parallel_enabled],HAVE_NETCDF4_MPI=yes,HAVE_NETCDF4_MPI=no)

# check for HDF5 and NetCDF-4
#if test "x$HAVE_HDF5" = xno && test "x$HAVE_NETCDF4" = xyes ; then
#    AC_MSG_ERROR([HDF5 is needed when NetCDF with NetCDF-4 capability is used])
#fi 

# check for curl when NetCDF-4 is selected
#if test "x$HAVE_NETCDF4" = xyes ; then
#    AC_CHECK_LIB([curl],[curl_version],,AC_MSG_NOTICE(curl library libcurl is not present, but may be needed for NetCDF4))
#fi

# substitution
AC_SUBST([CC_NOMPI])
AC_SUBST([FC_NOMPI])
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([SHLIBFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([MATH_LIBS])

# set installation prefix to current directory if none was given
if test "$prefix" = "NONE"; then
    prefix=`pwd`;
fi

# Define System in config.h
if test -z "$HOST" ; then HOST="unknown"; fi;
AC_DEFINE_UNQUOTED(CC_NOMPI, ["$CC_NOMPI"], [C Compiler without MPI Support])
AC_DEFINE_UNQUOTED(FC_NOMPI, ["$FC_NOMPI"], [Fortran 90/95 Compiler without MPI Support])
AC_DEFINE_UNQUOTED(CC, ["$CC"], [C Compiler])
AC_DEFINE_UNQUOTED(FC, ["$FC"], [Fortran 90/95 Compiler])
AC_DEFINE_UNQUOTED(OS, ["$host_os"], [Operating System Type])
AC_DEFINE_UNQUOTED(USER_NAME, ["$USERNAME"], [User name])
AC_DEFINE_UNQUOTED(HOST_NAME, ["$HOST"], [Host name])

if test "x$HAVE_MPI" = xyes; then
    AC_DEFINE(HAVE_MPI,1,Define 1 to run with the MPI parallelization)
else
    AC_DEFINE(HAVE_MPI,0,Define 1 to run with the MPI parallelization)
fi

if test "x$HAVE_NETCDF4" = xyes; then
    AC_DEFINE(HAVE_NETCDF4,1,Define 1 to run with NetCDF-4)
else
    AC_DEFINE(HAVE_NETCDF4,0,Define 1 to run with NetCDF-4)
    CPPFLAGS="$CPPFLAGS -DHIDE_NETCDF4"
fi

if test "x$HAVE_NETCDF4_MPI" = xyes; then
    AC_DEFINE(HAVE_NETCDF4_MPI,1,Define 1 to run with NetCDF-4 and MPI)
else
    AC_DEFINE(HAVE_NETCDF4_MPI,0,Define 1 to run with NetCDF-4 and MPI)
    CPPFLAGS="$CPPFLAGS -DHIDE_NETCDF4_MPI"
fi

AC_CONFIG_FILES([Makefile
		 src/Makefile
		 config/Makefile
                 scripts/Makefile
                 analysis/idl_lib/Makefile
                 analysis/idl/Makefile])
AC_OUTPUT

echo "

Pheno Analysis Configuration Summary:

Installation location:	${prefix}

Compiler:
MPI support: 		${HAVE_MPI}
C Compiler (No MPI):	${CC_NOMPI}
F90 Compiler (No MPI):	${FC_NOMPI}
C Compiler:    		${CC}
F90 Compiler:		${FC}
C Compiler Flags:	${CFLAGS}
F90 Compiler Flags:	${FCFLAGS}
Preprocessor Flags:	${CPPFLAGS}
Include Directives:	${INCLUDES}
Shared Library Flags:	${SHLIBFLAGS}
Linker Flags:		${LDFLAGS}
Libraries:		${LIBS}
Math Libraries:		${MATH_LIBS}

Data I/O Libraries:
NetCDF C:		${HAVE_NETCDF}
NetCDF Fortran 90:	${HAVE_NETCDF_F90}
NetCDF4:		${HAVE_NETCDF4}
NetCDF4-MPI:            ${HAVE_NETCDF4_MPI}
JPEG:			${HAVE_JPEG}
ZLIB:			${HAVE_ZLIB}
SZIP:			${HAVE_SZLIB}
HDF4:			${HAVE_HDF4}

Math Libraries:
VECLIB:			${HAVE_VECLIB}
LAPACK:	       		${HAVE_LAPACK}
BLAS:			${HAVE_BLAS}
ATLAS:			${HAVE_ATLAS}
ACML:			${HAVE_ACML}
MKL:			${HAVE_MKL}
ESSL:			${HAVE_ESSL}
MASS:			${HAVE_MASS}

Job Submission:
PBS:			${HAVE_PBS}
POE:			${HAVE_POE}
SLURM:			${HAVE_SLURM}
"