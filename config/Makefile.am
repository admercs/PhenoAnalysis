bin_PROGRAMS = phenoanalysis-config

phenoanalysis_config_SOURCES = \
../src/read_config.c \
phenoanalysis-config.c

# do not compile this code with MPI support since
# its shell execution happens in serial mode on
# front-end (login) nodes on HPC facilities
CC=@CC_NOMPI@
LIBS=
