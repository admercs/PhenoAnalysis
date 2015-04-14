# IDL scripts, compiled code
bin_SCRIPTS = maketopo.sav makepft.sav

# IDL scripts not to be installed but distributed
EXTRA_DIST = *.pro *.idl *.tbl *.txt

# IDL Compiler
IDL = idl

# IDL Compiler Flags
IDLFLAGS =

# Define new suffixes
SUFFIXES: .idl .sav

# Define code to be compiled by default
all: $(bin_SCRIPTS)

# Compile statement
.idl.sav:
	dir=`(cd $(srcdir); pwd)`; \
	IDL_PATH=$${IDL_PATH}:$${dir}; \
	echo $${IDL_PATH}; \
	$(IDL) $(IDLFLAGS) $<

# Define cleaning rules
CLEANFILES = $(bin_SCRIPTS)

# Hardcode dependencies since GNU automake cannot handle IDL dependencies automatically
maketopo.sav : maketopo.pro file_lock.pro read_cgiarsrtm.pro read_gtopo30.pro rebin_new.pro
makepft.sav : makepft.pro file_lock.pro read_avhrrvcf.pro read_cropdata.pro read_modistile.pro fill_dist.pro rebin_new.pro ll2rpol.pro rpol2ll.pro read_ecmwf.pro
