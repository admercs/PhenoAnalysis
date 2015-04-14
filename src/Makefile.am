bin_PROGRAMS = phenoanalysis quicksort_test reproject_test

reprojectsources = \
	reproject_module.F90 \
	proj_module.F90 \
	geosat_module.F90 \
	mfg_mviri_module.F90 \
	msg_seviri_module.F90 \
	modsin_module.F90 \
	rotpol_module.F90 \
	quicksort_module.F90 \
	bilinear_module.F90 \
	histogram_module.F90

phenoanalysis_SOURCES = \
	glob.c \
	read_config.c \
	parameter_module.F90 \
	file_module.F90 \
	driver_module.F90 \
	photoperiod_module.F90 \
	meteorology_average_module.F90 \
	phenology_model_module.F90 \
	time_module.F90 \
	phenology_module.F90 \
	topo_module.F90 \
	pft_module.F90 \
	restart_module.F90 \
	enkf_module_r8.F90 \
	hdf_module.F90 \
	modis_module.F90 \
	alma_module.F90 \
	ecmwf_module.F90 \
	cosmo_module.F90 \
	meteorology_module.F90 \
	netcdf_module.F90 \
	output_module.F90 \
	assimilation_module.F90 \
	init_module.F90 \
	exit_module.F90 \
	phenoanalysis.F90 \
	$(reprojectsources)

quicksort_test_SOURCES = quicksort_module.F90 quicksort_test.F90 

reproject_test_SOURCES = reproject_test.F90 $(reprojectsources)

noinst_HEADERS = read_config.h

CLEANFILES = *.mod

# hardcode dependencies (gnu automake cannot yet do F90 module interdependencies)
# use the script makedepend.sh (which calls makedepend.pl) in the scripts directory to 
# re-generate the dependencies below after you have possibly changed module use statements
# in the code or after you have added new modules to the code

# DO NOT EDIT BELOW THIS LINE
alma_module.o : alma_module.F90 parameter_module.o proj_module.o exit_module.o time_module.o reproject_module.o 
assimilation_module.o : assimilation_module.F90 file_module.o pft_module.o driver_module.o topo_module.o restart_module.o proj_module.o modis_module.o exit_module.o time_module.o enkf_module_r8.o 
bilinear_module.o : bilinear_module.F90  
cosmo_module.o : cosmo_module.F90 parameter_module.o proj_module.o exit_module.o time_module.o reproject_module.o 
driver_module.o : driver_module.F90 parameter_module.o 
ecmwf_module.o : ecmwf_module.F90 parameter_module.o proj_module.o exit_module.o time_module.o reproject_module.o 
enkf_module_r4.o : enkf_module_r4.F90  
enkf_module_r8.o : enkf_module_r8.F90  
exit_module.o : exit_module.F90 file_module.o 
file_module.o : file_module.F90  
geosat_module.o : geosat_module.F90 mfg_mviri_module.o msg_seviri_module.o 
hdf_module.o : hdf_module.F90  
histogram_module.o : histogram_module.F90 quicksort_module.o 
init_module.o : init_module.F90 file_module.o parameter_module.o meteorology_module.o pft_module.o driver_module.o assimilation_module.o phenology_module.o topo_module.o quicksort_module.o proj_module.o exit_module.o modis_module.o time_module.o 
meteorology_average_module.o : meteorology_average_module.F90  
meteorology_module.o : meteorology_module.F90 parameter_module.o cosmo_module.o driver_module.o phenology_module.o ecmwf_module.o exit_module.o meteorology_average_module.o time_module.o alma_module.o 
mfg_mviri_module.o : mfg_mviri_module.F90  
modis_module.o : modis_module.F90 file_module.o modsin_module.o parameter_module.o proj_module.o hdf_module.o 
modsin_module.o : modsin_module.F90  
msg_seviri_module.o : msg_seviri_module.F90  
netcdf_module.o : netcdf_module.F90 proj_module.o 
output_module.o : output_module.F90 netcdf_module.o driver_module.o exit_module.o time_module.o 
parameter_module.o : parameter_module.F90  
pft_module.o : pft_module.F90 file_module.o proj_module.o exit_module.o reproject_module.o 
phenoanalysis.o : phenoanalysis.F90 meteorology_module.o parameter_module.o init_module.o driver_module.o assimilation_module.o phenology_module.o restart_module.o exit_module.o time_module.o output_module.o 
phenology_model_module.o : phenology_model_module.F90  
phenology_module.o : phenology_module.F90 parameter_module.o driver_module.o phenology_model_module.o meteorology_average_module.o time_module.o 
photoperiod_module.o : photoperiod_module.F90  
proj_module.o : proj_module.F90 modsin_module.o rotpol_module.o geosat_module.o 
quicksort_module.o : quicksort_module.F90  
quicksort_test.o : quicksort_test.F90 quicksort_module.o 
reproject_module.o : reproject_module.F90 histogram_module.o proj_module.o bilinear_module.o 
reproject_test.o : reproject_test.F90 modsin_module.o histogram_module.o rotpol_module.o quicksort_module.o proj_module.o reproject_module.o 
restart_module.o : restart_module.F90 file_module.o pft_module.o driver_module.o exit_module.o 
rotpol_module.o : rotpol_module.F90  
time_module.o : time_module.F90  
topo_module.o : topo_module.F90 file_module.o proj_module.o exit_module.o reproject_module.o 
