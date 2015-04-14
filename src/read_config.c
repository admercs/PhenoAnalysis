/* Copyright 2013 EUMETSAT and MeteoSwiss
   Authors: Reto Stockli (MeteoSwiss) 
   
   This file is part of <INSERT_PACKAGE_NAME> and was created within the
   EUMETSAT Satellite Application Facility on Climate Monitoring (CM SAF)
   project by MeteoSwiss.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"

int read_config (int *length, char *string) {

  int retval = 0;

  if (!strncmp(string,"PACKAGE_VERSION",15)) {
    strncpy(string,PACKAGE_VERSION,*length);
  } else if (!strncmp(string,"PACKAGE_NAME",12)) {
    strncpy(string,PACKAGE_NAME,*length);
  } else if (!strncmp(string,"HAVE_PBS",8)) {
    memset(string,'\0',*length);
    snprintf(string,*length,"%d",HAVE_PBS);
  } else if (!strncmp(string,"HAVE_POE",8)) {
    memset(string,'\0',*length);
    snprintf(string,*length,"%d",HAVE_POE);
  } else if (!strncmp(string,"HAVE_SLURM",10)) {
    memset(string,'\0',*length);
    snprintf(string,*length,"%d",HAVE_SLURM);
  } else if (!strncmp(string,"ENABLE_DEBUG",12)) {
    memset(string,'\0',*length);
    snprintf(string,*length,"%d",ENABLE_DEBUG);
  } else if (!strncmp(string,"HAVE_MPI",8)) {
    memset(string,'\0',*length);
    snprintf(string,*length,"%d",HAVE_MPI);
  } else {
    strncpy(string,"N/A",*length);
    retval = -1;
  }

  return retval;
}

int read_config_wrapper(int argc, void *argv[])  
{  
  return read_config((int *) argv[0], (char *) argv[1]); 
}  
