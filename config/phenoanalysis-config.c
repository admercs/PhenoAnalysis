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
#include <getopt.h>

#include "../src/read_config.h"

/*
int main (int argc,char *argv[])
{

  int count;

  printf ("This program was called with \"%s\".\n",argv[0]);

  if (argc > 1)
    {
      for (count = 1; count < argc; count++)
	{
	  printf("argv[%d] = %s\n", count, argv[count]);
	}
    }
  else
    {
      printf("The command had no other arguments.\n");
    }

  return 0;

}
*/
     
int main (int argc, char **argv)
{
  int c;
  int retval;
  int length = 20;
     
  static int help_flag;
  char string[length];
  char *str = string;
     
  memset(str,'\0',length);

  if (argc == 1) {
    help_flag = 1;
  } else {
    help_flag = 0;
  }

  while (1)
    {
      static struct option long_options[] =
	{
	  /* These options set a flag. */
	  {"help",         no_argument, &help_flag, 1},
	  /* These options don't set a flag.
	     We distinguish them by their indices. */
	  {"all",          no_argument, 0, 'a'},
	  {"version",      no_argument, 0, 'b'},
	  {"name",         no_argument, 0, 'c'},
	  {"has-pbs",      no_argument, 0, 'd'},
	  {"has-poe",      no_argument, 0, 'e'},
	  {"has-slurm",    no_argument, 0, 'f'},
	  {"enable-debug", no_argument, 0, 'g'},
	  {"has-mpi",      no_argument, 0, 'h'},
	  {0, 0, 0, 0}
	};
      /* getopt_long stores the option index here. */
      int option_index = 0;
     
      c = getopt_long (argc, argv, "abcdefgh",
		       long_options, &option_index);
     
      /* Detect the end of the options. */
      if (c == -1)
	break;
     
      switch (c)
	{    
	case 'a':
	  printf ("option `--all' is not yet implemented \n");
	  break;
     
	case 'b':
	  strncpy(str,"PACKAGE_VERSION",length);
	  retval = read_config(&length,str);
	  printf ("%s \n",string);
	  break;
     
	case 'c':
	  strncpy(str,"PACKAGE_NAME",length);
	  retval = read_config(&length,str);
	  printf ("%s \n",string);
	  break;
     
	case 'd':
	  strncpy(str,"HAVE_PBS",length);
	  retval = read_config(&length,str);
	  printf ("%s \n",string);
	  break;
     
	case 'e':
	  strncpy(str,"HAVE_POE",length);
	  retval = read_config(&length,str);
	  printf ("%s \n",string);
	  break;
     
	case 'f':
	  strncpy(str,"HAVE_SLURM",length);
	  retval = read_config(&length,str);
	  printf ("%s \n",string);
	  break;
    
	case 'g':
	  strncpy(str,"ENABLE_DEBUG",length);
	  retval = read_config(&length,str);
	  printf ("%s \n",string);
	  break;
    
	case 'h':
	  strncpy(str,"HAVE_MPI",length);
	  retval = read_config(&length,str);
	  printf ("%s \n",string);
	  break;
    
	default:
	  /* getopt_long already prints error message for unrecognized option */
 	  help_flag = 1;
	  break;
	}
    }
     
  /* Instead of reporting ‘--help’
     and ‘--brief’ as they are encountered,
     we report the final status resulting from them. */
  if (help_flag) {
    printf ("Usage: %s [OPTION]\n\n",argv[0]);
    printf ("Available values for OPTION include:\n\n");
    printf ("  --help         display this help message and exit \n");
    printf ("  --all          display all options (INOP)\n");
    printf ("  --version      display package version number\n");
    printf ("  --name         display package name\n");
    printf ("  --has-pbs      whether package was compiled with PBS support\n");
    printf ("  --has-poe      whether package was compiled with POE support\n");
    printf ("  --has-slurm    whether package was compiled with SLURM support\n");
    printf ("  --enable-debug whether package was compiled with debugging support\n");
    printf ("  --has-mpi      whether package was compiled with mpi support\n");
    printf ("\n");
 
  }
    
  /* Print any remaining command line arguments (not options). */
  if (optind < argc)
    {
      while (optind < argc)
	printf ("Invalid non-option arguments provided: %s \n", argv[optind++]);
    }
     
  exit (0);
}


