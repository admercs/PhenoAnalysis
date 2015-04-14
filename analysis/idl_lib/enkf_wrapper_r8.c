#include <stdio.h>  
   
int enkf_wrapper_r8(int argc, void *argv[])  
{  
  extern void enkf_driver_r8_();	/* Fortran routine */  
  int ret;
        
  //  enkf_driver_(ndim, nrens, nrobs);
  enkf_driver_r8_((int *) argv[0],(int *) argv[1],(int *) argv[2],
	       (double *) argv[3], (double *) argv[4], (double *) argv[5], 
	       (double *) argv[6]);

  return ret;
}  
