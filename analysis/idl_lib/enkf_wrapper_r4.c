#include <stdio.h>  
   
int enkf_wrapper_r4(int argc, void *argv[])  
{  
  extern void enkf_driver_r4_();	/* Fortran routine */  
  int ret;
        
  //  enkf_driver_(ndim, nrens, nrobs);
  enkf_driver_r4_((int *) argv[0],(int *) argv[1],(int *) argv[2],
	       (float *) argv[3], (float *) argv[4], (float *) argv[5], 
	       (float *) argv[6]);

  return ret;
}  
