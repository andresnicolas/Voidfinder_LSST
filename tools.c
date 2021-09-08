
#include "global.h"

FILE *safe_open(char *filename, const char *mode)
{
  FILE *f;

  if (!(f = fopen(filename,mode))) {
     fprintf(stdout,"\nError. File '%s' not found. \n\n",filename);	  
     fflush(stdout);
     exit(EXIT_FAILURE);
  }   

  return f; 
}

float random_number()
{
   float x = (float) rand()/RAND_MAX;
   return x;
}
