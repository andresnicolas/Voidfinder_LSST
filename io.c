
#include "global.h"
#include "io.h"

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

void read_tracers(char *filename, vector <tracer> &tr)
{
  	
  FILE *f = safe_open(filename,"r");	
  int i = 0;
  float ra,dec,redshift,distance,x,y,z;

  while (fscanf(f,"%f %f %f %f %f %f %f",&ra,&dec,&redshift,&distance,&x,&y,&z) > 0) {
     tr.push_back(tracer());
     tr[i].coord.phi = ra*DEG2RAD;
     tr[i].coord.theta = (90.0 - dec)*DEG2RAD;
     tr[i].redshift = redshift;
     tr[i].weight = 1.0;
     i++;
  };

  fclose(f);

}

void read_randoms(char *filename, vector <tracer> &ran)
{
 
  FILE *f = safe_open(filename,"r");	
  int i = 0;
  float ra,dec,redshift,distance,x,y,z;

  while (fscanf(f,"%f %f %f %f %f %f %f",&ra,&dec,&redshift,&distance,&x,&y,&z) > 0) {
     ran.push_back(tracer());
     ran[i].coord.phi = ra*DEG2RAD;
     ran[i].coord.theta = (90.0 - dec)*DEG2RAD;
     ran[i].redshift = redshift;
     ran[i].weight = 1.0;
     i++;
  };

  fclose(f);
 
}

