
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
  float ra,dec,w;

  while (fscanf(f,"%f %f %f ",&ra,&dec,&w) > 0) {
     tr.push_back(tracer());
     tr[i].coord.phi = ra*DEG2RAD;
     tr[i].coord.theta = (90.0 - dec)*DEG2RAD;
     tr[i].weight = w;
     i++;
  };

  fclose(f);

}

void read_randoms(char *filename, vector <tracer> &ran)
{
 
  FILE *f = safe_open(filename,"r");	
  int i = 0;
  float ra,dec,w;

  while (fscanf(f,"%f %f %f",&ra,&dec,&w) > 0) {
     ran.push_back(tracer());
     ran[i].coord.phi = ra*DEG2RAD;
     ran[i].coord.theta = (90.0 - dec)*DEG2RAD;
     ran[i].weight = w;
     i++;
  };

  fclose(f);
 
}

void write_voids(char *filename, vector <voids> &v)
{

  FILE *f = safe_open(filename,"w");

  for (int i=0; i<v.size(); i++) {
      
      if (!v[i].tof) continue;

      float ra = v[i].coord.phi*RAD2DEG;
      float dec = 90.0 - v[i].coord.theta*RAD2DEG;  

      fprintf(f,"%8.6f %12.6f %12.6f %8.6f \n",v[i].radius,ra,dec,v[i].delta);

  }

}

