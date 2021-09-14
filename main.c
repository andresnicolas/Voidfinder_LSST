
#include "global.h"
#include "proto.h"

struct global_data G;

int main()
{
  	
  char filename[200];
  vector <tracer> tr;
  vector <tracer> ran;
  struct hpmap *map;

  sprintf(filename,"data/halos_shell.dat");
  read_tracers(filename,tr);

  sprintf(filename,"data/random_shell.dat");
  read_randoms(filename,ran);

  //float Reff = 2.0; // Pixel de R = 2.0 Mpc/h
  //float o = log(G.ShellDist / Reff / sqrt(3.0)) / log(2.0); 
  //int order = round(o);

  G.SurveyArea = 4.0 * M_PI / 8.0; // Area (angular) del survey: octante de esfera
  G.TracPerPix = 10;               // Tracers por pixel
  G.ShellDist = 950.0;             // Distancia media de la cascara [Mpc/h] 
  G.ShellThick = 100.0;            // Ancho de la cascara [Mpc/h]

  int order = round(0.5 * log(M_PI * (float) tr.size() / 3.0 / G.SurveyArea / G.TracPerPix) / log(2.0));  	  
  T_Healpix_Base<int> healpix(order,RING);
  map = (struct hpmap *) malloc(healpix.Npix()*sizeof(struct hpmap));
  create_map(tr,ran,healpix,map);

  FILE *f = fopen("mapa.dat","w");
  pointing ptg; 
  for (int ipix=0; ipix<healpix.Npix(); ipix++) {
      ptg = healpix.pix2ang(ipix);
      fprintf(f,"%f %f %f %f \n",ptg.phi,ptg.theta,map[ipix].delta,map[ipix].delta_smooth);      
  }
  fclose(f);

  tr.clear();
  ran.clear();
  for (int ipix=0; ipix<healpix.Npix(); ipix++) {
     map[ipix].tracer.clear();
     map[ipix].random.clear();
  }
  free(map);
}
