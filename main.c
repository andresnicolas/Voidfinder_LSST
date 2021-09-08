
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
  create_randoms(ran);

  int order = round((log(G.ShellDist)-0.5*log(3.0))/log(2.0)-1.0);
  T_Healpix_Base<int> healpix(order,RING);
  G.PixelArea = (4.0*M_PI*G.ShellDist*G.ShellDist) / (float)healpix.Npix();
  
  map = (struct hpmap *) malloc(healpix.Npix()*sizeof(struct hpmap));
  create_map(tr,ran,healpix,map);

  FILE *f = fopen("mapa.dat","w");
  pointing ptg; 
  for (int ipix=0; ipix<healpix.Npix(); ipix++) {
      ptg = healpix.pix2ang(ipix);
      fprintf(f,"%f %f %f \n",ptg.phi,ptg.theta,map[ipix].delta);      
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
