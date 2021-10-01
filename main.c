
#include "global.h"
#include "proto.h"

struct global_data G;

int main()
{
 	
  char filename[200];
  vector <tracer> tr;
  vector <tracer> ran;
  vector <voids> v;
  struct hpmap *map;

  sprintf(filename,"data/halos_shell.dat");
  read_tracers(filename,tr);

  sprintf(filename,"data/random_shell.dat");
  read_randoms(filename,ran);

  G.SurveyArea = 4.0 * M_PI / 8.0; // Area (angular) del survey: octante de esfera
  G.TracPerPix = 10;               // Tracers por pixel
  G.ShellDist = 950.0;             // Distancia media de la cascara [Mpc/h] 
  G.ShellThick = 100.0;            // Ancho de la cascara [Mpc/h]

  int order = round(0.5*log(M_PI*(float)tr.size()/3.0/G.SurveyArea/G.TracPerPix)/log(2.0));  
  T_Healpix_Base<int> healpix(order,RING);
  map = (struct hpmap *) malloc(healpix.Npix()*sizeof(struct hpmap));
  create_map(tr,ran,healpix,map);

  float delta_seed = -0.3;
  float delta_cut = -0.4;
  float tol = 0.0;
  find_centers(delta_seed,healpix,map,v);
  find_voids(delta_cut,healpix,map,tr,ran,v);
  clean_voids(tol,healpix,map,v);

  FILE *f1 = fopen("all.dat","w");
  FILE *f2 = fopen("clean.dat","w");
  for (int iv=0; iv<v.size(); iv++) {
      if (v[iv].radius > 0.0) fprintf(f1,"%f %f %f \n",v[iv].radius*G.ShellDist,v[iv].coord_init.phi*RAD2DEG,90.0-v[iv].coord_init.theta*RAD2DEG); 	  
      if (v[iv].tof) fprintf(f2,"%f %f %f \n",v[iv].radius*G.ShellDist,v[iv].coord_init.phi*RAD2DEG,90.0-v[iv].coord_init.theta*RAD2DEG);  
  }
  fclose(f1);
  fclose(f2);
   
 // FILE *f = fopen("mapa.dat","w");
 // pointing ptg; 
 // for (int ipix=0; ipix<healpix.Npix(); ipix++) {
 //     ptg = healpix.pix2ang(ipix);
 //     fprintf(f,"%f %f %f %f\n",ptg.phi,ptg.theta,map[ipix].delta,map[ipix].delta_smooth);      
 // }
 // fclose(f);

  tr.clear();
  ran.clear();
  for (int ipix=0; ipix<healpix.Npix(); ipix++) {
     free(map[ipix].tracer);
     free(map[ipix].random);
  }
  free(map);
}
