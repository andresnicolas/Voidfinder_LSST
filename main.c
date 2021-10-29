
#include "global.h"
#include "io.h"
#include "map.h"
#include "finder.h"

int main()
{
 	
  char filename[200];
  vector <tracer> tr;
  vector <tracer> ran;
  vector <voids> v;
  struct hpmap *map;

  int ncores = 8;
  omp_set_num_threads(ncores);

  sprintf(filename,"data/halos_shell.dat");
  read_tracers(filename,tr);

  sprintf(filename,"data/random_shell.dat");
  read_randoms(filename,ran);

  int TracPerPix = 10;                  // Tracers por pixel
  float SurveyArea = 4.0 * M_PI / 8.0;  // Area (angular) del survey: octante de esfera
  float ShellDist = 950.0;              // Distancia media de la cascara [Mpc/h] 
  float ShellThick = 100.0;             // Ancho de la cascara [Mpc/h]
  float delta_cut = -0.6;               // Delta de corte
  float tol = 0.0;                      // Tolerancia superposicion
  float rmax = 10.0 * DEG2RAD;          // Radio angular maximo  

  int order = round(0.5*log(M_PI*(float)tr.size()/3.0/SurveyArea/TracPerPix)/log(2.0));  

  T_Healpix_Base<int> healpix(order,RING);
  map = (struct hpmap *) malloc(healpix.Npix()*sizeof(struct hpmap));
 
  map_load_tracers(tr,healpix,map);
  map_load_randoms(ran,healpix,map);
  map_load_mask(healpix,map);
  map_compute_delta(tr,ran,healpix.Npix(),map);
  find_candidates(delta_cut,rmax,healpix,map,tr,ran,v);
  find_voids(delta_cut,healpix,map,tr,ran,v);
  map_load_voids(v,healpix,map);
  clean_voids(tol,healpix,map,v);

  sprintf(filename,"data/voids_2D.dat");
  write_voids(filename,v);

  tr.clear();
  ran.clear();
  v.clear();
  for (int ipix=0; ipix<healpix.Npix(); ipix++) {
     free(map[ipix].tracer);
     free(map[ipix].random);
  }
  free(map);
}
