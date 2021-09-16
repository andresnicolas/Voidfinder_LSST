
#include "global.h"
#include "proto.h"

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

void create_map(vector <tracer> &tr, vector <tracer> &ran, T_Healpix_Base<int> &hp, struct hpmap *map) 
{

  // Mascara angular
  for (int ipix=0; ipix<hp.Npix(); ipix++) {
      pointing p = hp.pix2ang(ipix);
      if (p.theta <= 0.5*M_PI && p.phi <= 0.5*M_PI)	
	 map[ipix].mask = true;
      else
	 map[ipix].mask = false;     
  }

  // Cargo tracers en mapa angular	
  for (int i=0; i<tr.size(); i++) {
     int ipix = hp.ang2pix(tr[i].coord);
     map[ipix].tracer.push_back(i);   
  }

  // Cargo randoms en mapa angular
  for (int i=0; i<ran.size(); i++) {
     int ipix = hp.ang2pix(ran[i].coord);
     map[ipix].random.push_back(i);   
  }

  // Calculo delta

  float wtrac = 0.0;
  float wrand = 0.0;
  for (int i=0; i<tr.size(); i++) 
      wtrac += tr[i].weight;	  
  for (int i=0; i<ran.size(); i++) 
      wrand += ran[i].weight;	
  float norm = wrand / wtrac;

  for (int ipix=0; ipix<hp.Npix(); ipix++) {
      map[ipix].delta = -1.0;
      if (!map[ipix].mask) continue;
      wtrac = 0.0;
      wrand = 0.0;
      for (int in=0; in<map[ipix].tracer.size(); in++) 
	  wtrac += tr[map[ipix].tracer[in]].weight;
      for (int in=0; in<map[ipix].random.size(); in++) 
	  wrand += ran[map[ipix].random[in]].weight;
      map[ipix].delta = (wtrac/wrand)*norm - 1.0;
  }

  // Calculo delta suavizada 
  float mean;
  int np;
  fix_arr <int,8> neigh;
  for (int ipix=0; ipix<hp.Npix(); ipix++) {
     if (!map[ipix].mask) {
	map[ipix].delta_smooth = map[ipix].delta;
        continue;	
     }	     
     mean = map[ipix].delta;
     np = 1;
     hp.neighbors(ipix,neigh);
     for (int in=0; in<8; in++) {
	 if (neigh[in] == -1) continue;    
	 if (!map[neigh[in]].mask) continue;
	 mean += map[neigh[in]].delta;
	 np++;
     }
     mean /= (float)np;
     map[ipix].delta_smooth = mean;
  }
  
}

