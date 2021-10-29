
#include "global.h"
#include "io.h"
#include "map.h"

void map_load_tracers(vector <tracer> &tr, T_Healpix_Base<int> &hp, struct hpmap *map)
{

  for (int ipix=0; ipix<hp.Npix(); ipix++) 
      map[ipix].ntrac = 0;

  for (int i=0; i<tr.size(); i++) {
     int ipix = hp.ang2pix(tr[i].coord);
     map[ipix].ntrac++;   
  }
  
  for (int ipix=0; ipix<hp.Npix(); ipix++) { 
      map[ipix].tracer = (int *) malloc(map[ipix].ntrac*sizeof(int));
      map[ipix].ntrac = 0;
  }

  for (int i=0; i<tr.size(); i++) {
     int ipix = hp.ang2pix(tr[i].coord);
     map[ipix].tracer[map[ipix].ntrac] = i;
     map[ipix].ntrac++;   
  }

}

void map_load_randoms(vector <tracer> &ran, T_Healpix_Base<int> &hp, struct hpmap *map) 
{

  for (int ipix=0; ipix<hp.Npix(); ipix++) 
      map[ipix].nrand = 0;      

  for (int i=0; i<ran.size(); i++) {
     int ipix = hp.ang2pix(ran[i].coord);
     map[ipix].nrand++;   
  }
  
  for (int ipix=0; ipix<hp.Npix(); ipix++) {
      map[ipix].random = (int *) malloc(map[ipix].nrand*sizeof(int));
      map[ipix].nrand = 0;
  }

  for (int i=0; i<ran.size(); i++) {
     int ipix = hp.ang2pix(ran[i].coord);
     map[ipix].random[map[ipix].nrand] = i;
     map[ipix].nrand++;   
  }

}

void map_load_voids(vector <voids> &v, T_Healpix_Base<int> &hp, struct hpmap *map)
{

  for (int ipix=0; ipix<hp.Npix(); ipix++)  
      map[ipix].nvoid = 0;	

  for (int i=0; i<v.size(); i++) {
      int ipix = hp.ang2pix(v[i].coord);
      map[ipix].nvoid++;      
  }

  for (int ipix=0; ipix<hp.Npix(); ipix++) { 
      map[ipix].voids = (int *) malloc(map[ipix].nvoid*sizeof(int));
      map[ipix].nvoid = 0;
  }

  for (int i=0; i<v.size(); i++) {
     int ipix = hp.ang2pix(v[i].coord);
     map[ipix].voids[map[ipix].nvoid] = i;
     map[ipix].nvoid++;   
  }

}

void map_compute_delta(vector <tracer> &tr, vector <tracer> &ran, int npix, struct hpmap *map)
{
  
  float wtrac = 0.0;
  float wrand = 0.0;
  for (int i=0; i<tr.size(); i++) 
      wtrac += tr[i].weight;	  
  for (int i=0; i<ran.size(); i++) 
      wrand += ran[i].weight;	

  const float norm = wrand / wtrac;

  for (int ipix=0; ipix<npix; ipix++) {
      map[ipix].delta = -1.0;
      if (!map[ipix].mask) continue;
      wtrac = 0.0;
      wrand = 0.0;
      for (int in=0; in<map[ipix].ntrac; in++) 
	  wtrac += tr[map[ipix].tracer[in]].weight;
      for (int in=0; in<map[ipix].nrand; in++) 
	  wrand += ran[map[ipix].random[in]].weight;
      map[ipix].delta = (wtrac/wrand)*norm - 1.0;
  }

}

void map_load_mask(T_Healpix_Base<int> &hp, struct hpmap *map)
{

  for (int ipix=0; ipix<hp.Npix(); ipix++) {
      pointing p = hp.pix2ang(ipix);
      if (p.theta <= 0.5*M_PI && p.phi <= 0.5*M_PI)	
	 map[ipix].mask = true;
      else
	 map[ipix].mask = false;     
  }

}

