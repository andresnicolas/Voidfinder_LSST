
#include "global.h"
#include "proto.h"

void find_centers(float delta_seed, T_Healpix_Base<int> hp, struct hpmap *map, vector <voids> v)
{

  int i = 0;	
  fix_arr <int,8> neigh;
  bool center;
  for (int ipix=0; ipix<hp.Npix(); ipix++) {
      if (!map[ipix].mask) continue;	  
      if (map[ipix].delta_smooth < delta_seed) {
	 center = true;     
         hp.neighbors(ipix,neigh);
	 for (int in=0; in<8; in++) {
             if (neigh[in] == -1) continue;
	     if (!map[neigh[in]].mask) continue;
             if (map[neigh[in]].delta_smooth < map[ipix].delta_smooth)
	        center = false;	     
	 }
	 if (!center) continue;     
     	 v.push_back(voids());
         v[i].coord_init = hp.pix2ang(ipix);
	 i++;
      }
  }

}

void find_voids(float delta_cut, T_Healpix_Base<int> hp, struct hpmap *map, vector <voids> v)
{

  float rmin = 5.0;
  float rmax = 100.0;
  float rbin = 5.0;  

  for (int iv; iv<v.size(); iv++) {

      	  
      	  
	  
  }	  

}
