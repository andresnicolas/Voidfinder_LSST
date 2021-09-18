
#include "global.h"
#include "proto.h"

void find_centers(float delta_seed, T_Healpix_Base<int> &hp, struct hpmap *map, vector <voids> &v)
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

void find_voids(float delta_cut, T_Healpix_Base<int> &hp, struct hpmap *map, 
		vector <tracer> &tr, vector <tracer> &ran, vector <voids> &v)
{

  float rmin = 2.0 * G.ShellDist / sqrt((double)hp.Npix());
  float rmax = 10.0 * rmin;
  float rbin = 2.0;  

  float norm = (float)ran.size() / (float)tr.size(); 

  for (int iv; iv<v.size(); iv++) {

      for (int ir=0; ir<10; ir++) {

	  float radius = (rmin + rbin * (float)ir)/G.ShellDist;    
          rangeset pix = hp.query_disc(v[iv].coord_init,radius);
          vector <int> pixels_inside = pix.toVector();

    	  int ntrac = 0;
          int nrand = 0;  
          for (int i=0; i<pixels_inside.size(); i++) {
	      int ipix = pixels_inside[i];	  
	      if (!map[ipix].mask) continue;	  
	      ntrac += map[ipix].tracer.size();
              nrand += map[ipix].random.size();	      
	  }

	  float delta = norm * ((float)ntrac / (float)nrand) - 1.0;

	  if (delta > delta_cut) {

	     if (ir > 0) { 	  

	        float radius_back = (rmin + rbin * (float)(ir-1))/G.ShellDist;	  
                rangeset pix_back = hp.query_disc(v[iv].coord_init,radius_back);
                vector <int> pixels_inside_back = pix_back.toVector();

		rangeset pix_ring = pix.op_andnot(pix_back);
		vector <int> pixels_ring = pix_ring.toVector();

		for (int i=0; i<pixels_ring.size(); i++) {
	            int ipix = pixels_ring[i]; 		
		    if (!map[ipix].mask) continue;	  
	            ntrac -= map[ipix].tracer.size();
                    nrand -= map[ipix].random.size();	 
		}
	 
	       	float delta = norm * ((float)ntrac / (float)nrand) - 1.0;

		fprintf(stdout," --- > d = %f | nt = %d | nr = %d \n",delta,ntrac,nrand);
		getchar();

	     } else {
	         // cruce la delta en el primer paso nomas   
	     }

             break;		  
	  }
      }
  }	  

}
