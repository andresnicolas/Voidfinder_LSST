
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

  struct sort *sort_gal,*sort_ran;	
  vector <float> dist_gal,dist_ran;

  float reff = 2.0 / sqrt((double)hp.Npix());
  float norm = (float)ran.size() / (float)tr.size(); 

  fprintf(stdout,"size = %d \n",v.size());
  getchar();

  for (int iv=0; iv<v.size(); iv++) {

      float theta1 = v[iv].coord_init.theta;	  
      float phi1 = v[iv].coord_init.phi;	  

      for (int ir=0; ir<10; ir++) { 

	  float radius = reff * (float)(ir+1);    
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

	     float r1 = reff * (float)ir;    
	     float r2 = reff * (float)(ir+1);
	     rangeset pix1 = hp.query_disc(v[iv].coord_init,r1);
	     rangeset pix2 = hp.query_disc_inclusive(v[iv].coord_init,r2,4);
	     rangeset pixd = pix2.op_andnot(pix1);
             vector <int> pixels_inside = pix1.toVector();
             vector <int> pixels_ring = pixd.toVector();
    	  
	     int ntrac = 0;
             int nrand = 0;  
             for (int i=0; i<pixels_inside.size(); i++) {
	         int ipix = pixels_inside[i];	  
	         if (!map[ipix].mask) continue;	  
	         ntrac += map[ipix].tracer.size();
                 nrand += map[ipix].random.size();	      
	     }

	     for (int i=0; i<pixels_ring.size(); i++) {
	         int ipix = pixels_ring[i];
	         if (!map[ipix].mask) continue;
	        
	         if (map[ipix].tracer.size() > 0) {
	            for (int j=0; j<map[ipix].tracer.size(); j++) {
	                int id = map[ipix].tracer[j];
	                float theta2 = tr[id].coord.theta;	  
                        float phi2 = tr[id].coord.phi;
	                float dist = acos(cos(theta1) * cos(theta2) 
	                           + sin(theta1) * sin(theta2) * cos(phi1 - phi2));
	                dist_gal.push_back(dist);
	            }
	         }

	         if (map[ipix].random.size() > 0) {
	            for (int j=0; j<map[ipix].random.size(); j++) {
	                int id = map[ipix].random[j];
	                float theta2 = ran[id].coord.theta;	  
                        float phi2 = ran[id].coord.phi;
	                float dist = acos(cos(theta1) * cos(theta2) 
	                           + sin(theta1) * sin(theta2) * cos(phi1 - phi2));
	                dist_ran.push_back(dist);
	            }
	         }
	     }

	     sort_gal = (struct sort *) malloc(dist_gal.size()*sizeof(struct sort));
	     sort_ran = (struct sort *) malloc(dist_ran.size()*sizeof(struct sort));
	     for (int k=0; k<dist_gal.size(); k++) {
	         sort_gal[k].val = dist_gal[k];
	         sort_gal[k].ord = k;	    
	     }
	     for (int k=0; k<dist_ran.size(); k++) {
	         sort_ran[k].val = dist_ran[k];
	         sort_ran[k].ord = k;	    
	     }
             QSort(sort_gal,0,dist_gal.size()-1);
             QSort(sort_ran,0,dist_ran.size()-1);

	     for (int k=0; k<dist_gal.size()-1; k++) {
	         float radius = 0.5*(sort_gal[k].val + sort_gal[k+1].val);
	         int kr = -1;
	         for (int kk=0; kk<dist_ran.size(); kk++) { 
	     	     if (sort_ran[kk].val < radius) 
	     	        kr = kk; 
	             else
	     	        break;
	         }	
	         if (kr == -1) continue;
	          
	         float delta = norm * ((float)(ntrac + k + 1) / (float)(nrand + kr + 1)) - 1.0;

	         if (delta < delta_cut) {
	            v[iv].radius = radius;
	            v[iv].delta = delta;	       
	            v[iv].tof = true;	    
	         } 
	     }

	     dist_gal.clear();
	     dist_ran.clear();
	     free(sort_gal);
	     free(sort_ran);

	     break;
	  }
      }
      
      if (v[iv].tof) fprintf(stdout,"%d %f %f \n",iv,v[iv].radius*G.ShellDist,v[iv].delta);
  }	  

}
