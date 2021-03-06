
#include "global.h"
#include "io.h"
#include "map.h"
#include "finder.h"
#include "qsort.h"

/*
 * Search void candidates by checking all pixels with negative density contrast. 
 * The routine estimates a radius by counting tracers of all pixels completely 
 * inside. 
 */

void find_candidates(float delta_cut, float rmax, T_Healpix_Base<int> &hp, struct hpmap *map, 
		     vector <tracer> &tr, vector <tracer> &ran, vector <voids> &v)
{

  /// Effective radius of a pixel
  const float reff = 2.0 / sqrt((double)hp.Npix());
  
  /// Number of spherical bins
  const int nr = round(rmax/reff);

  /// Normalization to compute density contrast
  const float norm = (float)ran.size() / (float)tr.size(); 
  
  int i = 0;	
  
  #pragma omp parallel for \
              schedule(dynamic) \
              default(none) \
              shared(tr,ran,v,norm,nr,reff,map,hp,delta_cut,i) 

  for (int ipix=0; ipix<hp.Npix(); ipix++) {
      
      if (!map[ipix].mask || map[ipix].delta > 0.0) continue; 

      for (int ir=nr; ir>=0; ir--) { 

          float radius = reff * (float)(ir+1);
          pointing coord = hp.pix2ang(ipix);	  
          rangeset pix = hp.query_disc(coord,radius);
          vector <int> pixels_inside = pix.toVector();

          int ntrac = 0;
          int nrand = 0;  
          for (int i=0; i<pixels_inside.size(); i++) {
              int ipix = pixels_inside[i];	  
              if (!map[ipix].mask) continue;	  
              ntrac += map[ipix].ntrac;
              nrand += map[ipix].nrand;	      
          }

	  /// density contrast considering only tracers in pixels within radius
          float delta = norm * ((float)ntrac / (float)nrand) - 1.0;

	  /// Void candidate found
	  if (delta < delta_cut) {
             #pragma omp critical 
	     {
      	     v.push_back(voids());
             v[i].coord.theta = coord.theta;
             v[i].coord.phi = coord.phi;
 	     v[i].radius = radius;
	     v[i].delta = delta;
             v[i].ntrac = ntrac;	     
             v[i].nrand = nrand;	     
	     v[i].tof = false;
             i++;
	     }	     
	     break;
		   
	  } 
      } /// end loop ir
  } /// end loop ipix

}

/*
 * This routine takes the void candidates and defines the true void radius by
 * counting tracers in all pixels which overlap the previous radius estimation.
 */

void find_voids(float delta_cut, T_Healpix_Base<int> &hp, struct hpmap *map, 
		vector <tracer> &tr, vector <tracer> &ran, vector <voids> &v)
{

  struct sort *sort_gal,*sort_ran;	
  vector <double> dist_gal,dist_ran;

  const float norm = (float)ran.size() / (float)tr.size(); 

#pragma omp parallel for \
        schedule(dynamic) \
        default(none) \
        shared(stdout,tr,ran,v,norm,map,hp,delta_cut) \
        private(dist_gal,dist_ran,sort_gal,sort_ran) 

  for (int iv=0; iv<v.size(); iv++) {

      double theta1 = v[iv].coord.theta;	  
      double phi1 = v[iv].coord.phi;
      int ntrac = v[iv].ntrac;
      int nrand = v[iv].nrand; 
    
      rangeset pix1 = hp.query_disc(v[iv].coord,v[iv].radius);
      rangeset pix2 = hp.query_disc_inclusive(v[iv].coord,v[iv].radius,4);
      rangeset pixd = pix2.op_andnot(pix1);
      vector <int> pixels_ring = pixd.toVector();
    	  
      for (int i=0; i<pixels_ring.size(); i++) {
          
	  int ipix = pixels_ring[i];
          if (!map[ipix].mask) continue;
        
          if (map[ipix].ntrac > 0) {
             for (int j=0; j<map[ipix].ntrac; j++) {
                 int id = map[ipix].tracer[j];
                 double theta2 = tr[id].coord.theta;	  
                 double phi2 = tr[id].coord.phi;
                 double dist = acos(cos(theta1) * cos(theta2) 
                             + sin(theta1) * sin(theta2) * cos(phi1 - phi2));
                 dist_gal.push_back(dist);
             }
          }

          if (map[ipix].nrand > 0) {
             for (int j=0; j<map[ipix].nrand; j++) {
                 int id = map[ipix].random[j];
                 double theta2 = ran[id].coord.theta;	  
                 double phi2 = ran[id].coord.phi;
                 double dist = acos(cos(theta1) * cos(theta2) 
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
	  
	  int kt = k + 1;
          int kr = 0;
          for (int kk=0; kk<dist_ran.size(); kk++) { 
      	      if (sort_ran[kk].val < radius) 
      	         kr = kk; 
              else
     	         break;
          }	
         
          float delta = norm * ((float)(ntrac+kt) / (float)(nrand+kr)) - 1.0;
    
          if (delta < delta_cut) {
             v[iv].radius = radius;
             v[iv].delta = delta;	       
             v[iv].ntrac = ntrac+kt;
             v[iv].nrand = nrand+kr;
             v[iv].tof = true;	    
          } 
      }


      dist_gal.clear();
      dist_ran.clear();
      free(sort_gal);
      free(sort_ran);

  }	  

}

/*
 * Clean voids by superposition. 
 * */

void clean_voids(float tol, T_Healpix_Base<int> &hp, struct hpmap *map, 
		 vector <voids> &v)
{

  int nv = 0;
  for (int i=0; i<v.size(); i++) 
      if (v[i].tof) nv++;

  struct sort *sort_void = (struct sort *) malloc(nv*sizeof(struct sort));

  nv = 0;
  for (int i=0; i<v.size(); i++) {
      if (v[i].tof) {
	 sort_void[nv].val = v[i].radius;
         sort_void[nv].ord = i;
         nv++;	 
      }
  }

  QSort(sort_void,0,nv-1);

  for (int i=nv-1; i>=0; i--) {

      int iv1 = sort_void[i].ord;	  
      if (!v[iv1].tof) continue;

      double theta1 = v[iv1].coord.theta;
      double phi1 = v[iv1].coord.phi;
      double rad1 = v[iv1].radius;
      rangeset pix = hp.query_disc_inclusive(v[iv1].coord,2.0*rad1,4);
      vector <int> pixels = pix.toVector();

      for (int ip=0; ip<pixels.size(); ip++) {
          int ipix = pixels[ip];
	  if (!map[ipix].mask) continue;
	         
	  if (map[ipix].nvoid > 0) {
	     for (int j=0; j<map[ipix].nvoid; j++) {
	         int iv2 = map[ipix].voids[j];
		 if (!v[iv2].tof || iv1 == iv2) continue;
	         double theta2 = v[iv2].coord.theta;	  
                 double phi2 = v[iv2].coord.phi;
		 double rad2 = v[iv2].radius;
	         double dist = acos(cos(theta1) * cos(theta2) 
	                     + sin(theta1) * sin(theta2) * cos(phi1 - phi2));
		 
		 if (dist < (1.0 - tol) * (rad1 + rad2)) v[iv2].tof = false;
	     }
	  }
      }
  }
     
  free(sort_void);

}

