
#include "global.h"
#include "proto.h"

void read_tracers(char *filename, vector <tracer> &tr)
{
  	
  FILE *f = safe_open(filename,"r");	
  int i = 0;
  float ra,dec,redshift,distance,x,y,z;

  while (fscanf(f,"%f %f %f %f %f %f %f",&ra,&dec,&redshift,&distance,&x,&y,&z) > 0) {
     tr.push_back(tracer());
     tr[i].phi = ra*DEG2RAD;
     tr[i].theta = (90.0 - dec)*DEG2RAD;
     tr[i].redshift = redshift;
     tr[i].weight = 1.0;
     i++;
  };

  fclose(f);

  G.NumTracer = tr.size();
  G.ShellDist = 950.0;
  G.ShellThick = 100.0;

}

void create_randoms(vector <tracer> &ran)
{
 
  float u;
  float R1 = G.ShellDist - 0.5*G.ShellThick; 
  float R2 = G.ShellDist + 0.5*G.ShellThick; 
  
  G.NumRandom = 30 * G.NumTracer;

  for (int i=0; i<G.NumRandom; i++) {
     ran.push_back(tracer());
     
     u = random_number();
     ran[i].phi = 2.0 * M_PI * u;

     u = random_number();
     ran[i].theta = acos(1.0 - 2.0 * u);
     
     u = random_number();
     ran[i].redshift = cbrt(u * pow(R2,3) + (1.0 - u) * pow(R1,3));

     ran[i].weight = 1.0;
  }  

}

void create_map(vector <tracer> &tr, vector <tracer> &ran, T_Healpix_Base<int> &hp, struct hpmap *map) 
{

  int i,ipix;
  pointing p;

  // Mascara angular
  for (ipix=0; ipix<hp.Npix(); ipix++) {
      p = hp.pix2ang(ipix);
      if (p.theta <= 0.5*M_PI && p.phi <= 0.5*M_PI)	
	 map[ipix].mask = true;
      else
	 map[ipix].mask = false;     
  }

  // Cargo tracers en mapa angular	
  for (i=0; i<tr.size(); i++) {
     p.theta = tr[i].theta;
     p.phi = tr[i].phi;
     ipix = hp.ang2pix(p);
     map[ipix].tracer.push_back(i);   
  }

  // Cargo randoms en mapa angular
  for (i=0; i<ran.size(); i++) {
     p.theta = ran[i].theta;
     p.phi = ran[i].phi;
     ipix = hp.ang2pix(p);
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

  for (ipix=0; ipix<hp.Npix(); ipix++) {
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
  for (ipix=0; ipix<hp.Npix(); ipix++) {
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
