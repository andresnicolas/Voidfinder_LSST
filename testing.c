
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "healpix_base.h"

#define deg2rad (M_PI/180.0)
#define rad2deg (180.0/M_PI)

using namespace std;

int main()
{
  struct tracer {
     float alpha;
     float delta;
     float redshift;
     float distance;
     float pos[3];
  }; 
  vector <tracer> H;
  int i,Ntracer;
  float alpha,delta,redshift,distance,x,y,z;
  FILE *fd; 
  char filename[200];

  // Read tracers

  sprintf(filename,"data/halos_shell.dat");
  fd = fopen(filename,"r");

  i = 0;
  while (fscanf(fd,"%f %f %f %f %f %f %f",&alpha,&delta,&redshift,&distance,&x,&y,&z) > 0) {
     H.push_back(tracer());
     H[i].alpha = alpha;
     H[i].delta = delta;
     H[i].redshift = redshift;
     i++;
  };
  Ntracer = H.size(); 
  //fprintf(stdout,"Number of tracers: %d \n",Ntracer);

  // HealPix map
  int nside = 128;
  pointing ptg;
  double radius;
  rangeset <int> ipix_1,ipix_2,ipix_diff;
  vector <int> ipix_in, ipix_ring;
  
  T_Healpix_Base pixmap(nside,RING,SET_NSIDE);
  //fprintf(stdout,"order = %ld | nside = %ld | npix = %ld \n",pixmap.Order(),pixmap.Nside(),pixmap.Npix());

  ptg.theta = deg2rad*20.0;
  ptg.phi = deg2rad*45.0;
  radius = deg2rad*5.0;

  pixmap.query_disc(ptg,radius,ipix_1);
  pixmap.query_disc_inclusive(ptg,radius,ipix_2,1);
  ipix_diff = ipix_2.op_andnot(ipix_1); 

  ipix_in = ipix_1.toVector();
  ipix_ring = ipix_diff.toVector();

  ipix_1.clear();
  ipix_2.clear();
  ipix_diff.clear();

  for (i=0; i<ipix_in.size(); i++) {
      int j=ipix_in[i];
      ptg = pixmap.pix2ang(j);
      x = cos(ptg.phi)*sin(ptg.theta);
      y = sin(ptg.phi)*sin(ptg.theta);
      z = cos(ptg.theta);           
      fprintf(stdout,"%f %f %f\n",x,y,z);
  }
//  for (i=0; i<ipix_ring.size(); i++) {
//      int j=ipix_ring[i];
//      ptg = pixmap.pix2ang(j);
//      x = cos(ptg.phi)*sin(ptg.theta);
//      y = sin(ptg.phi)*sin(ptg.theta);
//      z = cos(ptg.theta);           
//      fprintf(stdout,"%f %f %f\n",x,y,z);
//  }

}
