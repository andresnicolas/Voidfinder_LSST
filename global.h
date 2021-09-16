
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "healpix_base.h"

#define DEG2RAD (M_PI/180.0)
#define RAD2DEG (180.0/M_PI)

using namespace std;

struct tracer {
   pointing coord;
   float redshift; 
   float weight;
};

struct voids {
   pointing coord_init;	
   pointing coord;
   float radius;
   bool tof;
};

struct hpmap {
   vector <int> tracer; 
   vector <int> random; 
   float delta;
   float delta_smooth;   
   bool mask;
};

struct global_data {
   int TracPerPix;
   float ShellDist;
   float ShellThick;
   float SurveyArea;
};
extern struct global_data G;
