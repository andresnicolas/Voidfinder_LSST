
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "healpix_base.h"

#define DEG2RAD (M_PI/180.0)
#define RAD2DEG (180.0/M_PI)

using namespace std;

struct tracer {
   float phi;      /* RA en [rad] */
   float theta;    /* 90-DEC en [rad] */
   float redshift; 
   float weight;
};

struct hpmap {
   vector <int> tracer; 
   vector <int> random; 
   float delta;
   float delta_smooth;   
   bool mask;
};

struct global_data {
   int NumTracer;
   int NumRandom;   
   float ShellDist;
   float ShellThick;
   float PixelArea;
};
extern struct global_data G;
