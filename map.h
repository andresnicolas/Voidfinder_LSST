
struct hpmap {
   int ntrac;
   int nrand;
   int nvoid;
   int *tracer;
   int *random;
   int *voids;   
   float delta;
   bool mask;
};

void map_load_tracers(vector <tracer> &tr, T_Healpix_Base<int> &hp, struct hpmap *map);
void map_load_randoms(vector <tracer> &ran, T_Healpix_Base<int> &hp, struct hpmap *map); 
void map_load_voids(vector <voids> &v, T_Healpix_Base<int> &hp, struct hpmap *map);
void map_compute_delta(vector <tracer> &tr, vector <tracer> &ran, int npix, struct hpmap *map);
void map_load_mask(T_Healpix_Base<int> &hp, struct hpmap *map);

