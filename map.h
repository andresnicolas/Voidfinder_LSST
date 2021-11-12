
/// HealPix map 

struct hpmap {
   int ntrac;    /// number of tracers 
   int nrand;    /// number of randoms
   int nvoid;    /// number of void centers
   int *tracer;  /// tracers IDs
   int *random;  /// randoms IDs
   int *voids;   /// voids IDs
   float delta;  /// density contrast
   bool mask;    /// angular mask
};

/*
 * Prototypes for HealPix map routines 
 */

void map_load_tracers(vector <tracer> &tr, T_Healpix_Base<int> &hp, struct hpmap *map);
void map_load_randoms(vector <tracer> &ran, T_Healpix_Base<int> &hp, struct hpmap *map); 
void map_load_voids(vector <voids> &v, T_Healpix_Base<int> &hp, struct hpmap *map);
void map_compute_delta(vector <tracer> &tr, vector <tracer> &ran, int npix, struct hpmap *map);
void map_load_mask(T_Healpix_Base<int> &hp, struct hpmap *map);

