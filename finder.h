
/*
 * Finder functions prototypes 
 */

void find_candidates(float delta_cut, float rmax, T_Healpix_Base<int> &hp, struct hpmap *map, 
		     vector <tracer> &tr, vector <tracer> &ran, vector <voids> &v);
void find_voids(float delta_cut, T_Healpix_Base<int> &hp, struct hpmap *map, 
		vector <tracer> &tr, vector <tracer> &ran, vector <voids> &v);
void clean_voids(float tol, T_Healpix_Base<int> &hp, struct hpmap *map, 
		 vector <voids> &v);

