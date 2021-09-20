
/* tools */
FILE *safe_open(char *filename, const char *mode);
float random_number(void);

/* io */
void read_tracers(char *filename, vector <tracer> &tr);
void read_randoms(char *filename, vector <tracer> &ran);
void create_map(vector <tracer> &tr, vector <tracer> &ran, T_Healpix_Base<int> &hp, struct hpmap *map);

/* finder */
void find_centers(float delta_seed, T_Healpix_Base<int> &hp, struct hpmap *map, vector <voids> &v);
void find_voids(float delta_cut, T_Healpix_Base<int> &hp, struct hpmap *map, 
		vector <tracer> &tr, vector <tracer> &ran, vector <voids> &v);

/* sorting */
void QSort(struct sort *a, int start, int end);
