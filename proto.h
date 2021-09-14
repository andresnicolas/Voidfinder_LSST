
/* tools */
FILE *safe_open(char *filename, const char *mode);
float random_number(void);

/* io */
void read_tracers(char *filename, vector <tracer> &tr);
void read_randoms(char *filename, vector <tracer> &ran);
void create_map(vector <tracer> &tr, vector <tracer> &ran, T_Healpix_Base<int> &hp, struct hpmap *map);

