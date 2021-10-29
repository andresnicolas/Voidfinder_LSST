
struct tracer {
   pointing coord;
   float weight;
};

struct voids {
   pointing coord;
   float radius;
   float delta;
   bool tof;
};

FILE *safe_open(char *filename, const char *mode);
void read_tracers(char *filename, vector <tracer> &tr);
void read_randoms(char *filename, vector <tracer> &ran);
void write_voids(char *filename, vector <voids> &v);
