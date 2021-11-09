
/*
 * Structure with tracers/randoms information
 */

struct tracer {
   pointing coord; /// angular coordinated (phi,theta)
   float weight;   
};

/*
 * Structure with voids information
 */

struct voids {
   pointing coord; /// angular coordinated (phi,theta)
   float radius;   /// angular radius in [rad]
   float delta;    /// integrated density contrast at Rvoid
   bool tof;       /// internal flag to dismiss voids candidates
};

/*
 * Intpu/output functions prototypes
 */

FILE *safe_open(char *filename, const char *mode);
void read_tracers(char *filename, vector <tracer> &tr);
void read_randoms(char *filename, vector <tracer> &ran);
void write_voids(char *filename, vector <voids> &v);
