// partition.h

struct vertex;

struct vertex_params {
  double angle;
  double length;
  double position;
};

double calc_optimal_partition_offset(vertex* v);
double calc_offset(double l, double alpha, double delta);
double my_f(const gsl_vector* v, void* params);
void my_df(const gsl_vector* v, void* params, gsl_vector* df);
void my_fdf(const gsl_vector* v, void* params, double* f, gsl_vector* df);
void calc_new_lengths(double alpha, double delta, vertex_params const* p, 
                      double* d, double* alphas);
double calc_new_length(double alpha, double delta, double length);
void calc_dl_dalpha(double const* alphas, double delta, 
    vertex_params const* p, double const* d, double* dl_dalpha);
void calc_dl_ddelta(double const* alphas, double delta, 
    vertex_params const* p, double const* d, double* dl_ddelta);
void calc_vertex_params(vertex const* v, vertex_params* p);
void calc_sort_order(vertex const* v, geodesic* g);
