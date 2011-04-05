/* cricci.h
 */

typedef struct {
  double f_of_u;
  double u_norm;
  double K_norm;
  double step;
  int iterations;
} ricci_state;

typedef struct {
  mesh *m;
  ricci_state current_state;
  double *s;
  double ds;
  int verbose;
} ricci_solver;


void run_ricci_flow(mesh *m);
void initialize_ricci_solver(ricci_solver *r, mesh *m);
void deallocate_ricci_solver(ricci_solver *r);
void calc_initial_variables(ricci_solver *r);
void calc_flat_metric(ricci_solver *r);
int calc_hessian_product(double *x, double *y, int n,
    void *instance);

