/* cricci.h
 */

typedef enum {
  TOO_MANY_ITERATIONS = -100,
  LINE_SEARCH_FAILED = -99,
  RUNNING = 0,
  CONVERGENT = 1
} ricci_state;

/* Struct which encodes configuration information for the ricci
 * flow algorithm. 
 * relative_error, absolute_error are thresholds for convergence
 * wolfe_c1, wolfe_c2 are parameters for wolfe condition check
 * ds is the step size for the numerical integration
 * max_iterations is the maximum number of iterations of the algorithm
 * verbose: 0 for no messages, 1 for exit status, 2 for continual updates
 */
typedef struct {
  double relative_error;
  double absolute_error;
  double wolfe_c1;
  double wolfe_c2;
  int strong_wolfe;
  int integration_precision;
  int max_iterations;
  int max_line_steps;
  int cg_max_iterations;
  int cg_precon;
  double cg_tolerance;
  int verbose;
} ricci_config;


typedef struct {
  mesh *m;
  ricci_config *rc;
  double init_K_norm;
  double K_norm;
  double K_supnorm;
  double step_norm;
  double step_scale;
  int iteration;
  ricci_state status;
  double f;
  double f_next;
  double *s;
  double *s_next;
  double *K;
  double *K_next;
  double *step;
} ricci_solver;


void run_ricci_flow(mesh *m);
void initialize_ricci_solver(ricci_solver *r, mesh *m, ricci_config *rc);
void deallocate_ricci_solver(ricci_solver *r);
void calc_initial_variables(ricci_solver *r);
void calc_flat_metric(ricci_solver *r);
int calc_hessian_product(double *x, double *y, int n,
    void *instance);
int check_hessian_symmetry(ricci_solver *r);
ricci_state convergence_test(ricci_solver *r);
void update_hessian(ricci_solver *r);
void calc_next_step(ricci_solver *r);
void calc_line_search(ricci_solver *r);
void calc_next_s(ricci_solver *r);
int test_wolfe_conditions(ricci_solver *r);

void print_ricci_status(ricci_solver *r);

double sup_norm(const MV_Vector_double &x);
int vector_in_bounds(double *v, double min, double max, int n);
void clear_vector(double *v, int n);
void swap_vectors(double **v, double **w);
double dot_product(double *v, double *w, int n);
