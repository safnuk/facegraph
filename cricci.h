/* cricci.h
 */

typedef enum {
  TOO_MANY_ITERATIONS = -1,
  RUNNING = 0,
  CONVERGENT = 1,
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
  double ds;
  int max_iterations;
  int cg_max_iterations;
  double cg_tolerance;
  int verbose;
} ricci_config;


typedef struct {
  mesh *m;
  ricci_config *rc;
  double init_K_norm;
  double K_norm;
  double step_norm;
  double step_scale;
  int iteration;
  ricci_state status;
  double *s;
  double *K;
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
void print_ricci_status(ricci_solver *r);

double sup_norm(double *v, int n);

