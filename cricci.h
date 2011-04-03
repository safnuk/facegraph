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
  lbfgsfloatval_t *u;
  lbfgs_parameter_t param;
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

static lbfgsfloatval_t ricci_evaluate(
    void *instance,
    const lbfgsfloatval_t *u,
    lbfgsfloatval_t *K,
    const int n,
    const lbfgsfloatval_t step
    );

static int ricci_progress(
    void *instance,
    const lbfgsfloatval_t *u,
    const lbfgsfloatval_t *K,
    const lbfgsfloatval_t f_of_u,
    const lbfgsfloatval_t u_norm,
    const lbfgsfloatval_t K_norm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    );
