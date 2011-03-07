// ccirclepack.h

typedef struct {
  double error;
  double r_norm;
  double g_norm;
  double step;
  int iterations;
} state;

typedef struct {
  mesh *m;
  state current_state;
  lbfgsfloatval_t *ambient_lengths;
  lbfgsfloatval_t *r;
  lbfgs_parameter_t param;
  int verbose;
} optimizer;

double calc_circlepack_metric(mesh *m);
void calc_ambient_edge_lengths(optimizer *opt);
void calc_initial_radii(optimizer *opt);
void initialize_optimizer(optimizer *opt, mesh *m);
void deallocate_optimizer(optimizer *opt);
void calc_optimum_radii(optimizer *opt);
void record_output_in_mesh(optimizer *opt);
void copy_radii(optimizer *opt);
void calc_circle_angles(optimizer *opt);

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *r,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    );

static int progress(
    void *instance,
    const lbfgsfloatval_t *r,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t f_of_r,
    const lbfgsfloatval_t r_norm,
    const lbfgsfloatval_t g_norm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    );

static int quiet_progress(
    void *instance,
    const lbfgsfloatval_t *r,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t f_of_r,
    const lbfgsfloatval_t r_norm,
    const lbfgsfloatval_t g_norm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    );
