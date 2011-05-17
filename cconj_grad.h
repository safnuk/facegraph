/* cconj_grad.h
 */
 
void pccg_solve(int (*A)(double *, double *, int, void *), 
    int (*inv_diag)(double *, int, void *),
    double *x, double *b, double tolerance, int n, 
    int max_iterations, void *instance);
void initialize_preconditioner(double *precon, 
    int (*A)(double *, double *, int, void *),  
    int n, void *instance);
void cg_solve(int (*A)(double *, double *, int, void *), 
    double *x, double *b, double tolerance, int n, 
    int max_iterations, void *instance);
void scale_and_add(double *x, double *b, double alpha, double *c, int n);
double dot(double *x, double *y, int n);
void swap(double **r0, double **r1);
void copy_vector(double *x, double *y, int n);
