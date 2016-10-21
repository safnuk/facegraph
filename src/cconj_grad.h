/* cconj_grad.h
 */
#ifndef ConjGrad
#define ConjGrad
 
void pccg_solve(int (*A)(double *, double *, int, void *),
                double *x, double *b, double tolerance, 
                int n, int max_iterations,
                int verbose, void *instance);
void cg_solve(int (*A)(double *, double *, int, void *), 
    double *x, double *b, double tolerance, int n, 
    int max_iterations, int verbose, void *instance);
void scale_and_add(double *x, double *b, double alpha, double *c, int n);
double dot(double *x, double *y, int n);
void swap(double **r0, double **r1);
void copy_vector(double *x, double *y, int n);
void component_product(double *x, double *y, double *z, int n);

#endif
