/* cconj_grad.h
 */

void conjugate_gradient(int A(double *, double *, int, void *), 
    double *x, double *b, double max_error, int n, void *instance);
void scale_and_add(double *x, double *b, double alpha, double *c, int n);
double dot(double *x, double *y, int n);
void swap(double **r0, double **r1);
