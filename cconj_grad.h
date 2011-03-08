/* cconj_grad.h
 */

int conjugate_gradient(int A(double *, double *, int), 
    double *x, double *b, double max_error, int n);
