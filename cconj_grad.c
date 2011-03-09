/* cconj_grad.c
 * Conjugate gradient solving algorithm
 */

#include <stdlib.h>
#include <stdio.h>
#include "cconj_grad.h"

int B(double *x, double *y, int n, void* instance)
{
        int i, j;
        double M[2][2] = {
                2, 1,
                1, 3
        };
        for (i=0; i<n; i++) {
                y[i] = 0;
                for (j=0; j<n; j++) {
                        y[i] += M[i][j] * x[j];
                }
        }
        return 0;
}

int main(int argc, char *argv[])
{
        double x[2] = {0, 0};
        double b[2];

        if (argc != 3) {
                printf("  usage: %s x1 x2\n", argv[0]);
                return 0;
        }
        (double)sscanf(argv[1], "%f", &(x[0]));
        (double)sscanf(argv[2], "%f", &(x[1]));
        conjugate_gradient(&B, x, b, 0.0001, 2, NULL);
        return 0;
}
/* Uses conjugate gradient method to approximately solve the equation
 *      A x = b
 * Note that A is assumed to be symmetric and positive
 * definite.
 *
 * The answer is computed in-place, stored in x, and x is also used
 * as the initial guess at the answer - i.e. to be safe, it should
 * be initialized to the 0 vector before calling, but other initial
 * values may lead to faster convergence.
 *
 * The callback function A should be a linear operator 
 *      A: R^n -> R^n
 * with A(x, y, n, instance) computing A(x) and storing the result in y. 
 *
 * instance is a pointer to user-defined data which is used by A.
 */
int conjugate_gradient(int A(double *, double *, int, void *), 
                double *x, double *b, double max_error, 
                int n, void *instance)
{
        double *r0, *r1, *p, *Ap;
        int k = 0;
        double alpha, beta;
        error_squared = max_error * max_error;
        double r1_squared = error_squared + 1;
        double r0_squared;
        r0 = (double *)malloc(n * sizeof(double));
        r1 = (double *)malloc(n * sizeof(double));
        p = (double *)malloc(n * sizeof(double));
        Ap = (double *)malloc(n * sizeof(double));
        if(!r0 || !r1 || !p || !Ap) {
                printf("Memory allocation error.\n");
                exit(1);
        }
        
        A(x, Ap, n, instance); // Ap = A(x)
        scale_and_add(r0, b, -1, Ap, n); // r0 = b- A(x)
        copy_vector(p, r0, n); // p = r0
        A(p, Ap, n, instance);  // Ap = A(p);
        while (r1_squared > error_squared) {
                k++;
                r0_squared = dot(r0, r0, n); 
                alpha = r0_squared / dot(p, Ap, n);
                scale_and_add(x, x, alpha, p, n); // x = x + alpha p
                scale_and_add(r1, r0, -1*alpha, Ap, n); // r1=r0-alpha A(p)
                r1_squared = dot(r1, r1, n);
                if (r1_squared > max_error) {
                        beta = r1_squared / r0_squared;
                        scale_and_add(p, r1, beta, p, n); // p= r1 + beta p
                        swap(&r0, &r1); // r1 becomes r0, freeing 
                                      // r1 for next run
                        A(p, Ap, n, instance); // Ap = A(p)
                }
                printf("Iteration: %i  Error squared: %f\n", k, r1_squared);
        }
        free(Ap);
        free(p);
        free(r1);
        free(r0);
        return 0;
}

/* Calculates the vector sum
 *      b + alpha c
 * and stores the result in the array x. Note that x, b and c
 * need to be vectors of length n.
 */
void scale_and_add(double *x, double *b, double alpha, double *c, int n)
{
        int i;
        for (i=0; i<n; i++) {
                x[i] = b[i] + alpha * c[i];        
        }
}

/* Calculates the dot product of vectors x and y, and returns
 * the result. Note that x and y need to be arrays of length n.
 */
void dot(double *x, double *y, int n)
{
        int i;
        double sum = 0;
        for (i=0; i<n; i++) {
                sum += x[i] * y[i];
        }
        return sum;
}

/* Performs an in-place swap of the vectors r0 and r1.
 * Assumption is that r0 and r1 are arrays of doubles, i.e. 
 * pointers. Pointers-to-pointer addresses are passed, and after
 * being called r0 will point to the original r1 location, and
 * r1 will point to the original r0 location.
 */
void swap(double **ptr_r0, double **ptr_r1)
{
        double *tmp = *ptr_r0;
        *ptr_r0 = *ptr_r1;
        *ptr_r1 = tmp;
}
