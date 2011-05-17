/* cconj_grad.c
 * Conjugate gradient solving algorithm
 */

#include <stdlib.h>
#include <stdio.h>
#include "cconj_grad.h"

/*
int B(double *x, double *y, int n, void* instance)
{
        int i, j;
        double M[4][4] = {
                4, 1, -1, 0.5,
                1, 3, -0.25, 0.2,
                -1, -0.25, 5, 0,
                0.5, 0.2, 0, 2
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
        double x[4] = {0, 0, 0, 0};
        double b[4];
        float tmp;

        if (argc != 5) {
                printf("  usage: %s b1 b2 b3 b4\n", argv[0]);
                return 0;
        }
        (double)sscanf(argv[1], "%f", &tmp);
        b[0] = tmp;
        (double)sscanf(argv[2], "%f", &tmp);
        b[1] = tmp;
        (double)sscanf(argv[3], "%f", &tmp);
        b[2] = tmp;
        (double)sscanf(argv[4], "%f", &tmp);
        b[3] = tmp;
        cg_solve(&B, x, b, 0.0001, 4, 5, NULL);
        printf("x = [%f, %f, %f, %f]\n", x[0], x[1], x[2], x[3]);
        return 0;
}
*/

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
void cg_solve(int (*A)(double *, double *, int, void *), 
                double *x, double *b, double tolerance, 
                int n, int max_iterations, void *instance)
{
        int i;
        double *r0, *r1, *p, *Ap;
        int k = 0;
        double alpha, beta;
        double r1_squared, r0_squared, initial_residue;
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
        r0_squared = dot(r0, r0, n);
        initial_residue = r0_squared;
        while (r0_squared > tolerance * initial_residue && k < max_iterations) {
                k++;
                alpha = r0_squared / dot(p, Ap, n);
                scale_and_add(x, x, alpha, p, n); // x = x + alpha p
                scale_and_add(r1, r0, -1*alpha, Ap, n); // r1=r0-alpha A(p)
                r1_squared = dot(r1, r1, n);
                if (r1_squared > tolerance * initial_residue) {
                        beta = r1_squared / r0_squared;
                        scale_and_add(p, r1, beta, p, n); // p= r1 + beta p
                        swap(&r0, &r1); // r1 becomes r0, freeing 
                                      // r1 for next run
                        A(p, Ap, n, instance); // Ap = A(p)
                }
                r0_squared = r1_squared;
                printf("CG iteration: %i  Error squared: %f\n", k, r1_squared);
        }
        free(Ap);
        free(p);
        free(r1);
        free(r0);
}

/* Same as cg_solve, except it used the Jacobi
 * preconditioner (scales A by the inverse of its
 * diagonal elements).
 *
 * Callback function inv_diag should compute the vector
 * of the inverses of the diagonal elements of A.
 */
void pccg_solve(int (*A)(double *, double *, int, void *),
               int (*inv_diag)(double *, int, void *), 
                double *x, double *b, double tolerance, 
                int n, int max_iterations, void *instance)
{
        int i;
        double *r0, *r1, *p, *Ap, *precon;
        int k = 0;
        double alpha, beta;
        double r1_squared, r0_squared, initial_residue;
        r0 = (double *)malloc(n * sizeof(double));
        r1 = (double *)malloc(n * sizeof(double));
        p = (double *)malloc(n * sizeof(double));
        Ap = (double *)malloc(n * sizeof(double));
        precon = (double *)malloc(n * sizeof(double));
        if(!r0 || !r1 || !p || !Ap || !precon) {
                printf("Memory allocation error.\n");
                exit(1);
        }
        
        initialize_preconditioner(precon, A, n, instance);
        A(x, Ap, n, instance); // Ap = A(x)
        scale_and_add(r0, b, -1, Ap, n); // r0 = b- A(x)
        copy_vector(p, r0, n); // p = r0
        A(p, Ap, n, instance);  // Ap = A(p);
        r0_squared = dot(r0, r0, n);
        initial_residue = r0_squared;
        while (r0_squared > tolerance * initial_residue && k < max_iterations) {
                k++;
                alpha = r0_squared / dot(p, Ap, n);
                scale_and_add(x, x, alpha, p, n); // x = x + alpha p
                scale_and_add(r1, r0, -1*alpha, Ap, n); // r1=r0-alpha A(p)
                r1_squared = dot(r1, r1, n);
                if (r1_squared > tolerance * initial_residue) {
                        beta = r1_squared / r0_squared;
                        scale_and_add(p, r1, beta, p, n); // p= r1 + beta p
                        swap(&r0, &r1); // r1 becomes r0, freeing 
                                      // r1 for next run
                        A(p, Ap, n, instance); // Ap = A(p)
                }
                r0_squared = r1_squared;
                printf("Iteration: %i  Error squared: %f, x = [", k, r1_squared);
                for (i=0; i<n; i++) {
                        printf("%f, ", x[i]);
                }
                printf("]\n");
        }
        free(Ap);
        free(p);
        free(r1);
        free(r0);
}

void initialize_preconditioner(double *precon, int (*A)(double *, double *, int, void *),  
                int n, void *instance) 
{
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
double dot(double *x, double *y, int n)
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

/* Copies vector y into vector x. Assumes that both
 * x and y are arrays of length n.
 */
void copy_vector(double *x, double *y, int n)
{
        int i;
        for (i=0; i<n; i++) {
                x[i] = y[i];
        }
}
