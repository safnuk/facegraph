/* cconj_grad.c
 * Conjugate gradient solving algorithm
 */

#include <stdlib.h>
#include <stdio.h>
#include "cconj_grad.h"

/* Uses conjugate gradient method to solve the equation
 *      A x = b
 * Note that A is assumed to be symmetric and positive
 * definite.
 *
 * The answer is computed in place with x, and it is also used
 * as the initial guess at the answer - i.e. to be safe, it should
 * be initialized to 0 before calling. 
 */
int conjugate_gradient(int A(double *, double *, int), 
                double *x, double *b, double max_error, int n)
{
        double *r0, *r1, *p, *Ap;
        int k = 0;
        r0 = (double *)malloc(n * sizeof(double));
        r1 = (double *)malloc(n * sizeof(double));
        p = (double *)malloc(n * sizeof(double));
        Ap = (double *)malloc(n * sizeof(double));

        if(!r0 || !r1 || !p || !Ap) {
                printf("Memory allocation error.\n");
                exit(1);
        }

        free(Ap);
        free(p);
        free(r1);
        free(r0);
}
